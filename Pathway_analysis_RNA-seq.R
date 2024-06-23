setwd("/home/ubuntu/RNA-seq")

fora <- function(pathways, genes, universe, minSize=1, maxSize=length(universe)-1) {
  # Error if pathways is not a list
  if (!is.list(pathways)) {
    stop("pathways should be a list with each element containing genes from the universe")
  }
  
  # Warning message for duplicate gene names
  if (any(duplicated(universe))) {
    warning("There were duplicate genes in universe, they were collapsed")
    universe <- unique(universe)
  }
  
  minSize <- max(minSize, 1)
  
  pathwaysFiltered <- lapply(pathways, function(p) { unique(na.omit(match(p, universe))) })
  pathwaysSizes <- sapply(pathwaysFiltered, length)
  
  toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
  
  if (length(toKeep) == 0){
    return(data.table(pathway=character(),
                      pval=numeric(),
                      padj=numeric(),
                      overlap=integer(),
                      size=integer(),
                      overlapGenes=list()))
  }
  
  pathwaysFiltered <- pathwaysFiltered[toKeep]
  pathwaysSizes <- pathwaysSizes[toKeep]
  
  
  if (!all(genes %in% universe)) {
    warning("Not all of the input genes belong to the universe, such genes were removed")
  }
  genesFiltered <- unique(na.omit(match(genes, universe)))
  
  overlaps <- lapply(pathwaysFiltered, intersect, genesFiltered)
  
  overlapGenes <- lapply(overlaps, function(x) universe[x])
  
  overlapsT <- data.table(
    q=sapply(overlaps, length),
    m=sapply(pathwaysFiltered, length),
    n=length(universe)-sapply(pathwaysFiltered, length),
    k=length(genesFiltered))
  
  # q-1 because we want probability of having >=q white balls
  pathways.pvals <- with(overlapsT,
                         phyper(q-1, m, n, k, lower.tail = FALSE))
  
  res <- data.table(pathway=names(pathwaysFiltered),
                    pval=pathways.pvals,
                    padj=p.adjust(pathways.pvals, method="BH"),
                    overlap=overlapsT$q,
                    size=overlapsT$m,
                    overlapGenes=overlapGenes)
  
  res <- res[order(pval),]
  res
}

collapsePathwaysORA <- function(foraRes,
                                pathways,
                                genes,
                                universe,
                                pval.threshold=0.05) {
  
  pathways <- pathways[foraRes$pathway]
  pathways <- lapply(pathways, intersect, universe)
  
  parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))
  
  for (i in seq_along(pathways)) {
    p <- names(pathways)[i]
    if (!is.na(parentPathways[p])) {
      next
    }
    
    pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), p)
    
    if (length(pathwaysToCheck) == 0) {
      break
    }
    
    minPval <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)
    
    u1 <- setdiff(universe, pathways[[p]])
    foraRes1 <- fora(pathways = pathways[pathwaysToCheck],
                     genes=intersect(genes, u1),
                     universe=u1,
                     maxSize=length(u1)-1)
    minPval[foraRes1$pathway] <- pmin(minPval[foraRes1$pathway], foraRes1$pval)
    
    u2 <- pathways[[p]]
    foraRes2 <- fora(pathways = pathways[pathwaysToCheck],
                     genes=intersect(genes, u2),
                     universe=u2,
                     maxSize=length(u2)-1)
    minPval[foraRes2$pathway] <- pmin(minPval[foraRes2$pathway], foraRes2$pval)
    
    parentPathways[names(which(minPval > pval.threshold))] <- p
  }
  
  return(list(mainPathways=names(which(is.na(parentPathways))),
              parentPathways=parentPathways))
}

## Read deg data
deg <- read.csv("DEG_results_Primary_Tumor vs. Solid_Tissue_Normal.csv",
                 header = TRUE, sep = ",")

names(deg)[names(deg) == "X"] <- "genes"
# Connect to the Ensembl BioMart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query the database to convert gene symbols to Ensembl IDs
gene_conversion <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = deg$genes,
  mart = ensembl
)

deg_data <- merge(deg,gene_conversion, by.x = "genes", by.y="hgnc_symbol")
#------------------------------------------------Pathway Enrichment ------------------------------------------------
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(tidyr)

universe=deg_data %>% 
  filter(!is.na(adj.P.Val)) %>% pull(ensembl_gene_id)
mapping = AnnotationDbi::select(org.Hs.eg.db, universe, 'ENTREZID', 'ENSEMBL')

ranks_df <- deg_data %>% 
  filter(adj.P.Val < 0.05, !is.na(adj.P.Val)) %>% # only if enough DE genes  
  arrange((logFC)) %>% 
  left_join(mapping, by = c("ensembl_gene_id"="ENSEMBL")) %>% 
  filter(!is.na(ENTREZID) & !is.na(adj.P.Val))
ranks <- ranks_df$logFC
names(ranks) <- ranks_df$ENTREZID

all_in_life=list(
  msigdbr(species = "human", category = "H") %>% mutate(gs_subcat="Hallmark"),
  msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME"),
  msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG"),
  msigdbr(species = "human", category = "C2", subcategory = "CP:PID"),
  msigdbr(species = "human", category = "C5", subcategory = "GO:BP"),
  msigdbr(species = "human", category = "C5", subcategory = "GO:MF"),
  msigdbr(species = "human", category = "C5", subcategory = "HPO"),
  msigdbr(species = "human", category = "C3", subcategory = "TFT:GTRD"),
  msigdbr(species = "human", category = "C6") %>% mutate(gs_subcat="Oncogenic")
)


## ------------------------------------------------ gsea ------------------------------------------------
pathways_gsea_all = lapply(all_in_life, function(p){
  pathway = split(x = p$entrez_gene, f = p$gs_name)
  db_name = paste(p$gs_cat[1], p$gs_subcat[1],sep=":")
  respath <- fgsea(pathways = pathway, 
                   stats = ranks,
                   minSize  = 15,
                   maxSize  = 500)
  
  coll_respath = collapsePathways(respath[order(pval)][padj < 0.1], 
                                  pathway, ranks)
  as_tibble(respath[pathway %in% coll_respath$mainPathways])  %>% 
    mutate(database=db_name)
}) %>% bind_rows() %>% 
  mutate(analysis="GSEA")

gsea_tb=pathways_gsea_all %>% unnest(leadingEdge) %>%
  group_by(pathway) %>% 
  left_join(mapping, by =c("leadingEdge"="ENTREZID")) %>% 
  dplyr::select(pathway, padj, NES, ENSEMBL, analysis, database)

## ------------------------------------------------ ORA  ------------------------------------------------
# check cutoff
ora_input = deg_data %>% filter(!is.na(adj.P.Val), adj.P.Val<0.05, abs(logFC)>0.3) %>% pull(ensembl_gene_id)
input_entrezid <- AnnotationDbi::select(org.Hs.eg.db, ora_input, 'ENSEMBL', columns = c('ENTREZID', 'SYMBOL'))

total_deg=length(unique(ora_input))/length(unique(mapping$ENTREZID))

pathways_ora_all = lapply(all_in_life, function(p){
  pathway = split(x = p$entrez_gene, f = p$gs_name)
  db_name = paste(p$gs_cat[1], p$gs_subcat[1],sep=":")
  respath <- fora(pathways = pathway, 
                  genes = unique(input_entrezid$ENTREZID),
                  universe = unique(mapping$ENTREZID),
                  minSize  = 15,
                  maxSize  = 500)
  coll_respath = collapsePathwaysORA(respath[order(pval)][padj < 0.1], 
                                     pathway, unique(input_entrezid$ENTREZID), unique(mapping$ENTREZID))
  as_tibble(respath[pathway %in% coll_respath$mainPathways])  %>% 
    mutate(database=db_name, NES=(overlap/size)/(total_deg))
}) %>% bind_rows() %>% 
  mutate(analysis="ORA")

ora_tb = pathways_ora_all %>% unnest(overlapGenes) %>%
  group_by(pathway) %>% 
  left_join(mapping, by =c("overlapGenes"="ENTREZID")) %>% 
  dplyr::select(pathway, padj, NES, ENSEMBL, analysis,
                database)

pathways_long = bind_rows(gsea_tb, ora_tb)

###---------------------------annotate with gene_type --------------------------------------------------------
GTF <- "gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
gtf_df <- as.data.frame(rtracklayer::import(GTF))
gtf_df$gene_id <- sub("\\..*","\\", gtf_df$gene_id) # remove the ensembl numbering
tx2gene <- gtf_df %>% dplyr::select(transcript_id, gene_id) %>%
  filter(!is.na(transcript_id)) %>% 
  dplyr::distinct()

names(deg_data)[names(deg_data) == "ensembl_gene_id"] <- "gene_id"

res_ <- dplyr::inner_join(deg_data,gtf_df%>% dplyr::select(gene_id),by = "gene_id") %>% distinct()

pathways_long_ <- dplyr::inner_join(pathways_long,gtf_df%>% dplyr::select(gene_id,gene_name),by = c("ENSEMBL" = "gene_id")) %>% distinct()


## ----------------------------Match the gene panel ------------------------------
genes_PAAD <- sort(c(
  'TFAP2E', 'ASXL2', 'GSTM1', 'ONECUT2', 'TET1', 'HSPA5', 'SNRPF', 'PTGES2', 'CDA', 
  'SPARC', 'CES2', 'TSPAN1', 'ABCG2', 'BRCA1', 'CADM1', 'ABCB11', 'ABCB4', 'ABCC', 
  'ABCC10', 'ABCC3', 'ABCC5', 'SLC29A1', 'BCL2L1', 'ABCC2', 'FOXM1', 'TNFSF10', 
  'TM4SF1', 'DPYD', 'HMGA1', 'ABCB1', 'LCN2', 'SFN', 'RRM2', 'JAG1', 'NOTCH2', 
  'MUC5AC', 'SLC22A3', 'SLC29A3', 'MAP2', 'DCK', 'EGFR', 'MLH1', 'RRM1', 
  'SLC28A1', 'CHFR', 'MUC4', 'NRP1', 'SMARCA2', 'ABCC11', 'ARID1A', 'BNIP3', 
  'DKK3', 'IGFBP3', 'ISG15', 'LONRF2', 'MCL1', 'MUTYH', 'S100A4', 'SLC22A7', 
  'SLC28A2', 'SST', 'TYMS', 'NT5C1A', 'SLC22A2', 'SOX8', 'GSTM2', 'PROKR2', 
  'CTNNB1', 'ZEB1', 'VASH2', 'ABCB5', 'SLC28A3', 'MAP3K7', 'CXCL8', 'SLFN11', 
  'HNF1A', 'PYCARD',   
  'ADAMTS1', 'ACIN1', 'ALX4','BMP3', 'BNC1', 'CDKN2A', 'GSTP1', 'HIC1', 'MSX2', 
  'MEST', 'MGMT', 'MLH1', 'NEUROG1', 'NPTX2', 'PCDH10', 'PENK', 'POU4F1', 'SEMA5A', 
  'SEPTIN9', 'SFRP1', 'SFRP2', 'SPSB4', 'TFPI2', 'TNFRSF10C', 'WNT5A', 'UCHL1'))
## --------------------------- write down file for APP ------------------------------------------------------
pathway_genes <- pathways_long_[pathways_long_$gene_name %in% genes_PAAD, ]

write_csv(pathway_genes, "Pathway_genes_panel.csv")
write_csv(pathways_long_, "Pathway_genes_deg.csv")
