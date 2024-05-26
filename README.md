# Survival_Analysis

Set up R studio on AWS EC2 (https://docs.elmcomputing.io/ami/x86/rstudio-steps.html)

Data - TCGA Pancreatic Cancer RNA seq data (https://github.com/rmoffitt/pdacR/blob/master/R/parse_TCGA_PAAD.R)

Basics of Survival Analysis in R (http://www.sthda.com/english/wiki/survival-analysis-basics)

R shiny web app (http://bioinformatics.sdstate.edu/idep/)

Shiny code to refer - (https://www.youtube.com/watch?v=cOfKmgxBCso)

Day 1
- Create basic Kaplan M curve using (http://www.sthda.com/english/wiki/survival-analysis-basics)
- Create a shiny app (https://www.youtube.com/watch?v=cOfKmgxBCso) display title, any kind of button by clicking that button we should be able to see Kaplan M plot
- Upload the code on this github repo.

Day 2
- Try to implement (https://github.com/jeffreythy/survival-analysis-shiny-app/tree/master) repo using our lung dat and using 2 categories (status, sex)


RNA-seq + CHip-seq integrated analysis
- CHIP-seq data for PAAD cancer cell lines (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6122955/)
- RNA-seq data for PAAD cancer in also available in the same paper
- Integrated analysis (https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-020-07205-6/MediaObjects/12864_2020_7205_MOESM1_ESM.pdf)

- To install R studio and R from scratch select right ubuntu instance. 24 and 22 is not working. Check each step while installing these packages. 
https://jagg19.github.io/2019/08/aws-r/#long_way
https://posit.co/download/rstudio-server/


Day 4,5,6,7 (AAkash)
Now you have a bot of understanding of how R studio and R shiny works. A bit about backend, frontend and how to connect both of them. 
Let's dive deep into bigger stuff. Try to replicate (http://bioinformatics.sdstate.edu/idep96/) this web app. I want you to understand and figure out things by yourself. The art of coding is unserstanding each line of code, how to reuse codes from other people to build your own, which right questions to ask on google to solve your errors, understanding them and fixing it. I want you to learn it by yourself. If stuck just google or use chatgpt free version and try to learn.

Here is the code for the above website (https://github.com/iDEP-SDSU/idep/tree/master/shinyapps/idep96). From the entire code I just want you to try and replicate the first 2 tabs (load data and pre process). Don't just try to paste the code. Understand the code and build your own using their components. Google each and every term which you don't understand. I will leave you by yourself for 2 days and will meet on Tuesday. 

1) Understand the content in the first 2 tabs. Enough things are provided on the website for you to navigate and understand.
2) Try to understand the shiny code from their github account.
3) Now build your own shiny code just for the first 2 tabs.
4) YOu can use the demo data which they have provided as your data to work with.
