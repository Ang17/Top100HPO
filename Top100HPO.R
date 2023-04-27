library(shiny)
library(ontologySimilarity)
library(ontologyIndex)
library(ontologyPlot)
library(data.table)
library(tidyverse)

ui <- fluidPage(titlePanel("HPO Based Gene Priority"),
                
                sidebarLayout(
                  sidebarPanel(
                    helpText(
                      "Note: Files should be in .csv format
               and should not contain headers. The
               first column should include proband IDs
               and the second column should include HPO
               terms with a ; between subsequent terms,
               no space. EG: HP:0000001;HP:0000002"
                    ),
                    
                    fileInput("file", h3("File input")),
                    
                    actionButton(inputId = "submit_loc",
                                 label = "Submit"),
                    
                    downloadButton("downloadData", "Download Gene Matches"),
                    
                    downloadButton("downloadData1", "Download Disease Matches")
                    
                  ),
                  
                  mainPanel(
                    DT::dataTableOutput("matches"),
                    DT::dataTableOutput("matches1"))
                ))

server <- function(input, output) {
  values <- reactiveValues()
  observeEvent(eventExpr = input[["submit_loc"]],
               handlerExpr = {
                 print("PRESSED")
                 
                 #getting the HPO ontology dynamically; always the most updated ontology (every 2weeks)
                 #future look at other ontologies
                 ontology <-
                   get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags = "everything")
                 
                 
                 information_content <-
                   descendants_IC(ontology)
                 
                 #reference using method change to lin or resnick
                 term_sim_mat <- get_term_sim_mat(
                   ontology,
                   information_content = information_content,
                   method = "lin",
                   row_terms = names(information_content),
                   col_terms = names(information_content)
                 )
                 
                 gene2pheno <-
                   read.delim(
                     "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt",
                     header = T,
                     stringsAsFactors = F
                   )
                 
                 #loading the disease database
                 #header columns are good here
                 db2pheno <-
                   read.delim(
                     "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa",
                     skip = 4,
                     header = T,
                     stringsAsFactors = F
                   )
                 
                 #p terms are phenotype terms
                 # c terms are clinical terms
                 # i terms are inheritance terms
                 # letter designations used by database
                 # only want phenotype terms because other terms have hpo terms that will be included in analysis if you don't take them
                 # this is a research analysis; agnostic of what the inheritance mode is or just not assuming anything known about clinical onset kind of terms
                 nonp_terms <- db2pheno[db2pheno$aspect != "P", ]
                 
                 # getting HPO Ids of non phenotype terms to remove
                 nonp_terms1 <- unique(nonp_terms$hpo_id)
                 
                 #removing from gene database
                 gene2pheno <-
                   gene2pheno[!gene2pheno$hpo_id %in% nonp_terms1, ]
                 
                 #cleaning up disease database
                 #databases don't have same structure; disease db has more info and is easier to clean
                 #Qualifier column is not associations and makes clinical input a lot harder
                 #remove NOTs and only keep P terms for phenotypes
                 db2pheno <- db2pheno %>% filter(qualifier != "NOT",
                                                 aspect == "P")
                 
                 ## O&T Combined
                 
                 gn_ids <-
                   unique(gene2pheno$gene_symbol)
                 
                 #unique ensure no hpo term duplicates
                 gene_sets <-
                   lapply(gn_ids, function(x) {
                     unique(gene2pheno[x == gene2pheno$gene_symbol,][["hpo_id"]])
                   })
                 
                 names(gene_sets) <- 
                   gn_ids
                 
                 og_gene_sets <-
                   gene_sets
                 
                 #
                 
                 db2pheno <-
                   db2pheno %>% filter(!db2pheno$database_id %like% "DECIPHER")
                 
                 dz_id <-
                   unique(db2pheno$database_id)
                 
                 dz_sets <- lapply(dz_id, function(x) {
                   unique(db2pheno[x == db2pheno$database_id, ][["hpo_id"]])
                 })
                 
                 names(dz_sets) <-
                   dz_id
                 
                 og_dz_sets <-
                   dz_sets
                 
                 ## O&T Patient
                 #here put a reference
                 # preprocess file to match our input or change code to fit input file
                 #want a list of hpo terms
                 #importing ptdata
                 file1 <- input$file
                 ptdata <-
                   read.csv(file1$datapath,
                            stringsAsFactors = F,
                            header = F)
                 
                 #split hpo column based on semicolon
                 pt_term_sets <-
                   lapply(strsplit(ptdata$V2, ";"), unique)
                 
                 
                 #use names from column 1 in original pheno input to name pt_term_sets
                 names(pt_term_sets) <-
                   ptdata$V1
                 
                 #save original
                 og_pt_term_sets <- pt_term_sets
                 
                 # O&T Patient Combined
                 
                 for (a in 1:length(pt_term_sets)) {
                   gene_sets[length(gene_sets) + 1] <- pt_term_sets[a]
                 }
                 
                 names(gene_sets)[(length(og_gene_sets) + 1):(length(og_gene_sets) + length(pt_term_sets))] <-
                   names(pt_term_sets)
                 
                 for (a in 1:length(pt_term_sets)) {
                   dz_sets[length(dz_sets) + 1] <- pt_term_sets[a]
                 }
                 
                 names(dz_sets)[(length(og_dz_sets) + 1):(length(og_dz_sets) + length(pt_term_sets))] <-
                   names(pt_term_sets)
                 
                 ## Gene Analysis Combined
                 
                 # Make changes to pt_term_sets
                 pt_term_set3 <- pt_term_sets
                 
                 sim_mat_gn <- get_sim_grid(
                   term_sim_mat = term_sim_mat,
                   term_sets = gene_sets,
                   term_sets2 = pt_term_set3,
                   term_sim_method = "lin",
                   combine = "average"
                 )
                 
                 dist_mat_gn <-
                   max(sim_mat_gn) - sim_mat_gn
                 
                 sim_df_gn <-
                   as.data.frame(
                     sim_mat_gn,
                     row.names = rownames(sim_mat_gn),
                     col.names = colnames(sim_mat_gn)
                   )
                 
                 #
                 
                 sim_score_pt_gn <-
                   sim_df_gn[colnames(sim_df_gn) %in% names(pt_term_set3), ]
                 
                 sim_df_gn <-
                   sim_df_gn[!rownames(sim_df_gn) %in% names(pt_term_set3), ]
                 
                 sim_score_gn <-
                   sim_df_gn[colnames(sim_df_gn) %in% names(pt_term_set3)]
                 
                 #
                 
                 ord_sim_gn <- list()
                 
                 ord_sim_gn <- lapply(sim_score_gn, function(x) {
                   ord_sim_gn <-
                     x[order(x, decreasing = T)[1:100]]
                   names(ord_sim_gn) <-
                     rownames(sim_score_gn)[order(x, decreasing = T)[1:100]]
                   ord_sim_gn
                 })
                 
                 ptest_grps <- lapply(ord_sim_gn, names)
                 
                 names(ptest_grps) <-
                   NULL
                 
                 query <- names(pt_term_set3)
                 
                 #
                 
                 test1 <- ptest_grps
                 
                 names(test1) <- query
                 
                 values$test2 <- as.data.frame(test1)
                 
                 ## Disease Analysis Combined
                 
                 # Make changes to pt_term_sets
                 
                 sim_mat_dz <- get_sim_grid(
                   term_sim_mat = term_sim_mat,
                   term_sets = dz_sets,
                   term_sets2 = pt_term_set3,
                   term_sim_method = "lin",
                   combine = "average"
                 )
                 
                 dist_mat_dz <-
                   max(sim_mat_dz) - sim_mat_dz 
                 
                 sim_df_dz <-
                   as.data.frame(sim_mat_dz,
                                 row.names = rownames(sim_mat_dz),
                                 col.names = colnames(sim_mat_dz))
                 
                 #
                 
                 sim_score_pt_dz <-
                   sim_df_dz[colnames(sim_df_dz) %in% names(pt_term_set3), ]
                 
                 sim_df_dz <-
                   sim_df_dz[!rownames(sim_df_dz) %in% names(pt_term_set3), ]
                 
                 sim_score_dz <-
                   sim_df_dz[colnames(sim_df_dz) %in% names(pt_term_set3)]
                 
                 #
                 
                 ord_sim_dz <- lapply(sim_score_dz, function(x) {
                   ord_sim_dz <-
                     x[order(x, decreasing = T)[1:100]]
                   names(ord_sim_dz) <-
                     rownames(sim_score_dz)[order(x, decreasing = T)[1:100]]
                   ord_sim_dz
                 })
                 
                 ptest_grps_dz <- lapply(ord_sim_dz, names)
                 
                 names(ptest_grps_dz) <-
                   NULL 
                 
                 #
                 
                 test3 <- ptest_grps_dz
                 
                 names(test3) <- query
                 
                 values$test4 <- as.data.frame(test3)
                 
                 #
                 
                 output$matches <- DT::renderDataTable({
                   values$test2
                 })
                 
                 output$matches1 <- DT::renderDataTable({
                   values$test4
                 })
                 
                 output$downloadData <- downloadHandler(
                   filename = "Gene_Matches.csv",
                   content = function(file) {
                     write.csv(values$test2, file, row.names = F)
                   }
                 )
                 
                 output$downloadData1 <- downloadHandler(
                   filename = "Disease_Matches.csv",
                   content = function(file1) {
                     write.csv(values$test4, file1, row.names = F)
                   }
                 )
               })
  
}

shinyApp(ui, server)
