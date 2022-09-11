
setwd("~/Documents/main_files/AskExplain/Q3_2022/neurips_mental_data/shiny/")
load("../data/GTEX/GTEX_v8_brain.RDS")
load("../data/GTEX/GTEX_v8_brain_generated.RDS")

gtex_v8 <- read.delim(file="../data/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", skip=2)
gtex_v8_annotation <- read.delim(file="../data/GTEX/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
gtex_v8_TPM_genes <- read.delim(file="../data/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", skip=2)


gtex_v8_brain <- gtex_v8[,colnames(gtex_v8) %in% gsub(pattern = "-",replacement = ".",gtex_v8_annotation$SAMPID)[gtex_v8_annotation$SMTS=="Brain"]]
gtex_v8_brain_annotation <- gtex_v8_annotation[gsub(pattern = "-",replacement = ".",gtex_v8_annotation$SAMPID) %in% colnames(gtex_v8_brain),]
top_var_brain_genes <- order(apply(gtex_v8_brain,1,var),decreasing = T)

save(list = c("gtex_v8_brain","gtex_v8_brain_annotation","gtex_v8_TPM_genes","top_var_brain_genes"),file = "../data/GTEX/GTEX_v8_brain.RDS")




N_genes <- 1000
a <- Sys.time()
gtex_v8_brain_generated <- generate_new_samples(data_list = list(
  clean_gex(gtex_v8_brain[top_var_brain_genes[1:N_genes],]),
  do.call('cbind',lapply(unique(gtex_v8_brain_annotation$SMTSD),function(X){0+(gtex_v8_brain_annotation$SMTSD==X)}))
), N_new_samples = 1000, seed = 1)
b <- Sys.time()


gtex_v8_brain_generated_annotated <- t(apply(gtex_v8_brain_generated[[2]],1,function(X){1*(X==max(X))}))
save(list = c("gtex_v8_brain_generated","gtex_v8_brain_generated_annotated"),file = "../data/GTEX/GTEX_v8_brain_generated.RDS")









differential_test_two_tissues <- function(tissue_1,tissue_2,genes){
  
  all_de_genes <- c()
  for (i in 1:length(genes)){
    if (i%%10==1){
      print(i)
    }
    all_de_genes <- rbind(all_de_genes,
                          c(genes[i],t.test(
                            tissue_1[,i],
                            tissue_2[,i]
                            )$p.value
                          )
    )
  }
  return(all_de_genes)
}



full_cor <- c()
for (i in 1:13){
  row_cor <- c()
  for (j in 1:13){
    differential_test_generated <- differential_test_two_tissues(
      gtex_v8_brain_generated[[1]][gtex_v8_brain_generated_annotated[,i]==1,],
      gtex_v8_brain_generated[[1]][gtex_v8_brain_generated_annotated[,j]==1,],
      gtex_v8_TPM_genes[top_var_brain_genes[1:N_genes],2])
    
    differential_test_observed <- differential_test_two_tissues(
      t(gtex_v8_brain)[gtex_v8_brain_annotation$SMTSD==unique(gtex_v8_brain_annotation$SMTSD)[i],top_var_brain_genes[1:N_genes]],
      t(gtex_v8_brain)[gtex_v8_brain_annotation$SMTSD==unique(gtex_v8_brain_annotation$SMTSD)[j],top_var_brain_genes[1:N_genes]],
      gtex_v8_TPM_genes[top_var_brain_genes[1:N_genes],2])
    
    differential_test_generated <- data.frame(differential_test_generated[,1],as.numeric(differential_test_generated[,2]))
    differential_test_observed <- data.frame(differential_test_observed[,1],as.numeric(differential_test_observed[,2]))
    
    row_cor <- c(row_cor,cor(-log10(1e-300+as.numeric(differential_test_generated[,2])),-log10(1e-300+as.numeric(differential_test_observed[,2]))))
    plot(-log10(as.numeric(differential_test_generated[,2])),-log10(as.numeric(differential_test_observed[,2])))
    print(row_cor)
  }
  full_cor <- cbind(full_cor,row_cor)
}




colnames(full_cor) <- row.names(full_cor) <- unique(gtex_v8_brain_annotation$SMTSD)




save(list = c("full_cor"),file = "../data/GTEX/full_cor.RData")

write.table(full_cor,"../data/GTEX/full_cor.txt",quote = F,sep="\t")

write.table(apply(full_cor,2,function(X){X>0.85})*full_cor,"../data/GTEX/partial_cor.txt",quote = F,sep="\t")






