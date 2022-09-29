store_scores <- c()


# devtools::install_github("AskExplain/gcode_R",force=T)
library(mclust)

config <- gcode::extract_config(F)
config$init <- list(alpha="runif",beta="rnorm")
config$i_dim <- 70
config$j_dim <- 70
config$max_iter <- 100
config$tol <- 1e-5
config$dimension_reduction <- F
config$verbose <- F



library(mlbench)
library(pdfCluster)
name_data = (c("wine","oliveoil","Glass","thyroid","wdbc","PimaIndiansDiabetes","BostonHousing","Ionosphere","Shuttle","Satellite"))

for (name_data_ID in name_data){
  
  if (name_data_ID=="wine"){
    data(wine)
    id_remove <- c(1)
    main_data <- scale(wine[,-id_remove])
  }
  if (name_data_ID=="oliveoil"){
    data(oliveoil)
    id_remove <- c(1,2)
    main_data <- scale(oliveoil[,-id_remove])
    
  }
  if (name_data_ID=="Glass"){
    data(Glass)
    id_remove <- c(10)
    main_data <- scale(Glass[,-id_remove])
    
  }
  if (name_data_ID=="thyroid"){
    data(thyroid)
    id_remove <- c(1)
    main_data <- scale(thyroid[,-id_remove])
    
  }
  if (name_data_ID=="wdbc"){
    data(wdbc)
    id_remove <- c(1,2)
    main_data <- scale(wdbc[,-id_remove])
    
  }
  if (name_data_ID=="PimaIndiansDiabetes"){
    data(PimaIndiansDiabetes)
    id_remove <- c(1,9)
    main_data <- scale(PimaIndiansDiabetes[,-id_remove])
    
  }
  if (name_data_ID=="BostonHousing"){
    data(BostonHousing)
    id_remove <- c(4,9,10)
    main_data <- scale(BostonHousing[,-id_remove])
    
  }
  if (name_data_ID=="Ionosphere"){
    data(Ionosphere)
    id_remove <- c(1,2,35)
    main_data <- scale(Ionosphere[,-id_remove])
    
  }
  if (name_data_ID=="Shuttle"){
    data(Shuttle)
    id_remove <- c(10)
    main_data <- scale(Shuttle[,-id_remove])
    
  }
  if (name_data_ID=="Satellite"){
    data(Satellite)
    id_remove <- c(37)
    main_data <- scale(Satellite[,-id_remove])
    
  }
  
  
  
  
  
  for (id in c(1:dim(main_data)[2])){
    
    
    for (seed in 1:30){
      set.seed(seed)
      config$seed <- seed
      
      if (dim(main_data)[1]>500){
        main_data <- main_data[sample(c(1:dim(main_data)[1]),500),]
      }
      
      x <- as.matrix(main_data[,-c(id)])
      y <- (main_data[,id,drop=F])
      
      train_ids <- sample(c(1:dim(x)[1]),dim(x)[1]*0.5)
      test_ids <- c(1:dim(x)[1])[-train_ids]
      
      join <- gcode::extract_join_framework(F)
      join$complete <- lapply(c(join$complete),function(X){c(1)})
      join$covariance <- c(0)
      
      # Run generative encoding
      gcode.main <- gcode::gcode(data_list = list(cbind(main_data[train_ids,,drop=F])), config = config, join = join)
      
      # Extract alpha, the internal sample parameter in the gcode model
      sample_encoder <- gcode.main$main.parameters$alpha_signal[[1]]
      
      # Run a mixture model over alpha
      mclust_sample_encoder <- mclust::Mclust(t(sample_encoder))
      
      # Generate new samples of alpha, now called A
      sample_mclust <- function(mclust_sample_encoder,N_sample){
        
        return_vec <- c()
        for (i in 1:length(mclust_sample_encoder$parameters$pro)){
          set.seed(i)
          id_prob <- i
          
          return_vec <- rbind(return_vec,mvtnorm::rmvnorm(N_sample,mean = mclust_sample_encoder$parameters$mean[,id_prob],sigma = mclust_sample_encoder$parameters$variance$sigma[,,id_prob]))
          
        }
        
        return_vec
      }
      
      # Generate A
      extended_sample_decoder <- sample_mclust(mclust_sample_encoder = mclust_sample_encoder, N_sample = 1000)
      
      # From A and alpha, generate new samples
      extended_x <- (extended_sample_decoder)%*%MASS::ginv(t(extended_sample_decoder)%*%(extended_sample_decoder))%*%sample_encoder%*%as.matrix(x[train_ids,,drop=F])
      extended_y <- (extended_sample_decoder)%*%MASS::ginv(t(extended_sample_decoder)%*%(extended_sample_decoder))%*%sample_encoder%*%as.matrix(y[train_ids,,drop=F])
      
      # Run regression on observed data
      original_model <- glmnet::cv.glmnet(x=(x[train_ids,]), y=y[train_ids,])
      original_pred <- predict(object = original_model,newx = (x[test_ids,]))
      
      # Run regression on generated data
      extended_model <- glmnet::cv.glmnet(x=(extended_x), y=as.matrix(extended_y))
      extended_pred <- predict(object = extended_model,newx = (x[test_ids,]))
      
      store_scores <- rbind(store_scores,c(name_data_ID,
                                           id,mean(abs(y[test_ids,] - original_pred)),
                                           mean(abs(y[test_ids,] - extended_pred)),
                                           which(c(mean(abs(y[test_ids,] - original_pred)),
                                                   mean(abs(y[test_ids,] - extended_pred))
                                           )==min(c(mean(abs(y[test_ids,] - original_pred)),
                                                    mean(abs(y[test_ids,] - extended_pred))))))
      )
      
    }
    
  }
  
}




main_pvals <- lapply(c(1:length(unique(store_scores[,1]))),function(X){
  main_test <- t.test(as.numeric(store_scores[store_scores[,1]==unique(store_scores[,1])[X],3]),as.numeric(store_scores[store_scores[,1]==unique(store_scores[,1])[X],4]),paired = T)
})

save(list=c("store_scores","main_pvals"),file = "~/Documents/main_files/AskExplain/Q3_2022/neurips_mental_data/data/assorted_data/store_scores.RData")



write.table(do.call('c',lapply(main_pvals,function(X){X$p.value})),file = "~/Documents/main_files/AskExplain/Q3_2022/neurips_mental_data/data/assorted_data/pvalue_scores.txt",quote = F,col.names = F, row.names = F,sep="\t")





