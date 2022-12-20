devtools::install_github("AskExplain/gcode_R")
library(gcode)

install.packages("mclust")
library(mclust)
install.packages("mlbench")
library(mlbench)
install.packages("pdfCluster")
library(pdfCluster)


store_scores <- c()

config <- gcode::extract_config(F)
config$init <- list(alpha="runif",beta="runif")
config$i_dim <- 100
config$j_dim <- 100
config$max_iter <- 100
config$tol <- 1e-3
config$dimension_reduction <- F
config$verbose <- F
config$learn_rate <- 1e-2
config$batch_size <- 100
config$n.cores <- 6


name_data = (c("wine","oliveoil","Glass","thyroid","wdbc","PimaIndiansDiabetes","BostonHousing","Ionosphere"))

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

  
  
  
  
  for (id in c(1:dim(main_data)[2])){
    
    
    for (seed in 1:5){
      set.seed(seed)
      config$seed <- seed
      
      x <- as.matrix(main_data[,-c(id)])
      y <- (main_data[,id,drop=F])
      
      train_ids <- sample(c(1:dim(x)[1]),dim(x)[1]*0.8)
      test_ids <- c(1:dim(x)[1])[-train_ids]
      
      join <- gcode::extract_join_framework(F)
      join$complete <- lapply(c(join$complete),function(X){c(1)})
      join$covariance <- c(0)
      
      reference <- gcode::extract_references_framework(F)
      reference$data_list <- 1
      
      gcode.main <- gcode::gcode(data_list = list(main_data[train_ids,,drop=F]), config = config, join = join, references = reference)
      
      sample_encoder <- gcode.main$main.parameters$alpha_sample[[1]]
      extended_sample_decoder <- c()
      
      knn_cluster <- kmeans(t(sample_encoder),floor(length(train_ids)*0.99))
      knn_ids_to_cov <- FNN::knnx.index(data = t(sample_encoder), query = knn_cluster$centers, k = 5)
      
      knn_measure_properties <- lapply(c(1:dim(knn_ids_to_cov)[1]),function(X){
        internal_measure_Gaussian <- t(sample_encoder[,knn_ids_to_cov[X,],drop=F])
        return(list(mu=colMeans(internal_measure_Gaussian),cov = cov(internal_measure_Gaussian)))
      })
      
      N_sample <- 100
      
      extended_sample_decoder <- rbind(do.call('rbind',lapply(knn_measure_properties,function(X){
        mvtnorm::rmvnorm(N_sample, mean = X$mu, sigma = X$cov)
      })),t(sample_encoder))
      
      extended_x <- (extended_sample_decoder)%*%MASS::ginv((sample_encoder)%*%t(sample_encoder))%*%sample_encoder%*%as.matrix(x[train_ids,,drop=F] - gcode.main$main.parameters$intercept[[1]][-id]) + gcode.main$main.parameters$intercept[[1]][-id]
      extended_y <- (extended_sample_decoder)%*%MASS::ginv((sample_encoder)%*%t(sample_encoder))%*%sample_encoder%*%as.matrix(y[train_ids,,drop=F] - gcode.main$main.parameters$intercept[[1]][id]) + gcode.main$main.parameters$intercept[[1]][id]
      
      original_model.glmnet <- glmnet::cv.glmnet(y=y[train_ids,],x=x[train_ids,])
      original_pred.glmnet <- predict(object = original_model.glmnet, newx =x[test_ids,])
      extended_model.glmnet <- glmnet::cv.glmnet(y=extended_y,x=extended_x)
      extended_pred.glmnet <- predict(object = extended_model.glmnet, newx = x[test_ids,])
      
      original_model.xgboost <- xgboost::xgboost(data=(x[train_ids,]), label = (y[train_ids,]), nrounds = 100, verbose = F)
      original_pred.xgboost <- predict(object = original_model.xgboost,newdata= (x[test_ids,]))
      extended_model.xgboost <- xgboost::xgboost(data=(extended_x), label=(extended_y), nrounds = 100, verbose = F)
      extended_pred.xgboost <- predict(object = extended_model.xgboost,newdata = (x[test_ids,]))
      
      original_model.ranger <- ranger::ranger(x=(x[train_ids,]), y = (y[train_ids,]), verbose = F)
      original_pred.ranger <- predict(object = original_model.ranger,data = (x[test_ids,]))$predictions
      extended_model.ranger <- ranger::ranger(x=(extended_x), y=(extended_y), verbose = F)
      extended_pred.ranger <- predict(object = extended_model.ranger,data = (x[test_ids,]))$predictions
      
      
      
      
      store_scores <- rbind(store_scores,c(name_data_ID,
                                           id,
                                           cor(y[test_ids,] , original_pred.glmnet),
                                           cor(y[test_ids,] , extended_pred.glmnet),
                                           which(c(cor(y[test_ids,] , original_pred.glmnet),
                                                   cor(y[test_ids,] , extended_pred.glmnet)
                                           )==max(c(cor(y[test_ids,] , original_pred.glmnet),
                                                    cor(y[test_ids,] , extended_pred.glmnet)))),
                            cor(y[test_ids,] , original_pred.xgboost),
                            cor(y[test_ids,] , extended_pred.xgboost),
                            which(c(cor(y[test_ids,] , original_pred.xgboost),
                                    cor(y[test_ids,] , extended_pred.xgboost)
                            )==max(c(cor(y[test_ids,] , original_pred.xgboost),
                                     cor(y[test_ids,] , extended_pred.xgboost)))),
                            cor(y[test_ids,] , original_pred.ranger),
                            cor(y[test_ids,] , extended_pred.ranger),
                            which(c(cor(y[test_ids,] , original_pred.ranger),
                                    cor(y[test_ids,] , extended_pred.ranger)
                            )==max(c(cor(y[test_ids,] , original_pred.ranger),
                                     cor(y[test_ids,] , extended_pred.ranger))))
      ))
      print(store_scores)
      if (dim(store_scores)[1]>3){
      print("glmnet")
      print(t.test(as.numeric(store_scores[,3]),as.numeric(store_scores[,4])))
      print("xgboost")
      print(t.test(as.numeric(store_scores[,6]),as.numeric(store_scores[,7])))
      print("ranger")
      print(t.test(as.numeric(store_scores[,9]),as.numeric(store_scores[,10])))
      }
    }
    
  }
  
}




main_pvals <- lapply(c(1:length(unique(store_scores[,1]))),function(X){
  main_test <- t.test(as.numeric(store_scores[store_scores[,1]==unique(store_scores[,1])[X],3]),as.numeric(store_scores[store_scores[,1]==unique(store_scores[,1])[X],4]),paired = T)
})

save(list=c("store_scores"),file = "~/Documents/main_files/AskExplain/Q4_2022/generative_sampling/data/scores.RData")



write.table(do.call('c',lapply(main_pvals,function(X){X$p.value})),file = "~/Documents/main_files/AskExplain/Q3_2022/neurips_mental_data/data/assorted_data/pvalue_scores.txt",quote = F,col.names = F, row.names = F,sep="\t")





load("~/Documents/main_files/AskExplain/Q4_2022/generative_sampling/data/scores.RData")

png("~/Documents/main_files/AskExplain/Q4_2022/generative_sampling/figures/boxplot_glmnet.png", width = 400, height = 500)
boxplot(apply(store_scores[which(as.integer(store_scores[,5])>0),c(3,4)],2,as.numeric), ylab = "Pearson correlation", names = c("Original","***** \n ***** \n Synthetic"), main = "Linear model with regularisation")
dev.off()
png("~/Documents/main_files/AskExplain/Q4_2022/generative_sampling/figures/boxplot_xgboost.png", width = 400, height = 500)
boxplot(apply(store_scores[which(as.integer(store_scores[,5])>0),c(6,7)],2,as.numeric), ylab = "Pearson correlation", names = c("Original","***** \n ***** \n Synthetic"), main = "Boosting")
dev.off()
png("~/Documents/main_files/AskExplain/Q4_2022/generative_sampling/figures/boxplot_ranger.png", width = 400, height = 500)
boxplot(apply(store_scores[which(as.integer(store_scores[,5])>0),c(9,10)],2,as.numeric), ylab = "Pearson correlation", names = c("Original","***** \n ***** \n Synthetic"), main = "Random Forest")
dev.off()





