library(mclust)

generate_new_samples <- function(data_list,N_new_samples,dimension_reduction=70,seed=1){
  set.seed(seed)
  
  config <- gcode::extract_config(F)
  config$init <- list(alpha="runif",beta="rnorm")
  config$i_dim <- dimension_reduction
  config$j_dim <- dimension_reduction
  config$max_iter <- 5
  config$tol <- 1
  config$seed <- seed
  config$dimension_reduction <- F
  config$verbose <- T
  
  join <- gcode::extract_join_framework(F)
  join$complete <- lapply(c(join$complete),function(X){1:length(data_list)})
  join$complete$alpha_signal <- c(1,1)
  join$covariance <- rep(0,length(data_list))
  
  gcode.main <- gcode::gcode(data_list = data_list, config = config, join = join)
  
  sample_encoder <- gcode.main$main.parameters$alpha_signal[[1]]
  
  mclust_sample_encoder <- mclust::Mclust(t(sample_encoder))
  
  sample_mclust <- function(mclust_sample_encoder,N_sample){
    
    return_vec <- c()
    for (i in 1:length(mclust_sample_encoder$parameters$pro)){
      set.seed(i)
      id_prob <- i
      
      return_vec <- rbind(return_vec,mvtnorm::rmvnorm(N_sample,mean = mclust_sample_encoder$parameters$mean[,id_prob],sigma = mclust_sample_encoder$parameters$variance$sigma[,,id_prob]))
      
    }
    
    return_vec
  }
  
  extended_sample_decoder <- sample_mclust(mclust_sample_encoder = mclust_sample_encoder, N_sample = N_new_samples)
  
  extended <- lapply(c(1:length(data_list)),function(X){
    
    extended_x <- extended_sample_decoder%*%MASS::ginv(t(extended_sample_decoder)%*%(extended_sample_decoder))%*%sample_encoder%*%as.matrix(data_list[[X]])
    
    return(extended_x)
    
  })
  
  return(extended)
  
}






predict_comparison <- function(x,N_new_samples,seed){
  
  config <- gcode::extract_config(F)
  config$init <- list(alpha="rsvd",beta="rsvd")
  config$i_dim <- 10
  config$j_dim <- 10
  config$max_iter <- 500
  config$tol <- 1e-2
  config$dimension_reduction <- F
  config$verbose <- F
  
  store_scores <- c()
  
  top_var <- order(apply(x,2,var),decreasing = T)
  
  y <- x[,top_var[1],drop=F]
  x <- x[,top_var[1:101],drop=F]
  
  for (seed in 1:10){
    set.seed(seed)
    config$seed <- seed
    
    N = dim(x)[1]*0.6
    
    train_ids <- sample(c(1:dim(x)[1]),N)
    test_ids <- c(1:dim(x)[1])[-train_ids]
    
    
    join <- gcode::extract_join_framework(F)
    join$complete <- lapply(c(join$complete),function(X){c(1)})
    join$complete$alpha_sample<-c(1)
    join$complete$beta_sample<-c(1)
    join$complete$code<-c(1)
    join$complete$alpha_signal<-c(1)
    join$complete$beta_signal<-c(1)
    join$covariance <- c(0)
    
    gcode.main <- gcode::gcode(data_list = list(cbind(y,x)[train_ids,]), config = config, join = join)
    
    sample_encoder <- gcode.main$main.parameters$alpha_sample[[1]]
    
    extended_sample_decoder <- mvtnorm::rmvnorm(N_new_samples,mean = rowMeans(sample_encoder),sigma = cov(t(sample_encoder)))
    
    extended_x <- (extended_sample_decoder)%*%MASS::ginv(sample_encoder%*%t(sample_encoder))%*%sample_encoder%*%as.matrix(x[train_ids,])
    extended_y <- (extended_sample_decoder)%*%MASS::ginv(sample_encoder%*%t(sample_encoder))%*%sample_encoder%*%as.matrix(y[train_ids,])
    
    original_model <- glmnet::cv.glmnet(x=x[train_ids,], y=y[train_ids,],family="gaussian")
    original_pred <- predict(object = original_model,newx = x[test_ids,])
    
    extended_model <- glmnet::cv.glmnet(x=rbind(x[train_ids,],extended_x), y=rbind(y[train_ids,,drop=F],extended_y),family="gaussian")
    extended_pred <- predict(object = extended_model,newx = x[test_ids,])
    
    
    store_scores <- rbind(store_scores,c(mean(abs(y[test_ids,] - original_pred)),
                                         mean(abs(y[test_ids,] - extended_pred)),
                                         which(c(mean(abs(y[test_ids,] - original_pred)),
                                                 mean(abs(y[test_ids,] - extended_pred))
                                         )==min(c(mean(abs(y[test_ids,] - original_pred)),
                                                  mean(abs(y[test_ids,] - extended_pred))))))
    )
    
    print(store_scores)
    
  }
  
  store_scores
  
  
}



clean_gex <- function(main_gex){
  main_gex <- apply(main_gex,2,as.numeric)
  main_gex[is.null(main_gex)] <- 0
  main_gex[is.na(main_gex)] <- 0
  main_gex <- log10(1+main_gex/rowSums(main_gex+1))*1e6
  t(main_gex)
}

