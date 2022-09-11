
devtools::install_github("jlmelville/mnist")
library(mnist)

# fetch the data set from the MNIST website
mnist <- download_mnist()

# view the fifth digit
show_digit(mnist, 1)

# first 60,000 instances are the training set
mnist_train <- head(mnist, 60000)
# the remaining 10,000 are the test set
mnist_test <- tail(mnist, 10000)



test_mnist <- generate_new_samples(data_list = list(mnist[,-785],
                                                   as.matrix(mltools::one_hot(data.table::data.table(data.frame(mnist[,785])))))
                                  ,N_new_samples = 1000, dimension_reduction = 10, seed = 1)

save(list=c("new_mnist"),file="~/Documents/main_files/AskExplain/Q3_2022/neurips_mental_data/data/mnist/workflow/mnist_save_resampled.RData")




png("~/Documents/main_files/AskExplain/Q3_2022/neurips_mental_data/figures/mnist/all_mnist.png",width = 3000, height = 1400)

ids_main <- do.call('rbind',lapply(c(1:dim(new_mnist[[2]])[1]),function(Y){
  X <- new_mnist[[2]][Y,]
  sorted_X <- sort(X,decreasing = T)
  c(Y,which(X==max(X)),sorted_X[1] - sorted_X[2])
  }))


ids_main <-lapply(c(1:10),function(X){
  ids_main[ids_main[,2]==X,]
})

par(mfcol=c(5,10))
for (i in c(1:10)){
  set.seed(i)
  
  for (j in c(1:5)){
    a <- image(matrix(as.numeric(new_mnist[[1]][ids_main[[i]][order(ids_main[[i]][,3],decreasing = T)[j],1],]), nrow = 28)[, 28:1], 
               col = grDevices::gray(12:1 / 12),
               xaxt = "n", yaxt = "n", xlab="", ylab="")
  }

}
dev.off()




