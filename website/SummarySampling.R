library(shiny)
library(ggplot2)
library(bslib)

ui <- fluidPage(
  div(img(src="small_logo.png")),
  
  h1("AskExplain"),
  
  h2("Sampling new points from a known dataset"),
  
  fluidRow(
    column(2,wellPanel(
      fileInput(inputId = "upload_data", label = "CSV", accept = ".csv")
      )
    ),
    column(5,wellPanel(
      h3("Known dataset"),
      DT::dataTableOutput("data_original")
      )
    ),
    column(5,wellPanel(
      h3("New dataset"),
      DT::dataTableOutput("data_new")
      )
    )
  )
    
)


server <- function(input, output) {
  
  values <- reactiveValues(df_data = NULL)
  
  observeEvent(input$upload_data, {
    x <- values$df_data <- read.csv(input$upload_data$datapath)
    
    config <- gcode::extract_config(F)
    config$init <- list(alpha="rnorm",beta="rnorm")
    config$i_dim <- 3
    config$j_dim <- 3
    config$max_iter <- 100
    config$tol <- 1e-5
    config$dimension_reduction <- F
    config$verbose <- F
    
    join <- gcode::extract_join_framework(F)
    join$complete <- lapply(c(join$complete),function(X){c(1)})
    join$covariance <- c(0)
    
    reference <- gcode::extract_references_framework(F)
    reference$data_list <- 1
    
    extended_sample_decoder <- c()
    
    gcode.main <- gcode::gcode(data_list = list(x), config = config, join = join, references = reference)
    
    sample_encoder <- rbind(gcode.main$main.parameters$alpha_sample[[1]])
    extended_sample_decoder <- c()
    
    knn_cluster <- kmeans(t(sample_encoder),dim(x)[1]*0.95)
    
    knn_measure_properties <- c()
    for (i in c(3,9,15)){
      knn_ids_to_cov <- FNN::knnx.index(data = t(sample_encoder), query = knn_cluster$centers, k = i)
      
      knn_measure_properties <- c(knn_measure_properties,lapply(c(1:dim(knn_ids_to_cov)[1]),function(X){
        internal_measure_Gaussian <- t(sample_encoder[,knn_ids_to_cov[X,],drop=F])
        return(list(mu=colMeans(internal_measure_Gaussian),cov = cov(internal_measure_Gaussian)))
      }))
      
    }
    
    N_sample <- 5
    
    extended_sample_decoder <- rbind(do.call('rbind',lapply(knn_measure_properties,function(X){
      mvtnorm::rmvnorm(n = N_sample, mean = X$mu, sigma = lqmm::make.positive.definite(X$cov))
    })),extended_sample_decoder)
    
    
    sample_encoder <- rbind(1,sample_encoder)
    extended_sample_decoder <- cbind(1,extended_sample_decoder)
    
    values$extended_x <- round((extended_sample_decoder)%*%MASS::ginv((sample_encoder)%*%t(sample_encoder))%*%sample_encoder%*%as.matrix(x),3)
    
  })
  
  output$data_original <- DT::renderDataTable({
    
    DT::datatable(values$df_data[,-1], rownames = F)
    
  })
  
  output$data_new <- DT::renderDataTable({
    
    DT::datatable(values$extended_x[,-1], rownames = F)
    
  })
  
}



shinyApp(ui, server)


