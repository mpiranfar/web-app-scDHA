library(shiny)
library(ggplot2)
library(httr)
library(DT) 
library(scDHA)


# main analysis tool
data_analysis <- function(data_object) {
  
  data <- t(data_object$data)
  
  label <- as.character(data_object$label)
  
  data <- log2(data + 1)
  
  result <- scDHA(data, ncores = 2, seed = 1)
  
  cluster <- result$cluster
  
  mclust::adjustedRandIndex(cluster,label)
  
  result <- scDHA.pt(result, start.point = 1, ncores = 2, seed = 1)
  
  # Return a list containing both the processed data and the result
  return(list(data = data, result = result, label = label))
  
}


# Function to convert TSV to RDA
convert_tsv_to_rda <- function(input_file, output_file) {
  # Read TSV data using read.table
  tsv_data <- read.table(input_file, sep="\t", header=TRUE)
  
  # Save the data in RDA format using save function
  save(tsv_data, file = output_file)
  
  # Print a success message
  cat("TSV data converted and saved as RDA:", output_file, "\n")
}

# Create a scatter plot based on result
generate_scatter_plot <- function(data_plot, cell_colors) {
  
  if (!is.null(data_plot$cell_type)) {
  ggplot(data_plot, aes(x = data_plot$x, y = data_plot$y, color = data_plot$cell_type)) +
    geom_point() +
    scale_color_manual(values = cell_colors) +
    labs(title = "Scatter Plot of Trajectory Inference",
         x = "Sample Cells",
         y = "Pseudo Time",
         color = "Cell Type") +
    theme_minimal()
  } else {
    
    # If cell_type is NULL, color the points by pseudotime
    ggplot(data_plot, aes(x = data_plot$x, y = data_plot$y, color = data_plot$y)) +
      geom_point() +
      scale_color_viridis_c() + 
      labs(title = "Scatter Plot of Cell Representation",
           x = "Sample Cells",
           y = "Pseudo Time",
           color = "Pseudo Time") +
      theme_minimal()
  }
}


# UI
ui <- fluidPage(
  titlePanel("scDHA"),
  sidebarLayout(
    sidebarPanel(
      h4("Guideline Text"),
      radioButtons("choice", "If you want to upload your own data choose 'A'; otherwise if you want to select a built-in data set choose 'B':",
                   choices = c("A", "B")),
      
      conditionalPanel(
        condition = "input.choice == 'A'",
        fileInput("file1", "Upload a single-cell RNA-seq matrix in TSV format (required):"),
        fileInput("file2", "Upload a metadata file in TSV format (optional):")
      ),
      
      conditionalPanel(
        condition = "input.choice == 'B'",
        selectInput("data", "Choose a data source:",
                    choices = c("Yan Data", "Goolam Data", "Deng Data"))
      ),
      
      actionButton("submit", "Submit")
      
    ),
    
    mainPanel(
      downloadButton("downloadplot", "Download Plot as PNG"),
      downloadButton("download3", "Download Table as TSV"),
      downloadButton("download2", "Download Compressed Data as TSV"),
      uiOutput("unavailable_message"),
      plotOutput("scatter_plot"),
      tableOutput("my_table"),
      #tags$caption("Trajectory Inference")
    )
  )
)

# Server
server <- function(input, output) {
  # Reactive function to convert and display the table
  data_result <- reactive({
    if (input$submit > 0) {
      if (input$choice == "A") {
        # Get the file paths for the uploaded files
        input_tsv_file <- input$file1$datapath
        metadata_files <- input$file2$datapath
        
        # Check if file1 is uploaded
        if (!is.null(input_tsv_file) && !is.null(metadata_files)) {
          
          # combine data and metadata
          data_tsv <- cbind(input_tsv_file, metadata_files)
          
        } else if (!is.null(input_tsv_file) && is.null(metadata_files)){
          
          data_tsv <- input_tsv_file
          
        } else {
          
          # If no file is uploaded
          cat("Error: No uploaded files.")
          
          return(NULL)
        }
        
        # Define the local file path where the downloaded file is saved
        data_rda <- "output_file.rda"
        
        # Perform the conversion
        convert_tsv_to_rda(data_tsv, data_rda)
        nrow(data_tsv)
        
        # Load the RDA data
        rda_data <- load(data_rda)
        
        # Extract the data from the loaded object
        data_user <- get(rda_data)
        
        # Perform trajectory inference using the scDHA package
        list_data_user <- data_analysis(data_user)
        
        # Return the loaded data
        return(list_data_user)
        
      } else if (input$choice == "B") {
        # Process selected data source for option B
        selectedData <- input$data
        
        # Fetch online data based on the selected data source
        if (selectedData == "Yan Data") {
          
          output$unavailable_message <- renderUI({
            
            tags$p("Unfortunately, this dataset is currently unavailable.")
            
          })
          return(NULL)
        } else if (selectedData == "Goolam Data") {
          
          # Define the URL of the RDA file
          rda_file_url <- "https://raw.githubusercontent.com/duct317/scDHA/master/data/Goolam.rda"
          
          # Define the local file path where the downloaded file is saved
          local_file_path <- "Goolam.rda"
          
          # Download the RDA file from the URL to the local file path
          download.file(rda_file_url, destfile = local_file_path)
          
          # Load the data from the downloaded RDA file
          loaded_data <- load(local_file_path)
          
          # Extract the data from the loaded object
          data_to_save <- get(loaded_data)
          
          # Perform trajectory inference using the scDHA package
          Goolam_list <- data_analysis(data_to_save)
          
          # Return the loaded data
          return(Goolam_list)
          
        } else {
          cat("Error: Unable to download the file.")
        }
        
      } else if (selectedData == "Deng Data") {
        
        output$unavailable_message <- renderUI({
          
          tags$p("Unfortunately, this dataset is currently unavailable.")
          
        })
        
        return(NULL)
      } 
    } else {
      
      return(NULL)
      "Please select an option."
    }
  })
  
  # Render the table in the UI based on data
  output$my_table <- renderTable({
    
    data_reunion<- data_result()
    
    if (!is.null(data_reunion)) {
      data <- data_reunion$data
      result <- data_reunion$result
      label <- data_reunion$label
      
      # Create the table with cell ID and pseudotime columns
      table_data <- data.frame(Cell_ID = rownames(data), Pseudo_Time = result$pt)
      
      # If data is not NULL and metadata exists, add cell type column to the table
      if (nrow(table_data) > 0 && "label" %in% names(data_reunion)) {
        cell_types <- as.character(label)
        table_data$Cell_Type <- cell_types  # Insert the cell type column
      }
      
      return(table_data)
      
    } else {
      
      # Return an empty data frame if data is NULL
      return(data.frame())
    }
  })
  # Define the download handler
  output$download3 <- downloadHandler(
    
    # Specify the file name
    filename = function() {
      "table_data.tsv"
    },
    
    # Specify the file content
    content = function(file) {
      save_table <- data_result()
      
      # Create the table with cell ID and pseudotime columns
      table_data <- data.frame(Cell_ID = rownames(save_table$data), Pseudo_Time = save_table$result$pt)
      
      if (nrow(table_data) > 0 && !is.null(save_table$label)) {
        
        # If metadata exists, add cell type column to the table
        cell_types <- as.character(save_table$label)
        table_data$Cell_Type <- cell_types  # Insert the cell type column
      } 
      
      write.table(table_data, file, sep = "\t", row.names = FALSE)
    }
  )
  
  # Download button2 event handler
  output$download2 <- downloadHandler(
    
    # Specify the file name
    filename = function() {
      "latent_space.tsv"
    },
    
    # Specify the file content
    content = function(file) {
      save_latent <- data_result()
      
      # Write the data frame to the file in TSV format
      write.table(save_latent$result$latent, file, sep = "\t", row.names = FALSE)
    }
  )
  
  # Create a scatter plot based on result
  output$scatter_plot <- renderPlot({
    
    data_reunion <- data_result()
    
    if (!is.null(data_reunion)) {
      result <- data_reunion$result
      data <- data_reunion$data
      label <- data_reunion$label
      
      
      data_plot <- data.frame(
        x = data[, 1],
        y = result$pt,
        
        #cell_type = sample(c(unique(label)), nrow(data), replace = TRUE),
        pseudo_time = runif(nrow(data), min = 0, max = 5)
      )
      
      # If metadata exists, color the plot by cell type and add a legend
      if (!is.null(label)) {
        
        # Assuming the actual dataset has a "label" containing the correct cell types
        data_plot$cell_type <- factor(label, levels = c(unique(label)))
        
        # Define colors for each cell type
        cell_colors <- c(rainbow(length(unique(label))))
      } else{
        
        cell_colors <- NULL
      }
      
      # Call the function to generate and display the scatter plot
      generate_scatter_plot(data_plot, cell_colors)
    }
  })
  output$downloadplot <- downloadHandler (
    filename = function () {
      "plot.png"
    },
    content = function (file) {
      
      # Generate the data_plot and cell_colors based on the data_result()
      data_reunion <- data_result()
      
      if (!is.null(data_reunion)) {
        result <- data_reunion$result
        data <- data_reunion$data
        label <- data_reunion$label
        
        data_plot <- data.frame(
          x = data[, 1],
          y = result$pt,
          pseudo_time = runif(nrow(data), min = 0, max = 5)
        )
        
        if (!is.null(label)) {
          data_plot$cell_type <- factor(label, levels = c(unique(label)))
          cell_colors <- c(rainbow(unique(data)))
          
        } else {
          cell_colors <- NULL
        }
        
        # Generate the scatter plot using the function
        plot_to_save <- generate_scatter_plot(data_plot, cell_colors)
        
        # Save the plot as PNG
        png(file)
        print(plot_to_save)
        dev.off()
      }
    }
  )
}

# Run the app
shinyApp(ui, server)
=======
library(shiny)
library(ggplot2)
library(httr)
library(DT) 
library(scDHA)


# main analysis tool
data_analysis <- function(data_object) {
  
  data <- t(data_object$data)
  
  label <- as.character(data_object$label)
  
  data <- log2(data + 1)
  
  result <- scDHA(data, ncores = 2, seed = 1)
  
  cluster <- result$cluster
  
  mclust::adjustedRandIndex(cluster,label)
  
  result <- scDHA.pt(result, start.point = 1, ncores = 2, seed = 1)
  
  # Return a list containing both the processed data and the result
  return(list(data = data, result = result, label = label))
  
}


# Function to convert TSV to RDA
convert_tsv_to_rda <- function(input_file, output_file) {
  # Read TSV data using read.table
  tsv_data <- read.table(input_file, sep="\t", header=TRUE)
  
  # Save the data in RDA format using save function
  save(tsv_data, file = output_file)
  
  # Print a success message
  cat("TSV data converted and saved as RDA:", output_file, "\n")
}

# Create a scatter plot based on result
generate_scatter_plot <- function(data_plot, cell_colors) {
  
  if (!is.null(data_plot$cell_type)) {
  ggplot(data_plot, aes(x = data_plot$x, y = data_plot$y, color = data_plot$cell_type)) +
    geom_point() +
    scale_color_manual(values = cell_colors) +
    labs(title = "Scatter Plot of Trajectory Inference",
         x = "Sample Cells",
         y = "Pseudo Time",
         color = "Cell Type") +
    theme_minimal()
  } else {
    
    # If cell_type is NULL, color the points by pseudotime
    ggplot(data_plot, aes(x = data_plot$x, y = data_plot$y, color = data_plot$y)) +
      geom_point() +
      scale_color_viridis_c() + 
      labs(title = "Scatter Plot of Cell Representation",
           x = "Sample Cells",
           y = "Pseudo Time",
           color = "Pseudo Time") +
      theme_minimal()
  }
}


# UI
ui <- fluidPage(
  titlePanel("scDHA"),
  sidebarLayout(
    sidebarPanel(
      h4("Guideline Text"),
      radioButtons("choice", "If you want to upload your own data choose 'A'; otherwise if you want to select a built-in data set choose 'B':",
                   choices = c("A", "B")),
      
      conditionalPanel(
        condition = "input.choice == 'A'",
        fileInput("file1", "Upload a single-cell RNA-seq matrix in TSV format (required):"),
        fileInput("file2", "Upload a metadata file in TSV format (optional):")
      ),
      
      conditionalPanel(
        condition = "input.choice == 'B'",
        selectInput("data", "Choose a data source:",
                    choices = c("Yan Data", "Goolam Data", "Deng Data"))
      ),
      
      actionButton("submit", "Submit")
      
    ),
    
    mainPanel(
      downloadButton("downloadplot", "Download Plot as PNG"),
      downloadButton("download3", "Download Table as TSV"),
      downloadButton("download2", "Download Compressed Data as TSV"),
      uiOutput("unavailable_message"),
      plotOutput("scatter_plot"),
      tableOutput("my_table"),
      #tags$caption("Trajectory Inference")
    )
  )
)

# Server
server <- function(input, output) {
  # Reactive function to convert and display the table
  data_result <- reactive({
    if (input$submit > 0) {
      if (input$choice == "A") {
        # Get the file paths for the uploaded files
        input_tsv_file <- input$file1$datapath
        metadata_files <- input$file2$datapath
        
        # Check if file1 is uploaded
        if (!is.null(input_tsv_file) && !is.null(metadata_files)) {
          
          # combine data and metadata
          data_tsv <- cbind(input_tsv_file, metadata_files)
          
        } else if (!is.null(input_tsv_file) && is.null(metadata_files)){
          
          data_tsv <- input_tsv_file
          
        } else {
          
          # If no file is uploaded
          cat("Error: No uploaded files.")
          
          return(NULL)
        }
        
        # Define the local file path where the downloaded file is saved
        data_rda <- "output_file.rda"
        
        # Perform the conversion
        convert_tsv_to_rda(data_tsv, data_rda)
        nrow(data_tsv)
        
        # Load the RDA data
        rda_data <- load(data_rda)
        
        # Extract the data from the loaded object
        data_user <- get(rda_data)
        
        # Perform trajectory inference using the scDHA package
        list_data_user <- data_analysis(data_user)
        
        # Return the loaded data
        return(list_data_user)
        
      } else if (input$choice == "B") {
        # Process selected data source for option B
        selectedData <- input$data
        
        # Fetch online data based on the selected data source
        if (selectedData == "Yan Data") {
          
          output$unavailable_message <- renderUI({
            
            tags$p("Unfortunately, this dataset is currently unavailable.")
            
          })
          return(NULL)
        } else if (selectedData == "Goolam Data") {
          
          # Define the URL of the RDA file
          rda_file_url <- "https://raw.githubusercontent.com/duct317/scDHA/master/data/Goolam.rda"
          
          # Define the local file path where the downloaded file is saved
          local_file_path <- "Goolam.rda"
          
          # Download the RDA file from the URL to the local file path
          download.file(rda_file_url, destfile = local_file_path)
          
          # Load the data from the downloaded RDA file
          loaded_data <- load(local_file_path)
          
          # Extract the data from the loaded object
          data_to_save <- get(loaded_data)
          
          # Perform trajectory inference using the scDHA package
          Goolam_list <- data_analysis(data_to_save)
          
          # Return the loaded data
          return(Goolam_list)
          
        } else {
          cat("Error: Unable to download the file.")
        }
        
      } else if (selectedData == "Deng Data") {
        
        output$unavailable_message <- renderUI({
          
          tags$p("Unfortunately, this dataset is currently unavailable.")
          
        })
        
        return(NULL)
      } 
    } else {
      
      return(NULL)
      "Please select an option."
    }
  })
  
  # Render the table in the UI based on data
  output$my_table <- renderTable({
    
    data_reunion<- data_result()
    
    if (!is.null(data_reunion)) {
      data <- data_reunion$data
      result <- data_reunion$result
      label <- data_reunion$label
      
      # Create the table with cell ID and pseudotime columns
      table_data <- data.frame(Cell_ID = rownames(data), Pseudo_Time = result$pt)
      
      # If data is not NULL and metadata exists, add cell type column to the table
      if (nrow(table_data) > 0 && "label" %in% names(data_reunion)) {
        cell_types <- as.character(label)
        table_data$Cell_Type <- cell_types  # Insert the cell type column
      }
      
      return(table_data)
      
    } else {
      
      # Return an empty data frame if data is NULL
      return(data.frame())
    }
  })
  # Define the download handler
  output$download3 <- downloadHandler(
    
    # Specify the file name
    filename = function() {
      "table_data.tsv"
    },
    
    # Specify the file content
    content = function(file) {
      save_table <- data_result()
      
      # Create the table with cell ID and pseudotime columns
      table_data <- data.frame(Cell_ID = rownames(save_table$data), Pseudo_Time = save_table$result$pt)
      
      if (nrow(table_data) > 0 && !is.null(save_table$label)) {
        
        # If metadata exists, add cell type column to the table
        cell_types <- as.character(save_table$label)
        table_data$Cell_Type <- cell_types  # Insert the cell type column
      } 
      
      write.table(table_data, file, sep = "\t", row.names = FALSE)
    }
  )
  
  # Download button2 event handler
  output$download2 <- downloadHandler(
    
    # Specify the file name
    filename = function() {
      "latent_space.tsv"
    },
    
    # Specify the file content
    content = function(file) {
      save_latent <- data_result()
      
      # Write the data frame to the file in TSV format
      write.table(save_latent$result$latent, file, sep = "\t", row.names = FALSE)
    }
  )
  
  # Create a scatter plot based on result
  output$scatter_plot <- renderPlot({
    
    data_reunion <- data_result()
    
    if (!is.null(data_reunion)) {
      result <- data_reunion$result
      data <- data_reunion$data
      label <- data_reunion$label
      
      
      data_plot <- data.frame(
        x = data[, 1],
        y = result$pt,
        
        #cell_type = sample(c(unique(label)), nrow(data), replace = TRUE),
        pseudo_time = runif(nrow(data), min = 0, max = 5)
      )
      
      # If metadata exists, color the plot by cell type and add a legend
      if (!is.null(label)) {
        
        # Assuming the actual dataset has a "label" containing the correct cell types
        data_plot$cell_type <- factor(label, levels = c(unique(label)))
        
        # Define colors for each cell type
        cell_colors <- c(rainbow(length(unique(label))))
      } else{
        
        cell_colors <- NULL
      }
      
      # Call the function to generate and display the scatter plot
      generate_scatter_plot(data_plot, cell_colors)
    }
  })
  output$downloadplot <- downloadHandler (
    filename = function () {
      "plot.png"
    },
    content = function (file) {
      
      # Generate the data_plot and cell_colors based on the data_result()
      data_reunion <- data_result()
      
      if (!is.null(data_reunion)) {
        result <- data_reunion$result
        data <- data_reunion$data
        label <- data_reunion$label
        
        data_plot <- data.frame(
          x = data[, 1],
          y = result$pt,
          pseudo_time = runif(nrow(data), min = 0, max = 5)
        )
        
        if (!is.null(label)) {
          data_plot$cell_type <- factor(label, levels = c(unique(label)))
          cell_colors <- c(rainbow(unique(data)))
          
        } else {
          cell_colors <- NULL
        }
        
        # Generate the scatter plot using the function
        plot_to_save <- generate_scatter_plot(data_plot, cell_colors)
        
        # Save the plot as PNG
        png(file)
        print(plot_to_save)
        dev.off()
      }
    }
  )
}

# Run the app
shinyApp(ui, server)
>>>>>>> 0816ad43f62475fb3ea41a68168edf52be9eb0a0
