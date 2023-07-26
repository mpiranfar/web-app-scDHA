# Web App scDHA 

This is a Shiny web application that performs single-cell RNA-seq data analysis using the scDHA package.



## How to Run the App

1. Make sure you have R and Torch package and libtorch library installed on your computer.
libtorch can be installed using: torch::install_torch()

2. Install the required R packages by running the following in RStudio:

   install.packages(c("scDHA", "shiny", "ggplot2", "httr", "DT", "mclust", "viridis"))

3. Open the web-app.R file in RStudio.

4. Load the script by selecting all the code and clicking "Run."

5. The app will launch in your default web browser.

6. Use the app to analyze your single-cell RNA-seq data utilizing the scDHA package. You can either upload your own data in TSV format or choose a built-in dataset.

7. Click the "Submit" button to perform the analysis.

8. The results will be displayed as a scatter plot and a table.

9. Click the "Download Table as TSV" button to download the table data. 

10. Click the "Download Compressed Data as TSV" button to download the latent space data. 

11. Click the "Download Plot as PNG" button to download the plot as a png file. 


Built-in Datasets
The app provides the following built-in datasets:

Yan Data
Goolam Data
Deng Data

Note: Some datasets may not be currently available.