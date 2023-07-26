# web-app-scDHA
This is a Shiny web application that performs single-cell RNA-seq data analysis using the scDHA package.

## How to Run the App
1. Make sure you have R and the scDHA package installed on your computer.
To see how to install the scDHA package see: https://github.com/duct317/scDHA

2. Install the other required R packages by running the following in RStudio:
3. 
   install.packages(c("shiny", "ggplot2", "httr", "DT", "mclust", "viridis"))

4. Open the web-app.R file in RStudio.

5. Load the script by selecting all the code and clicking "Run."

6. The app will launch in your default web browser.

7. Use the app to analyze your single-cell RNA-seq data utilizing the scDHA package. You can either upload your own data in TSV format or choose a built-in dataset.

8. Click the "Submit" button to perform the analysis.

9. The results will be displayed as a scatter plot and a table.

10. Click the "Download Table as TSV" button to download the table data. 

11. Click the "Download Compressed Data as TSV" button to download the latent space data. 

12. Click the "Download Plot as PNG" button to download the plot as a png file.


# Built-in Datasets
The app provides the following built-in datasets:

* Yan Data

* Goolam Data

* Deng Data

Note: Some datasets may not be currently available.
