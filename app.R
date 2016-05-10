# Emilia F Gan
# Student Number: 1138607
# Email: efgan@uw.edu
# Copyright E Gan 2016
# Submitted in partial fulfillment of the requirements
# of the degree Master of Science

# Box Plot Application

# Opening a CSV File and reading into R
# Modeled after: http://shiny.rstudio.com/gallery/upload-file.html
# Downloading plots learned from: https://www.youtube.com/watch?v=LSnWGmVkB6A

library(shiny)
library(ggplot2) # for plotting color plots with legend
library(edgeR) # for DGE  
library(limma)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 9*1024^2)

#--------------------------------------------------- UI CODE --------------------------------------------------------------
ui <- fluidPage(
    titlePanel("Boxplot Analysis"),
    hr(),
    tabsetPanel(
        tabPanel("Boxplots",
                 sidebarLayout(
                     # These elements will display in the sidebar panel at the L of the page.
                     sidebarPanel(
                         # Allows user to select data file from local file system.
                         fileInput('file1', 'Choose file to upload',
                                   accept = c(
                                       'text/csv',
                                       'text/comma-separated-values',
                                       'text/tab-separated-values',
                                       'text/plain',
                                       '.csv',
                                       '.tsv'
                                   )
                         ),
                         tags$hr(),
                         
                         # Indicates whether or not the data file has column headers.
                         # Note: Plot will only display color-coded groups if the header "Type" is present
                         # when data file is formatted with samples displayed in rows.
                         radioButtons('sep', 'Separator:',
                                      c(Comma=',',
                                        Semicolon=';',
                                        Tab='\t'),
                                      ','),
                         tags$hr(),
                         
                         # Collects sample formatting data for proper analysis.
                         # For PCA function, samples should be listed one per row.
                         radioButtons('dir', 'Genes are listed in:',
                                      c(Rows = 'R',
                                        Columns = 'C'),
                                      'R'),
                         
                         tags$hr(),
                         
                         # Determines whether counts file represents normalized or unnormalized data.
                         # Normalized data are not processed further prior to PCA.
                         # Unnormalized data are normalized by the app using the TMM method as recommended
                         # by Peter Linsley of BRI using function written by Kristen Dang.
                         radioButtons('norm', 'Apply TMM Normalization:',
                                      c(No = 'No_Norm',
                                        Yes = 'Norm'),
                                      'No_Norm'),
                         tags$hr(),
                         
                         # Allows user to input a custom title for the plot.
                         # If left blank, no title will be displayed.
                         textInput('title', "Plot Title", ''),
                         
                         # Suppresses error messages about NULL while awaiting user selection.
                         # Credit to: Joe Cheng[R Studio] at
                         # https://groups.google.com/forum/#!topic/shiny-discuss/FyMGa2R_Mgs
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: hidden; }",
                                    ".shiny-output-error:before { visibility: hidden; }"
                         ),
                         width = 3
                     ),
                     
                     # These elements will display in the main panel (R side of page).
                     mainPanel(
                         uiOutput('outSelections'),
                         tags$hr(),
                         plotOutput('boxPlot'),
                         downloadButton('download1', 'Download', class = "buttCol"),
                         tags$hr(),
                         fluidRow(
                             column(4,
                                    radioButtons('customY', "Set Custom Y Range",  c(No = 'NotCustom', Yes = 'Custom'), 'NotCustom', inline = TRUE)
                             ),
                             column(7, offset = 1,
                                    conditionalPanel(
                                    condition = "input.customY == 'Custom'",
                                    uiOutput('yMaxSelector')
                                    )
                                  
                             )
                         ),  # end fluidRow
                         conditionalPanel(
                             condition = "input.customY == 'Custom'",               
                             plotOutput('boxPlotCustomY'),
                             downloadButton('download2', 'Download', class = "buttCol")
                         ),
                         tags$head(tags$style(".buttCol{background-color:#faf0e6;} .buttCol{color: black;}")),
                         tags$hr(),
                         tableOutput('table'),
                         tags$hr()
                     ) # end mainPanel
                 ) # end sidebarLayout
                 
        ), # end tabPanel
        tabPanel("Help",
                 tags$h1("File Formatting Issues", style = "color:darkblue"),
                 tags$p(" It is anticipated that most issues that arise while trying to use this application will be related
                        to file formatting. A very specific file format is required. This format has several variations, 
                        depending on the type of data included in the file. The figures below illustrate accepatable file formats, 
                        as they would appear in a spreadsheet program (such as Excel) before being coverted into one of the accepted
                        file formats (.csv, .tsv, or .txt).", style="font-size:20px"),
                 tags$img(height = 443, width = 904, src = "Permitted.PNG"),
                 tags$p(" In the figure above, genes are arranged by column and samples by row. Transposing this arrangement,
                        so that genes are in rows, while samples are in columns is permitted and should not cause any issues 
                        in the functioning of the application.", style="font-size:20px")
                 )  # end tabPanel
                 )  # end tabsetPanel
                 )  #end fluidPage


#----------------------------------------------- SERVER CODE ----------------------------------------
server <- function(input,output){
    # Function to read data from a file, after checking file in not NULL.
    getData <- function(dataFile){
        inFile <- dataFile
        if (is.null(inFile)) {
            return(NULL)
        }
        else {
            read.csv(inFile$datapath, header = TRUE, sep = input$sep, quote = input$quote, row.names = 1, stringsAsFactors = FALSE)
        }
    }
    
    # Read in the data from the input file.
    myData <- reactive({getData(input$file1)})
    
    # Get data dimensions as 2 item list [#rows, rcols].
    dataDimensions <- reactive(dim(myData()))
    
    # Get information on whether sample data are in ROWS or COLUMNS.
    direction <- reactive({input$dir})
    
    #-----------------------------------------------------------------------------------------
    #           TASK ONE: FORMAT THE DATA
    #-----------------------------------------------------------------------------------------
    
    # Get the data into a "tidier" format with row names identified as such and
    # Sample/Gene identifiers separated from the numerical data table elements.
    
    # This section of code checks whether data are arranged with genes or samples in rows and
    # separates any type data present from the counts data.
    
    # 1) Check that ROW TYPE DATA is present 
    # i.e. Last COLUMN contains TYPE data
    hasGroupDataInLastCol <- function(df){
        if(!is.null(df$Type)){
            return(TRUE)
        }
        return(FALSE)
    }
    
    # 2) Check that COLUMN TYPE DATA is present 
    # i.e. Last ROW contains TYPE data
    hasGroupDataInLastRow <- function(df){
        numRows = dataDimensions()[1]
        if(row.names(df)[numRows] == "Type"){
            return(TRUE)
        }
        return(FALSE)
    }
    
    # 3) Now, we're ready to do the actual reformatting:
    # Remove the first column from the data frame -- it should contain text (gene or sample names).
    # If it contains data, the input file is NOT PROPERLY FORMATTED!
    # At the end of this code block, only numerical values should be in the body of the data frame.
    # The data frame's column and row names should now be appropriately recognized as such by R.
    formattedData <- function(){
        
        dimensions <- dim(myData())
        numRows <- dimensions[1]
        numCols <- dimensions[2]
        
        preformatted <- myData()
        # If there is GROUP data in the last column, it needs to be removed.
        if(isTRUE(hasGroupDataInLastCol(preformatted)) == TRUE) {
            numCols = numCols - 1
        }
        # If there is group data in the last row, it needs to be removed.
        if(isTRUE(hasGroupDataInLastRow(preformatted)) == TRUE){
            numRows = numRows - 1
        }
        # reformat
        preformatted <- preformatted[1:numRows, 1:numCols]
        # get the row names and store them as "r_names"
        r_names <- row.names(preformatted)
        
        preformatted <- lapply(preformatted, as.numeric)
        preformatted <- as.data.frame(preformatted)
        
        # Add back the row name info as formal row names since they got lost 
        # in the two steps above.
        row.names(preformatted) <- r_names
        if(input$dir == 'C'){
            preformatted <- t(preformatted)
        }
        return(preformatted)
    }
    
    #-----------------------------------------------------------------------------------------
    #           TASK TWO: BE READY TO NORMALIZE DATA IF THIS OPTION IS SELECTED
    #-----------------------------------------------------------------------------------------
    
    # Normalize the counts data using TMM, if user selects this option. 
    # TMM normalization requires that genes be in ROWS -- This is very important for correct functioning
    # of the TMM algorithm. In addition to incorrect results, improper orientation could cause the program
    # to crash if a gene has zero counts across all samples.
    normalized <- function(){
        
        counts_m <- formattedData()
        
        dge <- DGEList(counts=counts_m)
        # can now use TMM to normalize the counts data
        # this function Calculates normalization factors to scale the raw library sizes
        # using the weighted trimmed mean of M-values (to the reference), 
        # where the weights are from the delta method on Binomial data
        data.dge = calcNormFactors(dge, method=c("TMM"), refColumn = NULL,logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
        #### NEXT: Use this function from Kristen Dang ####
        # This function uses the library size from the datafile and the
        # normalization factors calculated to return a matrix of normalized counts
        # See DGEList documentation: default for lib.size = colSums(counts)
        # default norm.factors = rep(1, ncol(counts)) (i.e. a vector of '1's)
        return_TMM_matrix_noRound=function(inDGEObject){
            eff.libsize = inDGEObject$samples$lib.size * inDGEObject$samples$norm.factors
            tndata = 1e6*t(inDGEObject$counts)/eff.libsize
            return(t(tndata))
        }
        
        # Call the function, then convert matrix obtained to a data frame:
        data.tmm = return_TMM_matrix_noRound(data.dge) #### A function from Kristen Dang
        data.tmm = as.data.frame(data.tmm)

        return(data.tmm)
    }
    
    #-----------------------------------------------------------------------------------------
    #           TASK THREE: PREPARE GENE NAMES AND SLIDER INPUT MAX VALUE
    #-----------------------------------------------------------------------------------------
    
    # Get gene names for drop down menu from file read in to the app.
    output$outSelections <- renderUI({
        genes <- row.names(formattedData())
        selectInput("genes", "Please click the box below for drop-down menu to select genes to plot:", choices = genes, multiple=TRUE, selectize=TRUE)
        
    })
    
    # Get maximum count value in data to use as y-max for slider input widget.
    output$yMaxSelector <- renderUI({
        if(input$norm == 'Norm'){
            maxNorm <- round(max(normalized()), 0)
            #print("Printing maxNorm:")
            #print(maxNorm)
            sliderInput('yRange', "To change range of y-axis, enter a new maximum value", 0, maxNorm, value = maxNorm)
        }
        else {
            maxForm <- round(max(formattedData()), 0)
            #print("Printing maxForm")
            #print(maxForm)
            sliderInput('yRange', "To change range of y-axis, enter a new maximum value", 0, maxForm, value = maxForm)
        }
    })
    
    #-----------------------------------------------------------------------------------------
    #           TASK FOUR: CREATE BOXPLOT INPUTS
    #-----------------------------------------------------------------------------------------
    
    # Create the boxplots
    boxPlotInput <- function(){
        
        # Build the y-axis label
        if(input$norm == 'Norm'){
            normalBool = "Normalized"
            dataForPlot <- normalized()
        }
        else {
            normalBool = ""
            dataForPlot <- formattedData()
        }

        # Combine y-axis label components
        yaxisLabel = paste(normalBool," Gene Expression Counts ", sep="")
        
        dataForPlot <- t(dataForPlot)
        dataToPlot <- dataForPlot[,input$genes]

        boxplot(dataToPlot, col="lightblue", pars = list(boxwex = 0.6), main=input$title, xlab="Genes", ylab = yaxisLabel) 
    }
    
    # Creates boxplot with customized y-axis range.
    boxPlotCustomYInput <- function(){
        # store the data to plot
        dataForPlot <- formattedData()
        # Build the y-axis label
        if(input$norm == 'Norm'){
            normalBool = "Normalized"
            dataForPlot <- normalized()
        }
        else {
            normalBool = ""
            dataForPlot <- formattedData()
        }
        
        # Combine y-axis label components
        yaxisLabel = paste(normalBool," Gene Expression Counts ", sep="")
        dataForPlot <- t(dataForPlot)
        dataToPlot <- dataForPlot[,input$genes]

        boxplot(dataToPlot, col="lightblue", pars = list(boxwex = 0.6), main=input$title, xlab="Genes", ylab = yaxisLabel, ylim = c(0, input$yRange))
    }
    
    #-----------------------------------------------------------------------------------------
    #           TASK FIVE: DISPLAY BOXPLOTS
    #-----------------------------------------------------------------------------------------
    
    output$boxPlot <- renderPlot({
        boxPlotInput()
    })
    
    output$boxPlotCustomY <- renderPlot({
        boxPlotCustomYInput()
    })
    
    #-----------------------------------------------------------------------------------------
    #           TASK SIX: DOWNLOAD PLOTS
    #-----------------------------------------------------------------------------------------
    
    # Download handler for boxplot plot no scaling (top plot).
    output$download1 <- downloadHandler(           
        filename = 'downloadedPlot.pdf',
        content = function(file) {
            pdf(file)
            dataForPlot <- formattedData()
            # Build the y-axis label
            if(input$norm == 'Norm'){
                normalBool = "Normalized"
                dataForPlot <- normalized()
            }
            else {
                normalBool = ""
                dataForPlot <- formattedData()
            }

            # Combine y-axis label components
            yaxisLabel = paste(normalBool," Gene Expression Counts ", sep="")
            dataForPlot <- t(dataForPlot)
            dataToPlot <- dataForPlot[,input$genes]

            boxplot(dataToPlot, col="lightblue", main=input$title, xlab="Genes", ylab = yaxisLabel)
            dev.off()
        }
    )
    
    # Download handler for boxplot plot with scaling (bottom plot).
    output$download2 <- downloadHandler( 
        filename = 'downloadedPlot.pdf',
        content = function(file) {
            pdf(file)
            dataForPlot <- formattedData()
            # Build the y-axis label
            if(input$norm == 'Norm'){
                normalBool = "Normalized"
                dataForPlot <- normalized()
            }
            else {
                normalBool = ""
                dataForPlot <- formattedData()
            }
            
            # Combine y-axis label components
            yaxisLabel = paste(normalBool," Gene Expression Counts ", sep="")
            dataForPlot <- t(dataForPlot)
            dataToPlot <- dataForPlot[,input$genes]

            boxplot(dataToPlot, col="lightblue", main=input$title, xlab="Genes", ylab = yaxisLabel, ylim = c(0, input$yRange))
            dev.off()
            
        }
    )
    
    #-----------------------------------------------------------------------------------------
    #           TASK SEVEN: DISPLAY DATA FILE
    #-----------------------------------------------------------------------------------------
    
    # Display data table below plot.
    # Data are displayed in the same format as they are in the input file.
    output$table <- renderTable(
        if(input$dir == 'C'){
            t(formattedData())
        }
        else {
            formattedData()
        }
    )
}

shinyApp(ui = ui, server = server)