



##-----------------
##  [01] ADD EXCEL
##-----------------
## FUNCTIONS
## writes data to excel worksheet
## wb=existing wb object or if NULL then a new wb is created using createWorkbook()
## sheet="mysheetname" or sheet=1, if NULL then Sheet1,Sheet2 etc.
## x=data to write, rownames are included
## open=logical to open the wb after data is written (this is a temp file) T wb is opened F=wb is not
## wb object is returned if wb is NULL then use mywb <- add_excel() to store results
## in mywb rather than wb which would overwrite existing wb in GlobalEnv
## USE: wb<-createWorkbook()
##      add_excel(wb=wb,sheet="mysheet", x=mydf, open=F)
##      wb2<-add_excel(wb=NULL,sheet=NULL,x=mydf,open=F)
##      saveWorkbook(wb2, file="workbookname.xlsx",overwrite=TRUE,returnValue=TRUE)
## "http://www.ensembl.org/id/"                           ## ensembl
## "http://www.ncbi.nlm.nih.gov/gene/"                    ## entrez gene
## "https://www.uniprot.org/uniprot/"                     ## uniprot
## "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=" ## Accession/mature_Acc
## "http://www.mirbase.org/cgi-bin/mirna_summary.pl?fam=" ## FamilyAccession
## "http://www.mirbase.org/textsearch.shtml?q="           ## search all
## "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="   ## Precursor_Acc
# ens.url ="http://www.ensembl.org/id/"
# enz.url ="http://www.ncbi.nlm.nih.gov/gene/"
# uni.url = "https://www.uniprot.org/uniprot/"
# mb.mat.acc.url = "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc="
# mb.fam.acc.url = "http://www.mirbase.org/cgi-bin/mirna_summary.pl?fam="
# mb.pre.acc.url = "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="
# mb.sch.all.url = "http://www.mirbase.org/textsearch.shtml?q="
# pub.url = "https://pubmed.ncbi.nlm.nih.gov/"
#' Add data to Excel worksheet
#'
#' Internal function to add data to the cells of an Excel worksheet. The function
#' can also be used to add hyperlinks to designated column cells of the Excel output of the  to a column in a data frame to
#' be exported to an excel file. Implicitly, works only on one column.
#'
#'Excel hyperlink formulas will be added to the values in each column a Hyperlinkdefined in the url argument.
#' make_excel_hyperlinks
#' 
#' 
#' @param data The data such as a data.frame, matrix, vector, or string of text 
#'  to be written to an Excel worksheet.(see sheet argument).
#' @param url.cols a vector of column names in data that you want to add hyperlinks
#'  to. Here data should be a data.frame or matrix. The values in a particular 
#'  column are typically an ensembl/entrez gene id, uniprot id, miRBase id etc.
#' @param urls a character vector of identical length to url.cols naming the
#'  resource url for each url.cols value. The urls and url.cols values provided
#'  are used to create Excel hyperlink formulas that will be added to the workbook   
#'  See some examples below for ensembl gene_id, NCBI's entrez id, UniProt, miRBase,
#'  and PubMed 
#' Ensembl: "http://www.ensembl.org/id/"                           
#' NCBI GeneID: "http://www.ncbi.nlm.nih.gov/gene/" 
#' UniProt:  "https://www.uniprot.org/uniprot/"  
#' miRBase microRNA Precursor_Acc: "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc="
#' miRBase microRNA Accession/mature_Acc: "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=" 
#' miRBase microRNA FamilyAccession: "http://www.mirbase.org/cgi-bin/mirna_summary.pl?fam=" 
#' miRBase microRNA general search: "http://www.mirbase.org/textsearch.shtml?q=" 
#' PubMed: "https://pubmed.ncbi.nlm.nih.gov/" 
#' @param wb A workbook object. The default is NULL which creates a new workbook
#'  object wb in the global environment. If wb already exists it is overwritten. 
#' @param sheet character string of lenght 1 specifying the name of a new or existing
#'  worksheet in or to be added to to the wb workbook object. The default is NULL,
#'  which adds a new worksheet with Excel's default naming conventions
#'  e.g. Sheet1. If Sheet1 exists Sheet2 is added and so on.
#' @param rowNames TRUE/FALSE logical. If TRUE, row names of data are written to
#'  the sheet as well. Default is FALSE
#' @param tabColour character string of length 1 specifiying the worksheet tab 
#' color hex code. The default is gray 
#' i.e. "#D8D8D8"
#' @param fgFill character string of length 1 specifying the column heading color 
#' i.e. the background color (hex code) for the column cells in row 1 of the data 
#' written to the worksheet. The default is gray i.e. "#D8D8D8",
#' @param startCol numeric value of length 1 specifying the starting column to 
#' begin writting the data. Default is column 1.
#' @param startRow numeric value of length 1 specifying the starting row to begin 
#' writting the data. Default is row 1.
#' @param withFilter TRUE/FALSE logical or NA. If TRUE or NA Excel's builtin filter
#'  will be added to the columns in row 1.
#' @param save TRUE/FALSE logical. If TRUE after writing data the wb workbook object
#'  is saved as an .xlsx file, specified
#'  in the file argument.
#' @param file A character string of length 1 naming the path/to/file.xlsx file. 
#' The default is NULL. If file is NULL and save = TRUE then a file named 
#' "workbook.xlsx" is saved in the project directory based on the overwrite parameter.
#' If save = FALSE this parameter is ignored.
#' @param overwrite TRUE/FALSE logical. If TRUE overwrite any existing file.
#' Default is FALSE
#' @param open TRUE/FALSE logical. If TRUE the wb workbook object is opened in 
#' Excel (temporary file does not require file to be saved and does not save the file)
add_excel <- function(data,
                      url.cols = NULL,
                      urls = NULL,
                      wb = NULL,
                      sheet = NULL,
                      tabColour = "#d8d8d8",
                      fgFill = "#d8d8d8",
                      startCol = 1,
                      startRow = 1,
                      rowNames = FALSE,
                      withFilter = FALSE,
                      save = FALSE,
                      file = NULL,
                      overwrite = FALSE,
                      open = FALSE){
    
    # library(openxlsx)
    ## add hyperlink formulas to columns of input data.frame
    ## the columns of the output data.frame returned contains excel hyperlink formula text
    add.hyper <- FALSE
    if(all(!is.null(url.cols) & !is.null(urls))){
        if(all(url.cols %in% colnames(data))){
            if(length(url.cols)==length(urls)){
                for(i in 1:length(url.cols)){
                    data <- make_excel_hyperlinks(data=data, url.col=url.cols[i], url=urls[i])
                    print(paste(url.cols[i],":", urls[i], "added ..."))
                }
                add.hyper <- TRUE
            }}
    }
    
    
    ## if wb doesn't exist a new workbook object is created. The sheet argument is then used to determine
    ## the name of the worksheet to add. If sheet is NULL the name of the worksheet is defined as Sheet#,
    ## where # = number worksheets + 1. If sheet index is supp
    ## If a worksheet with input sheet name does not exist
    ## a new worksheet to the wb and write data to worksheet using the startCol and startRow
    ## values
    ## truncate if worksheet name is > 30 characters
    ## if
    # library(openxlsx)
    if(is.null(wb)){ wb <- openxlsx::createWorkbook() }
    if(is.null(sheet)){ sheet <- paste("Sheet",(length(names(wb))+1),sep="") };sheet
    if(is.numeric(sheet)){ sheet <- ifelse(sheet %in% names(wb), paste0("Sheet",sheet), sheet) };sheet
    sheet <- ifelse(nchar(sheet) > 30, substring(sheet,1,30), sheet);sheet
    if(sheet %in% names(wb) == FALSE){ openxlsx::addWorksheet(wb, sheetName=sheet, tabColour=tabColour) }
    names(wb)
    
    ## add column header style (grey background, filter, with increased row height)
    colStyle <- openxlsx::createStyle(fgFill = fgFill, halign="center",
                                      valign="center", border="TopBottomLeftRight",
                                      textDecoration="bold", fontColour="black", 
                                      fontSize=11, numFmt="TEXT")
    
    openxlsx::writeData(wb = wb, 
                        sheet = sheet, 
                        x = data,
                        startCol = startCol,
                        startRow = startRow,
                        colNames = TRUE,
                        rowNames = rowNames,
                        borders = "all",
                        keepNA = TRUE, 
                        na.string = "NA",
                        headerStyle = colStyle,
                        withFilter = withFilter,
                        sep = "\t")
    
    if(add.hyper==TRUE){
        if(rowNames == TRUE){
            tmp=colnames(data);names(tmp)<-1:length(colnames(data))
            hyper.cols <- as.numeric(names(tmp)[tmp%in%url.cols]) + startCol;hyper.cols
        }
        if(rowNames==FALSE){
            tmp=colnames(data);names(tmp)<-1:length(colnames(data))
            hyper.cols <- as.numeric(names(tmp)[tmp%in%url.cols]) + (startCol-1);hyper.cols
        }
        
        hyperlinkStyle <- openxlsx::createStyle(fontColour="#0000FF",halign="left",
                                                valign="center",textDecoration="underline")
        openxlsx::addStyle(wb=wb, sheet=sheet, style=hyperlinkStyle, cols=hyper.cols, 
                           rows=(startRow+1):(nrow(data)+startRow),
                           gridExpand=TRUE, stack=TRUE)
    }
    openxlsx::setRowHeights(wb=wb, sheet=sheet, rows = 1, heights = 25)
    
    openxlsx::activeSheet(wb = wb) <- sheet
    cli::cli_inform("data was added to worksheet: {.val {sheet}}")
    if(save==TRUE){
        if(is.null(file)){ file <- "workbook.xlsx"}
        openxlsx::saveWorkbook(wb=wb, file=file, overwrite=overwrite, returnValue = TRUE)
        cli::cli_inform("The {.file {file}} was saved.")
    }
    if(open==TRUE){ openxlsx::openXL(wb) }
    
    
    return(invisible(wb))
    
} ## ADD EXCEL




#' Make hyperlinks for an excel column
#'
#' Internal functions to add hyperlinks to a column in a data frame to
#' be exported to an excel file. Implicitly, works only on one column.
#'
#' @param data The data frame in which to add hyperlinks
#' @param url.col The name of the column to add hyperlinks to
#' @param url The url that will be prepended to the info in the column
#'
#' @return The original dataframe, now with hyperlinks in the desired column
#'
#' @examples
#' # No examples yet
make_excel_hyperlinks <- function(data, url.col, url) {
    
    ids <- data[,url.col]
    tmp <- is.na(ids)
    url2 <- paste0(url, ids)
    ids2 <- paste0("HYPERLINK(\"", url2, "\", \"", ids, "\")")
    ids2[tmp] <- NA
    data[,url.col] <- ids2
    class(data[,url.col]) <- "formula"
    
    return(data)
}

