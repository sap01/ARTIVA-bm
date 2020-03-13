#!/usr/bin/env Rscript
## Goal of this program: Run ARTIVA algorithm [1, 2] using the ARTIVA package
## on the given datasets.
##
## Author: Saptarshi Pyne (saptarshipyne01@gmail.com)
## Last modified on: Nov 2, 2017
##
## How to execute this script:
## Let us assume that this script is inside directory '/home/saptarshi/R/R-3.3.2/projects/repoARTIVA' and
## the Rscript file is inside directory '/home/saptarshi/R/R-3.3.2/bin'.
## Then, execute this script using the following commands:
## $ cd /home/saptarshi/R/R-3.3.2/projects/repoARTIVA/  
## $ nohup time /home/saptarshi/R/R-3.3.2/bin/Rscript ARTIVA.R input.json &
## where 'asset/input.json' contains the user-defined parameters. A file 
## named 'nohup.out' will be generated inside 
## '/home/saptarshi/R/R-3.3.2/projects/repoARTIVA/'.
##
## Input: A time series gene expression dataset with multiple time series.
##
## Output: Time-varying Gene Regulatory Networks and a corresponding rolled up network.
##
## Note: When ARTIVA package is installed via the R command 'install.packages('ARTIVA')' using India (IIT Madras aka IITM) 
## CRAN mirror, ARTIVA::ARTIVAnet() failed for some values of its 'dataDescription' parameter. But then the 
## GitHub ARTIVA repo (https://github.com/cran/ARTIVA) is downloaded as a ZIP file (saved at asset/ARTIVA-master.zip) 
## and ARTIVA package is installed from
## this ZIP file, using the following steps:
## Step 1: If 'ARTIVA' package is already installed, uninstall it. For
## that purpose, open a R prompt from '/home/saptarshi/R/R-3.3.2/projects/repoARTIVA'
## and issue the following command:
## > remove.packages('ARTIVA')
## Step 2: Unzip the ZIP file.  A directory named 'ARTIVA-master' will be created.
## Step 3: Enter the following commands in a R prompt:
## > library(devtools) ## If package 'devtools' is not already installed, please install it using 
## install.packages('devtools')
## > source <- devtools:::source_pkg('asset/ARTIVA-master')
## Please note that there are triple colons in the above mentioned command.
## > install(source)
## After completing the steps, the ARTIVA package is installed. ARTIVA::ARTIVAnet() does not fail anymore
## for those values of its 'dataDescription' parameter for which it was earlier failing.
## One of the probable reason is that at least one of the files among 'init.R' and 'main.R' differs in 
## the GitHub CRAN repo from the IITM CRAN mirror.

## Remove all objects in the current workspace
rm(list = ls())

##------------------------------------------------------------
## Begin: Load the Required Packages
##------------------------------------------------------------
## For reading from and writing to '.json' files
library(rjson)

library(ARTIVA)
##------------------------------------------------------------
## End: Load the Required Packages
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Read User-defined input Params
##------------------------------------------------------------
input.args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1)
{
  stop("Exactly one input file must be supplied.", call.=FALSE)
}

input.params <- rjson::fromJSON(file = paste(getwd(), 'asset', input.args, sep = '/'))
rm(input.args)

## Input file for time-series gene expression data
input.data.filename <- input.params$input.data.filename
input.data.filename <- paste(getwd(), 'asset', input.data.filename, sep = '/')

## Number of time points (T)
num.timepts <- input.params$num.timepts

## Number of time series (S)
num.time.series <- input.params$num.time.series

## True rolled network file.
## If 'true.net.filename' is an empty string, then the true rolled network
## is not known a prior. Otherwise, it is known a prior and would be used
## to evaluate performance metrics of the predicted rolled network.
true.net.filename <- input.params$true.net.filename
if (true.net.filename != '')
{
  true.net.filename <- paste(getwd(), 'asset', true.net.filename, sep = '/')
}

## Maximum num of regulators for a gene in a time segment. 
## Please see the 'maxPred' input parameter of the 
## ARTIVA::ARTIVAnet() function.
max.fanin <- input.params$max.fanin

## Num of RJMCMC iterations.
## Please see the 'niter' input param of the 
## ARTIVA::ARTIVAnet() function.
num.itns <- input.params$num.itns

## Set maximum number of change points.
## Please see the 'maxCP' input param of the 
## ARTIVA::ARTIVAnet() function.
## "For long temporal courses (more than 20 time points), we advise - for computational 
## reasons - to limit the maximal number of CP to 15" (p. 4, 'maxCP')[3].
max.num.cp <- input.params$max.num.cp

## Posterior probability threshold for the selection of edges.
## Only the edges with posterior probability >= edge.post.prob.threshold,
## will be selected to generate the final network.
## Please see the 'edgesThreshold' input param of the
## ARTIVA::ARTIVAnet() function.
edge.post.prob.threshold <- input.params$edge.post.prob.threshold

rm(input.params)
##------------------------------------------------------------
## End: Read User-defined input Params
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Create the output directory
##------------------------------------------------------------
init.path <- getwd()

## Output directory name
output.dirname <- paste('asset/output', format(Sys.time(), "%Y%m%d%H%M%S"), sep = '')
output.dirname <- paste(init.path, output.dirname, sep = '/')

if(.Platform$OS.type == "unix") {
  if(! output.dirname %in% system("ls" ,intern=TRUE))
  {
    system(paste('mkdir ', output.dirname, sep = ''))
  }
} else{# if(.Platform$OS.type == "unix"){
  shell(paste('mkdir ', output.dirname, sep = ''), intern = TRUE, mustWork =NA)
}
##------------------------------------------------------------
## End: Create the output directory
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Load the Required External Functions
##------------------------------------------------------------
source(paste(init.path, 'calcPerfDiNet.R', sep = '/'))
source(paste(init.path, 'RDataToCytoscape.R', sep = '/'))
##------------------------------------------------------------
## End: Load the Required External Functions
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Main program
##------------------------------------------------------------

## Print the output dir name in 'nohup.out'
print('The output directory name is:')
print(output.dirname)
print('') ## to append a blank line

## Save console output in a file named 'output.txt' inside the output directory.
output.filename <- paste(output.dirname, 'output.txt', sep = '/')
output.file.conn <- file(output.filename, open = "wt")
sink(output.file.conn)

##------------------------------------------------------------
## Begin: Read input data file
##------------------------------------------------------------

## Begin: Find file extension of the input data file. Only '.tsv' and '.RData'
## are allowed.
## Split the string at every '.' and consider the last substring as the 
## file extension.
input.data.filename.ext <- unlist(strsplit(input.data.filename, '[.]'))
## End: Find file extension of the input data file. Only '.tsv' and '.RData'
## are allowed.

if (input.data.filename.ext[length(input.data.filename.ext)] == 'tsv')
{
  input.data <- read.table(input.data.filename, header = TRUE, sep="\t")
  
  ## Remove first col i.e. the time point names
  input.data <- input.data[, -1]
  
} else if (input.data.filename.ext[length(input.data.filename.ext)] == 'RData')
{
  ## Loads an object named input.data
  load(input.data.filename)
}

num.nodes <- ncol(input.data)
##------------------------------------------------------------
## End: Read input data
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Learn the ARTIVA model
##------------------------------------------------------------
start.time <- proc.time() # start the timer

## In targetData and parentData, rows = nodes, cols = time pts.
targetData <- t(input.data)
parentData <- targetData

node.names <- colnames(input.data)

## Ordering of the time points in the dataset.
## Please see input param 'dataDescription' of
## function ARTIVA::ARTIVAnet().
time.pt.desc <- rep(1:num.timepts, num.time.series)

## Please see input param 'outputPath' of
## function ARTIVA::ARTIVAnet().
output.path <- paste(output.dirname, 'ARTIVAnet', sep = '/')

artiva.out <- ARTIVA::ARTIVAnet(targetData, parentData, 
                        targetNames = node.names, parentNames = node.names, 
                        dataDescription = time.pt.desc,
                        saveEstimations = TRUE, 
                        saveIterations = FALSE,
                        savePictures = TRUE, 
                        outputPath = output.path, 
                        dyn = 1, 
                        segMinLength = 2, 
                        maxCP = max.num.cp, 
                        maxPred = max.fanin, 
                        nbCPinit = NULL, 
                        CPinit = NULL, 
                        niter = num.itns, 
                        burn_in = NULL, 
                        PSRFactor = FALSE, PSRF_thres = 1.05, segmentAnalysis = TRUE, 
                        edgesThreshold = edge.post.prob.threshold, 
                        layout = "fruchterman.reingold", 
                        cCP = 0.5, cEdges = 0.5, alphaCP = 1, betaCP = 0.5, 
                        alphaEdges = 1, betaEdges = 0.5, v0 = 1, gamma0 = 0.1, 
                        alphad2 = 2, betad2 = 0.2, 
                        silent = FALSE)

## At this point, 
## > artiva.out
## returns a set of time-varying edges (which is equivalent to a set of 
## time-varying network structures) in the form of a table. In this 
## table, each row represents an edge. Columns 'CPstart' and 'CPend' 
## defines the time segment when the
## corresponding edge is found. Column 'PostProb' represents the 
## posterior probability of that edge during that time segment.
## An example is given below.
## Parent Target CPstart CPend PostProb CoeffMean edgesThreshold
## 1 CG9536 CG9536       2    58   0.7845   0.0032            0.5
## 1 CG9536 CG9536       59    67   0.8828   0.33166            0.5
## 2 CG4920 CG9536       2    67   0.5312   0.03910            0.5
## 3 CG9536 CG4920       2    67   0.2677   0.00000            0.5
## 4 CG4920 CG4920       2    67   1.0000   0.82820            0.5
##
## The edges can be shortlisted based on their PostProb values, as shown below.
## > artiva.out[artiva.out$PostProb >= 0.5, ]
## Parent Target CPstart CPend PostProb CoeffMean edgesThreshold
## 1 CG9536 CG9536       2    58   0.7845   0.0032            0.5
## 2 CG9536 CG9536       59    67   0.8828   0.33166            0.5
## 3 CG4920 CG9536       2    67   0.5312   0.03910            0.5
## 5 CG4920 CG4920       2    67   1.0000   0.82820            0.5

## Retain only the edges with posterior probability >= edge.post.prob.threshold
## in the unrolled time-varying network structures
save(artiva.out, file = paste(output.dirname, 'artiva.out.RData', sep = '/'))
if (nrow(artiva.out) > 0)
{
  artiva.out <- artiva.out[artiva.out$PostProb >= edge.post.prob.threshold, ]
}
rm(edge.post.prob.threshold)
unrolled.dbn.adj.list <- artiva.out
save(unrolled.dbn.adj.list, file = paste(output.dirname, 'unrolled.dbn.adj.list.RData', sep = '/'))
rm(artiva.out)
##------------------------------------------------------------
## End: Learn the ARTIVA model
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Generate unrolled DBN adjacency matrices from unrolled
## DBN adjacency lists
##------------------------------------------------------------

## Initialize the list of unrolled DBN adjacency matrices
unrolled.dbn.adj.matrix.list <- vector(mode = 'list', length = num.timepts)

## Initialize each time pt specific adjacency matrix.
## Rows = source nodes, cols =  target nodes.
time.pt.spec.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes)
rownames(time.pt.spec.adj.matrix) <- node.names
colnames(time.pt.spec.adj.matrix) <- node.names
for (time.pt.idx in 1:length(unrolled.dbn.adj.matrix.list)) {
  unrolled.dbn.adj.matrix.list[[time.pt.idx]] <- time.pt.spec.adj.matrix
}
rm(time.pt.idx)
rm(time.pt.spec.adj.matrix)

num.edges <- nrow(unrolled.dbn.adj.list)
if (num.edges > 0) {
  for (edge.idx in 1:num.edges) {
    
    ## Using 'as.character()' is crucial. Otherwise, they 
    ## are considered as factors and incorrect outputs
    ## are observed.
    src.node.name <- as.character(unrolled.dbn.adj.list[edge.idx, 'Parent'])
    tgt.node.name <- as.character(unrolled.dbn.adj.list[edge.idx, 'Target'])
    
    start.time.pt <- unrolled.dbn.adj.list[edge.idx, 'CPstart']
    end.time.pt <- unrolled.dbn.adj.list[edge.idx, 'CPend']
    
    for (time.pt.idx in start.time.pt:end.time.pt) {
      adj.matrix <- unrolled.dbn.adj.matrix.list[[time.pt.idx]]
      adj.matrix[src.node.name, tgt.node.name] <- 1
      unrolled.dbn.adj.matrix.list[[time.pt.idx]] <- adj.matrix
    }
    rm(time.pt.idx)
    
  }
  rm(edge.idx)
}

save(unrolled.dbn.adj.matrix.list, file = paste(output.dirname, 'unrolled.dbn.adj.matrix.list.RData', sep = '/'))
##------------------------------------------------------------
## End: Generate unrolled DBN adjacency matrices from unrolled
## DBN adjacency lists
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Roll up the unrolled DBN to generate a rolled DBN adjacency matrix
##------------------------------------------------------------

## Initialize the rolled DBN adjacency matrix.
## Rows = source nodes, cols =  target nodes.
rolled.dbn.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes)
rownames(rolled.dbn.adj.matrix) <- node.names
colnames(rolled.dbn.adj.matrix) <- node.names  

num.edges <- nrow(unrolled.dbn.adj.list)
if (num.edges > 0)
{
  for (edge.idx in 1:num.edges)
  {
    ## Using 'as.character()' is crucial. Otherwise, they 
    ## are considered as factors and incorrect outputs
    ## are observed.
    src.node.name <- as.character(unrolled.dbn.adj.list[edge.idx, 'Parent'])
    tgt.node.name <- as.character(unrolled.dbn.adj.list[edge.idx, 'Target'])
    rolled.dbn.adj.matrix[src.node.name, tgt.node.name] <- 1
    
  }
  rm(edge.idx)
}

save(rolled.dbn.adj.matrix, file = paste(output.dirname, 'rolled.dbn.adj.matrix.RData', sep = '/'))
rm(unrolled.dbn.adj.list)
##------------------------------------------------------------
## End: Roll up the unrolled DBN to generate a rolled DBN adjacency matrix
##------------------------------------------------------------

## Create an '.sif' file equivalent to the rolled DBN adjacency matrix
## that is readable in Cytoscape.
adjmxToSif(rolled.dbn.adj.matrix, output.dirname)

## If the true rolled network is known a prior, then evaluate the performance
## metrics of the predicted rolled network.
if (true.net.filename != '') {
  
  ## Loads R obj 'true.net.adj.matrix'
  true.net.adj.matrix <- NULL
  load(true.net.filename)
  
  ## Begin: Create the format for result
  Result <- matrix(0, nrow = 1, ncol = 11)
  colnames(Result) <- list('TP', 'TN', 'FP', 'FN', 'TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F1')
  # ## End: Create the format for result
  
  if (is.matrix(true.net.adj.matrix)) {
    ## True net is time-invariant. Therefore, 
    ## 'true.net.adj.matrix' is a single matrix.
    
    predicted.net.adj.matrix <- rolled.dbn.adj.matrix
    
    ResultVsTrue <- calcPerfDiNet(predicted.net.adj.matrix, true.net.adj.matrix, Result, num.nodes)
    rm(predicted.net.adj.matrix)
    writeLines('Result ARTIVA vs True = \n')
    print(ResultVsTrue)
    rm(ResultVsTrue)
    
  } else if (is.list(true.net.adj.matrix)) {
    ## True nets are time-varying. Therefore, 
    ## 'true.net.adj.matrix' is a list of matrices.
    
    for (net.idx in 1:length(unrolled.dbn.adj.matrix.list)) {
      
      predicted.net.adj.matrix <- unrolled.dbn.adj.matrix.list[[net.idx]]
      
      ResultVsTrue <- calcPerfDiNet(predicted.net.adj.matrix, true.net.adj.matrix[[net.idx]], Result, num.nodes)
      rm(predicted.net.adj.matrix)
      Result <- rbind(Result, matrix(ResultVsTrue[1, ], nrow = 1, ncol = ncol(Result)))
      
      # rm(ResultVsTrue)
    }
    rm(net.idx)
    
    ## Print mean performance averaged over all time-varying networks
    ResultVsTrue <- colMeans(Result)
    ResultVsTrue <- matrix(colMeans(Result), nrow = 1, ncol = ncol(Result))
    colnames(ResultVsTrue) <- colnames(Result)
    writeLines('Result ARTIVA vs True = \n')
    print(ResultVsTrue)
    rm(ResultVsTrue)
  }
  
  save(Result, file = paste(output.dirname, 'Result.RData', sep = '/'))
  rm(Result)
}

rm(rolled.dbn.adj.matrix)

elapsed.time <- (proc.time() - start.time) # Stop the timer
writeLines('elapsed.time = \n')
print(elapsed.time)

sink()
close(output.file.conn)

##------------------------------------------------------------
## Begin: Save R session info in a File
##------------------------------------------------------------
sink(paste(output.dirname, 'sessionInfo.txt', sep = '/'))
sessionInfo()
sink()
##------------------------------------------------------------
## End: Save R session info in a File
##------------------------------------------------------------

##------------------------------------------------------------
## End: Main program
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: References
##------------------------------------------------------------
## 1. Lebre, Sophie, et al. "Statistical inference of the time-varying structure of gene-regulation networks." BMC systems biology 4.1 (2010): 130.
## 2. The 'ARTIVA' package in R: https://cran.r-project.org/package=ARTIVA
## 3. The 'ARTIVA' package manual version 1.2.3: https://cran.r-project.org/web/packages/ARTIVA/ARTIVA.pdf
##------------------------------------------------------------
## End: References
##------------------------------------------------------------
