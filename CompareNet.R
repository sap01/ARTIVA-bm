## Returns 1 if 'di.net.adj.matrix' = 'cmi.net.adj.matrix' else returns 0
##
CompareNet <-function(di.net.adj.matrix,cmi.net.adj.matrix,num.node)
{
  for(i in 1:num.node)
  {
    for(j in 1:num.node)
    {
      if(di.net.adj.matrix[i,j] != cmi.net.adj.matrix[i,j])
      {
        return(0)
      }
    }
  }
  return (1)
}

################################################################################################3
## Goal: Given two di nets' adjacency matrices (must of of same dim, same rownames, same colnames
## where rows = source nodes, cols = tgt nodes, 0 and 1 represent absence and presence of 
## an edge, resp.), this function prints the common edges in an output text file.
##
print.common.di.edges <- function(di.net.adj.matrix1, di.net.adj.matrix2)
{
  setwd('/home/saptarshi/R/R-3.3.2/projects/repoprojldbn')
  
  if ((dim(di.net.adj.matrix1)[1] != dim(di.net.adj.matrix2)[1]) | 
      (dim(di.net.adj.matrix1)[2] != dim(di.net.adj.matrix2)[2]))
  {
    stop('Dimensions are not comparable!')
  }
  
  output.filename <- paste(getwd(), 'asset/common.di.edgelist.txt', sep = '/')
  output.file <- file(output.filename, 'w')
  rm(output.filename)
  
  cat('source', 'target', '\n', file = output.file, sep = '\t')
  
  edge.cnt <- 0
  for (row.idx in 1:nrow(di.net.adj.matrix1))
  {
    for (col.idx in 1:ncol(di.net.adj.matrix1))
    {
      if ((di.net.adj.matrix1[row.idx, col.idx] == 1) & (di.net.adj.matrix2[row.idx, col.idx] == 1))
      {
        cat(rownames(di.net.adj.matrix1)[row.idx], colnames(di.net.adj.matrix1)[col.idx], '\n', 
            file = output.file, sep = '\t')
        
        edge.cnt <- edge.cnt + 1
      }
    }
  }
  
  close(output.file)
  
  print(edge.cnt)
}