# Goal: Calculate performance metrics of the directed net 'inferredNet' w.r.t. the directed net 'targetNet'.
# n = Number of nodes in each of the nets.
calcPerfDiNet <-function(inferredNet, targetNet, Result, n) {
  TrPos = 0
  TrNeg = 0
  FlPos = 0
  FlNeg = 0
  
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      if(inferredNet[i,j] == 1 )
      {
        if(inferredNet[i,j] == targetNet[i,j])
        {
          #True Positive
          TrPos = TrPos + 1	 	
        }
        else
        {
          FlPos = FlPos +1
        }
      }	   
      if(inferredNet[i,j] == 0 )
      {
        if (inferredNet[i,j] == targetNet[i,j])
        {
          # True Negative
          TrNeg = TrNeg + 1
        }
        else
        {
          FlNeg = FlNeg +1
        }
        
      }
    }	
  }
  
  #------------------------------------------------------------
  # Begin: Calculate Performance Metrics
  #------------------------------------------------------------
  TPR <- TrPos/(TrPos + FlNeg)
  FPR <- FlPos/(FlPos + TrNeg)
  
  ## Calculate FDR
  FDR <- NULL
  if ((FlPos == 0) & (TrPos == 0)) {
    FDR <- 0
  } else {
    FDR <- FlPos/(FlPos + TrPos)
  }
  
  ## Calculate PPV
  PPV <- NULL
  if ((FlPos == 0) & (TrPos == 0)) {
    PPV <- 0
  } else {
    PPV <- TrPos/(TrPos + FlPos)
  }

  ACC <- (TrPos + TrNeg)/(TrPos + FlPos + TrNeg + FlNeg)
  
  ## Calculate F1-score
  F1 <- NULL
  if ((PPV == 0) & (TPR == 0)) {

    F1 <- 0
    ## Ref: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
    
  } else {
    F1 <- 2 * PPV * TPR / (PPV + TPR)
  }
  
  ## Calculate MCC.
  ## '((TrPos == 0) & (FlPos == 0))' => Null graph.
  ## '((TrNeg == 0) & (FlNeg == 0))' => Complete graph.
  MCC <- NULL
  if (((TrPos == 0) & (FlPos == 0)) | ((TrNeg == 0) & (FlNeg == 0))) {
    
    MCC <- 0
    ## Ref: https://lettier.github.io/posts/2016-08-05-matthews-correlation-coefficient.html
    
  } else {
    MCC <- ((TrPos * TrNeg) - (FlNeg * FlPos)) / sqrt((TrPos + FlPos) * (TrPos + FlNeg) * (TrNeg + FlPos) * (TrNeg+FlNeg))
  }
  
  
  ## Calculate AUC under ROC
  # table <- minet::validate(inferredNet, targetNet)
  # AUC <- minet::auc.roc(table)
  #------------------------------------------------------------
  # End: Calculate Performance Metrics
  #------------------------------------------------------------
  
  Result[1, 1] <- TrPos
  Result[1, 2] <- TrNeg
  Result[1, 3] <- FlPos
  Result[1, 4] <- FlNeg
  Result[1, 5] <- TPR
  Result[1, 6] <- FPR
  Result[1, 7] <- FDR
  Result[1, 8] <- PPV
  Result[1, 9] <- ACC
  Result[1, 10] <- MCC
  Result[1, 11] <- F1
  # Result[1,8] <- AUC
  
  return (Result)
}

################################################################################################
## Goal: Given a predicted directed gene reuglatory network adjacency matrix, evaluate its 
## accuracy w.r.t. known gene-gene interactions as done in [Table 1, 1].
## Reference(s):
## [1] Song, Le, Mladen Kolar, and Eric P. Xing. "KELLER: estimating time-varying interactions 
## between genes." Bioinformatics 25.12 (2009): i128-i136.
##
eval.wrt.known.gene.ias <- function(di.net.adj.matrix)
{
  di.net.adj.matrix['CycE', 'CycA']
  di.net.adj.matrix['CycA', 'CycE']
  
  di.net.adj.matrix['CycE', 'Rca1']
  di.net.adj.matrix['Rca1', 'CycE']
  
  di.net.adj.matrix['CycE', 'Dp']
  di.net.adj.matrix['Dp', 'CycE']
  
  # di.net.adj.matrix['Gint3', 'Got2']
  # di.net.adj.matrix['Got2', 'Gint3']
  di.net.adj.matrix['G.ialpha65A', 'G.oalpha47A']
  di.net.adj.matrix['G.oalpha47A', 'G.ialpha65A']
  
  di.net.adj.matrix['Hem', 'blow']
  di.net.adj.matrix['blow', 'Hem']
  
  di.net.adj.matrix['Ice', 'Ark']
  di.net.adj.matrix['Ark', 'Ice']
  
  di.net.adj.matrix['Jra', 'dnc']
  di.net.adj.matrix['dnc', 'Jra']
  
  di.net.adj.matrix['Nf1', 'dnc']
  di.net.adj.matrix['dnc', 'Nf1']
  
  di.net.adj.matrix['Pak', 'trio']
  di.net.adj.matrix['trio', 'Pak']
  
  di.net.adj.matrix['Sb', 'Rho1']
  di.net.adj.matrix['Rho1', 'Sb']
  
  di.net.adj.matrix['Snap', 'Syx1A']
  di.net.adj.matrix['Syx1A', 'Snap']
  
  di.net.adj.matrix['Src42A', 'ksr']
  di.net.adj.matrix['ksr', 'Src42A']
  
  di.net.adj.matrix['w', 'nej']
  di.net.adj.matrix['nej', 'w']
  
  di.net.adj.matrix['brk', 'tkv']
  di.net.adj.matrix['tkv', 'brk']
  
  di.net.adj.matrix['brm', 'N']
  di.net.adj.matrix['N', 'brm']
  
  di.net.adj.matrix['brm', 'shg']
  di.net.adj.matrix['shg', 'brm']
  
  di.net.adj.matrix['btl', 'stumps']
  di.net.adj.matrix['stumps', 'btl']
  
  di.net.adj.matrix['cact', 'dl']
  di.net.adj.matrix['dl', 'cact']
  
  di.net.adj.matrix['caps', 'chi']
  di.net.adj.matrix['chi', 'caps']
  
  di.net.adj.matrix['da', 'Dl']
  di.net.adj.matrix['Dl', 'da']
  
  di.net.adj.matrix['dally', 'sgl']
  di.net.adj.matrix['sgl', 'dally']
  
  di.net.adj.matrix['dl', 'Dif']
  di.net.adj.matrix['Dif', 'dl']
  
  di.net.adj.matrix['dom', 'E(z)']
  di.net.adj.matrix['E(z)', 'dom']
  
  di.net.adj.matrix['ea', 'Tl']
  di.net.adj.matrix['Tl', 'ea']
  
  di.net.adj.matrix['emc', 'bs']
  di.net.adj.matrix['bs', 'emc']
  
  di.net.adj.matrix['esc', 'E(z)']
  di.net.adj.matrix['E(z)', 'esc']
  
  di.net.adj.matrix['gl', 'peb']
  di.net.adj.matrix['peb', 'gl']
  
  di.net.adj.matrix['hep', 'p53']
  di.net.adj.matrix['p53', 'hep']
  
  di.net.adj.matrix['mam', 'wg']
  di.net.adj.matrix['wg', 'mam']
  
  di.net.adj.matrix['msn', 'Nrt']
  di.net.adj.matrix['Nrt', 'msn']
  
  di.net.adj.matrix['msn', 'dock']
  di.net.adj.matrix['dock', 'msn']
  
  di.net.adj.matrix['nej', 'th']
  di.net.adj.matrix['th', 'nej']
  
  di.net.adj.matrix['numb', 'Rca1']
  di.net.adj.matrix['Rca1', 'numb']
  
  di.net.adj.matrix['pbl', 'CycE']
  di.net.adj.matrix['CycE', 'pbl']
  
  di.net.adj.matrix['pbl', 'Src64B']
  di.net.adj.matrix['Src64B', 'pbl']
  
  di.net.adj.matrix['pbl', 'dl']
  di.net.adj.matrix['dl', 'pbl']
  
  di.net.adj.matrix['pbl', 'tum']
  di.net.adj.matrix['tum', 'pbl']
  
  di.net.adj.matrix['pnr', 'svr']
  di.net.adj.matrix['svr', 'pnr']
  
  di.net.adj.matrix['pros', 'Abl']
  di.net.adj.matrix['Abl', 'pros']
  
  di.net.adj.matrix['pros', 'pnt']
  di.net.adj.matrix['pnt', 'pros']
  
  di.net.adj.matrix['sdt', 'baz']
  di.net.adj.matrix['baz', 'sdt']
  
  di.net.adj.matrix['sno', 'Dl']
  di.net.adj.matrix['Dl', 'sno']
  
  di.net.adj.matrix['spen', 'ksr']
  di.net.adj.matrix['ksr', 'spen']
  
  di.net.adj.matrix['tsh', 'wg']
  di.net.adj.matrix['wg', 'tsh']
  
  di.net.adj.matrix['up', 'Mhc']
  di.net.adj.matrix['Mhc', 'up']
}