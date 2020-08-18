

#' SIMILARITY LEARNING FOR RAFSIL
#'
#' @param F matrix Feature matrix cells x genes
#' @param nrep integer The number of replicates to average over for RAFSIL1
#' @param oob logical Should out of bag proximity be used?
#' @param method string "RAFSIL1" or "RAFSIL2"
#' @return matrix of dissimilarities between cells
#' @importFrom randomForest randomForest
#' @importFrom fpc pamk
rafsilRF <- function(F,nrep=50,oob=TRUE,method="RAFSIL1"){
#==========================================================

  dis = NULL

  #- RAFSIL1
  #---------
  if(method=="RAFSIL1"){

    sims = replicate(nrep,randomForest(x=F,y=NULL,proximity=TRUE,oob.prox=oob)$proximity)
    sim  = apply(sims,c(1,2),mean)
    dis  = (1-sim)

  #- RAFSIL 2
  #----------
  } else if( method == "RAFSIL2"){

    featFu <- function(i){
      cl   = pamk(F[,i])
      labs = as.factor(cl$pamobject$clustering)
      return(randomForest(x=F[,-i], y=labs, proximity=TRUE,oob.prox=oob)$proximity)
    }

    sims = vapply(1:ncol(F),featFu,array(0,dim=c(nrow(F),nrow(F))))
    sim  = apply(sims,c(1,2),mean)
    dis  = (1-sim)

  } else {
    stop("method not implemented")
  }

  return(dis)

}
