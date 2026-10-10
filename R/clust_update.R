clust_update <- function(m,ClustIds,Ly,Lt,K,ref_hatyc,optns,robust,LOO,
                         mfpca_scatter = c("sign", "moment")){
  J <- length(Ly)
  # m: index of curve (m-th curve )
  # ClustIds m,ClustIds,Ly,Lt,K,ref_hatyc,optns,robust
  # packages and functions for parallel ----
  # suppressMessages(require(MASS))
  # suppressMessages(require(fdapace))
  # suppressMessages(require(KFPCA))
  # start of cluster updating----
  ith_yy <- ith_tt <- list()
  for(j in 1:length(Ly)){
    ith_yy[[j]] <- Ly[[j]][[m]]
    ith_tt[[j]] <- Lt[[j]][[m]]
  }
  hatyc <- multi_eigen <- multi_score <- multi_eigen_val <- list()
  for(u in 1:K){	## index for cluster group
    uu <- ClustIds[[u]]  ## index of curves in cluster u
    if(LOO){
      if( sum(m == uu)== 1  ){ ## If m-th curve is in the group u

        yy <- list()
        tt <- list()
        tart <- list()
        for(j in 1:length(Ly)){
          yy[[j]] <- lapply( uu,  function(x) Ly[[j]][[x]])
          yy[[j]][[ which(m==uu)  ]] <- NULL
          tt[[j]] <- lapply( uu,  function(x) Lt[[j]][[x]])
          tart[[j]] <- tt[[j]][[ which(m==uu)  ]]
          tt[[j]][[ which(m==uu)  ]] <- NULL
        }
        if(robust){
          temp_r <- lapply(1:J, function(x)
            KFPCA(Lt = tt[[x]],  Ly = yy[[x]], interval = optns$interval,
                          nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                          nRegGrid = optns$nRegGrid, fdParobj = optns$basis) )
        }else{
          temp_r <-  lapply(1:length(Ly), function(x)
            fdapace::FPCA(yy[[x]], tt[[x]], optns =  optns ))
        }
        # todo here
        yhatc_pred <- predicted_yhat(temp_r, ith_yy,  ith_tt, robust, optns, mfpca_scatter)

      }

      if( sum(m == uu)== 0  ){## If m-th curve is not in the group u
        yhatc_pred <- predicted_yhat(ref_hatyc[[u]], ith_yy,  ith_tt, robust, optns, mfpca_scatter)
      }
    }else{
      yhatc_pred <- predicted_yhat(ref_hatyc[[u]], ith_yy,  ith_tt, robust, optns, mfpca_scatter)
    }

    hatyc[[u]] <- yhatc_pred$fittedX
    multi_eigen[[u]] <- yhatc_pred$multi_eigen
    multi_score[[u]] <- yhatc_pred$multi_score
  }
  # calculate predict error
  err <- list()
  seCost <- vector()
  for(u in 1:K){
    y <- list()
    uu <- ClustIds[[u]]
    for(j in 1:length(Ly)){
      if(sum(m == uu)==1){
        err[[j]] <- hatyc[[u]][[j]]-  (  ith_yy[[j]] )
      }
      # y[[j]] <-  (hatyc[[u]][[j]] - ith_yy[[j]] )^2 / (ith_yy[[j]])^2
      y[[j]] <-  (hatyc[[u]][[j]] - ith_yy[[j]] )^2
    }
    seCost[u] <- sum(unlist(y) ) #   trapzRcpp(X = fpcaObjY$workGrid, Y = y)
  }
  clust_i <- which.min(seCost)

  return(list(clust_i = clust_i, err = err,
              hatyc = hatyc, seCost = seCost))
}
