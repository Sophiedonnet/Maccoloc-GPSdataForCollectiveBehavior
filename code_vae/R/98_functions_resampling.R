library(graphclust)
library(gtools)
########################################
### get best order to plot 
get_order <- function(res_iidcolsbm){
  
  gmin <- which.min(rowSums(as.matrix(res_iidcolsbm$alpha)))
  Q <- res_iidcolsbm$Q
  
  alpha <- res_iidcolsbm$alpha
  
  commu<- order(res_iidcolsbm$pi[[1]],decreasing = TRUE)[1:2]
  g1 <- commu[1]
  g2 <- commu[2]
  u1 <- alpha[g1,]
  u1[g1] <- -10
  g11 <- which.max(u1)
  u2 <- alpha[g2,]
  u2[c(g2,g11)] <- -10
  g22 <- which.max(u2)
  
  if(mean(alpha[g22,])<= mean(alpha[g11,])){ 
    best_order <- c(g22,g2,g11,g1,(1:Q)[-c(gmin,g1,g11,g2,g22)], gmin)
  }else{
    best_order <- c(g11,g1,g22,g2,(1:Q)[-c(gmin,g1,g11,g2,g22)], gmin)
  }
  return(best_order)
}

########################################
### get best order to plot 

get_order_withref <- function(fit, fit_ref,order_ref){
  
  theta_ref <-  list(pi=fit_ref$pi[[1]][order_ref],gamma=fit_ref$alpha[order_ref,order_ref])
  theta <- list(pi=fit$pi[[1]],gamma=fit$alpha)

  Q = length(theta$pi)
  List_ord <- permutations(n = Q , r = Q)
  
  best_order <- 1:Q
  best_dist <- Inf
  for (p in 1:nrow(List_ord)){
    
    perm <- List_ord[p, ]
    theta_perm <- list(pi=theta$pi[perm],gamma=theta$gamma[perm,perm])
    d_perm <- graphonL2norm(theta_ref, theta_perm)
    if (d_perm<best_dist){
      best_dist <- d_perm
      best_order <- perm 
    }
  }
  return(best_order)
}
  

#==================================================
#  get Z_trim from estim under iidcolsbm
#=================================================
get_Ztrim <- function(res_iidcolsbm){
  
  Z_trim <- sapply(seq(20), function(x) {
    y <- rep(NA, 33)
    names_boats <- names(res_iidcolsbm$Z[[x]])
    if(length(names_boats)==0){
      names(res_iidcolsbm$Z[[x]]) <- rownames(res_iidcolsbm$A[[x]])
    }
    y[(seq(33) %in% names(res_iidcolsbm$Z[[x]]))] <- res_iidcolsbm$Z[[x]]
    return(y)
  })
  Z_trim
  return(Z_trim)
}

#==================================================
#  function to delabel switching Z_trim from estim under iidcolsbm
#=================================================

get_Ztrim_relabel <- function(res_iidcolsbm){
  
  Z_trim <- sapply(seq(20), function(x) {
    y <- rep(NA, 33)
    y[(seq(33) %in% names(res_iidcolsbm$Z[[x]]))] <- res_iidcolsbm$Z[[x]]
    return(y)
  })
  commu<- order(res_iidcolsbm$pi[[1]],decreasing = TRUE)[1:2]
  g1 <- commu[1]
  g2 <- commu[2]
  u1 <- res_iidcolsbm$alpha[g1,]
  u1[g1] <- -10
  g11 <- which.max(u1)
  u2 <- res_iidcolsbm$alpha[g2,]
  u2[c(g2,g11)] <- -10
  g22 <- which.max(u2)
  map_from <- c(g1, g11, g2, g22,(1:res_iidcolsbm$Q)[-c(g1, g11, g2, g22)])
  map_to   <- c(g2, g22, g1, g11,(1:res_iidcolsbm$Q)[-c(g1, g11, g2, g22)])
  
  Z_trim_nolabelsw <- Z_trim
  for(t in 1:20){
    freq_g1_ind1 = mean(Z_trim[1,]==g1,na.rm = TRUE)
    freq_g2_ind1 = mean(Z_trim[1,]==g2,na.rm = TRUE)
    change_g <- (freq_g1_ind1<freq_g2_ind1)*g1 + (freq_g1_ind1>=freq_g2_ind1)*g2
    if(!is.na(Z_trim[1,t])){
      if(Z_trim[1,t]==change_g){ #on affecte bateau 1 au block 1 donc on echange label 2 3 pour 1 4
        Z_trim_nolabelsw[,t] <- map_to[ match(Z_trim[,t], map_from) ]
      }}
    
  }
  Z_trim <- Z_trim_nolabelsw
  return(Z_trim)
}
#=======================================
# == Get Zmap from Z trim 
#==================================== 

relabel_Ztrim <- function(Ztrim,best_order){
  
  Q <- max(best_order)
  map_to  <- order(best_order)
  Ztrim_relabel <- Ztrim
  for(t in 1:20){
    Ztrim_relabel[,t] <- map_to[Ztrim[,t]]
  }
  return(Ztrim_relabel)
}



#=======================================
# == Get Zmap from Z trim 
#==================================== 

get_Ztrim_map <- function(Ztrim){
  Q <- max(c(Ztrim),na.rm = TRUE)
  Z_trim_map <- (outer(Ztrim, seq(Q), "==") |> 
                 apply(c(1,3), "mean", na.rm = TRUE)) |> 
  apply(1, "which.max")
  return(Z_trim_map)}





#=================================================
# =============== function to delabel switching and fusionnate into two big communities
#===================================================

get_Ztrim_bigcommu <- function(res_iidcolsbm){
  
  #-------------------------------------------- 
  Z_trim <- sapply(seq(20), function(x) {
    y <- rep(NA, 33)
    y[(seq(33) %in% names(res_iidcolsbm$Z[[x]]))] <- res_iidcolsbm$Z[[x]]
    return(y)
  })
  commu<- order(res_iidcolsbm$pi[[1]],decreasing = TRUE)[1:2]
  g1 <- commu[1]
  g2 <- commu[2]
  u1 <- res_iidcolsbm$alpha[g1,]
  u1[g1] <- -10
  g11 <- which.max(u1)
  u2 <- res_iidcolsbm$alpha[g2,]
  u2[c(g2,g11)] <- -10
  g22 <- which.max(u2)
  map_from <- c(g1, g11, g2, g22,(1:res_iidcolsbm$Q)[-c(g1, g11, g2, g22)])
  map_to   <- c(g2, g22, g1, g11,(1:res_iidcolsbm$Q)[-c(g1, g11, g2, g22)])
  
  Z_trim_new <- Z_trim
  for(t in 1:20){
    freq_g1_ind1 = mean(Z_trim[1,]==g1,na.rm = TRUE)
    freq_g2_ind1 = mean(Z_trim[1,]==g2,na.rm = TRUE)
    change_g <- (freq_g1_ind1<freq_g2_ind1)*g1 + (freq_g1_ind1>=freq_g2_ind1)*g2
    if(!is.na(Z_trim[1,t])){
      if(Z_trim[1,t]==change_g){ #on affecte bateau 1 au block 1 donc on echange label 2 3 pour 1 4
        Z_trim_new[,t] <- map_to[ match(Z_trim[,t], map_from) ]
      }}
    
  }
  
  map_from_2 <- c(g11, g22,(1:res_iidcolsbm$Q)[-c(g11,g22)])
  map_to_2   <- c(g1, g2, (1:res_iidcolsbm$Q)[-c(g11, g22)])
  Z_trim_new_new <- Z_trim_new
  for(t in 1:20){
    Z_trim_new_new[,t] <- map_to_2[ match(Z_trim_new[,t], map_from_2) ]
  }
  return(Z_trim_new_new)
}



