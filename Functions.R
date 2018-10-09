library(ks)
library(rootSolve)
library(parallel)
library(MBESS)
library(MASS)
library(extraDistr)
library(microbenchmark)
library(Surrogate)
library(Matrix)
library(reshape)
library(Cairo)
library(ggplot2)
library(pbapply)
library(gridExtra)

##### 1. CREATING CORRELATION MATRIX
### 1.1 Generation of random values for Theta
#  Denstiy function of Thera
f.fun = function(theta,k,r){
  # theta: Theta values (vector)
  # k: parameter k > 0. For a uniform distribution k=0
  # r: dimension of the correlation matrix
  sin(theta)^(2*k+r)
}
Vectorize(f.fun)
# Cumulative function of Theta
g.fun = function(theta,k,r){
  # theta: Theta values (vector)
  # k: parameter k > 0. For a uniform distribution k=0
  # r: dimension of the correlation matrix
  sin(theta)^(2*k+r)/integrate(f.fun,k=k,r=r,lower=0,upper=pi)$value
}
Vectorize(g.fun)
G.fun = function(theta,k,r){
  # theta: Theta values (vector)
  # k: parameter k > 0. For a uniform distribution k=0
  # r: dimension of the correlation matrix
  integrate(g.fun,k=k,r=r,lower=0,upper=theta)$value
}
Vectorize(G.fun)
# inverse Cumulative function of Theta
G.inv <- function(U,k,r){
  # U: Cumulative probability 0 <= U <= 1
  # k: parameter k > 0. For a uniform distribution k=0
  # r: dimension of the correlation matrix
  uniroot(function(x){G.fun(x,k=k,r=r)-U},interval=c(0,pi))$root
}
Vectorize(G.inv)
# random values for Theta
r.G.fun = function(n,k,r){
  # n: number of random theta vectors
  # k: parameter k > 0. For a uniform distribution k=0
  # r: dimension of the correlation matrix
  
  u = runif(n)
  theta =sapply(u,G.inv,k=k,r=r,simplify=TRUE)
  return(theta)
}

# random single correlation based on partial correlations
r.random = function(j,k,R){
  d = ncol(R)
  r1 = R[j,(j+1):(j+k-1)]
  r3 = R[(j+1):(j+k-1),j+k]
  R2 = solve(R[(j+1):(j+k-1),(j+1):(j+k-1)])
  r.c = rnsbeta(1,1+0.5*(d-1-k),1+0.5*(d-1-k),-1,1)
  r = tcrossprod(crossprod(r1,R2),r3) + r.c*sqrt( (1 - tcrossprod(crossprod(r1,R2),r1) )*(1 - tcrossprod(crossprod(r3,R2),r3) )  )
  return(r)
}

is.PD = function(X,Vech=FALSE){
  if(Vech==TRUE){
    X= invvech(vechX)
  }
  min.lambda = min(eigen(X)$values)
  return(min.lambda > 0)
}

### 1.2 Generation of a random correlation matrix (d x d) 
# based on the Hyperspherical parameterization of the Cholesky factor
GenCorrMatrix.HP <- function(d,k){
  # p: dimension of the matrix
  # k: parameter k>=0 (for uniform correlation matrix k = 0)
  ThetaMat = matrix(0,d,d)
  for(i in 1:(d-1)){
    ThetaMat[,i] = c(rep(0,i),r.G.fun(d-i,k=k,r=d-i))
  }
  U = matrix(0,d,d)
  for(i in 1:d){
    for(j in 1:d){
      if(j >= 2 & j <= i-1){U[i,j]=prod(sin(ThetaMat[i,1:(j-1)]))*cos(ThetaMat[i,j])}
      if(i==j){U[i,j] = prod(sin(ThetaMat[i,1:(j-1)])) }
    }
  }
  U[,1] = cos(ThetaMat[,1])
  R = U%*%t(U)
  R = as.matrix(forceSymmetric(R))
  R
  return(list(R,ThetaMat))
}
# based on partial correlations
GenCorrMatrix.PC = function(d){
  #Rf = R
  #Rf[is.na(Rf)] = 0
  R = diag(d)
  Rf = R
  j.ind = do.call('c',lapply(1:(d-1),function(x){x:1}))
  k.ind = do.call('c',lapply(1:(d-1),function(x){1:x}))
  for(i in 1:(length(j.ind))){
    j = j.ind[i]
    k = k.ind[i]
      if(k==1){
        R[j,j+k] = R[j+k,j] = rnsbeta(1,d/2,d/2,-1,1)
      }else{
        R[j,j+k] = R[j+k,j]  = r.random(j,k,R)
      }
    
  }
  if(min(eigen(R)$values) < 0.0001){
    alpha  = uniroot(Det.fun,c(0,1),R=R,Rfixed=Rf)$root
    R = alpha*R + (1-alpha)*Rf
  }
  return(R)
}



#### Generation of a random correlation matrix with constraints
##### 2.1 SCALING METHOD
# Generation of the correlation matrix using the theta values
ThetaMatFun = function(d,Thetaval){
  # p: dimension the correlation matrix
  # Thetaval: vector of the Theta matix (Vec Theta)
  ThetaMat = matrix(0,d,d)
  ThetaMat[lower.tri(ThetaMat)] =Thetaval
  U = matrix(0,d,d)
  for(i in 1:d){
    for(j in 1:d){
      if(j >= 2 & j <= i-1){U[i,j]=prod(sin(ThetaMat[i,1:(j-1)]))*cos(ThetaMat[i,j])}
      if(i==j){U[i,j] = prod(sin(ThetaMat[i,1:(j-1)])) }
      #if(i> j){U[i,j]=0}
    }
  }
  U[,1] = cos(ThetaMat[,1])
  R = U%*%t(U)
  return(R)
}
# Target function to minimize 
TargetFun = function(Thetaval,R,W){
  # Thetaval: vector of the Theta matix (Vec Theta)
  # p: dimension of the correlation matrix
  # R: pseudo d x d correlation matrix
  # W: weight d x d matrix
  d=dim(R)[1]
  W.vec = vech(W)
  R.adj = ThetaMatFun(d,Thetaval)
  R.adj.vec = vech(R.adj)
  R.vec = vech(R)
  S = sum(W.vec*(R.adj.vec-R.vec)^2)
  return(S)
}
# Generation of random correlation matrix with fixed values (using scaling method). Result in vech form
GenCor.Scaling <- function(k,Rfixed,Vech = FALSE){
  # k: parameter k >= 0 (for uniform distribution, k = 0)
  # Rfixed: fixed values (not equal to 0)
  # Vech: result in vech form or not
  d = dim(Rfixed)[1]
  Rinfo = GenCorrMatrix.HP(d,k)
  R = Rinfo[[1]]
  Theta = Rinfo[[2]][lower.tri(Rinfo[[2]])]
  R[!is.na(Rfixed)] = Rfixed[!is.na(Rfixed)]  
  Rf = Rfixed
  Rf[is.na(Rf)] = 0
  is.pd = sum(eigen(R)$values < 0)  == 0
  if(is.pd == FALSE){
    W = matrix(0.001,d,d)
    W[!is.na(Rfixed)] = (!is.na(Rfixed))[!is.na(Rfixed)]*1e+100
    diag(W) = rep(1,d)
    #Sol = optim(par=Theta,fn=TargetFun,R=R,W=W) #### change function to nlm()
    Sol = optim(par=Theta,fn=TargetFun,R=R,W=W) #### change function to nlm()
    Ralt = ThetaMatFun(d,Sol$par)
    is.pd2 = sum(eigen(Ralt)$values < 0.0001)  == 0
    if(is.pd2 == FALSE){
      alpha  = uniroot(Det.fun,c(0,1),R=Ralt,Rfixed=Rf)$root
      u = alpha
      Ralt = u*Ralt + (1-u)*Rf
    }
    
  }else{
    Ralt=R
    is.pd2 = sum(eigen(Ralt)$values < 0.0001)  == 0
    if(is.pd2 == FALSE){
      alpha  = uniroot(Det.fun,c(0,1),R=Ralt,Rfixed=Rf)$root
      u = alpha
      Ralt = u*Ralt + (1-u)*Rf
    }
  }
  if(Vech == TRUE){
    diffMax = max(abs(Ralt - Rfixed),na.rm = TRUE)
    return(c(vech(Ralt),is.pd,diffMax))
  }else{
    return(Ralt)
  }
  # Result in vech form + if the first pseudo matrix is pd + shrinkage parameter
}
##### 2.2 SHRINKAGE METHOD
# function to compute the min eigen value of the shrinkage correlation matrix. 
Det.fun = function(alpha,R,Rfixed){
  # alpha: shrinkage parameter
  # R: pseudo correlation matrix
  # Rfixed: correlation matrix
  if(anyNA(Rfixed)){
    Rfixed[is.na(Rfixed)] = 0
  }
  f = alpha*R + (1-alpha)*Rfixed
  min(eigen(f)$values) - 0.0001
}
# Generation of random correlation matrix wiht fixed values (using shrinkage method)
GenCor.shrinkage = function(k,Rfixed,Vech=FALSE){
  # k: parameter k >= 0 (for uniform distribution, k = 0)
  # Rfixed: fixed values (not equal to 0)
  # Vech: result in vech form or not
  d = dim(Rfixed)[1]
  R = GenCorrMatrix.HP(d,k)[[1]]
  R[!is.na(Rfixed)] = Rfixed[!is.na(Rfixed)]  
#  R = as.matrix(forceSymmetric(R))
  is.pd = sum(eigen(R)$values < 0.0001)  == 0
  if(is.pd == FALSE){
    alpha  = uniroot(Det.fun,c(0,1),R=R,Rfixed=Rfixed)$root
    u = alpha*runif(1,0,1)
    Rf = Rfixed
    Rf[is.na(Rf)] = 0
    Ralt = u*R + (1-u)*Rf
  }else{
    Ralt=R
    alpha = NA
  }

  if(Vech == TRUE){
    return(c(vech(Ralt),is.pd,alpha))
  }else{
    return(Ralt)
    }
  # Result in vech form + if the first pseudo matrix is pd + shrinkage parameter
}
##### 2.3. Rejection sampling algorithm
# function to compute a correlation matrix with restriction
GenCor.NA = function(k,Rfixed,Vech=FALSE){
  # k: parameter k >= 0 (for uniform distribution, k = 0)
  # Rfixed: fixed values (not equal to 0)
  # Vech: result in vech form or not
  d = dim(Rfixed)[1]
  is.pd = FALSE
  num.cor = d*(d+1)/2 - d
  count = 0
  while(is.pd == FALSE){
    R = diag(d)
    R[upper.tri(R)] =runif(num.cor,-1,1)
    R= as.matrix(forceSymmetric(R))
    R[!is.na(Rfixed)] = Rfixed[!is.na(Rfixed)]  
    is.pd = sum(eigen(R)$values < 0.0001)  == 0
    count = count + 1
  }

  if(Vech == TRUE){
    return(c(vech(R),count))
  }else{
    return(R)
  }
}

GenCorr.gradualRS = function (Rfixed, G = seq(from = -1, to = 1, by = 1e-05)) {
  num_elem <- dim(Rfixed)[1]^2
  Rfixed.orig <- Rfixed
  size <- row_num <- col_num_ind <- 3
  total_size <- dim(Rfixed)[1]
  here <- Rfixed
  count = 0
  while (size <= total_size) {
    here <- Rfixed[(1:size), (1:size)]
    here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), 
                                replace = TRUE)
    here[upper.tri(here)] = t(here)[upper.tri(here)]
    while (det(here) < 0) {
      here <- Rfixed[(1:size), (1:size)]
      here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), 
                                  replace = TRUE)
      here[upper.tri(here)] = t(here)[upper.tri(here)]
      count = count + 1
    }
    print(paste('submatrix',size,'x',size))
    flush.console()
    
    Rfixed[1:row_num, 1:col_num_ind] <- here
    row_num <- row_num + 1
    col_num_ind <- col_num_ind + 1
    size <- size + 1
  }
  Rfixed <- here
  return(Rfixed)
}

#### 2.4 Algorithm based on partial correlation
r.random = function(j,k,R,Range=c(-1,1)){
  d = ncol(R)
  r1 = R[j,(j+1):(j+k-1)]
  r3 = R[(j+1):(j+k-1),j+k]
  R2 = solve(R[(j+1):(j+k-1),(j+1):(j+k-1)])
  r.c = rnsbeta(1,1+0.5*(d-1-k),1+0.5*(d-1-k),Range[1],Range[2])
  r = tcrossprod(crossprod(r1,R2),r3) + r.c*sqrt( (1 - tcrossprod(crossprod(r1,R2),r1) )*(1 - tcrossprod(crossprod(r3,R2),r3) )  )
  return(r)
}
RCorr.fun = function(R,Range=c(-1,1)){
  Rf = R
  Rf[is.na(Rf)] = 0
  d = ncol(R)
  j.ind = do.call('c',lapply(1:(d-1),function(x){x:1}))
  k.ind = do.call('c',lapply(1:(d-1),function(x){1:x}))
  for(i in 1:(length(j.ind))){
    j = j.ind[i]
    k = k.ind[i]
    if(is.na(R[j,j+k])){
      if(k==1){
        R[j,j+k] = R[j+k,j] = rnsbeta(1,d/2,d/2,Range[1],Range[2])
      }else{
        R[j,j+k] = R[j+k,j]  = r.random(j,k,R,Range)
      }
    }
  }
  if(min(eigen(R)$values) < 0){
    alpha  = uniroot(Det.fun,c(0,1),R=R,Rfixed=Rf)$root
    R = alpha*R + (1-alpha)*Rf
  }
  return(R)
}

####
GenerateRandCor = function(d,M,rho,rpc=TRUE,grs=TRUE,sca=TRUE,shr=TRUE){
  Rfixed = matrix(NA,d,d)
  R1 = matrix(rho,d/2,d/2);diag(R1) = 1
  Rfixed[1:(d/2),1:(d/2)] = R1
  Rfixed[d/2 + 1:(d/2),d/2 + 1:(d/2)] = R1
  
  Comb = combn(rep(1:d,each=2),2)
  LABEL = paste('r',unique(apply(Comb,2,function(x){paste(x,collapse = '')})),sep='')
  FILE = paste('d',d,'rho',rho,sep='.')
  if(rpc==TRUE){
    RPC = mapply( function(x){
      start_time <- Sys.time()
      R = RCorr.fun(Rfixed)
      end_time <- Sys.time()
      time = end_time - start_time
      c(vech(R),time)
    },1:M)
    RPC = t(RPC)
    colnames(RPC) = c(LABEL,'time')
    write.csv(RPC,paste('simulations/RPC',FILE,'csv',sep = '.'),row.names = FALSE)
  }
  if(grs==TRUE){
    gRS = mapply( function(x){
      start_time <- Sys.time()
      R = GenCorr.gradualRS(Rfixed)
      end_time <- Sys.time()
      time = end_time - start_time
      c(vech(R),time)
    },1:M)
    gRS = t(gRS)
    colnames(gRS) = c(LABEL,'time')
    write.csv(RPC,paste('simulations/GRS',FILE,'csv',sep = '.'),row.names = FALSE)
  }
  if(sca == TRUE){
    SCA = mapply( function(x){
      start_time <- Sys.time()
      R = GenCor.Scaling(0,Rfixed,TRUE)
      end_time <- Sys.time()
      time = end_time - start_time
      c(R,time)
    },1:M)
    SCA = t(SCA)
    colnames(SCA) = c(LABEL,'iniPD','maxdiff','time')
    write.csv(SCA,paste('simulations/SCA',FILE,'csv',sep = '.'),row.names = FALSE)
  }
  if(shr == TRUE){
    SHR = mapply( function(x){
      start_time <- Sys.time()
      R = GenCor.shrinkage(0,Rfixed,TRUE)
      end_time <- Sys.time()
      time = end_time - start_time
      c(R,time)
    },1:M)
    SHR = t(SHR)
    colnames(SHR) = c(LABEL,'iniPD','alpha','time')
    write.csv(SHR,paste('simulations/SHR',FILE,'csv',sep = '.'),row.names = FALSE)
  }
}

##### 3.INDIVIDUAL CAUSAL EVALUATION (ICA) - SINGLE TRIAL - SENSITIVITY ANALYSIS. 
ICA.fun = function(R,Sigma){
  
  rho.delta = (sqrt(Sigma[1]*Sigma[3])*R[1,3] + sqrt(Sigma[2]*Sigma[4])*R[2,4] - sqrt(Sigma[2]*Sigma[3])*R[2,3] - sqrt(Sigma[1]*Sigma[4])*R[1,4])/
    sqrt( (Sigma[1]+Sigma[2] -2*sqrt(Sigma[1]*Sigma[2])*R[1,2])*(Sigma[3] + Sigma[4] - 2*sqrt(Sigma[3]*Sigma[4])*R[3,4]) )
  
  sigma.delta.T = Sigma[1] + Sigma[2] - 2*sqrt(Sigma[1]*Sigma[2])*R[1,2]
  delta = sigma.delta.T*(1-rho.delta^2)
  return(c(rho.delta,delta))
} 
### multivariate ICA
MutivarICA.fun = function(R,Sigma,N){
  d = nrow(R)
  sdMat = diag(sqrt(Sigma))
  rtn = sdMat %*% R %*% t(sdMat)
  var_diff <- function(cov_mat) {
    cov_val <- cov_mat[1, 1] + cov_mat[2, 2] - (2 * 
                                                  cov_mat[1, 2])
    fit <- c(cov_val)
    fit
  }
  cov_2_diffs <- function(cov_mat) {
    cov_val <- (cov_mat[2, 2] - cov_mat[1, 2]) - 
      (cov_mat[2, 1] - cov_mat[1, 1])
    fit <- c(cov_val)
    fit
  }
  A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)
  B <- NULL
  Aantal <- (dim(R)[1] - 2)/2
  rtn_part <- rtn[c(3:dim(rtn)[1]), c(1, 2)]
  for (z in 1:Aantal) {
    cov_mat_here <- rtn_part[c((z * 2) - 1, z * 2), 
                             c(1:2)]
    B <- rbind(B, cov_2_diffs(cov_mat_here))
  }
  Dim <- dim(R)[1]
  Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]
  C <- matrix(NA, nrow = Aantal, ncol = Aantal)
  for (l in 1:Aantal) {
    for (k in 1:Aantal) {
      Sub_mat_var_hier <- Sub_mat_var[c((k * 2) - 
                                          1, k * 2), c((l * 2) - 1, l * 2)]
      C[k, l] <- cov_2_diffs(cov_mat = Sub_mat_var_hier)
    }
  }
  Delta <- cbind(rbind(A, B), rbind(t(B), C))
  ICA <- (t(B) %*% solve(C) %*% B)/A
  Adj.ICA <- 1 - (1 - ICA) * ((N - 1)/(N - Aantal - 
                                         1))
  return(c(ICA, Adj.ICA))
}

##### 4. COMPUTATION TIME EVALUATION
####### computation time
Computation.time <- function(d,rho,Times,grs=T,pc=T,sca=T,shr=T){
  ###d = 22; rho=0.5
  Rfixed = matrix(NA,d,d)
  R1 = matrix(rho,d/2,d/2);diag(R1) = 1
  Rfixed[1:(d/2),1:(d/2)] = R1
  Rfixed[d/2 + 1:(d/2),d/2 + 1:(d/2)] = R1
  
  if(pc==TRUE){PC = microbenchmark('PC' = {RCorr.fun(Rfixed)},times=Times)}else{PC = NULL}
  if(grs==TRUE){GRS = microbenchmark('gRS' = {GenCorr.gradualRS(Rfixed)},times=Times)}else{GRS = NULL}
  if(sca==TRUE){SCA = microbenchmark('SCA' = {GenCor.Scaling(0,Rfixed,FALSE)},times=Times)}else{SCA = NULL}
  if(shr==TRUE){SHR = microbenchmark('SHR' = {GenCor.shrinkage(0,Rfixed,FALSE)},times=Times)}else{SHR = NULL}
  
  return(list(PC,GRS,SCA,SHR))
}

###############
ICAvech.fun = function(Rvech,Sigma,N){
  R = invvech(Rvech)
  d=ncol(R)
  IND =  vec(matrix(1:d,ncol=2),byrow = TRUE) 
  R = R[IND,IND]

  if(d == 4){
    ICA = ICA.fun(R,Sigma)[1]
  }
  if(d>4){
    ICA = MutivarICA.fun(R,Sigma,N)[1]
  }
  return(ICA)
}

ICA.table = function(d,rho){
# d = 12; rho=0.2
  file = paste('d',d,'rho',rho,sep = '.')
  Sigma = rep(100,d)
  R.RPC = tryCatch({read.csv(paste('Simulations/RPC',file,'csv',sep='.'))},error =function(e){NULL})
  R.SCA = tryCatch({read.csv(paste('Simulations/SCA',file,'csv',sep='.'))},error =function(e){NULL})
  R.SHR = tryCatch({read.csv(paste('Simulations/SHR',file,'csv',sep='.'))},error =function(e){NULL})
  R.GRS = tryCatch({read.csv(paste('Simulations/GRS',file,'csv',sep='.'))},error =function(e){NULL},warning =function(w){NULL})
  
  vech.IND = d*(d+1)/2
  ICA.RPC = tryCatch({apply(as.matrix(R.RPC[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SCA = tryCatch({apply(as.matrix(R.SCA[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SHR = tryCatch({apply(as.matrix(R.SHR[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.GRS = tryCatch({apply(as.matrix(R.GRS[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  
  Results = rbind(c(mean(ICA.GRS),sd(ICA.GRS),quantile(ICA.GRS,c(0.05,0.95))),
  c(mean(ICA.RPC),sd(ICA.RPC),quantile(ICA.RPC,c(0.05,0.95))),
  c(mean(ICA.SHR),sd(ICA.SHR),quantile(ICA.SHR,c(0.05,0.95))),
  c(mean(ICA.SCA),sd(ICA.SCA),quantile(ICA.SCA,c(0.05,0.95))))
  colnames(Results) = c('mean','sd','p0.05','p0.95')
  rownames(Results) = c('gRS','PC','SHR','SCA')
  Time = c(tryCatch({median(R.GRS$time)},warning=function(w){NA}),median(R.RPC$time),median(R.SHR$time),median(R.SCA$time))
  names(Time) = c('gRS','PC','SHR','SCA')
  R.SHR$alpha[is.na(R.SHR$alpha)] = 1
  Other = c(mean(R.SHR$iniPD),mean(R.SHR$alpha),sd(R.SHR$alpha),mean(R.SCA$maxdiff),sd(R.SCA$maxdiff),max(R.SCA$maxdiff))
  names(Other) = c('iniPD','mean.alpha','sd.alpha','mean.diff','sd.diff','max.diff')
  Res = list(ICA = Results,time = Time,others = Other)
  return(Res)
}

addres.SHR.SCA = function(d,rho){
  file = paste('d',d,'rho',rho,sep = '.')  
  R.SCA = tryCatch({read.csv(paste('Simulations 1/SCA',file,'csv',sep='.'))},error =function(e){NULL})
  R.SHR = tryCatch({read.csv(paste('Simulations 1/SHR',file,'csv',sep='.'))},error =function(e){NULL})
  R.SHR$alpha[is.na(R.SHR$alpha)] =1  
  
  c(mean(R.SHR$iniPD),mean(R.SCA$maxdiff),
    mean(R.SHR$alpha))
}

ICA.plots = function(d,rho,Ylim=NULL,MAIN=NULL){
  # d = 8; rho=0.8
  Rfixed = matrix(NA,d,d)
  R1 = matrix(rho,d/2,d/2);diag(R1) = 1
  Rfixed[1:(d/2),1:(d/2)] = R1
  Rfixed[d/2 + 1:(d/2),d/2 + 1:(d/2)] = R1
  Rfixed[is.na(Rfixed)] = 0
  file = paste('d',d,'rho',rho,sep = '.')
  Sigma = rep(100,d)
  R.RPC = tryCatch({read.csv(paste('Simulations 1/RPC',file,'csv',sep='.'))},error =function(e){NULL})
  R.SCA = tryCatch({read.csv(paste('Simulations 1/SCA',file,'csv',sep='.'))},error =function(e){NULL})
  R.SHR = tryCatch({read.csv(paste('Simulations 1/SHR',file,'csv',sep='.'))},error =function(e){NULL})
  R.GRS = tryCatch({read.csv(paste('Simulations 1/GRS',file,'csv',sep='.'))},error =function(e){NULL},warning =function(w){NULL})
  
  vech.IND = d*(d+1)/2
  ICA.RPC = tryCatch({apply(as.matrix(R.RPC[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SCA = tryCatch({apply(as.matrix(R.SCA[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SHR = tryCatch({apply(as.matrix(R.SHR[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.GRS = tryCatch({apply(as.matrix(R.GRS[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  
  
  ICAs = melt(as.data.frame(cbind(ICA.RPC,ICA.SHR,ICA.SCA,ICA.GRS)))
  ICA0 = ICAvech.fun(vech(Rfixed),Sigma,N=20)
  PLOT = ggplot(aes(x=value, group=variable,linetype=variable,col=variable), data=ICAs)  + geom_density(show.legend = F,size=1) + 
        theme_bw() + labs(x=expression(R[H]^2),size=5) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=20),
              axis.text = element_text(size = 16))   +
              ylim(Ylim) + xlim(c(0,1)) + ggtitle(MAIN) + theme(plot.title = element_text(hjust = 0.5,size=40)) +
              scale_linetype_manual(values=c('dashed','solid',"dashed", "solid")) +  scale_colour_manual( values = c("black","lightgrey",'lightgrey','black'))
  return(PLOT)
}


Corr.plots = function(d,rho,r=10,Ylim=NULL,MAINTITLE=c(' ',' ',' ',' ')){
  # d = 8; rho=0.2
  Rfixed = matrix(NA,d,d)
  R1 = matrix(rho,d/2,d/2);diag(R1) = 1
  Rfixed[1:(d/2),1:(d/2)] = R1
  Rfixed[d/2 + 1:(d/2),d/2 + 1:(d/2)] = R1
  IND = vech(Rfixed)
  Rfixed[is.na(Rfixed)] = 0
  file = paste('d',d,'rho',rho,sep = '.')
  Sigma = rep(100,d)
  R.RPC = tryCatch({read.csv(paste('Simulations 1/RPC',file,'csv',sep='.'))},error =function(e){NULL})
  R.SCA = tryCatch({read.csv(paste('Simulations 1/SCA',file,'csv',sep='.'))},error =function(e){NULL})
  R.SHR = tryCatch({read.csv(paste('Simulations 1/SHR',file,'csv',sep='.'))},error =function(e){NULL})
  R.GRS = tryCatch({read.csv(paste('Simulations 1/GRS',file,'csv',sep='.'))},error =function(e){NULL},warning =function(w){NULL})
  
  Rind = (1:length(IND))[is.na(IND)]
  ncor = sample(Rind,r)
  
  vech.IND = d*(d+1)/2
  ICA.RPC = tryCatch({apply(as.matrix(R.RPC[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SCA = tryCatch({apply(as.matrix(R.SCA[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SHR = tryCatch({apply(as.matrix(R.SHR[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.GRS = tryCatch({apply(as.matrix(R.GRS[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  
  RPC = melt(R.RPC[,ncor])
  RPC.plot = ggplot(aes(x=value, col=variable), data=RPC) + geom_density(show.legend = F,size=1) + theme_bw() + labs(x='correlation') +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=20),
          axis.text = element_text(size = 16))  +
    xlim(c(-1,1)) + ggtitle(MAINTITLE[2]) + theme(plot.title = element_text(hjust = 0.5,size=40)) + ylim(Ylim)
  
  SCA = melt(R.SCA[,ncor])
  SCA.plot = ggplot(aes(x=value, col=variable), data=SCA) + geom_density(show.legend = F,size=1) + theme_bw() + labs(x='correlation') +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=20),
          axis.text = element_text(size = 16))  +
    xlim(c(-1,1)) + ggtitle(MAINTITLE[4]) + theme(plot.title = element_text(hjust = 0.5,size=40)) + ylim(Ylim)
  
  SHR = melt(R.SHR[,ncor])
  SHR.plot = ggplot(aes(x=value, col=variable), data=SHR) + geom_density(show.legend = F,size=1) + theme_bw() + labs(x='correlation') +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=20),
          axis.text = element_text(size = 16))  +
    xlim(c(-1,1)) + ggtitle(MAINTITLE[3]) + theme(plot.title = element_text(hjust = 0.5,size=40)) + ylim(Ylim)
  
  GRS = melt(R.GRS[,ncor])
  GRS.plot = ggplot(aes(x=value, col=variable), data=GRS) + geom_density(show.legend = F,size=1) + theme_bw() + labs(x='correlation') +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=20),
          axis.text = element_text(size = 16))  +
    xlim(c(-1,1)) + ggtitle(MAINTITLE[1]) + theme(plot.title = element_text(hjust = 0.5,size=40)) + ylim(Ylim)

    return(list(GRS.plot,RPC.plot,SHR.plot,SCA.plot))
}

## ICA computation using several cores
ICA.ContCont.MultS.PC.MC =function (M = 1000, N, Sigma, Seed = 123, Show.Progress = FALSE){
  is.PD = function(X, tol = 1e-08) {
    X[is.na(X)] = 0
    min.lambda = min(eigen(X, only.values = T, symmetric = T)$values)
    return(min.lambda > tol)
  }
  r.random = function(j, k, R) {
    d = ncol(R)
    r1 = R[j, (j + 1):(j + k - 1)]
    r3 = R[(j + 1):(j + k - 1), j + k]
    R2.inv = R[(j + 1):(j + k - 1), (j + 1):(j + k - 1)]
    R2 = solve(R2.inv)
    r.c = extraDistr::rnsbeta(1, 1 + 0.5 * (d - 1 - k), 1 + 
                                0.5 * (d - 1 - k), -1, 1)
    D2 = (1 - tcrossprod(crossprod(r1, R2), r1)) * (1 - tcrossprod(crossprod(r3, 
                                                                             R2), r3))
    if (D2 < 0 & D2 > -1e-08) {
      D2 = 0
    }
    r = tcrossprod(crossprod(r1, R2), r3) + r.c * sqrt(D2)
    return(r)
  }
  Correlation.matrix.PC = function(R, Range = c(-1, 1)) {
    Rf = R
    Rf[is.na(Rf)] = 0
    d = ncol(R)
    j.ind = do.call("c", lapply(1:(d - 1), function(x) {
      x:1
    }))
    k.ind = do.call("c", lapply(1:(d - 1), function(x) {
      1:x
    }))
    for (i in 1:(length(j.ind))) {
      j = j.ind[i]
      k = k.ind[i]
      if (is.na(R[j, j + k])) {
        if (k == 1) {
          R[j, j + k] = R[j + k, j] = extraDistr::rnsbeta(1, 
                                                          d/2, d/2, Range[1], Range[2])
        }
        else {
          R[j, j + k] = R[j + k, j] = r.random(j, k, 
                                               R)
        }
        if (is.nan(R[j, j + k])) {
          stop("error")
        }
      }
    }
    return(R)
  }
  MutivarICA.fun = function(R, Sigma, N) {
    d = nrow(R)
    sdMat = diag(sqrt(Sigma))
    rtn = sdMat %*% R %*% t(sdMat)
    var_diff <- function(cov_mat) {
      cov_val <- cov_mat[1, 1] + cov_mat[2, 2] - (2 * cov_mat[1, 
                                                              2])
      fit <- c(cov_val)
      fit
    }
    cov_2_diffs <- function(cov_mat) {
      cov_val <- (cov_mat[2, 2] - cov_mat[1, 2]) - (cov_mat[2, 
                                                            1] - cov_mat[1, 1])
      fit <- c(cov_val)
      fit
    }
    A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)
    B <- NULL
    Aantal <- (dim(R)[1] - 2)/2
    rtn_part <- rtn[c(3:dim(rtn)[1]), c(1, 2)]
    for (z in 1:Aantal) {
      cov_mat_here <- rtn_part[c((z * 2) - 1, z * 2), c(1:2)]
      B <- rbind(B, cov_2_diffs(cov_mat_here))
    }
    Dim <- dim(R)[1]
    Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]
    C <- matrix(NA, nrow = Aantal, ncol = Aantal)
    for (l in 1:Aantal) {
      for (k in 1:Aantal) {
        Sub_mat_var_hier <- Sub_mat_var[c((k * 2) - 1, 
                                          k * 2), c((l * 2) - 1, l * 2)]
        C[k, l] <- cov_2_diffs(cov_mat = Sub_mat_var_hier)
      }
    }
    Delta <- cbind(rbind(A, B), rbind(t(B), C))
    ICA <- (t(B) %*% solve(C) %*% B)/A
    Adj.ICA <- 1 - (1 - ICA) * ((N - 1)/(N - Aantal - 1))
    return(c(ICA, Adj.ICA))
  }
  #set.seed(Seed)
  d = nrow(Sigma)[1]
  Vars = diag(Sigma)
  IND = ks::vec(matrix(1:d, ncol = 2, byrow = T), byrow = F)
  Sigma = Sigma[IND, IND]
  R = cov2cor(Sigma)
  if (!is.PD(R)) {
    alpha = uniroot(function(alpha, R, Rfixed, tol = 1e-04) {
      if (anyNA(R)) {
        R[is.na(R)] = 0
      }
      f = alpha * R + (1 - alpha) * Rfixed
      min(eigen(f)$values) - tol
    }, c(0, 1), R = R, Rfixed = diag(d), tol = 1e-08)$root
    R = R * alpha + (1 - alpha) * diag(d)
    warning(paste("The initial correlation matrix is not PD. TThe matrix was shrunk by a factor alpha=", 
                  alpha, " for correction", sep = ""))
  }
  cl <- makeCluster(getOption("cl.cores", 3))
  clusterExport(cl=cl, varlist=c('R','d',"Vars",'N', 'is.PD', 'r.random','Correlation.matrix.PC','MutivarICA.fun'), envir=environment())
  clusterSetRNGStream(cl,Seed)
  
  Results = pblapply(X=1:M,function(X){
    R.random = Correlation.matrix.PC(R)
    IND = ks::vec(matrix(1:d, ncol = 2), byrow = TRUE)
    R.random = R.random[IND, IND]
    ICA = MutivarICA.fun(R.random, Vars, N)
    return(c(ICA, R.random[lower.tri(R.random)]))
  },cl=cl)
  stopCluster(cl) 
  Results = do.call('rbind',Results)
  R2_H = Results[,1 ]
  Corr.R2_H = Results[,2 ]
  Outcome = list(R2_H = R2_H, Corr.R2_H = Corr.R2_H, Lower.Dig.Corrs.All = t(Results[-c(1:2), 
                                                                                     ]))
  class(Outcome) <- "ICA.ContCont.MultS"
  return(Outcome)
}

#### Variable selection function
Var.selection.ICA = function(True.Sur.Selected,Sur.candidates,treat, M=500){
  N = length(treat)
  p = ncol(True.Sur.Selected)
  if(is.null(p)){p=1}
  step = p
  q = ncol(Sur.candidates)
  Complete.Data = cbind(True.Sur.Selected,Sur.candidates)
  Sigma.control = cov(Complete.Data[treat==-1,])
  Sigma.treat = cov(Complete.Data[treat==1,])  
  Sigma.NA = matrix(NA,p+q,p+q)
  Sigma = rbind(cbind(Sigma.control,Sigma.NA),cbind(Sigma.NA,Sigma.treat))
  m = p+q
  d = (p+1)*2
  IND =  ks::vec(matrix(1:d,ncol=2),byrow = TRUE) 
  
  ICA.add.one.sur = matrix(NA,q,5)
  Time = NULL
  for(i in 1:q){
    print(paste('step',step,'surrogate',i,'out of',q,sep = ' '))
    sur =  i + p
    Sigma.eval = Sigma[c(1:p,sur,(1:p)+m,m+sur),c(1:p,sur,(1:p)+m,m+sur)]
    Sigma.eval = Sigma.eval[IND,IND]
    tic()
    Fit = ICA.ContCont.MultS.PC.MC(M,N,Sigma.eval,Seed=123+i)
    exectime <- toc()
    Time[i] <- exectime$toc - exectime$tic
    
    ICA.add.one.sur[i,] = c(median(Fit$Corr.R2_H),mean(Fit$Corr.R2_H),sd(Fit$Corr.R2_H),range(Fit$Corr.R2_H,Fit$Corr.R2_H))
    
  }
  
  ICA.add.one.sur = cbind(ICA.add.one.sur,Time)
  rownames(ICA.add.one.sur) = colnames(Sur.candidates)
  colnames(ICA.add.one.sur) = c('median','mean','sd','min','max','time')
  label.max = (1:q)[ICA.add.one.sur[,1]==max(ICA.add.one.sur[,1])]
  sur.name = rownames(ICA.add.one.sur)[label.max]
  max.ICA = c(ICA.add.one.sur[label.max,])
  write.csv(ICA.add.one.sur,paste('Surrogate Selection/Step',step,'txt',sep='.'))
  results = list(ICA.add.one.sur=ICA.add.one.sur, max.ICA = max.ICA, sur.name=sur.name)
  return(results)
}
