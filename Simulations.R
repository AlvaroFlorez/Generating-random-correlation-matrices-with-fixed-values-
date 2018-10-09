source('Functions.R')
#### SIMULATION FOR SECTION 5.2
## Computational time (Results of Table 1)
Timep1r0.2 = Computation.time(4,0.2,100,grs=T,pc=T,sca=T,shr=T)
Timep1r0.5 = Computation.time(4,0.5,100,grs=T,pc=T,sca=T,shr=T)
Timep1r0.8 = Computation.time(4,0.8,100,grs=T,pc=T,sca=T,shr=T)
Timep3r0.2 = Computation.time(8,0.2,100,grs=T,pc=T,sca=T,shr=T)
Timep3r0.5 = Computation.time(8,0.5,100,grs=T,pc=T,sca=T,shr=T)
Timep3r0.8 = Computation.time(8,0.8,100,grs=T,pc=T,sca=T,shr=T)
Timep5r0.2 = Computation.time(12,0.2,100,grs=T,pc=T,sca=T,shr=T)
Timep5r0.5 = Computation.time(12,0.5,100,grs=T,pc=T,sca=T,shr=T)
Timep5r0.8 = Computation.time(12,0.8,100,grs=T,pc=T,sca=T,shr=T)
Timep10r0.2 = Computation.time(22,0.2,100,grs=F,pc=T,sca=T,shr=T)
Timep10r0.5 = Computation.time(22,0.5,100,grs=F,pc=T,sca=T,shr=T)
Timep10r0.8 = Computation.time(22,0.8,100,grs=F,pc=T,sca=T,shr=T)

Table.1 = rbind(
c(mapply(function(x){summary(x,unit='s')$median},Timep1r0.2),mapply(function(x){summary(x,unit='s')$median},Timep1r0.5),
  mapply(function(x){summary(x,unit='s')$median},Timep1r0.8)),
c(mapply(function(x){summary(x,unit='s')$median},Timep3r0.2),mapply(function(x){summary(x,unit='s')$median},Timep3r0.5),
  mapply(function(x){summary(x,unit='s')$median},Timep3r0.8)),
c(mapply(function(x){summary(x,unit='s')$median},Timep5r0.2),mapply(function(x){summary(x,unit='s')$median},Timep5r0.5),
  mapply(function(x){summary(x,unit='s')$median},Timep5r0.8)),
c(mapply(function(x){summary(x,unit='s')$median},Timep10r0.2),mapply(function(x){summary(x,unit='s')$median},Timep10r0.5),
  mapply(function(x){summary(x,unit='s')$median},Timep10r0.8)))

### Generating and storing random correlation matrices
rho= c(0.2,0.5,0.8)
##  one surrogate (4 x 4 matrices)
mapply(GenerateRandCor,rho=rho, MoreArgs = list(d=4,M=1000,rpc=T,grs=T,sca=T,shr=T))
##  three surrogate (8 x 8 matrices)
mapply(GenerateRandCor,rho=rho, MoreArgs = list(d=8,M=1000,rpc=T,grs=T,sca=T,shr=T))
##  five surrogate (12 x 12 matrices)
grs = c(T,T,F) # excluding GRS algorithm for the d=12 rho=0.8 case (too time-consuming)
mapply(GenerateRandCor,rho=rho,grs=grs, MoreArgs = list(d=12,M=1000,rpc=T,sca=T,shr=T))
##  10 surrogate (4 x 4 matrices)
mapply(GenerateRandCor,rho=rho, MoreArgs = list(d=22,M=1000,rpc=T,grs=F,sca=T,shr=T)) # # excluding GRS algorithm for all cases (time-consuming)


### ICA densities (Figure 1)
#### FIGURES: density of (multivariate) ICA
d4rl = ICA.plots(4,0.2,c(0,11),expression(paste('p=1  ', rho[ij],'=0.2')))
d4rm = ICA.plots(4,0.5,c(0,11),expression(paste('p=1  ', rho[ij],'=0.5')))
d4rh = ICA.plots(4,0.8,c(0,11),expression(paste('p=1  ', rho[ij],'=0.8')))

d8rl = ICA.plots(8,0.2,c(0,20),expression(paste('p=3  ', rho[ij],'=0.2')))
d8rm = ICA.plots(8,0.5,c(0,20),expression(paste('p=3  ', rho[ij],'=0.5')))
d8rh = ICA.plots(8,0.8,c(0,20),expression(paste('p=3  ', rho[ij],'=0.8')))

d12rl = ICA.plots(12,0.2,c(0,30),expression(paste('p=5 ', rho[ij],'=0.2')))
d12rm = ICA.plots(12,0.5,c(0,30),expression(paste('p=5 ', rho[ij],'=0.5')))
d12rh = ICA.plots(12,0.8,c(0,30),expression(paste('p=5 ', rho[ij],'=0.8')))

d22rl = ICA.plots(22,0.2,c(0,100),expression(paste('p=10  ', rho[ij],'=0.2')))
d22rm = ICA.plots(22,0.5,c(0,100),expression(paste('p=10  ', rho[ij],'=0.5')))
d22rh = ICA.plots(22,0.8,c(0,100),expression(paste('p=10  ', rho[ij],'=0.8')))
# Figure 1
grid.arrange(d4rl, d4rm,d4rh,d8rl, d8rm,d8rh,d12rl, d12rm,d12rh,d22rl, d22rm,d22rh, ncol=3)

### correlations density
d8r2 = Corr.plots(8,0.2,r=10,Ylim=c(0,4),MAINTITLE=c('GRS','PC','SHR','SCA'))
d8r5 = Corr.plots(8,0.5,r=10,Ylim=c(0,7),MAINTITLE=c('GRS','PC','SHR','SCA'))
d8r8 = Corr.plots(8,0.8,r=10,Ylim=c(0,18),MAINTITLE=c('GRS','PC','SHR','SCA'))
# Figure A.1
grid.arrange(d8r2[[1]], d8r2[[2]], d8r2[[3]], d8r2[[4]], d8r5[[1]], d8r5[[2]], d8r5[[3]], d8r5[[4]], 
             d8r8[[1]], d8r8[[2]], d8r8[[3]], d8r8[[4]],ncol=4)

addres.d4 = mapply(addres.SHR.SCA,rho=c(0.2,0.5,0.8),MoreArgs = list(d=4))
addres.d8 = mapply(addres.SHR.SCA,rho=c(0.2,0.5,0.8),MoreArgs = list(d=8))
addres.d12 = mapply(addres.SHR.SCA,rho=c(0.2,0.5,0.8),MoreArgs = list(d=12))
addres.d22 = mapply(addres.SHR.SCA,rho=c(0.2,0.5,0.8),MoreArgs = list(d=22))

# Table A.1
cbind(addres.d4,addres.d8,addres.d12,addres.d22)

#### SIMULATION FOR SECTION 5.3
## Covariance matrix for two surrogates (based on variable selection on the transPAT data)
Sigma = matrix( c(1.767791e+02, NA, 9.362633e-02, NA, -9.616815e-04, NA, NA, 11.8330574635, NA, -2.432785e-04, NA, -3.566809e-04,
                  9.362633e-02, NA, 6.212567e-05, NA, -2.708876e-07, NA, NA, -0.0002432785, NA, 2.060978e-07, NA, 3.415315e-08,
                  -9.616815e-04, NA,-2.708876e-07, NA, 1.831039e-08, NA, NA, -0.0003566809, NA, 3.415315e-08, NA, 2.471260e-08), 6,6)
## Generating parameter space
Fit.Total = ICA.contcont.MultS(10000,15,Sigma,Method='GRS',Seed=12345)
# Removing possible NA values
for(i in 1:10000){
  if(length(Fit.Total[[i]])==37){Fit.Total[[i]][38]=NA}
}
write.csv(Fit.Total,'Parameter.space.csv')

Parameter.space = read.csv('Parameter.space.csv')
Parameter.space = Parameter.space[,-1]

### determining 'true' correlation matrix
R.low = Parameter.space[round(Parameter.space[,d*(d+1)/2 +1],5) == 0.7499,1:(d*(d+1)/2)]
R.medium = Parameter.space[round(Parameter.space[,d*(d+1)/2 +1],4) == 0.85,1:(d*(d+1)/2)]
R.high = Parameter.space[round(Parameter.space[,d*(d+1)/2 +1],5) == 0.9499,1:(d*(d+1)/2)]

# means and variances to generate random datasets
Means =  c(2.499532e+01, 6.930726e-03, 3.488782e-04, 1.176323e+01, 1.858426e-04, 3.305136e-04)
Vars = diag(Sigma)
### defining the 'true' covariance matrices
Cov.low = cor2cov(invvech(unlist(R.low)),sqrt(Vars))
Cov.medium = cor2cov(invvech(unlist(R.medium)),sqrt(Vars))
Cov.high = cor2cov(invvech(unlist(R.high)),sqrt(Vars))
# coverage for a sample size of 100
treat = c(rep(-1,50),rep(1,50))
Coverage.low.N100 = CoverageFun(Means,treat, Cov.low,0.75,Method='PC',Seed=123)
Coverage.low.N100[[1]]
# coverage for a sample size of 200
treat = c(rep(-1,100),rep(1,100))
Coverage.low.N200 = CoverageFun(Means,treat, Cov.low,0.75,Method='PC',Seed=123)
Coverage.low.N200[[1]]
# coverage for a sample size of 500
treat = c(rep(-1,250),rep(1,250))
Coverage.low.N500 = CoverageFun(Means,treat, Cov.low,0.75,Method='PC',Seed=123)
Coverage.low.N500[[1]]
load("CoverageData2.RData")

## Table 2.
rbind(Coverage.low.N100[[1]][,4],Coverage.low.N200[[1]][,4],Coverage.low.N500[[1]][,4],
                Coverage.medium.N100[[1]][,4],Coverage.medium.N200[[1]][,4],Coverage.medium.N500[[1]][,4],
                Coverage.high.N100[[1]][,4],Coverage.high.N200[[1]][,4],Coverage.high.N500[[1]][,4])

#### RESULTS FOR SECTION 6.
# reading data
Data.micro = read.csv('Microbiome_data.csv')

Data = Data.micro$True
names(Data) = 'true'
Sur.candidates = Data.micro[,-c(1:3)]
treat = Data.micro$Treat
# variable selection procedure
for(i in 1:67){
  res = Var.selection.ICA(Data,Sur.candidates,treat,M=10000)  
  new.sur = Sur.candidates[,colnames(Sur.candidates)==res$sur.name]
  Sur.candidates = Sur.candidates[,colnames(Sur.candidates)!=res$sur.name]
  Data = cbind(Data,new.sur)
  colnames(Data)[i+1] = res$sur.name
  write.csv(colnames(Data),'Surrogate Selection/SurrogateOrder.txt')
}

### identification of the first five
surr.names = colnames(Data)[2:6]
OTUs.names = colnames(Data.micro[,-c(1:3)])
pos = NULL
for(i in 1:5){
  pos[i] = (1:67)[OTUs.names == surr.names[i]]
}
pos
Data.micro2 = load('otu_full_with_col_column.rda')

surr.names[[1]]
otu[otu$OTU == '221429',][1,]
surr.names[[2]]
otu[otu$OTU == '593043',][1,]$Family
surr.names[[3]]
otu[otu$OTU == 'New.ReferenceOTU334',][1,]$Family
surr.names[[4]]
otu[otu$OTU == '178213',][1,]$Family
surr.names[[5]]
otu[otu$OTU == '135956',][1,]$Family

#### computation time

Time =NULL
NumSurr = c(1,(1:12)*5,67)
for(i in 1:length(NumSurr)){
  labels = paste('Surrogate Selection/Step',NumSurr[[i]],'txt',sep='.')
  res = read.csv(labels)  
  Time[i] = median(res$time)
}
# Figure 2
plot(NumSurr,Time,type='l',xlab='Number of surrogates',ylab='time (in minutes)',lwd=2,cex.lab=1.2)

## Median ICA and range plot
# Figure 3.
par(mfrow=c(1,2))
res.max = as.data.frame(matrix(0,38,7))
plot(NULL,xlim=c(0,67),ylim=c(0,1),xlab='number of surrogates',ylab='',cex.lab=1.5)
for(i in 1:67){
  labels = paste('Surrogate Selection/Step',i,'txt',sep='.')
  res = read.csv(labels)  
  res.max[i,] = as.data.frame(matrix(0,38,7))
  
  points(rep(i,nrow(res)),res$median,col='grey')
  points(i,max(res$median),pch=19)
}
points(10:67,rep(1,58),pch=19)
title(ylab=expression(paste('Median  ',R[H]^2)), line=2.2, cex.lab=1.5)

Range = res.max[,5:6]
plot(NULL,xlim=c(0,67),ylim=c(0,1),xlab='number of surrogates',ylab='',cex.lab=1.5)
for(i in 1:167){lines(c(i,i),Range[i,],lwd=3)}
title(ylab=expression(paste('Range ',R[H]^2)), line=2.2, cex.lab=1.5)

## Additional simulations

## correlation matrix for case: p=10 (large correlations)
d = 22
Rfixed.p10.high = matrix(NA,d,d)
R1 = round(GenCorrMatrix.PC(11,0.7,0.8,9),2)
R2 = round(GenCorrMatrix.PC(11,0.7,0.8,13),2)
IND1 = order(R1[1,],decreasing = T)
IND2 = order(R2[1,],decreasing = T)
R1 = R1[IND1,IND1]
R2 = R2[IND2,IND2]
# high correlations
Rfixed.p10.high[1:(d/2),1:(d/2)] = R1
Rfixed.p10.high[d/2 + 1:(d/2),d/2 + 1:(d/2)] = R2
# medium correlations
Rfixed.d10.medium = Rfixed.d10.high*0.7
diag(Rfixed.d10.medium) = 1
# low correlations
Rfixed.d10.low = Rfixed.d10.high*0.3
diag(Rfixed.d10.low) = 1

## p = 5
Rfixed.d5.high = Rfixed.d10.high[c(1:6,12:17),c(1:6,12:17)]
Rfixed.d5.medium = Rfixed.d10.medium[c(1:6,12:17),c(1:6,12:17)]
Rfixed.d5.low = Rfixed.d10.low[c(1:6,12:17),c(1:6,12:17)]
## p = 3
Rfixed.d3.high = Rfixed.d10.high[c(1:4,12:15),c(1:4,12:15)]
Rfixed.d3.medium = Rfixed.d10.medium[c(1:4,12:15),c(1:4,12:15)]
Rfixed.d3.low = Rfixed.d10.low[c(1:4,12:15),c(1:4,12:15)]
## p = 1
Rfixed.d1.high = Rfixed.d10.high[c(1:2,12:13),c(1:2,12:13)]
Rfixed.d1.medium = Rfixed.d10.medium[c(1:2,12:13),c(1:2,12:13)]
Rfixed.d1.low = Rfixed.d10.low[c(1:2,12:13),c(1:2,12:13)]

GenerateRandCor.Sim2 = function(M,Rfixed,rho,rpc=TRUE,grs=TRUE,sca=TRUE,shr=TRUE){
  d=ncol(Rfixed)
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
    write.csv(RPC,paste('simulations 2/RPC',FILE,'csv',sep = '.'),row.names = FALSE)
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
    write.csv(gRS,paste('simulations 2/GRS',FILE,'csv',sep = '.'),row.names = FALSE)
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
    write.csv(SCA,paste('simulations 2/SCA',FILE,'csv',sep = '.'),row.names = FALSE)
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
    write.csv(SHR,paste('simulations 2/SHR',FILE,'csv',sep = '.'),row.names = FALSE)
  }
}
ICA.plots.Sim2 = function(d,rho,Ylim=NULL,MAIN=NULL){
  # d = 8; rho=0.8
  Rfixed = matrix(NA,d,d)
  R1 = matrix(rho,d/2,d/2);diag(R1) = 1
  Rfixed[1:(d/2),1:(d/2)] = R1
  Rfixed[d/2 + 1:(d/2),d/2 + 1:(d/2)] = R1
  Rfixed[is.na(Rfixed)] = 0
  if(rho == 0.2){rho.lab='low'}
  if(rho == 0.5){rho.lab='medium'}
  if(rho == 0.8){rho.lab='high'}
  #rho.lab=rho
  file = paste('d',d,'rho',rho.lab,sep = '.')
  Sigma = rep(100,d)
  R.RPC = tryCatch({read.csv(paste('Simulations 2/RPC',file,'csv',sep='.'))},error =function(e){NULL})
  R.SCA = tryCatch({read.csv(paste('Simulations 2/SCA',file,'csv',sep='.'))},error =function(e){NULL})
  R.SHR = tryCatch({read.csv(paste('Simulations 2/SHR',file,'csv',sep='.'))},error =function(e){NULL})
  R.GRS = tryCatch({read.csv(paste('Simulations 2/GRS',file,'csv',sep='.'))},error =function(e){NULL},warning =function(w){NULL})
  
  vech.IND = d*(d+1)/2
  ICA.RPC = tryCatch({apply(as.matrix(R.RPC[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SCA = tryCatch({apply(as.matrix(R.SCA[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.SHR = tryCatch({apply(as.matrix(R.SHR[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  ICA.GRS = tryCatch({apply(as.matrix(R.GRS[,1:vech.IND]),1,ICAvech.fun,Sigma=Sigma,N=20)},error =function(e){NULL})
  
  
  ICAs = melt(as.data.frame(cbind(ICA.RPC,ICA.SHR,ICA.SCA,ICA.GRS)))
  PLOT = ggplot(aes(x=value, group=variable,linetype=variable,col=variable), data=ICAs)  + geom_density(show.legend = F,size=1) + 
    theme_bw() + labs(x=expression(R[H]^2),size=5) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=20),
          axis.text = element_text(size = 16))   +
    ylim(Ylim) + xlim(c(0,1)) + ggtitle(MAIN) + theme(plot.title = element_text(hjust = 0.5,size=40)) +
    scale_linetype_manual(values=c('dashed','solid',"dashed", "solid")) +  scale_colour_manual( values = c("black","lightgrey",'lightgrey','black'))
  return(PLOT)
}

#### generating random correlations with fixed values
## One surrogate
GenerateRandCor.Sim2(1000,Rfixed.d1.high,'high',rpc=F,grs=TRUE,sca=F,shr=F)
GenerateRandCor.Sim2(1000,Rfixed.d1.medium,'medium',rpc=F,grs=TRUE,sca=F,shr=F)
GenerateRandCor.Sim2(1000,Rfixed.d1.low,'low',rpc=F,grs=TRUE,sca=F,shr=F)
## Three surrogate
GenerateRandCor.Sim2(1000,Rfixed.d3.high,'high',rpc=F,grs=TRUE,sca=F,shr=F)
GenerateRandCor.Sim2(1000,Rfixed.d3.medium,'medium',rpc=F,grs=TRUE,sca=F,shr=F)
GenerateRandCor.Sim2(1000,Rfixed.d3.low,'low',rpc=F,grs=TRUE,sca=F,shr=F)
## Five surrogate
GenerateRandCor.Sim2(1000,Rfixed.d5.high,'high',grs=FALSE)
GenerateRandCor.Sim2(1000,Rfixed.d5.medium,'medium',rpc=F,grs=TRUE,sca=F,shr=F)
GenerateRandCor.Sim2(1000,Rfixed.d5.low,'low',rpc=F,grs=TRUE,sca=F,shr=F)
## 10 surrogate
GenerateRandCor.Sim2(1000,Rfixed.d10.high,'high',grs=FALSE)
GenerateRandCor.Sim2(1000,Rfixed.d10.medium,'medium',grs=FALSE)
GenerateRandCor.Sim2(1000,Rfixed.d10.low,'low',grs=FALSE)

#### FIGURES: density of (multivariate) ICA
d4rl = ICA.plots(4,0.2,c(0,11),expression(paste('p=1  low ', rho[ij])))
d4rm = ICA.plots(4,0.5,c(0,11),expression(paste('p=1  medium ', rho[ij])))
d4rh = ICA.plots(4,0.8,c(0,11),expression(paste('p=1  high ', rho[ij])))


d8rl = ICA.plots(8,0.2,c(0,20),expression(paste('p=3  low ', rho[ij])))
d8rm = ICA.plots(8,0.5,c(0,20),expression(paste('p=3  medium ', rho[ij])))
d8rh = ICA.plots(8,0.8,c(0,20),expression(paste('p=3  high ', rho[ij])))

d12rl = ICA.plots(12,0.2,c(0,30),expression(paste('p=5  low ', rho[ij])))
d12rm = ICA.plots(12,0.5,c(0,30),expression(paste('p=5  medium ', rho[ij])))
d12rh = ICA.plots(12,0.8,c(0,30),expression(paste('p=5  high ', rho[ij])))


d22rl = ICA.plots(22,0.2,c(0,100),expression(paste('p=10  low ', rho[ij],'=0.2')))
d22rm = ICA.plots(22,0.5,c(0,100),expression(paste('p=10  medium ', rho[ij],'=0.5')))
d22rh = ICA.plots(22,0.8,c(0,100),expression(paste('p=10  high ', rho[ij],'=0.8')))

# Figure A.2
grid.arrange(d4rl, d4rm,d4rh,d8rl, d8rm,d8rh,d12rl, d12rm,d12rh,d22rl, d22rm,d22rh, ncol=3)

## Computational time (Results of Table 2A)
Timep1r0.2 = Computation.time(4,0.2,100,grs=T,pc=T,sca=T,shr=T)
Timep1r0.5 = Computation.time(4,0.5,100,grs=T,pc=T,sca=T,shr=T)
Timep1r0.8 = Computation.time(4,0.8,100,grs=T,pc=T,sca=T,shr=T)
Timep3r0.2 = Computation.time(8,0.2,100,grs=T,pc=T,sca=T,shr=T)
Timep3r0.5 = Computation.time(8,0.5,100,grs=T,pc=T,sca=T,shr=T)
Timep3r0.8 = Computation.time(8,0.8,100,grs=T,pc=T,sca=T,shr=T)
Timep5r0.2 = Computation.time(12,0.2,100,grs=T,pc=T,sca=T,shr=T)
Timep5r0.5 = Computation.time(12,0.5,100,grs=T,pc=T,sca=T,shr=T)
Timep5r0.8 = Computation.time(12,0.8,100,grs=T,pc=T,sca=T,shr=T)
Timep10r0.2 = Computation.time(22,0.2,100,grs=F,pc=T,sca=T,shr=T)
Timep10r0.5 = Computation.time(22,0.5,100,grs=F,pc=T,sca=T,shr=T)
Timep10r0.8 = Computation.time(22,0.8,100,grs=F,pc=T,sca=T,shr=T)
#Table A.2
Table = rbind(
  c(mapply(function(x){summary(x,unit='s')$median},Timep1r0.2),mapply(function(x){summary(x,unit='s')$median},Timep1r0.5),
    mapply(function(x){summary(x,unit='s')$median},Timep1r0.8)),
  c(mapply(function(x){summary(x,unit='s')$median},Timep3r0.2),mapply(function(x){summary(x,unit='s')$median},Timep3r0.5),
    mapply(function(x){summary(x,unit='s')$median},Timep3r0.8)),
  c(mapply(function(x){summary(x,unit='s')$median},Timep5r0.2),mapply(function(x){summary(x,unit='s')$median},Timep5r0.5),
    mapply(function(x){summary(x,unit='s')$median},Timep5r0.8)),
  c(mapply(function(x){summary(x,unit='s')$median},Timep10r0.2),mapply(function(x){summary(x,unit='s')$median},Timep10r0.5),
    mapply(function(x){summary(x,unit='s')$median},Timep10r0.8)))

