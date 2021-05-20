#!/bin/RScript
#HSD
#AML

rm(list=ls())
wd<-"./" #set working directory
setwd(wd)

#Packages----
library(ggplot2)

#ML framework----
source("HSD-AML-20210305-ML-v2.R")

#Global Variables----

CAGE1.PATH<-"../data/raw_cage1.txt" #empirical data
CAGE2.PATH<-"../data/raw_cage2.txt"
CAGE3.PATH<-"../data/raw_cage3.txt"

DR.a<-0.767 #drive conversion rate, cage 1 & 2
DRR.a<-0.22 #germline resistance formation rate, cage 1 & 2
ER.a<-0.522 #embryo resistance formation rate, cage 1 & 2 
H.a<-1 #ratio of drive heterozygotes in generation 0, cage 1 & 2

DR.b<-0 #drive conversion rate, cage 3
DRR.b<-0 #germline resistance formation rate, cage 3
ER.b<-0 #embryo resistance formation rate, cage 3
H.b<-0.59 #ratio of drive heterozygotes in generation 0, cage 3 

X.ab <- 1 #germline cut rate, off-target, all cages 
E.ab <- 0 #embryo cut rate, off-target, all cages
T.ab <- 0.5 #already cut in generation 0, off-target, all cages 

MODES.a<-c(0,1,2,3,4,5,6,10) #tested modes, cage 1 & 2 
MODES.b<-c(0,1,3,10) #tested modes, cage 3 

D.ab<-rep("c",3) #fitness modes (mate choice, fecundity, viability) drive, all cages 
O.ab<-rep("c",3) #fitness modes (mate choice, fecundity, viability) off-target, all cages

DIST.ab<-0.01 #convergence distance for CI determination, all cages
A.ab<-0.05 #type 1 error, all cages 
RM.START<-1 #number of generations to remove from the beginning of an experiment
D.LIM<-0 #add to 1 for upper limit for drive fitness estimate in ML
O.LIM<-0 #add to  1 for upper limit for off-target fitness estimate in ML

TOSTORE.ab<-"./" #store ML output, all cages


#Additional Functions----

#read_emp_data
#x ... path to empirical cage data (.txt)
#return(y) ... data frame formatted for ML approach
read_emp_data<-function(x){
  y<-read.table(x)
  colnames(y)<-c("g","r","wt","N")
  y$N<-NULL
  return(y)
  }

#plot_emp_data
#x ... output from read_emp_data()
#return(y) ... drive carrier frequency plot over time 
plot_emp_data<-function(x){
  y<-x
  y$r<-x$r/(x$r+x$wt)
  y$wt<-x$wt/(x$wt+x$r)
  y$g<-as.factor(y$g)
  g<-ggplot(data=y,aes(x=g,y=r))+geom_line(col="red",group=1)+geom_point(col="red")
  g<-g+theme_minimal()+coord_cartesian(ylim = c(0,1),xlim=c(0,25))
  g<-g+scale_y_continuous("f(r)",breaks = c(0,0.25,0.5,0.75,1))
  g<-g+xlab("Generation")+theme(text=element_text(size = 15))
  g<-g+scale_x_discrete("Generation",breaks=seq(0,25,5))
  return(g)
}

#Load & EDA Data----

(cage1<-read_emp_data(CAGE1.PATH))
(cage2<-read_emp_data(CAGE2.PATH))
(cage3<-read_emp_data(CAGE3.PATH))


#plot_emp_data(cage1)
#plot_emp_data(cage2)
#plot_emp_data(cage3)

if(RM.START>0){ #remove generations at the beginning of experiment (optional)
  cage1<-cage1[-c(1:RM.START),]
  cage1$g<-cage1$g-RM.START
  cage2<-cage2[-c(1:RM.START),]
  cage2$g<-cage2$g-RM.START
  cage3<-cage3[-c(1:RM.START),]
  cage3$g<-cage3$g-RM.START
}

#ML cage 1&2 (joint analysis)----
obs.data<-list()
obs.data[[1]]<-cage1
obs.data[[2]]<-cage2


for(mode.iter in MODES.a){
  print(paste("Mode:",mode.iter)) #determine limits
  my.res.iter<-list() #set storage for individual runs 
  toStore<-paste(TOSTORE.ab,"ML-cage12-mode",mode.iter,".rds",sep="")
  my.lim.iter<-set_limits(fmode = mode.iter,fd.lim = D.LIM,fo.lim = O.LIM) #determine parameters to estimate
  print("Limits:")
  print(my.lim.iter)
  print(paste("Storage:",toStore))
  
  if(mode.iter==10) { #reset D.ab for leaky somatic expression
    D.ab<-c("n","l","n")
  } else {
    D.ab<-rep("c",3)
  }

  #ML estimate: 
  ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
            lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
            lN=obs.data,ld=D.ab,lo=O.ab,lt=rep(T.ab,2),lmpath=M.path,
            lg=G,lx=X.ab,le=E.ab,ldr=DR.a,ldrr=DRR.a,ler=ER.a,lh=H.a,lmode=mode.iter)
  print(ml)
  
  #AICc calculation:
  n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
  print(paste("Number of transistions:",n.trans))
  aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
  print(paste("AICc:",aicc))
  
  #CI estimate:
  ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                   fd = D.ab,fo = O.ab,fg = G,ft = rep(T.ab,2),fmpath = M.path,
                   fx = X.ab,fe = E.ab,fdr = DR.a,fdrr = DRR.a,fer = ER.a,
                   fh = H.a,fmode = mode.iter)
  print("CI interval:")
  print(ci)
  
  my.res.iter[[1]]<-ml
  my.res.iter[[2]]<-aicc
  my.res.iter[[3]]<-ci
  
  saveRDS(my.res.iter,toStore)
  rm(my.res.iter,ci,aicc,n.trans,ml,toStore,my.lim.iter)
  print("------------------")
}

#ML cage 3----
obs.data<-list()
obs.data[[1]]<-cage3

for(mode.iter in MODES.b){
  
  print(paste("Mode:",mode.iter)) #determine limits
  my.res.iter<-list() #set storage for individual runs 
  toStore<-paste(TOSTORE.ab,"ML-cage3-mode",mode.iter,".rds",sep="")
  my.lim.iter<-set_limits(fmode = mode.iter,fd.lim = D.LIM,fo.lim = O.LIM) #determine parameters to estimate
  print("Limits:")
  print(my.lim.iter)
  print(paste("Storage:", toStore))
  
  
  if(mode.iter==10) { #reset D.ab for leaky somatic expression
    D.ab<-c("n","l","n")
  } else {
    D.ab<-rep("c",3)
  }
  
  #ML estimate: 
  ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
            lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
            lN=obs.data,ld=D.ab,lo=O.ab,lt=T.ab,lmpath=M.path,
            lg=G,lx=X.ab,le=E.ab,ldr=DR.b,ldrr=DRR.b,ler=ER.b,lh=H.b,lmode=mode.iter)
  print(ml)
  
  #AICc calculation:
  n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
  print(paste("Number of transitions:",n.trans))
  aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
  print(paste("AICc:",aicc))
  
  #CI estimate:
  ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                   fd = D.ab,fo = O.ab,fg = G,ft = T.ab,fmpath = M.path,
                   fx = X.ab,fe = E.ab,fdr = DR.b,fdrr = DRR.b,fer = ER.b,
                   fh = H.b,fmode = mode.iter)
  print("CI interval:")
  print(ci)
  
  my.res.iter[[1]]<-ml
  my.res.iter[[2]]<-aicc
  my.res.iter[[3]]<-ci
  
  saveRDS(my.res.iter,toStore)
  rm(my.res.iter,ci,aicc,n.trans,ml,toStore,my.lim.iter)
  print("------------------")
  }

#ML cage 2 (separate)----
obs.data<-list()
obs.data[[1]]<-cage2


for(mode.iter in MODES.a){
  print(paste("Mode:",mode.iter)) #determine limits
  my.res.iter<-list() #set storage for individual runs 
  toStore<-paste(TOSTORE.ab,"ML-cage2-mode",mode.iter,".rds",sep="")
  my.lim.iter<-set_limits(fmode = mode.iter,fd.lim = D.LIM,fo.lim = O.LIM) #determine parameters to estimate
  print("Limits:")
  print(my.lim.iter)
  print(paste("Storage:",toStore))
  
  
  if(mode.iter==10) { #reset D.ab for leaky somatic expression
    D.ab<-c("n","l","n")
  } else {
    D.ab<-rep("c",3)
  }
  
  #ML estimate: 
  ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
            lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
            lN=obs.data,ld=D.ab,lo=O.ab,lt=T.ab,lmpath=M.path,
            lg=G,lx=X.ab,le=E.ab,ldr=DR.a,ldrr=DRR.a,ler=ER.a,lh=H.a,lmode=mode.iter)
  print(ml)
  
  #AICc calculation:
  n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
  print(paste("Number of transistions:",n.trans))
  aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
  print(paste("AICc:",aicc))
  
  #CI estimate:
  ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                   fd = D.ab,fo = O.ab,fg = G,ft = T.ab,fmpath = M.path,
                   fx = X.ab,fe = E.ab,fdr = DR.a,fdrr = DRR.a,fer = ER.a,
                   fh = H.a,fmode = mode.iter)
  print("CI interval:")
  print(ci)
  
  my.res.iter[[1]]<-ml
  my.res.iter[[2]]<-aicc
  my.res.iter[[3]]<-ci
  
  saveRDS(my.res.iter,toStore)
  rm(my.res.iter,ci,aicc,n.trans,ml,toStore,my.lim.iter)
  print("------------------")
}

#ML cage 1 (separate analysis)----
obs.data<-list()
obs.data[[1]]<-cage1


for(mode.iter in MODES.a){
  print(paste("Mode:",mode.iter)) #determine limits
  my.res.iter<-list() #set storage for individual runs 
  toStore<-paste(TOSTORE.ab,"ML-cage1-mode",mode.iter,".rds",sep="")
  my.lim.iter<-set_limits(fmode = mode.iter,fd.lim = D.LIM,fo.lim = O.LIM) #determine parameters to estimate
  print("Limits:")
  print(my.lim.iter)
  print(paste("Storage:",toStore))
  
  
  if(mode.iter==10) { #reset D.ab for leaky somatic expression
    D.ab<-c("n","l","n")
  } else {
    D.ab<-rep("c",3)
  }
  
  #ML estimate: 
  ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
            lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
            lN=obs.data,ld=D.ab,lo=O.ab,lt=T.ab,lmpath=M.path,
            lg=G,lx=X.ab,le=E.ab,ldr=DR.a,ldrr=DRR.a,ler=ER.a,lh=H.a,lmode=mode.iter)
  print(ml)
  
  #AICc calculation:
  n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
  print(paste("Number of transistions:",n.trans))
  aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
  print(paste("AICc:",aicc))
  
  #CI estimate:
  ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                   fd = D.ab,fo = O.ab,fg = G,ft = T.ab,fmpath = M.path,
                   fx = X.ab,fe = E.ab,fdr = DR.a,fdrr = DRR.a,fer = ER.a,
                   fh = H.a,fmode = mode.iter)
  print("CI interval:")
  print(ci)
  
  my.res.iter[[1]]<-ml
  my.res.iter[[2]]<-aicc
  my.res.iter[[3]]<-ci
  
  saveRDS(my.res.iter,toStore)
  rm(my.res.iter,ci,aicc,n.trans,ml,toStore,my.lim.iter)
  print("------------------")
}
