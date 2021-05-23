#!/bin/RScript
#HSD
#AML

rm(list=ls())
wd<-"./" #choose your working directory
setwd(wd)
getwd()

#Packages----
library(ggplot2)

#ML framework----
source("HSD-AML-20210519-ML-v5.R")

#Global Variables----

CAGE1.PATH<-"../data/raw_cage1.txt" #empirical data homing drives
CAGE2.PATH<-"../data/raw_cage2.txt"

DR.a<-0.767 #drive conversion rate, cage 1 & 2
DRR.a<-0.22 #germline resistance formation rate, cage 1 & 2
ER.a<-0.522 #embryo resistance formation rate, cage 1 & 2 
H.a<-1 #ratio of drive heterozygotes in generation 0, cage 1 & 2

MODES.cage1<-c(12) #tested mode, cage 1
MODES.cage2<-c(11) #tested mode, cage 2

D.cage1<-c("n","l","n") #fitness modes (mate choice, fecundity, viability) cage 1
D.cage2<-rep("c",3) #fitness modes (mate choice, fecundity, viability) cage 2
R1.ab<-0 #r1 resistance formation rate, all cages (will be estimated)

DIST.ab<-0.001 #convergence distance for CI determination, all cages
A.ab<-0.05 #type 1 error, all cages 
RM.START<-1 #number of generations to remove from the beginning of an experiment
D.LIM<-0 #add to 1 for upper limit for drive fitness estimate in ML

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


plot_emp_data(cage1)
plot_emp_data(cage2)

if(RM.START>0){ #remove generations at the beginning of experiment (optional)
  cage1<-cage1[-c(1:RM.START),]
  cage1$g<-cage1$g-RM.START
  cage2<-cage2[-c(1:RM.START),]
  cage2$g<-cage2$g-RM.START
}

#ML cage 2 (separate)----
obs.data<-list()
obs.data[[1]]<-cage2


for(mode.iter in MODES.cage2){
  print(paste("Mode:",mode.iter)) #determine limits
  my.res.iter<-list() #set storage for individual runs 
  toStore<-paste(TOSTORE.ab,"ML-cage2-mode",mode.iter,".rds",sep="")
  my.lim.iter<-set_limits(fmode = mode.iter,fd.lim = D.LIM) #determine parameters to estimate
  print("Limits:")
  print(my.lim.iter)
  print(paste("Storage:",toStore))
  
  #ML estimate: 
  ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
            lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
            lN=obs.data,ld=D.cage2,lmpath=M.path,
            lg=G,ldr=DR.a,ldrr=DRR.a,ler=ER.a,lh=H.a,lr1=R1.ab,lmode=mode.iter)
  print(ml)
  
  #AICc calculation:
  n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
  print(paste("Number of transistions:",n.trans))
  aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
  print(paste("AICc:",aicc))
  
  #CI estimate:
  ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                   fd = D.cage2,fg = G,fmpath = M.path,
                   fdr = DR.a,fdrr = DRR.a,fer = ER.a,
                   fh = H.a,fr1 = R1.ab,fmode = mode.iter)
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


for(mode.iter in MODES.cage1){
  print(paste("Mode:",mode.iter)) #determine limits
  my.res.iter<-list() #set storage for individual runs 
  toStore<-paste(TOSTORE.ab,"ML-cage1-mode",mode.iter,".rds",sep="")
  my.lim.iter<-set_limits(fmode = mode.iter,fd.lim = D.LIM) #determine parameters to estimate
  print("Limits:")
  print(my.lim.iter)
  print(paste("Storage:",toStore))
  
  #ML estimate: 
  ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
            lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
            lN=obs.data,ld=D.cage1,lmpath=M.path,
            lg=G,ldr=DR.a,ldrr=DRR.a,ler=ER.a,lh=H.a,lr=R1.ab,lmode=mode.iter)
  print(ml)
  
  #AICc calculation:
  n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
  print(paste("Number of transistions:",n.trans))
  aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
  print(paste("AICc:",aicc))
  
  #CI estimate:
  ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                   fd = D.cage1,fg = G,fmpath = M.path,
                   fdr = DR.a,fdrr = DRR.a,fer = ER.a,
                   fh = H.a,fr1 = R1.ab,fmode = mode.iter)
  print("CI interval:")
  print(ci)
  
  my.res.iter[[1]]<-ml
  my.res.iter[[2]]<-aicc
  my.res.iter[[3]]<-ci
  
  saveRDS(my.res.iter,toStore)
  rm(my.res.iter,ci,aicc,n.trans,ml,toStore,my.lim.iter)
  print("------------------")
}
