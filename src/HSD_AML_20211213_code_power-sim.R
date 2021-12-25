#!/bin/RScript
#HSD
#AML
#Cornell University

#Set I/O----
rm(list = ls())
#Packages----

#Functions----
source("HSD-AML-20211213-ML-v7.R") #load functions

start_sim<-function(x){ #start simulation, generate phenotype frequency vector
  y<-c(x,1-x)
  names(y)<-c("r","wt")
  return(y)
}

#calculate CI based on FI
#https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
calculate.ci<-function(my.hessian){
  if(det(my.hessian)==0){
    return(NA)
  }
  fisher<-solve(c(-1)*my.hessian)
  sigma_estimate<-sqrt(diag(fisher))
  return(sigma_estimate*1.96)
}

#call_sim ... simulates data using sim_data() for <nrep> replicates and a <ncencus> census size
#p ... parental phenotype counts
#g ... names of genotypes
#fne ... effective population size
#f .... fitness params (named!)
#fd... gene drive fitness cost
#n ... number of generations to simulate
#fmpath ... path to transition matrix
#fh ... ratio of heterozygous flies at gen. 0
#fdr ...gene drive conversion rate
#fdrr ... gene drive resistance rate
#fer ... embryo resistance rate
#fr1 ... rate of r1 resistance mutations
#nrep ... number or replicates
#ncensus ... census size
#fmode ... mode for fitness estimates
#fmodesim ... id for type of simulated data
#fnsim ... simulation iteration number
call_sim<-function(p,g,fne,f,fd,n,fmpath,fh,fdr,fdrr,fer,fr1,nrep,ncensus,fmode,fmodesim,fnsim){
  N<-replicate(nrep,sim_data(p = p,g = g,fne = fne,f = f,fd = fd,n = n,fmpath = fmpath,fh = fh,fdr = fdr,fdrr = fdrr,fer = fer,fr1 = fr1),simplify = F) #generate simulated data
  N<-lapply(N,sim_geno_to_pheno) #calculate phenotype

  N<-lapply(N,function(k) { #calculate counts
    k[["r"]]<-k[["r"]]*ncensus
    k[["wt"]]<-k[["wt"]]*ncensus
    return(k)})
  #print(N)
  my_lim<-set_limits(fmode = fmode,fd.lim = 0.1) #set limits for ML
  #print(my_lim)
  ml_out<-NULL

  ml_out<-tryCatch( #trial of ML run
    {
      ml<-optim(my_lim[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1), #ML estimation
                lower=my_lim[[1]],upper = my_lim[[3]],lN=N,ld=fd,
                lmpath=fmpath,lg=g,ldr=fdr,ldrr=fdrr,ler=fer,lh=fh,lr1=fr1,lmode=fmode,hessian = T)
    },
    error=function(cond)  #if error, return NA
    {
      message(cond)
      return(NULL)
    })


  if(!is.null(ml_out)) { #if valid output, store results

    ml_CI<-calculate.ci(ml$hessian)
    CIu<-ml$par+ml_CI
    CId<-ml$par-ml_CI

    res<-data.frame("sim"=fnsim,"est.ne"=ml$par[1],"est.p1"=ml$par[2],
                  "ci.ne.upper"=CIu[1], "ci.ne.lower"=CId[1],
                  "ci.p1.upper"=CIu[2], "ci.p1.lower"=CId[2],
                  "L"=ml$value, "conv"=ml$convergence,
                  "mode"=fmode,"data"=fmodesim,"ngen"=sum(unlist(lapply(N,nrow)))-length(N))
  }

  if (is.null(ml_out)) { #if non valid, return df with NA
    res<-data.frame("sim"=fnsim, "est.ne"=NA, "est.p1"=NA,
                    "ci.ne.upper"=NA, "ci.ne.lower"=NA,
                    "ci.p1.upper"=NA, "ci.p1.lower"=NA,
                    "L"=NA, "conv"=NA,"mode"=fmode,"data"=fmodesim,"ngen"=sum(unlist(lapply(N,nrow)))-length(N))
  }
  return(res)
}

#launch_sim ... call call_sim() <n.sim> times
#p ... parental phenotype counts
#g ... names of genotypes
#fne ... effective population size
#p1 ... input parameter 1 (viability, fecundity)
#n ... number of generations to simulate
#fmpath ... path to transition matrix
#fh ... ratio of heterozygous flies at gen. 0
#fdr ...gene drive conversion rate
#fdrr ... gene drive resistance rate
#fer ... embryo resistance rate
#fr1 ... rate of r1 resistance mutations
#nrep ... number or replicates
#ncensus ... census size
#fmode ... mode for fitness estimates
#fmodesim ... id for type of simulated data
#n.sim ... number of simulations

launch_sim<-function(p,g,fne,p1,n,fmpath,fh,fdr,fdrr,fer,fr1,nrep,ncensus,fmode,fmodesim,n.sim) {

  f<-c(fne,rep(1,3)) #determine fitness vector for simulations
  names(f)<-c("ne","fdm","fdf","fdv")

  if(fmodesim==3){
    D<-rep("c",3) #determine fitness cost vector
    f["fdv"]<-p1 #vary viability
  }

  if(fmodesim==10){
    D<-c("n","l","n")
    f["fdf"]<-p1 #vary fecundity , l = leaky somatic expression
  }

  runs<-sapply(c(1:n.sim), function(j) call_sim(p = p,g = g,fne = fne,f = f,fd = D ,n = n,
                                                fmpath = fmpath,fh = fh,fdr = fdr,fdrr = fdrr,
                                                fer = fer,fr1 = fr1,nrep = nrep,ncensus = ncensus,fmode = fmode,fmodesim = fmodesim,fnsim = j),simplify = F)
  runs<-do.call("rbind",runs)
  return(runs)
}


#Global Variables----
DR<-0.767 #drive conversion rate
DRR<-0.22 #germline resistance rate
ER<-0.522 #embryo resistance rate
H<-1 #ratio of drive heterozygotes in generation 0
P0<-0.1 #drive frequency, generation 0
R1<-0 #r1 formation rate 

#Variable Params----
args = commandArgs(trailingOnly=TRUE)
params <- read.table(file = args[1],header = F,sep = ",")
colnames(params)<-c("ncensus","ne","nrep","ngen","nsim","p1","mode") #1st entry: parameters
iter <- as.integer(args[2]) #2nd entry: line to run
params<-params[iter,]

#Simulations----
p<-start_sim(P0) #drive frequency

my.sim<-launch_sim(p = p,g = G,fne = params[["ne"]],p1 = params[["p1"]],n = params[["ngen"]], #simulation
           fmpath = M.path,fh = H,fdr = DR,fdrr = DRR,fer = ER,fr1 = R1, nrep = params[["nrep"]],
           ncensus = params[["ncensus"]],fmode = params[["mode"]],fmodesim = params[["mode"]],n.sim = params[["nsim"]])

my.sim$ne<-rep(params[["ne"]],nrow(my.sim)) #add parameters to result
my.sim$p1<-rep(params[["p1"]],nrow(my.sim))
my.sim$nreps<-rep(params[["nrep"]],nrow(my.sim))
my.sim$n<-rep(params[["ncensus"]],nrow(my.sim))
my.sim$mode<-rep(params[["mode"]],nrow(my.sim))

#print(head(my.sim))

#Output----
tosave<-paste("p1",params$p1,"ne",params$ne,"n",params$ncensus,
              "reps",params$nrep,"gens",params$ngen,"mode",params$mode,sep="-")
tosave<-paste(tosave,".rds",sep = "")
print(tosave)

saveRDS(my.sim, tosave)
