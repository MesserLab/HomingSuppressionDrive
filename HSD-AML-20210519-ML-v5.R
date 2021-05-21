#!/bin/RScript
#AML
#HSD
#ML framework
#v4 -> v5: change way of r1 calculation 

#Packages----
library(stringr)
library(ggplot2)

#Global Variables----
M.path<-"../data/HSD_AML_20210518_data_trans-matrix-FILLED.csv" 

#w=wildtype allele
#d=drive allele
#r1=resistance allele, gene functioning
#r2=resistance allele, gene not functioning 

G<-c("ww","wd","wr1","wr2","dd","dr1","dr2","r1r1","r1r2","r2r2")

#Drive Functions----

#init_m
#x ... path to raw transition matrix (counts)
#y ... names of count columns
#return(z) ... expected offspring frequencies
init_m<-function(x,y){
  z<-read.csv(x,header = T,sep=";")
  z<-cbind(z[,1:2],z[,c(3:12)]/rowSums(z[,c(3:12)]))
  colnames(z)<-c("m","f",y)
  return(z)
}

#init_genotype
#x ... names of genotypes
#return(y) ... empty named genotype vector
init_genotype<-function(x){
  y<-rep(0,10)
  names(y)<-x
  return(y)
}


#germ_drive
#v ... original offspring frequencies (length=12)
#x ... gene drive conversion rate 
#y ... total germline resistance formation rate
#z... r1 resistance formation rate (r1 rate = z*y, r2 rate = (1-z)*y)
germ_drive<-function(v,x,y,z){
  w<-v
  r1<-z*y #r1 resistance formation rate 
  r2<-(1-z)*y #r2 resistance formation rate 
  w["wd"]<-v["wd"]*(1-x-r1-r2)
  w["dr1"]<-v["dr1"]+v["wd"]*r1
  w["dr2"]<-v["dr2"]+v["wd"]*r2
  w["dd"]<-v["dd"]+v["wd"]*x
  if(all.equal(sum(w),1,tolerance=1e-5)==F) print("Normalization problem!")
  return(w)
}

#embryo_drive
#v ... original offspring frequencies (length=10)
#e ... total embryo resistance formation rate given mom has at least one drive allele
#z ... r1 resistance formation rate (r1 rate = z*e, r2 rate = (1-z)*y)
embryo_drive<-function(v,e,z){
  w<-v
  r1<-z*e #r1 formation rate = fraction of r2 formation rate 
  r2<-e #r2 formation rate experimentally provided
  e.total<-r1+r2
  w["ww"]<-v["ww"]*(1-e.total)^2
  w["wd"]<-v["wd"]*(1-e.total)
  w["wr1"]<-v["wr1"]*(1-e.total)+v["ww"]*2*(1-e.total)*r1
  w["wr2"]<-v["wr2"]*(1-e.total)+v["ww"]*2*(1-e.total)*r2
  w["r1r1"]<-v["r1r1"]+v["ww"]*r1^2+v["wr1"]*r1
  w["r2r2"]<-v["r2r2"]+v["ww"]*r2^2+v["wr2"]*r2
  w["r1r2"]<-v["r1r2"]+v["wr1"]*r2+v["wr2"]*r1+v["ww"]*r1*r2*2
  w["dr1"]<-v["dr1"]+v["wd"]*r1
  w["dr2"]<-v["dr2"]+v["wd"]*r2
  if(all.equal(sum(w),1,tolerance=1e-5)==F) print("Normalization problem!")
  return(w)
}

#mod_transition_m
#m ... transition matrix 
#g ... names of genotypes
#dr ... drive conversion rate
#drr ... drive resistance formation rate
#er ... embryo resistance formation rate 
#r1 .... r1 resistance formation rate : r1rate= r1*drr, r1rate=r1*er, r2rate=(1-r1)*drr, r2rate=(1-r1)*er
#return(y) ... transition matrix with adapted offspring frequencies after germline + embryo cutting
mod_transition_m<-function(m,g,dr,drr,er,r1){
  y<-m
  #Step 1: Germline cuts 
  for(i in c(1:nrow(y))){
    fy_iter<-y[i,] #store raw  transition probabilities + cross
    fmale_genotype<-init_genotype(x = g) #empty male genotype vector
    ffemale_genotype<-init_genotype(x = g) #empty female genotype vector
    fmale_genotype[fy_iter[["m"]]]<-1 #set male genotype
    ffemale_genotype[fy_iter[["f"]]]<-1 #set female genotype
    fmale_genotype<-germ_drive(fmale_genotype,x = dr ,y = drr,z = r1) #drive conversion + resistance
    ffemale_genotype<-germ_drive(ffemale_genotype,x = dr, y = drr,z = r1) #drive conversion + resistance
    fnew_mates<-as.vector(outer(ffemale_genotype,fmale_genotype))  #all possible mating probs
    fnew_progeny<-fnew_mates*m[,(3:12)] #determine new weights
    y[i,(3:12)]<-colSums(fnew_progeny)/sum(fnew_progeny) #replace marginal probs for offspring
  }
  
  #Stemp 2: Embryo cuts
  fd_carrier<-grep("d",colnames(y))-2
  for(i in c(1:nrow(y))){
    if(y$f[i]%in%fd_carrier){
      y[i,3:12]<-embryo_drive(v = y[i,3:12],e = er,z = r1)

    }
  }
  return(y)
}

#determine_fitness_coefficients
#m... determines the fitness costs applied 
#fg ... names of genotypes
#a ... determine allele for which fitness costs apply
#p... fitness cost parameters
#return(x) ... fitness vector for genotypes (n=10); rel weights
determine_fitness_coefficients<-function(m,fg,a,p){
  x<-rep(1,10)
  names(x)<-fg
  fa_count<-str_count(names(x),a)
  if(m=="n"){ #neutral -> nothing
  }
  if(m=="d"){ #dominant
    x[which(fa_count>=1)]<-p[1]
  } 
  if(m=="c"){ #co-dominant
    x[which(fa_count==2)]<-p[1]
    x[which(fa_count==1)]<-sqrt(p[1])
  }
  if(m=="r"){ #recessive
    x[which(fa_count==2)]<-p[1]
  }
  if(m=="m") { #multiplicative
    x[which(fa_count==1)]<-p[1]
    x[which(fa_count==2)]<-p[1]^2
  }
  if(m=="l") { #leaky somatic expression
    #w|a are affected
    wa<-paste("w",a,sep = "")
    fa_count<-str_count(names(x),wa)
    x[which(fa_count==1)]<-p[1]
  }
  return(x)
}

#female_fertility
#fg ... names of genotypes
#return(x) ... fitness vector for genotypes (n=10); rel weights 
#any genotype being a combination of drive + r2 resistance allele is sterile 
female_fertility<-function(fg){
  x<-rep(1,10)
  names(x)<-fg
  sterile<-sapply(c("dd","dr2","r2r2"), function(k) grep(k,names(x)))
  sterile<-as.vector(sterile)
  x[sterile]<-0
  return(x)
}

#geno_to_pheno
#x ... named vector of genotypes (m+f,length=20)
#return= vector of phenotype frequencies | x (r,wt)
geno_to_pheno<-function(x){
  fd_id<-grep("d",names(x))
  r<-sum(x[fd_id])
  fwt<-sum(x[-fd_id])
  y<-c(r,fwt)
  y<-y/sum(y)
  names(y)<-c("r","wt")
  return(y)
}

#pheno_to_geno
#g ... named vector of genotypes (m+f, length=20)
#p ... observed phenotypes (r,wt)
#return(fgcorr) ... vector of genotype frequencies, corrected by observed phenotypes (m+f,length=20)
pheno_to_geno<-function(g,p){
  fpobs<-p/sum(p)
  fpraw<-geno_to_pheno(g) #marginal sums of genotypes giving the same phenotype 
  fgcorr<-numeric(length = 20)
  names(fgcorr)<-names(g)
  fd_id<-grep("d",names(fgcorr))
  fgcorr[fd_id]<-g[fd_id]/fpraw["r"]*fpobs["r"]   #standardize 
  fgcorr[-fd_id]<-g[-fd_id]/fpraw["wt"]*fpobs["wt"]
  return(fgcorr)
}

#start_genotype .... function to start a genotype 
#x ... counts (r|wt)
#fg ... genotype names
#h... ratio of heterozygous flies 
#return(y)... genotype frequencies m+f (length=20)
start_genotype<-function(x,fg,h){
  x<-unlist(x)
  names(x)<-c("r","wt")
  x<-x/sum(x)
  y<-numeric(length = 10)
  names(y)<-fg
  y["ww"]<-x["wt"]/2
  y["wd"]<-x["r"]*h/2
  y["dd"]<-x["r"]*(1-h)/2
  y<-rep(y,2)
  return(y)
}

#propagate_genotypes
#p ... parental genotypes (20 (m+f))
#m ... transition matrix (modified)
#g ... names of genotypes
#fd ... gen drive fitness cost mode (mate choice, fecundity, viability)
#f... fitness parameters (named!)
#return=expected genotype frequencies next generation (36 (m+f))
propagate_genotypes<-function(p,m,g,fd,f){
  
  fmale<-p[1:10]/sum(p[1:10]) #standardize male freqs
  ffemale<-p[11:20]/sum(p[11:20]) #standardize female freqs
  
  #mate choice selection 
  fmate_cas9<-determine_fitness_coefficients(m = fd[1],fg = g,a = "d",p = f[c("fdm")])
  fmale<-fmale*fmate_cas9/sum(fmale*fmate_cas9)
  fpairs<-as.vector(outer(ffemale,fmale))
  fmmate<-fpairs*m[,3:12]
  
  #fecundity selection
  ffec_cas9<-determine_fitness_coefficients(m = fd[2],fg = g,a = "d",p = f[c("fdf")])
  ffec_cas9<-rep(ffec_cas9,length(fmale))
  ffec_sterile<-female_fertility(g) #sterility caused by drive
  ffec_sterile<-rep(ffec_sterile,length(fmale))
  ffec<-ffec_cas9*ffec_sterile
  fmfec<-ffec*fmmate
  
  #viability selection
  fvia_cas9<-determine_fitness_coefficients(m = fd[3],fg = g,a = "d",p = f[c("fdv")])
  fmvia<-t(fvia_cas9*t(fmfec))
  
  #extend to males + females again 
  g<-colSums(fmvia)/sum(fmvia)
  g<-rep(g,2)
  g<-g/sum(g)
  return(g)
  
}

#Sim Functions ----
#sim_data
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

sim_data<-function(p,g,fne,f,fd,n,fmpath,fh,fdr,fdrr,fer,fr1){
  r<-p/sum(p)
  r<-start_genotype(x = r,fg = g,h = fh)
  y<-r
  m<-init_m(x = fmpath,y = g) #raw transition matrix
  fmmod<-mod_transition_m(m = m, g = g,dr = fdr,drr = fdrr ,er = fer,r1 = fr1) #drive conversion + resistance
  for(i in c(1:n)){ #for each generation 
    fp1<-propagate_genotypes(p = r,m = fmmod,g = g,fd = fd,f = f)
    fpm<-fp1[1:10]/sum(fp1[1:10])
    fpf<-fp1[11:20]/sum(fp1[11:20])
    if(fne>0){
      gm<-t(rmultinom(1,fne/2,fpm))
      gm<-gm/sum(gm)
      gf<-t(rmultinom(1,fne/2,fpf))
      gf<-gf/sum(gf)
    } else {
      gm<-fpm
      gf<-fpf
    }
    r<-c(gm,gf)
    r<-r/sum(r)
    y<-rbind(y,r)
  }
  rownames(y)<-c(0:n)
  colnames(y)<-rep(g,2)
  return(y)
}

#sim_geno_to_pheno 
#x ... output from sim_data
#return(y) ...data frame with phenotypes (r+wt)
sim_geno_to_pheno<-function(x){
  g<-colnames(x)
  y<-data.frame(g=numeric(),r=numeric(),wt=numeric())
  for(i in c(1:nrow(x))){
    xi<-x[i,]
    names(xi)<-g
    yi<-geno_to_pheno(xi)
    y<-rbind(y,data.frame(g=i,r=yi[["r"]],wt=yi[["wt"]]))
  }
  return(y)
}

#sim_help_calc_allele_freq
#x ... output from sim_data
#y ... allele of interest (w,d,r1,r2)
#return(y) ... data frame with freq for y
sim_help_calc_allele_freq<-function(x,y){
  z<-apply(x,1,function(k) sum(k[y==1])/2)
  z2<-apply(x,1,function(k) sum(k[y==2]))
  z3<-z+z2
  return(z3)
  }

#sim_calc_allele_freq
#x ... output form sim_data
#return(z) ... data frame with drive, wt, resistance freq
sim_calc_allele_freq<-function(x){
  y<-x[,1:10]
  y<-y*2
  id<-str_count(colnames(y),"d")
  d_freq<-sim_help_calc_allele_freq(x=y,y=id)
  id<-str_count(colnames(y),"w")
  w_freq<-sim_help_calc_allele_freq(x=y,y=id)
  id<-str_count(colnames(y),"r1")
  r1_freq<-sim_help_calc_allele_freq(x = y,y=id)
  id<-str_count(colnames(y),"r2")
  r2_freq<-sim_help_calc_allele_freq(x=y,y=id)
  z<-data.frame(g=c(1:nrow(x)),w=w_freq,d=d_freq,r1=r1_freq,r2=r2_freq)
  return(z)
}

#sim_calc_genotype_freq
#x ... output from sim_data
#return(z) ... matrix with genotype freqs for drive locus
sim_calc_genotype_freq<-function(x){
  y<-x[,1:10]
  y<-y*2
  ww<-apply(y[,grep("ww",colnames(y))],1,function(k) sum(k))
  wd<-apply(y[,grep("wd",colnames(y))],1,function(k) sum(k))
  wr1<-apply(y[,grep("wr1",colnames(y))],1,function(k) sum(k))
  wr2<-apply(y[,grep("wr2",colnames(y))],1,function(k) sum(k))
  dd<-apply(y[,grep("dd",colnames(y))],1,function(k) sum(k))
  dr1<-apply(y[,grep("dr1",colnames(y))],1,function(k) sum(k))
  dr2<-apply(y[,grep("dr2",colnames(y))],1,function(k) sum(k))
  r1r1<-apply(y[,grep("r1r1",colnames(y))],1,function(k) sum(k))
  r1r2<-apply(y[,grep("r1r2",colnames(y))],1,function(k) sum(k))
  r2r2<-apply(y[,grep("r2r2",colnames(y))],1,function(k) sum(k))

  z<-cbind(ww,wd,wr1,wr2,dd,dr1,dr2,r1r1,r1r2,r2r2)
  colnames(z)<-c("ww","wd","wr1","wr2","dd","dr1","dr2","r1r1","r1r2","r2r2")
  return(z)
}


#sim_plot
#x ... output from sim_geno_to_pheno (data.frame)
#return(g) ... returns frequency plot of drive carriers
sim_plot<-function(x){
  g<-ggplot(data=x,aes(x=g,y=r))+geom_line(col="red")+geom_point(col="red")
  g<-g+theme_minimal()+coord_cartesian(ylim = c(0,1))
  g<-g+scale_y_continuous("f(r)",breaks = c(0,0.25,0.5,0.75,1))
  g<-g+xlab("Generation")+theme(text=element_text(size = 15))
  return(g)
}


#sim_plot_alleles
#x ... output from sim_calc_allele_freq
#return(g) ... allele frequency plot 
sim_plot_alleles<-function(x){
  y<-reshape2::melt(x,id.vars="g")
  colnames(y)<-c("g","allele","freq")
  g<-ggplot(data=y,aes(x=g,y=freq,col=allele))
  g<-g+geom_line()+geom_point()
  g<-g+theme_minimal()
  g<-g+theme(text=element_text(size=15))
  g<-g+scale_y_continuous("Allele Frequencies",breaks = seq(0,1,0.2))
  g<-g+scale_x_continuous("Generation",breaks = seq(0,40,5))
  g<-g+scale_color_manual("",values = c("w"="blue","d"="red","r1"="purple","r2"="orchid"))
  g<-g+theme(legend.position = "top")
  return(g)
}

#ML Functions----

#calc_rho
#fp ... predicted frequencies (r+wt)
#fo ... observed frequencies (r+wt)
#fn ... effective population size
#return(r) ...  logL calculation
calc_rho<-function(fp,fo,fn){
  fp<-fp/sum(fp)
  fo<-fo/sum(fo)
  fo.count<-fo*fn #observed counts 
  log.a<-log(fn)-log(1-(fp[1]^fn+fp[2]^fn)/2)
  r<-log.a+lgamma(fn+1)-lgamma(fo.count[1]+1)-lgamma(fo.count[2]+1)
  if(fp[1]>0) r<-r+fo.count[1]*log(fp[1])
  if(fp[2]>0) r<-r+fo.count[2]*log(fp[2])
  return(r)
  
}


#logL
#lf ... parameter to optimize
#lN ... list of cages to analyze (phenotype counts)
#ld ... gene drive fitness mode (mate choice, fecundity, viability)
#lmpath ... path to transition matrix
#lg ... genotypes
#ldr ... drive conversion rate
#ldrr ... drive resistance formation rate
#ler ... embryo resistance formation rate
#lh ... ratio of heterozygotes in generation 0 
#lr1 ... r1 resistance allele formation rate 
#lmode ... mode of setting parameters 
#return(l) ... log-Likelihood of data in lN 

logL<-function(lf,lN,ld,lmpath,lg,ldr,ldrr,ler,lh,lr1,lmode){
  
  #neutral set up:
  l.fdm<-1 #drive, mate-choice
  l.fdf<-1 #drive, fecundity
  l.fdv<-1 #drive, viability
  l.ne<-1 #effective population size 
  l.r1<-lr1 #r1 resistance allele formation rate 
  
  
  l.m<-init_m(lmpath,lg)
  if(lmode==0)  l.ne<-lf["ne"]
  
  if(lmode==1){ #fitness cost: drive; mate choice = fecundity
    l.ne<-lf["ne"]
    l.fdm<-lf["fdm"]
    l.fdf<-lf["fdm"]
  }
  
  if(lmode==3){ #fitness cost: drive; viability
    l.ne<-lf["ne"]
    l.fdv<-lf["fdv"]
  }
  
  if(lmode==7){ #fitness cost: drive; mate choice = fecundity; viability
    l.ne<-lf["ne"]
    l.fdm<-lf["fdm"]
    l.fdf<-lf["fdm"]
    l.fdv<-lf["fdv"]
  }
  
  if(lmode==10) { #fitness cost: drive; fecundity
    l.ne<-lf["ne"]
    l.fdf<-lf["fdf"]
  }
  
  if(lmode==11) { #estimate r1 formation rate for fixed direct fitness costs, viability
    l.r1<-lf["fr1"]
    l.fdv<-0.7969326
    l.ne<-round(240.9729031)
  }
  
  if(lmode==12){ #estimate r1 formation rate for fixed leaky somatic expression
    l.r1<-lf["fr1"]
    l.fdf<-0.4250285
    l.ne<-round(228.7897653)
  }
  
  l.params<-c(l.ne,l.fdm,l.fdf,l.fdv)
  names(l.params)<-c("ne","fdm","fdf","fdv")
  #print(l.r1)
  l.mmod<-mod_transition_m(m = l.m,g = lg,dr = ldr, drr = ldrr,er = ler,r1 = l.r1) #modified transition matrix
  #print(l.params)
  
  l<-0 #start logL calculation
  l.n.sets<-length(lN)
  for(j in c(1:l.n.sets)){
    l.Nj<-lN[[j]] #store counts
    l.Nfreq<-data.frame(r=l.Nj$r,wt=l.Nj$wt) #normalize counts
    l.Nfreq$n<-l.Nfreq$r+l.Nfreq$wt
    l.Nfreq$r<-l.Nfreq$r/l.Nfreq$n
    l.Nfreq$wt<-l.Nfreq$wt/l.Nfreq$n
    l.Nfreq$n<-NULL
    l.obs.geno<-start_genotype(x = l.Nfreq[1,],h = lh,fg = lg) #initialize
    #print(l.obs.geno)
    for(i in c(2:nrow(l.Nj))){
      l.obs.pheno<-unlist(l.Nfreq[i,]) #observed phenotype
      #print(l.obs.pheno)
      l.pred.geno<-propagate_genotypes(p = l.obs.geno,m = l.mmod,g = lg,fd = ld,f = l.params) #expected genotype
      #print(l.pred.geno)
      l.pred.pheno<-geno_to_pheno(l.pred.geno) #expected phenotype|expected genotype
      #print(l.pred.pheno)
      l<-l+calc_rho(fp = l.pred.pheno,fo = l.obs.pheno,fn = l.ne) #logL one transition step
      l.obs.geno<-pheno_to_geno(p = l.obs.pheno,g = l.pred.geno) #offspring in j  = parents in j+1
    }
  }
  return(l)
}

#set_limits
#fmode ... mode of selection applied 
#fd.lim ... step to upper drive parameter limit (optional)

set_limits<-function(fmode,fd.lim=1){
  
  if(fmode==0){
    y<-100
    names(y)<-c("ne")
    y.up<-y+49900
    y.down<-y-75
  }
  
  if(fmode==1){
    y<-c(100,1)
    names(y)<-c("ne","fdm")
    y.up<-y+c(49900,fd.lim)
    y.down<-y-c(75,0.99)
  }
  
  if(fmode==3){
    y<-c(100,1)
    names(y)<-c("ne","fdv")
    y.up<-y+c(49900,fd.lim)
    y.down<-y-c(75,0.99)
  }

  
  if(fmode==7){
    y<-c(100,1,1)
    names(y)<-c("ne","fdm","fdv")
    y.up<-y+c(49900,fd.lim,fd.lim)
    y.down<-y-c(75,0.99,0.99)
  }
  
  if(fmode==10){
    y<-c(100,1)
    names(y)<-c("ne","fdf")
    y.up<-y+c(49900,fd.lim)
    y.down<-y-c(75,0.99)
  }
  
  if(fmode%in%c(11,12)){
    y<-c(0.005)
    names(y)<-c("fr1")
    y.up<-1
    y.down<-0
  }
  
return(list(y.down,y,y.up))
}

#calculate_AICc
#flnL ... log likelihood
#p ... number of params
#n ... number of data points
#return(aicc)... corrected AIC for small sample size
calculate_AICc<-function(flnL,p,n){
  aic=2*p-2*flnL+c(2*p*p+2*p)/c(n-p-1)
  return(aic)
}


#calculate_CI
#fml ... ML result (optim(logL))
#fN  ... observed data 
#fdist ... min. dist for CI search
#fa ... type 1 error
#fd ... gene drive fitness mode (mate choice, fecundity, viability)
#fg ... genotypes
#fmpath ... path to transition matrix
#fdr ... drive conversion rate
#fdrr ... drive resistance formation rate
#fer ... embryo resistance formation rate
#fr1 ... r1 allele resistance formation rate
#fh ... ratio of heterozygotes in generation 0 
#fmode ... mode of setting parameters 
#return(y) ... data frame with CI based on fa type 1 error rate 

calculate_CI<-function(fml,fN,fdist,fa,fd,fg,fmpath,fdr,fdrr,fer,fr1,fh,fmode){
  fdist_i<-fdist
  fsig<-qchisq(1-fa,1) #calculate Chi^2 stat
  fLmax<-fml$value 
  y<-data.frame(par=fml$par)
  y$low_ci<-c(-1)
  y$up_ci<-c(-1)
  
  for(i in c(1:length(fml$par))){
    if(fmode%in%c(11,12)==F) fdist<-ifelse(i>1,fdist_i,1) #ne: other step limit
    #keep step limit for r1 formation rates
    fpara_dynamic<-fml$par #init
    fpara_test<-fpara_dynamic[i]
    fstep<-fpara_test*0.99-fdist #initialize step size
    frun<-T
    if (ml$par==0) frun<-F #avoid endless loop, if estimate==0
    while(frun){
      fpara_test<-max(fdist,fpara_test-fstep)
      if(fpara_test==fdist) frun<-F #avoid endless loop
      fpara_dynamic[i]<-fpara_test #shift parameter 
      fl<-logL(lf = fpara_dynamic,lN = fN,ld = fd,lmpath = fmpath,lg = fg,ldr = fdr,ldrr = fdrr,ler = fer,lh = fh,lr1=fr1,lmode = fmode)
      print(paste(fpara_test,fl,sep=":"))
      fstat<-2*(fLmax-fl) #calculate Chi^2
      if(fstat>fsig){ #if significantly different
        if(2*fstep<fdist) {
          frun<-F
          y$low_ci[i]<-fpara_test
        } else {
          fpara_test<-fpara_test+fstep #restore parameter
          fstep<-fstep/2 #lower step size
        }
      }
    }
    
    fpara_dynamic<-fml$par #init upstream 
    fpara_test<-fpara_dynamic[i]
    fstep<-fpara_test*0.99 #initialize step size
    if(fpara_test==0) fstep<-0.001 #initialize step size if estimate==0
    frun<-T
    
    while(frun){
      fpara_test<-fpara_test+fstep #shift parameter
      if((fpara_test>=50000 & i==1) | (fpara_test>=5 & i>1)) run<-F #avoid endless loop 
      fpara_dynamic[i]<-fpara_test
      fl<-logL(lf = fpara_dynamic,lN = fN,ld = fd,lmpath = fmpath,lg = fg,ldr = fdr,ldrr = fdrr,ler = fer,lh = fh,lr1=fr1,lmode = fmode)
      print(paste(fpara_test,fl,sep=":"))
      fstat<-2*(fLmax-fl) #calculate Chi^2
      if(fstat>fsig){ #if significantly different
        if(2*fstep<fdist) {
          frun<-F
          y$up_ci[i]<-fpara_test
        } else {
          fpara_test<-fpara_test-fstep #restore parameter
          fstep<-fstep/2 #lower step size
        }
      }
    }
  }
  return(y)
}
