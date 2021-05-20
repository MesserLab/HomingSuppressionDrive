#!/bin/RScript
#AML
#HSD
#ML framework
#v1 --> v2: include fitness scenario of leaky somatic expression 

#Packages----
library(stringr)
library(ggplot2)

#Global Variables----
M.path<-"../data/HSD_AML_20210201_data_trans-matrix-FILLED.csv"

#w=wildtype, drive locus
#d=drive, drive locus
#r=resistance, drive locus
#u=uncut, off-target locus
#c=cut, off-target locus
G<-c("wwuu","wwuc","wwcc","wduu","wduc","wdcc","wruu","wruc","wrcc",
     "rduu","rduc","rdcc","dduu","dduc","ddcc","rruu","rruc","rrcc")

#Drive Functions----

#init_m
#x ... path to raw transition matrix (counts)
#y ... names of count columns
#return(z) ... expected offspring frequencies
init_m<-function(x,y){
  z<-read.csv(x,header = T)
  z<-cbind(z[,1:2],z[,c(3:20)]/rowSums(z[,c(3:20)]))
  colnames(z)<-c("m","f",y)
  return(z)
}

#init_genotype
#x ... names of genotypes
#return(y) ... empty named genotype vector
init_genotype<-function(x){
  y<-rep(0,18)
  names(y)<-x
  return(y)
}


#germ_cut_d 
#v ... original offspring frequencies (length = 18)
#x ... cut rate u -> c in germline; on/two drive alleles act the same 
#return(w) ... adapted offspring frequencies
germ_cut_d<-function(v,x){
  w<-v
  w["wduu"]<-v["wduu"]*(1-x)^2
  w["dduu"]<-v["dduu"]*(1-x)^2
  w["rduu"]<-v["rduu"]*(1-x)^2
  w["wduc"]<-v["wduc"]*(1-x)+v["wduu"]*2*x*(1-x)
  w["dduc"]<-v["dduc"]*(1-x)+v["dduu"]*2*x*(1-x)
  w["rduc"]<-v["rduc"]*(1-x)+v["rduu"]*2*x*(1-x)
  w["wdcc"]<-v["wdcc"]+v["wduc"]*x+v["wduu"]*x^2
  w["ddcc"]<-v["ddcc"]+v["dduc"]*x+v["dduu"]*x^2
  w["rdcc"]<-v["rdcc"]+v["rduc"]*x+v["rduu"]*x^2
  return(w)
}

#embryo_cut_d
#v ... original offspring frequencies (length = 18)
#e .. cut rate u -> c if mom had at least one drive allele 
#return(w) ... adapted offspring frequencies 
embryo_cut_d<-function(v,e){
  w<-v
  w["wwuu"]<-v["wwuu"]*(1-e)^2
  w["wduu"]<-v["wduu"]*(1-e)^2
  w["wruu"]<-v["wruu"]*(1-e)^2
  w["rduu"]<-v["rduu"]*(1-e)^2
  w["dduu"]<-v["dduu"]*(1-e)^2
  w["rruu"]<-v["rruu"]*(1-e)^2
  w["wwuc"]<-v["wwuu"]*2*e*(1-e)+v["wwuc"]*(1-e)
  w["wduc"]<-v["wduu"]*2*e*(1-e)+v["wduc"]*(1-e)
  w["wruc"]<-v["wruu"]*2*e*(1-e)+v["wruc"]*(1-e)
  w["rduc"]<-v["rduu"]*2*e*(1-e)+v["rduc"]*(1-e)
  w["dduc"]<-v["dduu"]*2*e*(1-e)+v["dduc"]*(1-e)
  w["rruc"]<-v["rruu"]*2*e*(1-e)+v["rruc"]*(1-e)
  w["wwcc"]<-v["wwcc"]+v["wwuc"]*e+v["wwuu"]*e^2
  w["wdcc"]<-v["wdcc"]+v["wduc"]*e+v["wduu"]*e^2
  w["wrcc"]<-v["wrcc"]+v["wruc"]*e+v["wruu"]*e^2
  w["rdcc"]<-v["rdcc"]+v["rduc"]*e+v["rduu"]*e^2
  w["ddcc"]<-v["ddcc"]+v["dduc"]*e+v["dduu"]*e^2
  w["rrcc"]<-v["rrcc"]+v["rruc"]*e+v["rruu"]*e^2
  return(w)
}

#germ_drive
#v ... original offspring frequencies (length=18)
#x ... gene drive conversion rate 
#y ... germ line resistance formation rate 
germ_drive<-function(v,x,y){
  w<-v
  w["wduu"]<-v["wduu"]*(1-x-y)
  w["wduc"]<-v["wduc"]*(1-x-y)
  w["wdcc"]<-v["wdcc"]*(1-x-y)
  w["dduu"]<-v["dduu"]+v["wduu"]*x
  w["dduc"]<-v["dduc"]+v["wduc"]*x
  w["ddcc"]<-v["ddcc"]+v["wdcc"]*x
  w["rduu"]<-v["rduu"]+v["wduu"]*y
  w["rduc"]<-v["rduc"]+v["wduc"]*y
  w["rdcc"]<-v["rdcc"]+v["wdcc"]*y
  return(w)
}

#embryo_drive
#v ... original offspring frequencies (length=18)
#e ... embryo resistance formation rate given mom has at least one drive allele
embryo_drive<-function(v,e){
  w<-v
  w["wwuu"]<-v["wwuu"]*(1-e)^2
  w["wwuc"]<-v["wwuc"]*(1-e)^2
  w["wwcc"]<-v["wwcc"]*(1-e)^2
  w["wduu"]<-v["wduu"]*(1-e)
  w["wduc"]<-v["wduc"]*(1-e)
  w["wdcc"]<-v["wdcc"]*(1-e)
  w["wruu"]<-v["wwuu"]*2*e*(1-e)+v["wruu"]*(1-e)
  w["wruc"]<-v["wwuc"]*2*e*(1-e)+v["wruc"]*(1-e)
  w["wrcc"]<-v["wwcc"]*2*e*(1-e)+v["wrcc"]*(1-e)
  w["rduu"]<-v["rduu"]+v["wduu"]*e
  w["rduc"]<-v["rduc"]+v["wduc"]*e
  w["rdcc"]<-v["rdcc"]+v["wdcc"]*e
  w["rruu"]<-v["rruu"]+v["wruu"]*e+v["wwuu"]*e^2
  w["rruc"]<-v["rruc"]+v["wruc"]*e+v["wwuc"]*e^2
  w["rrcc"]<-v["rrcc"]+v["wrcc"]*e+v["wwcc"]*e^2
  return(w)
}

#mod_transition_m
#m ... transition matrix 
#x ... germline cut rate (off-target)
#e ... embryo cut rate (off-target)
#g ... names of genotypes
#dr ... drive conversion rate
#drr ... drive resistance formation rate
#er ... embryo resistance formation rate 
#return(y) ... transition matrix with adapted offspring frequencies after germline + embryo cutting
mod_transition_m<-function(m,x,e,g,dr,drr,er){
  y<-m
  #Step 1: Germline cuts 
  for(i in c(1:nrow(y))){
    fy_iter<-y[i,] #store raw  transition probabilities + cross
    fmale_genotype<-init_genotype(x = g) #empty male genotype vector
    ffemale_genotype<-init_genotype(x = g) #empty female genotype vector
    fmale_genotype[fy_iter[["m"]]]<-1 #set male genotype
    ffemale_genotype[fy_iter[["f"]]]<-1 #set female genotype
    fmale_genotype<-germ_drive(fmale_genotype,x = dr ,y = drr) #drive conversion + resistance
    ffemale_genotype<-germ_drive(ffemale_genotype,x = dr, y = drr) #drive conversion + resistance
    fmale_genotype<-germ_cut_d(fmale_genotype,x)  #apply germline cuts in males
    ffemale_genotype<-germ_cut_d(ffemale_genotype,x) #apply germline cuts in females
    fnew_mates<-as.vector(outer(ffemale_genotype,fmale_genotype))  #all possible mating probs
    fnew_progeny<-fnew_mates*m[,(3:20)] #determine new weights
    y[i,(3:20)]<-colSums(fnew_progeny)/sum(fnew_progeny) #replace marginal probs for offspring
  }
  
  #Stemp 2: Embryo cuts
  fd_carrier<-grep("d",colnames(y))-2
  for(i in c(1:nrow(y))){
    if(y$f[i]%in%fd_carrier){
      y[i,3:20]<-embryo_drive(y[i,3:20],er)
      y[i,3:20]<-embryo_cut_d(y[i,3:20],e)
    }
  }
  return(y)
}

#determine_fitness_coefficients
#m... determines the fitness costs applied 
#fg ... names of genotypes
#a ... determine allele for which fitness costs apply
#p... fitness cost parameters
#return(x) ... fitness vector for genotypes (n=18); rel weights
determine_fitness_coefficients<-function(m,fg,a,p){
  x<-rep(1,18)
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
#return(x) ... fitness vector for genotypes (n=18); rel weights 
#any genotype being a combination of drive + resistance allele is sterile 
female_fertility<-function(fg){
  x<-rep(1,18)
  names(x)<-fg
  sterile<-sapply(c("rd","dd","rr"), function(k) grep(k,names(x)))
  sterile<-as.vector(sterile)
  x[sterile]<-0
  return(x)
}

#geno_to_pheno
#x ... named vector of genotypes (m+f,length=36)
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
#g ... named vector of genotypes (m+f, length=36)
#p ... observed phenotypes (r,wt)
#return(fgcorr) ... vector of genotype frequencies, corrected by observed phenotypes (m+f,length=36)
pheno_to_geno<-function(g,p){
  fpobs<-p/sum(p)
  fpraw<-geno_to_pheno(g) #marginal sums of genotypes giving the same phenotype 
  fgcorr<-numeric(length = 36)
  names(fgcorr)<-names(g)
  fd_id<-grep("d",names(fgcorr))
  fgcorr[fd_id]<-g[fd_id]/fpraw["r"]*fpobs["r"]   #standardize 
  fgcorr[-fd_id]<-g[-fd_id]/fpraw["wt"]*fpobs["wt"]
  return(fgcorr)
}

#start_genotype .... function to start a genotype 
#x ... counts (r|wt)
#fg ... genotype names
#t ... already cut rate
#h... ratio of heterozygous flies 
#return(y)... genotype frequencies m+f (length=36)
start_genotype<-function(x,fg,t,h){
  x<-unlist(x)
  names(x)<-c("r","wt")
  x<-x/sum(x)
  y<-numeric(length = 18)
  names(y)<-fg
  y["wwuu"]<-x["wt"]/2 #no cut in without drive
  y["wduu"]<-x["r"]*h*(1-t)^2/2 #cut in heterozygotes
  y["wduc"]<-x["r"]*h*2*t*(1-t)/2
  y["wdcc"]<-x["r"]*h*t^2/2
  y["dduu"]<-x["r"]*(1-h)*(1-t)^2/2 #cut in homozygotes
  y["dduc"]<-x["r"]*(1-h)*2*t*(1-t)/2
  y["ddcc"]<-x["r"]*(1-h)*t^2/2
  y<-rep(y,2)
  return(y)
}

#propagate_genotypes
#p ... parental genotypes (36 (m+f))
#m ... transition matrix (modified)
#g ... names of genotypes
#fd ... gen drive fitness cost mode (mate choice, fecundity, viability)
#fo ... off target fitness cost (mate choice, fecundity, viability)
#f... fitness parameters (named!)
#return=expected genotype frequencies next generation (36 (m+f))
propagate_genotypes<-function(p,m,g,fd,fo,f){
  fmale<-p[1:18]/sum(p[1:18]) #standardize male freqs
  ffemale<-p[19:36]/sum(p[19:36]) #standardize female freqs
  
  #mate choice selection 
  fmate_cas9<-determine_fitness_coefficients(m = fd[1],fg = g,a = "d",p = f[c("fdm")])
  fmate_ot<-determine_fitness_coefficients(m = fo[1],fg = g,a = "c",p = f[c("fom")])
  fmate<-fmate_cas9*fmate_ot
  fmale<-fmale*fmate/sum(fmale*fmate)
  fpairs<-as.vector(outer(ffemale,fmale))
  fmmate<-fpairs*m[,3:20]
  
  #fecundity selection
  ffec_cas9<-determine_fitness_coefficients(m = fd[2],fg = g,a = "d",p = f[c("fdf")])
  ffec_cas9<-rep(ffec_cas9,length(fmale))
  ffec_ot<-determine_fitness_coefficients (m = fo[2],fg = g,a = "c",p = f[c("fof")])
  ffec_ot<-rep(ffec_ot,length(fmale))
  ffec_sterile<-female_fertility(g) #sterility caused by drive
  ffec_sterile<-rep(ffec_sterile,length(fmale))
  ffec<-ffec_cas9*ffec_ot
  ffec<-ffec*ffec_sterile
  fmfec<-ffec*fmmate
  
  #viability selection
  fvia_cas9<-determine_fitness_coefficients(m = fd[3],fg = g,a = "d",p = f[c("fdv")])
  fvia_ot<-determine_fitness_coefficients(m=fo[3],fg = g,a = "c",p = f[c("fov")])
  fvia<-fvia_cas9*fvia_ot
  fmvia<-t(fvia*t(fmfec))
  
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
#fo ... off-target fitness cost
#e ... embryo cut rate
#x ... germline cut rate
#n ... number of generations to simulate
#fmpath ... path to transition matrix 
#ft ... already cut parameter
#fh ... ratio of heterozygous flies at gen. 0
#fdr ...gene drive conversion rate
#fdrr ... gene drive resistance rate
#fer ... embryo resistance rate 
sim_data<-function(p,g,fne,f,fd,fo,e,x,n,fmpath,ft,fh,fdr,fdrr,fer){
  r<-p/sum(p)
  r<-start_genotype(x = r,fg = g,t = ft,h = fh)
  y<-r
  m<-init_m(x = fmpath,y = g) #raw transition matrix
  fmmod<-mod_transition_m(m = m, x = x, e = e, g = g,dr = fdr,drr = fdrr ,er = fer ) #germline/embryo cut + drive conversion + resistance
  for(i in c(1:n)){ #for each generation 
    fp1<-propagate_genotypes(p = r,m = fmmod,g = g,fd = fd,fo = fo,f = f)
    fpm<-fp1[1:18]/sum(fp1[1:18])
    fpf<-fp1[19:36]/sum(fp1[19:36])
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
#y ... allele of interest (r,w,d)
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
  y<-x[,1:18]
  y<-y*2
  id<-str_count(colnames(y),"d")
  d_freq<-sim_help_calc_allele_freq(x=y,y=id)
  id<-str_count(colnames(y),"w")
  w_freq<-sim_help_calc_allele_freq(x=y,y=id)
  id<-str_count(colnames(y),"r")
  r_freq<-sim_help_calc_allele_freq(x = y,y=id)
  z<-data.frame(g=c(1:nrow(x)),w=w_freq,d=d_freq,r=r_freq)
  return(z)
}

#sim_calc_genotype_freq
#x ... output from sim_data
#return(z) ... matrix with genotype freqs for drive locus
sim_calc_genotype_freq<-function(x){
  y<-x[,1:18]
  y<-y*2
  ww<-apply(y[,grep("ww",colnames(y))],1,function(k) sum(k))
  wd<-apply(y[,grep("wd",colnames(y))],1,function(k) sum(k))
  dd<-apply(y[,grep("dd",colnames(y))],1,function(k) sum(k))
  wr<-apply(y[,grep("wr",colnames(y))],1,function(k) sum(k))
  dr<-apply(y[,grep("rd",colnames(y))],1,function(k) sum(k))
  rr<-apply(y[,grep("rr",colnames(y))],1,function(k) sum(k))
  z<-cbind(ww,wd,wr,dd,dr,rr)
  colnames(z)<-c("ww","wd","wr","dd","dr","rr")
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
  g<-g+scale_color_manual("",values = c("w"="blue","d"="red","r"="purple"))
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
#lo ... off-target fitness mode (mate choice, fecundity, viability)
#lt ... already cut rate (length = length lN)
#lmpath ... path to transition matrix
#lg ... genotypes
#lx ... germline cut rate
#le ... embryo cut-rate 
#ldr ... drive conversion rate
#ldrr ... drive resistance formation rate
#ler ... embryo resistance formation rate
#lh ... ratio of heterozygotes in generation 0 
#lmode ... mode of setting parameters 
#return(l) ... log-Likelihood of data in lN 

logL<-function(lf,lN,ld,lo,lt,lmpath,lg,lx,le,ldr,ldrr,ler,lh,lmode){
  
  #neutral set up:
  l.fdm<-1 #drive, mate-choice
  l.fdf<-1 #drive, fecundity
  l.fdv<-1 #drive, viability
  l.fom<-1 #off-target, mate-choice
  l.fof<-1 #off-target, fecundity
  l.fov<-1 #off-target, viability
  l.ne<-1 #effective population size 
  
  l.ne<-lf["ne"] #always estimate Ne
  
  l.m<-init_m(lmpath,lg)
  l.mmod<-mod_transition_m(m = l.m,x = lx,e = le,g = lg,dr = ldr, drr = ldrr,er = ler) #modified transition matrix
  
  if(lmode==1){ #fitness cost: drive; mate choice = fecundity
    l.fdm<-lf["fdm"]
    l.fdf<-lf["fdm"]
  }
  
  if(lmode==2){ #fitness cost: off-target; mate choice = fecundity
    l.fom<-lf["fom"]
    l.fof<-lf["fom"]
  }
  
  if(lmode==3){ #fitness cost: drive; viability
    l.fdv<-lf["fdv"]
  }
  
  if(lmode==4){ #fitness cost: off-target; viability
    l.fov<-lf["fov"]
  }
  
  if(lmode==5){ #fitness cost: drive + off-target; mate choice = fecundity
    l.fdm<-lf["fdm"]
    l.fdf<-lf["fdm"]
    l.fom<-lf["fom"]
    l.fof<-lf["fom"]
  }
  
  if(lmode==6){ #fitness cost: drive + off-target; viability
    l.fdv<-lf["fdv"]
    l.fov<-lf["fov"]
  }
  
  if(lmode==7){ #fitness cost: drive; mate choice = fecundity; viability
    l.fdm<-lf["fdm"]
    l.fdf<-lf["fdm"]
    l.fdv<-lf["fdv"]
  }
  
  if(lmode==8){ #fitness cost: off-target; mate choice = fecundity; viability
    l.fom<-lf["fom"]
    l.fof<-lf["fom"]
    l.fov<-lf["fov"]
  }
  
  if(lmode==9) { #fitness cost: drive + off-target; 
    #mate choice = fecundity; viability
    l.fdm<-lf["fdm"]
    l.fdf<-lf["fdm"]
    l.fdv<-lf["fdv"]
    l.fom<-lf["fom"]
    l.fof<-lf["fom"]
    l.fov<-lf["fov"]
  }
  
  if(lmode==10) { #fitness cost: drive; fecundity
    l.fdf<-lf["fdf"]
  }
  l.params<-c(l.ne,l.fdm,l.fdf,l.fdv,l.fom,l.fof,l.fov)
  names(l.params)<-c("ne","fdm","fdf","fdv","fom","fof","fov")
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
    l.obs.geno<-start_genotype(x = l.Nfreq[1,],h = lh,fg = lg,t = lt[j]) #initialize
    #print(l.obs.geno)
    for(i in c(2:nrow(l.Nj))){
      l.obs.pheno<-unlist(l.Nfreq[i,]) #observed phenotype
      #print(l.obs.pheno)
      l.pred.geno<-propagate_genotypes(p = l.obs.geno,m = l.mmod,g = lg,fd = ld,fo = lo,f = l.params) #expected genotype
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
#fo.lim ... step upper off-target parameter limit (optional)

set_limits<-function(fmode,fd.lim=1,fo.lim=1){
  
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
  
  if(fmode==2){
    y<-c(100,1)
    names(y)<-c("ne","fom")
    y.up<-y+c(49900,fo.lim)
    y.down<-y-c(75,0.99)
  }
  
  if(fmode==3){
    y<-c(100,1)
    names(y)<-c("ne","fdv")
    y.up<-y+c(49900,fd.lim)
    y.down<-y-c(75,0.99)
  }
  
  if(fmode==4){
    y<-c(100,1)
    names(y)<-c("ne","fov")
    y.up<-y+c(49900,fo.lim)
    y.down<-y-c(75,0.99)
  }
  
  if(fmode==5){
    y<-c(100,1,1)
    names(y)<-c("ne","fdm","fom")
    y.up<-y+c(49900,fd.lim,fo.lim)
    y.down<-y-c(75,0.99,0.99)
  }
  
  if(fmode==6){
    y<-c(100,1,1)
    names(y)<-c("ne","fdv","fov")
    y.up<-y+c(49900,fd.lim,fo.lim)
    y.down<-y-c(75,0.99,0.99)
  }
  
  if(fmode==7){
    y<-c(100,1,1)
    names(y)<-c("ne","fdm","fdv")
    y.up<-y+c(49900,fd.lim,fd.lim)
    y.down<-y-c(75,0.99,0.99)
  }
  
  if(fmode==8){
    y<-c(100,1,1)
    names(y)<-c("ne","fom","fov")
    y.up<-y+c(49900,fo.lim,fo.lim)
    y.down<-y-c(75,0.99,0.99)
  }
  
  if(fmode==9){
    y<-c(100,1,1,1,1)
    names(y)<-c("ne","fdm","fdv","fom","fov")
    y.up<-y+c(49900,fd.lim,fd.lim,fo.lim,fo.lim)
    y.down<-y-c(75,rep(0.99,4))
  }
  
  if(fmode==10){
    y<-c(100,1)
    names(y)<-c("ne","fdf")
    y.up<-y+c(49900,fd.lim)
    y.down<-y-c(75,0.99)
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
#fo ... off-target fitness mode (mate choice, fecundity, viability)
#fg ... genotypes
#ft .... already cut rate
#fmpath ... path to transition matrix
#fx ... germline cut rate
#fe ... embryo cut-rate 
#fdr ... drive conversion rate
#fdrr ... drive resistance formation rate
#fer ... embryo resistance formation rate
#fh ... ratio of heterozygotes in generation 0 
#fmode ... mode of setting parameters 
#return(y) ... data frame with CI based on fa type 1 error rate 

calculate_CI<-function(fml,fN,fdist,fa,fd,fo,fg,ft,fmpath,fx,fe,fdr,fdrr,fer,fh,fmode){
  fdist_i<-fdist
  fsig<-qchisq(1-fa,1) #calculate Chi^2 stat
  fLmax<-fml$value 
  y<-data.frame(par=fml$par)
  y$low_ci<-c(-1)
  y$up_ci<-c(-1)
  
  for(i in c(1:length(fml$par))){
    fdist<-ifelse(i>1,fdist_i,1) #ne: other step limit
    fpara_dynamic<-fml$par #init
    fpara_test<-fpara_dynamic[i]
    fstep<-fpara_test*0.99-fdist #initialize step stize
    frun<-T
    while(frun){
      fpara_test<-max(fdist,fpara_test-fstep)
      if(fpara_test==fdist) frun<-F #avoid endless loop
      fpara_dynamic[i]<-fpara_test #shift parameter 
      fl<-logL(lf = fpara_dynamic,lN = fN,ld = fd,lo = fo,lt=ft,lmpath = fmpath,lg = fg,lx = fx,le = fe,ldr = fdr,ldrr = fdrr,ler = fer,lh = fh,lmode = fmode)
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
    frun<-T
    
    while(frun){
      fpara_test<-fpara_test+fstep #shift parameter
      if((fpara_test>=50000 & i==1) | (fpara_test>=5 & i>1)) run<-F #avoid endless loop 
      fpara_dynamic[i]<-fpara_test
      fl<-fl<-logL(lf = fpara_dynamic,lN = fN,ld = fd,lo = fo,lt=ft,lmpath = fmpath,lg = fg,lx = fx,le = fe,ldr = fdr,ldrr = fdrr,ler = fer,lh = fh,lmode = fmode)
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
