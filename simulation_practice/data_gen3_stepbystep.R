
# GENERALIZATIONS FROM SMALL SAMPLES

#0) Load packages.

library(pacman)
p_load("MatchIt","sae","tmle","plyr","rstan","gam","dbarts",
       "BayesTree","generalize")
library(generalize)
rstan::rstan_options(auto_write = TRUE)
rstan::rstan_options(javascript = FALSE)
options(mc.cores = parallel::detectCores())
devAskNewPage(ask = FALSE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')



#1) Input data.

dsim <- read.csv("raw/ind_sim1.csv", sep = ",", header = TRUE)

size_gen <- 1000

num.strata <- 5

num.samp <- round(0.03*nrow(dsim))

chains <- 4

iter_hb <- 2000

#Covariate Set

covars <- c("z_s09_ELA", "z_s09_math", "Attendance_0809", "FTE08",
            "MEMBER08","PUPTCH08","county.pop","TITLEI08",
            "STITLI08","male","white","sped","frpl",
            "lep")

ns <- c(3,4,5)

est_smallsamp <- c("sub", "odds_ipw", "overlap_ipw", "bart", "drwlsl",
                   "drwls2","drwls3","tmle","eblup","hb")



### INITIALIZATION OF OUTPUT ARRAYS ###

# ALLOCATION OF SCHOOLS #

allocation_AAA<-array(NA,dim=c(num.strata,2,size_gen,length(ns)))
dimnames(allocation_AAA)<-list(paste("Stratum",1:num.strata,sep=""),c("Control","Trt"),paste("sim",1:size_gen,sep=""),paste("SS",ns,sep=""))

# MAX NO. OF STRATA #

maxstrata_AAA<-array(NA,dim=c(size_gen,1,length(ns)))
dimnames(maxstrata_AAA)<-list(NULL,"Max",paste("SS",ns,sep=""))

mse_est<-array(NA,dim=c(length(est_smallsamp),1,size_gen,length(ns)))
dimnames(mse_est)<-list(c(est_smallsamp),NULL,paste("sim",1:size_gen,sep=""),paste(as.character(ns),"Strata",sep=""))

bias_est<-array(NA,dim=c(length(est_smallsamp),1,size_gen,length(ns)))
dimnames(bias_est)<-list(c(est_smallsamp),NULL,paste("sim",1:size_gen,sep=""),paste(as.character(ns),"Strata",sep=""))

estimators<-array(NA,dim=c(length(est_smallsamp),1,size_gen,length(ns)))
dimnames(estimators)<-list(c(est_smallsamp),NULL,paste("sim",1:size_gen,sep=""),paste(as.character(ns),"Strata",sep=""))

samp_ID<-matrix(NA,(num.samp),size_gen)



### set seed

## instead of for k in each, just do for one number, 100


#2) Use a PPS sampling structure to select sample schools.
k=100
sel_samp <- sample(dsim$SchoolNo,
                   size = num.samp,
                   prob = dsim$prob,
                   replace = FALSE)
sel_samp

samp_ID[,100] <- sel_samp

dsim$study<-ifelse((dsim$SchoolNo %in% sel_samp),1,0)

D <- subset(dsim,dsim$study==1)[,"SchoolNo"]

sel_samp2 <- sample(D,size=round(0.5*num.samp),replace=F)
sel_samp2

dsim$trt <- ifelse((dsim$SchoolNo %in% sel_samp2),1,0)

#3) Get covariates for the simulated data.

covars_sim <- as.matrix(scale(data.frame(dsim[c(covars)])))

d1 <- data.frame(cbind(dsim$SchoolNo,data.frame(scale(dsim[c(covars)])),dsim$study,dsim$trt))



pos <- match(4,ns) 

mod <- matchit(dsim.study~covars_sim, data = data.frame(d1), method = "subclass", subclass = 4, estimand = "ATE", distance = "glm")

matched<- match.data(mod)

dsim_new<-data.frame(cbind(dsim,matched[c("subclass","distance")]))

sample.dsim<-subset(dsim_new,study==1)



d_use <- dsim_new[c("SchoolNo",covars,"study","subclass","distance","trt")]

beta <- as.matrix(rep(1,10))


for (f in 1:nrow(d_use)){
  d_use$v[f] <- rnorm(1)
  d_use$yt[f] <- 1 + sum(beta[2:5]*d_use[f,c("z_s09_ELA","z_s09_math","frpl","lep")])+ beta[6]*d_use$v[f]+beta[7]*1 + beta[8]*1*d_use$z_s09_ELA[f] + beta[9]*1*d_use$z_s09_math[f]+beta[10]*((d_use$frpl[f])^2)
  d_use$yc[f] <- 1 + sum(beta[2:5]*d_use[f,c("z_s09_ELA","z_s09_math","frpl","lep")])+ beta[6]*d_use$v[f]
}

samp_use<-data.frame(subset(d_use,study==1))

for (g in 1:nrow(samp_use)){
  samp_use$y[g] <-ifelse(samp_use$trt[g]==1,samp_use$yt[g], samp_use$yc[g])
}

d_non <- data.frame(subset(d_use,study==0)) # Non-sampled population

d_non$y <- d_non$trt <- NA


d_combined <- data.frame(rbind(samp_use,d_non)) # Stacked data frame for the non-sampled and sampled units

truth <- mean(d_use$yt) - mean(d_use$yc) # True PATE
truth



# SUBCLASSIFICATION ESTIMATOR

for(ss in 1:nrow(samp_use)){
  for(tt in 1:4){
    samp_use[ss,paste("sub",tt,sep="")] <- ifelse(samp_use$subclass[ss]==tt,1,0)
    samp_use[ss,paste(paste("sub",tt,sep=""),"T",sep="")] <- samp_use[ss,paste("sub",tt,sep="")]*samp_use[ss,"trt"]
  } # Creates indicators for each subclass and subclass x trt interaction in the sample
}

for(ss in 1:nrow(d_non)){
  for(tt in 1:4){
    d_non[ss,paste("sub",tt,sep="")] <- ifelse(d_non$subclass[ss]==tt,1,0)
    d_non[ss,paste(paste("sub",tt,sep=""),"T",sep="")] <- d_non[ss,paste("sub",tt,sep="")]*d_non[ss,"trt"]
  } # Creates indicators for each subclass and subclass x trt interaction in the non-sampled pop
}

mod_sub <- lm(as.formula(paste(paste(paste("y~",paste(paste(c(paste("sub",c(1:4),sep=""),paste(paste("sub",c(1:4),sep=""),"T",sep="")),collapse="+"),sep=""))),"-1",sep="")),data=samp_use)

estimators["sub",1,k,pos] <- mean(mod_sub$coefficients[(4+1):(2*4)])
estimators



# WEIGHTING BY ODDS ESTIMATOR

gen_object <- generalize(outcome="y",treatment = "trt",trial = "study",selection_covariates = covars,data=d_combined,is_data_disjoint = T)

estimators["odds_ipw",1,k,pos] <- gen_object$TATE$estimate

# WEIGHTING BY OVERLAP ESTIMATOR

estimators["overlap_ipw",1,k,pos] <- (sum(samp_use$y[samp_use$trt==1]*(samp_use$distance[samp_use$trt==1]))/(sum(samp_use$trt[samp_use$trt==1]*samp_use$distance[samp_use$trt==1]))) - 
  (sum(samp_use$y[samp_use$trt==0]*(samp_use$distance[samp_use$trt==0]))/(sum(1-samp_use$trt[samp_use$trt==0]*samp_use$distance[samp_use$trt==0])))


estimators[1]
