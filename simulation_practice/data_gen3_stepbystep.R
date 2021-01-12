
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

samp_ID[,k] <- sel_samp

dsim$study<-ifelse((dsim$SchoolNo %in% sel_samp),1,0)

D <- subset(dsim,dsim$study==1)[,"SchoolNo"]

sel_samp2 <- sample(D,size=round(0.5*num.samp),replace=F)

dsim$trt <- ifelse((dsim$SchoolNo %in% sel_samp2),1,0)

#3) Get covariates for the simulated data.

covars_sim <- as.matrix(scale(data.frame(dsim[c(covars)])))

d1 <- data.frame(cbind(dsim$SchoolNo,data.frame(scale(dsim[c(covars)])),dsim$study,dsim$trt))
