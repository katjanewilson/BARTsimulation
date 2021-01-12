
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

dsim <- read.csv("ind_sim1.csv",sep=",",header=TRUE)

size_gen <- 1000

num.strata <- 5

num.samp <- round(0.03*nrow(dsim))

chains <- 4

iter_hb <- 2000

# Covariate Set

covars <- c("z_s09_ELA","z_s09_math","Attendance_0809","FTE08","MEMBER08","PUPTCH08","county.pop","TITLEI08","STITLI08","male","white","sped","frpl","lep")

ns <- c(3,4,5)

est_smallsamp <- c("sub","odds_ipw","overlap_ipw","bart","drwls1","drwls2","drwls3","tmle","eblup","hb")

# source(paste(paste("stan_start",sep=""),".R",sep=""))

# INITIATE STAN CODE #

# load(file=paste(paste("samp_use",(num.samp),sep=""),".Rdata",sep=""))
# 
# covars2 <- c(covars,paste("sub",c(1:num.strata),sep=""),paste(paste("sub",c(1:num.strata),sep=""),"T",sep=""))
# 
# data_stan = list(
#   N=nrow(samp_use),
#   num_strata=num.strata,
#   p=length(covars2),
#   subclass=as.numeric(samp_use$subclass),
#   y_out=samp_use$y
# )
# 
# for(cc in 1:data_stan$p){
#   data_stan[[paste("x",cc,sep="")]] <- samp_use[,covars2[cc]]
# }
# 
# stanmodel <- stan(file="sdl_5.stan", data=data_stan, chains=0)

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

##################################################################################################################################
set.seed(888)

for (k in 1:size_gen){
  # browser()
  
  #2) Use a PPS sampling structure to select sample schools.
  
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
  
  for (ff in (3:num.strata)){
    # browser()
    
    pos <- match(ff,ns) # Match the number of strata with the index in the arrays
    
    #4)  Run the propensity score model.
    
    mod <- matchit(dsim.study~covars_sim, data = data.frame(d1), method = "subclass", subclass = ff, estimand = "ATE", distance = "glm")
    
    matched<- match.data(mod)
    
    dsim_new<-data.frame(cbind(dsim,matched[c("subclass","distance")]))
    
    sample.dsim<-subset(dsim_new,study==1)
    
    #5) Double check the number of schools in each stratum. Make sure there are no strata with < 1 control or treatment schools.
    
    t2_samp <- as.matrix(table(factor(sample.dsim$subclass,levels=1:ff),sample.dsim$trt))
    
    colnames(t2_samp)<-c("Control","Trt")
    
    counter <- 1
    
    # The following code is to prevent a subclass from having zero control or zero treatment units
    
    while((any(t2_samp[,1] == 0) | any(t2_samp[,2] == 0))){ 
      
      sel_samp <- sample(dsim$SchoolNo,
                         size = num.samp,
                         prob = dsim$prob,
                         replace = FALSE)
      
      samp_ID[,k] <- sel_samp
      
      dsim$study<-ifelse((dsim$SchoolNo %in% sel_samp),1,0)
      
      ssamp2 <-subset(dsim,dsim$study==1)[,"SchoolNo"]
      
      ssamp2_trt <-sample(ssamp2,size=round(0.5*num.samp),replace=F)
      
      dsim$trt<-ifelse((dsim$SchoolNo %in% ssamp2_trt),1,0)
      
      covars_sim <- as.matrix(scale(data.frame(dsim[c(covars)]))) # Scaled covariates
      
      d1 <- data.frame(cbind(dsim$SchoolNo,data.frame(scale(dsim[c(covars)])),dsim$study,dsim$trt))
      
      mod2 <- matchit(dsim.study~covars_sim, data = data.frame(d1), method = "subclass", subclass = ff, estimand = "ATE", distance = "glm")
      
      matched<- match.data(mod2)
      
      dsim_new<-data.frame(cbind(dsim,matched[c("subclass","distance")]))
      
      sample.dsim<-subset(dsim_new,study==1)
      
      t2 <- as.matrix(table(factor(dsim_new$subclass,levels=1:ff),dsim_new$study))
      
      colnames(t2) <- c("Pop","Sample")
      
      t2_samp <- as.matrix(table(factor(sample.dsim$subclass,levels=1:ff),sample.dsim$trt))
      
      colnames(t2_samp) <- c("Control","Trt")
      
      counter <- counter + 1
      
      if(counter > 200){
        break
      }
    }
    
    #6) Collect the covariates and propensity score logits in a new data frame.
    
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
    
    ###########################################################################################################
    
    #7) ESTIMATORS
    
    d_combined <- data.frame(rbind(samp_use,d_non)) # Stacked data frame for the non-sampled and sampled units
    
    truth <- mean(d_use$yt) - mean(d_use$yc) # True PATE
    
    # SUBCLASSIFICATION ESTIMATOR
    
    for(ss in 1:nrow(samp_use)){
      for(tt in 1:ff){
        samp_use[ss,paste("sub",tt,sep="")] <- ifelse(samp_use$subclass[ss]==tt,1,0)
        samp_use[ss,paste(paste("sub",tt,sep=""),"T",sep="")] <- samp_use[ss,paste("sub",tt,sep="")]*samp_use[ss,"trt"]
      } # Creates indicators for each subclass and subclass x trt interaction in the sample
    }
    
    for(ss in 1:nrow(d_non)){
      for(tt in 1:ff){
        d_non[ss,paste("sub",tt,sep="")] <- ifelse(d_non$subclass[ss]==tt,1,0)
        d_non[ss,paste(paste("sub",tt,sep=""),"T",sep="")] <- d_non[ss,paste("sub",tt,sep="")]*d_non[ss,"trt"]
      } # Creates indicators for each subclass and subclass x trt interaction in the non-sampled pop
    }
    
    mod_sub <- lm(as.formula(paste(paste(paste("y~",paste(paste(c(paste("sub",c(1:ff),sep=""),paste(paste("sub",c(1:ff),sep=""),"T",sep="")),collapse="+"),sep=""))),"-1",sep="")),data=samp_use)
    
    estimators["sub",1,k,pos] <- mean(mod_sub$coefficients[(ff+1):(2*ff)])
    
    # WEIGHTING BY ODDS ESTIMATOR
    
    gen_object <- generalize(outcome="y",treatment = "trt",trial = "study",selection_covariates = covars,data=d_combined,is_data_disjoint = T)
    
    estimators["odds_ipw",1,k,pos] <- gen_object$TATE$estimate
    
    # WEIGHTING BY OVERLAP ESTIMATOR
    
    estimators["overlap_ipw",1,k,pos] <- (sum(samp_use$y[samp_use$trt==1]*(samp_use$distance[samp_use$trt==1]))/(sum(samp_use$trt[samp_use$trt==1]*samp_use$distance[samp_use$trt==1]))) - 
      (sum(samp_use$y[samp_use$trt==0]*(samp_use$distance[samp_use$trt==0]))/(sum(1-samp_use$trt[samp_use$trt==0]*samp_use$distance[samp_use$trt==0])))
    
    # BART
    
    gen_bart <- generalize_bart(outcome="y",treatment = "trt",trial = "study",selection_covariates = covars,data=d_combined,is_data_disjoint = T)
    
    estimators["bart",1,k,pos] <- gen_bart$TATE$estimate
    
    # DOUBLY ROBUST WEIGHTING ESTIMATORS (Kern et al, 2016; Dahabreh et al, 2019) - There are three types of DR estimators. DR1 and DR2 use the conditional expectation of the outcome and the propensity score. DR3 uses a weighted outcome regression
    
    dr_weights <- generate_weights(Smod=as.formula(paste("study",paste(covars,collapse="+"),sep="~")),data=d_combined)
    
    d_combined2 <- dr_weights$dat
    
    d_combined2$y[is.na(d_combined2$y)] <-0
    
    d_combined2$w[is.na(d_combined2$w)] <-0
    
    d_combined2$trt[is.na(d_combined2$trt)] <-0
    
    OM<-OM_est(data=d_combined2,covars = covars,outcome = "y")
    
    d_combined2$p1<-OM$p1
    d_combined2$p0<-OM$p0
    
    estimators["drwls1",1,k,pos] <- dr1_est(d_combined2)$DR1
    
    estimators["drwls2",1,k,pos] <- dr2_est(d_combined2)$DR2
    
    estimators["drwls3",1,k,pos] <- dr3_est(d_combined2,covars = covars,outcome = "y")$DR3
    
    # TMLE ESTIMATOR
    
    gen_tmle <- generalize_tmle2(outcome="y", treatment="trt", trial="study", selection_covariates=covars, 
                                 data=d_combined,is_data_disjoint = TRUE)
    
    estimators["tmle",1,k,pos] <- gen_tmle$TATE$estimate
    
    # MODEL-BASED ESTIMATORS USING SMALL AREA ESTIMATION
    
    popn <- matrix(NA,ff,2)
    colnames(popn) <- c("subclass","popsize")
    
    for(g in 1:nrow(popn)){
      popn[g,1] <- g
      popn[g,2] <- nrow(d_use[d_use$subclass==g,])
    }
    
    #8) EBLUP
    
    cov_eblup_pop <- ddply(.data = d_non[,c("subclass",c(covars,paste("sub",c(1:ff),sep=""),paste(paste("sub",c(1:ff),sep=""),"T",sep="")))],.variables = "subclass",colwise(mean)) # pop means of covars_use
    
    cov_eblup_pop[,c(paste(paste("sub",c(1:ff),sep=""),"T",sep=""))] <- 0
    
    mod_eblup <- eblupBHF(as.formula(paste(paste(paste("y~",paste(paste(c(paste("sub",c(1:ff),sep=""),paste(paste("sub",c(1:ff),sep=""),"T",sep="")),collapse="+"),paste("+",paste(covars,collapse="+"),collapse = "+"),sep=""))),"-1",sep="")),
                          dom = subclass,meanxpop = cov_eblup_pop,popnsize = popn,method="REML",data=samp_use)
    
    estimators["eblup",1,k,pos] <- mean(mod_eblup$fit$fixed[c((ff+1):(2*ff))])
    
    #9) HB ESTIMATOR
    data_stan = list(
      N=nrow(samp_use),
      num_strata=ff,
      p=ncol(cov_eblup_pop)-1,
      subclass=as.numeric(samp_use$subclass),
      y_out=samp_use$y
    )
    
    for(cc in 1:(ncol(cov_eblup_pop)-1)){
      data_stan[[paste("x",cc,sep="")]] <- samp_use[,colnames(cov_eblup_pop)[cc+1]]
    }
    
    # RUN WITH THE CAUCHY/t PRIOR #
    
    stanmodel <- stan(file=paste(paste("sdl_",ff,sep=""),".stan",sep=""),data=data_stan, chains=0)
    
    stan.mod=stan(fit=stanmodel,data=data_stan,chains=chains,iter=iter_hb,control = list(adapt_delta=0.85))
    
    par_out<-array(NA,dim=c(nrow(stan.mod),ncol(cov_eblup_pop)+ff-1,chains)) # Get the parameters
    dimnames(par_out)<-list(NULL,c(paste("alpha",1:ff,sep=""),paste("beta",1:(ncol(cov_eblup_pop)-1),sep="")),paste("chain",1:chains,sep="")) # nrow(stan.mod) corresponds to half of iter_hb
    
    for(f in 1:chains){
      for(g in 1:(ncol(par_out)-ff)){
        for(j in 1:ff){
          # browser()
          par_out[,j,f]<-as.matrix(rstan::extract(stan.mod,paste(paste("alpha[",j,sep=""),"]",sep=""),permuted=FALSE,inc_warmup=FALSE)[,f,1])
          par_out[,g+ff,f]<-as.matrix(rstan::extract(stan.mod,paste(paste("beta[",g,sep=""),"]",sep=""),permuted=FALSE,inc_warmup=FALSE)[,f,1])
        }
      }
    }
    
    est_trt_hb <- matrix(NA,ff,chains)
    rownames(est_trt_hb) <- c(paste("Stratum",1:ff,sep=""))
    colnames(est_trt_hb) <- c(paste("chain",1:chains,sep=""))
    
    for(m in 1:ff){
      for(j in 1:ncol(est_trt_hb)){
        est_trt_hb[m,j] <- apply(data.frame(par_out[,(ncol(par_out)-ff+m),j]),2,mean)
      }
    }
    est_hb <- matrix(c(c(1:ff),apply(est_trt_hb,1,mean)),nrow=ff,ncol=2)
    
    estimators["hb",1,k,pos] <- mean(est_hb)
    
    ##############################################################################################################################################
    
    row_add <- num.strata - dim(t2_samp)[1]
    
    allocation_AAA[,,k,pos]<-as.matrix(rbind(t2_samp[1:(dim(t2_samp)[1]),1:2],matrix(NA,row_add,2)))
    
    maxstrata_AAA[k,,pos]<-ifelse(((any(allocation_AAA[,,k,pos][,1]==1,na.rm=TRUE)==TRUE | any(allocation_AAA[,,k,pos][,1]==2,na.rm=TRUE)==TRUE)
                                   && (any(allocation_AAA[,,k,pos][,2]==1,na.rm=TRUE)==TRUE|any(allocation_AAA[,,k,pos][,2]==2,na.rm=TRUE)==TRUE)),
                                  1,0)
    
    #10) OUTCOMES TO ASSESS PERFORMANCE OF EACH ESTIMATOR
    
    # MEAN SQUARED ERROR
    
    mse_est[1,1,k,pos]<- (truth - estimators[1,1,k,pos])^2
    mse_est[2,1,k,pos]<- (truth - estimators[2,1,k,pos])^2
    mse_est[3,1,k,pos]<- (truth - estimators[3,1,k,pos])^2
    mse_est[4,1,k,pos]<- (truth - estimators[4,1,k,pos])^2
    mse_est[5,1,k,pos]<- (truth - estimators[5,1,k,pos])^2
    mse_est[6,1,k,pos]<- (truth - estimators[6,1,k,pos])^2
    mse_est[7,1,k,pos]<- (truth - estimators[7,1,k,pos])^2
    mse_est[8,1,k,pos]<- (truth - estimators[8,1,k,pos])^2
    mse_est[9,1,k,pos]<- (truth - estimators[9,1,k,pos])^2
    mse_est[10,1,k,pos]<- (truth - estimators[10,1,k,pos])^2
    
    # BIAS
    
    bias_est[1,1,k,pos]<- (truth - estimators[1,1,k,pos])
    bias_est[2,1,k,pos]<- (truth - estimators[2,1,k,pos])
    bias_est[3,1,k,pos]<- (truth - estimators[3,1,k,pos])
    bias_est[4,1,k,pos]<- (truth - estimators[4,1,k,pos])
    bias_est[5,1,k,pos]<- (truth - estimators[5,1,k,pos])
    bias_est[6,1,k,pos]<- (truth - estimators[6,1,k,pos])
    bias_est[7,1,k,pos]<- (truth - estimators[7,1,k,pos])
    bias_est[8,1,k,pos]<- (truth - estimators[8,1,k,pos])
    bias_est[9,1,k,pos]<- (truth - estimators[9,1,k,pos])
    bias_est[10,1,k,pos]<- (truth - estimators[10,1,k,pos])
    
  }
}

save.image("ind_sim3.RData")