# CODE FOR GENERATING SELECTION PROBABILITIES (INDIANA CRT)

#0) Load packages

library(foreign)
library(MatchIt)
library(plyr)
library(xtable)
library(boot)

#1) Read in RCT data.
library(readr)
data <- read_csv("simulation_practice/indiana_data.csv")

#2) Covariates
covars_list<-c("z_s09_ELA","z_s09_math","z_f08_ELA","z_f08_math","z_05_ELA","z_05_math","z_06_ELA","z_06_math","z_07_ELA","z_07_math","avg_teach_salary_08","Attendance_0809","Attendance_0708","FTE08","MEMBER08","PUPTCH08","urban08","suburb08","county.pop","TITLEI08","STITLI08","CHARTR08","male","white","sped","frpl","lep") # All covars

###### With county.pop and no CHARTR08, Attendance_0708, avg_teach_salary_08

covars_red<-c("z_s09_ELA","z_s09_math","Attendance_0809","FTE08","MEMBER08","PUPTCH08","county.pop","TITLEI08","STITLI08","male","white","sped","frpl","lep") 

#3) Create variables for sample and population schools
data$N.EXPT <- ifelse((data$study==0 && data$trt==0),1,0)

attach(data)

data <- cbind(SchoolNo,data.frame(scale(data[c(covars_red)])),study,trt,N.EXPT)

#4) Generate selection probabilities

alpha <- c(-2,-2,rep(0.00001,5),rep(0.00005,3))

for(i in 1:nrow(data)){
  
  data[i,"prob"] <- inv.logit(sum(alpha[1:2]*c(data[i,"z_s09_ELA"]^2,data[i,"z_s09_math"]^2)) + sum(alpha[3:7]*data[i,c("FTE08","MEMBER08","PUPTCH08","county.pop","TITLEI08")]) + sum(alpha[8:10]*c(data[i,"sped"]^2,data[i,"frpl"]^2,data[i,"lep"]^2)))
  
}

summary(data[,"prob"])

