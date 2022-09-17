rm(list=ls())

setwd("~/wanglijie/Bay")
library("igraph");library("qtl");library(BayesNetBP);library(bnlearn)
library("dplyr")
library(riskRegression)
library(cmprsk)
library(survival)
# #########################################################
# factor_names;numeric_names
# #########################################################
load("data_bayesian.RData")
type_data_bayesian <- as.data.frame(unlist(sapply(data_bayesian, class)))
names(type_data_bayesian) <- "TYPE"

load("data_cohort_smo_no.RData")
data_cohort_cox <- subset(data_cohort_sub1, Tab=="data_training")
var_need <- c(rownames(type_data_bayesian),"Cancer_target","Time_of_follow_up_duration")
data_cohort_cox1 <- data_cohort_cox[,var_need]
type_data_bayesian <- as.data.frame(unlist(sapply(data_cohort_cox1, class)))
names(type_data_bayesian) <- "TYPE"

factor_names <- rownames(subset(type_data_bayesian, TYPE=="factor"))
factor_names <- setdiff(factor_names,"Cancer_target")
for (i in factor_names){
  data_cohort_cox1[,i] <- as.factor(as.character(as.integer(data_cohort_cox1[,i])))
  print (i)
}

quantcut_new <- function(x, q=5, na.rm=T)
{
  # x <- a[218:249]
  # q <- 4
  if(length(q)==1)
    q <- seq(0,1, length.out=q+1)
  
  quant <- quantile(x, q, na.rm=T)
  dups <- duplicated(quant)
  if(any(dups))
  {
    flag <- x %in% unique(quant[dups])
    retval <- ifelse(flag,
                     paste("[",
                           as.character(x),
                           "]",
                           sep=''),
                     NA)
    uniqs <- unique(quant)
    
    # move cut points over a bit...
    reposition <- function(cut)
    {
      flag <- (x>=cut)
      if(sum(flag, na.rm=T)==0) #
        return(cut)
      else
        return(min(x[flag], na.rm=T))
    }
    
    newquant <- sapply(uniqs, reposition)
    retval[!flag] <- as.character(cut(x[!flag],
                                      breaks=newquant,
                                      include.lowest=TRUE)) #,...
    
    levs <- unique(retval[order(x)]) # ensure factor levels are
    levs <- as.character(na.omit(levs))
    # properly ordered
    retval <- factor(retval, levels=levs)
    
    ## determine open/closed interval ends
    mkpairs <- function(x) # make table of lower, upper
    {sapply(x,
            function(y) if(length(y)==2) y[c(2,2)] else y[2:3])}
    pairs <- mkpairs(strsplit(levs, '[^0-9+\\.\\-]+'))
    rownames(pairs) <- c("lower.bound","upper.bound")
    colnames(pairs) <- levs
    
    closed.lower <- rep(F,ncol(pairs)) # default lower is open
    closed.upper <- rep(T,ncol(pairs)) # default upper is closed
    closed.lower[1] <- TRUE             # lowest interval is always closed
    
    for(i in 2:ncol(pairs))            # open lower interval if above singlet
      if(pairs[1,i]==pairs[1,i-1] && pairs[1,i]==pairs[2,i-1])
        closed.lower[i] <- FALSE
    
    for(i in 1:(ncol(pairs)-1))        # open upper interval if below singlet
      if(pairs[2,i]==pairs[1,i+1] && pairs[2,i]==pairs[2,i+1])
        closed.upper[i] <- FALSE 
    levs <- ifelse(pairs[1,]==pairs[2,],
                   pairs[1,],
                   paste(ifelse(closed.lower,"[","("),
                         pairs[1,],
                         ",",
                         pairs[2,],
                         ifelse(closed.upper,"]",")"),
                         sep='')) 
    levels(retval) <- levs
  }
  else
    retval <- cut(x, quant, include.lowest=TRUE)
  return(retval)
}
numeric_names <- rownames(subset(type_data_bayesian, TYPE=="numeric"))
numeric_names <- setdiff(numeric_names,"Time_of_follow_up_duration")
for (i in numeric_names){
  data_cohort_cox1[,i] <- quantcut_new(data_cohort_cox1[,i]) 
  print (i)
}
for (i in numeric_names){
  data_cohort_cox1[,i] <- as.factor(as.character(as.integer(data_cohort_cox1[,i])))
  print (i)
}

# #########################################################
# Imputed
# #########################################################
load("Bay_01.RData")
data_cohort_imputed <- as.data.frame(lung@imputed.data)
type_data_bayesian <- as.data.frame(unlist(sapply(data_cohort_imputed, class)))
names(type_data_bayesian) <- "TYPE"
for (i in factor_names){
  data_cohort_imputed[,i] <- as.factor(as.character(data_cohort_imputed[,i]))
  print (i)
}
for (i in numeric_names){
  data_cohort_imputed[,i] <- quantcut_new(data_cohort_imputed[,i]) 
  print (i)
}
for (i in numeric_names){
  data_cohort_imputed[,i] <- as.factor(as.character(as.integer(data_cohort_imputed[,i])))
  print (i)
}
# #########################################################
# Cox
# #########################################################
library(rms);library(pec)
var_final <- c(factor_names,numeric_names)
var_end <- c("Time_of_follow_up_duration","Cancer_target")
data_result <- data_cohort_cox1[,var_end]
data_cox <- cbind(data_cohort_imputed,data_result)        

data_cox$Cancer_target <- ifelse(data_cox$Cancer_target=="Yes",1,0)
BaSurv <- Surv(data_cox$Time_of_follow_up_duration, data_cox$Cancer_target)
FML <- as.formula(paste0("BaSurv~",paste0(var_final, collapse = "+")))

dd <- datadist(data_cox)
options(datadist='dd')
cox.back <- selectCox(FML,data = data_cox, rule ="aic")
GSum <- cox.back$fit

var_final <- c("Age_group","R59","J60_J65","J12_J18","NEUT","Ca")

#All possible combination
data_event <- expand.grid(levels(data_cox$Age_group), levels(data_cox$R59),levels(data_cox$J60_J65),
                          levels(data_cox$J12_J18),levels(data_cox$NEUT),levels(data_cox$Ca))
names(data_event) <- c('Age_group',"R59","J60_J65","J12_J18","NEUT","Ca")
data_event <- data.frame(data_event)
data_event <- cbind(data_event,data_event)
names(data_event)[7:12] <- c('Age_group_score',"R59_score","J60_J65_score","J12_J18_score","NEUT_score","Ca_score")
for(i in 1:ncol(data_event)){
  data_event[,i] <- as.character(data_event[,i])
}
# #########################################################
# Risk
# #########################################################
# Cox Ft 
cox.fit1  <- cph(Surv(time,status)~X1+X2,data=train,surv=TRUE)
#Model prediction
# C-index-internal
set.seed(20191139)
p <-(survConcordance(Surv(trainset_na$time, trainset_na$status) ~  predict(cox.fit1)))$`concordance`
# C-index-extenal
q <- 1-(rcorrcens(Surv(time, status) ~ predict(cox.fit1, newdata=testset_na), data = testset_na))[1]
1-predictSurvProb(cox.fit1,newdata=test,times=c(1000,2000,3000))

#Cause-specific hazard regression Ft #status0 1 2 numeric factor
CSH<-CSC(Hist(time,status)~sex+age+invasion,data=Melanoma)
summary(CSH)
#Model prediction
pec:::predictEventProb(CSH,cause=1,newdata=test,time=c(1000,2000,3000))
predictRisk(CSH,cause=1,newdata=test,time=c(1000,2000,3000))

#Subdistribution hazards (SHs) model
SH <- FGR(Hist(time,status)~sex+age+invasion,data=Melanoma)
#Model prediction
predictRisk(SH,newdata,times=c(3500),cause=1)

# C-index for Cause-specific hazard regression and SHs 
library(party)
f <- pecCforest(Surv(time,status) ~age+ sim_patho_diagnosis+single_tumor
                +karnofsky_performance_score+midline_pzt +degree_of_resection,
                data=data,controls = cforest_unbiased(mtry = 2))

forest <- 1-predictSurvProb(f,newdata=data,times=365)#
forest <- predictRisk(CSH,newdata=data,times=365)
data$end1 <- ifelse(data$time>=365,0,data$status)
data$sftime1 <- ifelse(data$time>=365,365,data$time)

cfor=1-rcorr.cens(forest,Surv(data$time,data$status==1))
cforq=1-rcorr.cens(forest,Surv(data$sftime1,data$end1))
cfor
cforq
forestp <- 1-predictSurvProb(f,newdata=testset,times=365)

testset$end1 <- ifelse(testset$time>=365,0,testset$status)
testset$sftime1 <- ifelse(testset$time>=365,365,testset$time)

cfor=1-rcorr.cens(forestp,Surv(testset$time,testset$status))
cforq=1-rcorr.cens(forestp,Surv(testset$sftime1,testset$end1))
cfor
cforq

data_cohort_cox2 <- data_cox[,c("Time_of_follow_up_duration","Cancer_target",var_final)]

BaSurv <- Surv(data_cohort_cox2$Time_of_follow_up_duration, data_cohort_cox2$Cancer_target)
FML <- as.formula(paste0("BaSurv~",paste0(var_final, collapse = "+")))

res.coxs<-coxph(FML,data=data_cohort_cox2)

# #########################################################
# Combine with Bay
# #########################################################
var_cox <- c("Age_group","R59","J60_J65","J12_J18","NEUT","Ca")
var_all <- c(factor_names,numeric_names)
data <- data_cox[,var_all]

colnames(data) <- paste0("a",1:ncol(data))
rownames(data) <- 1:nrow(data)

lable <- cbind(colnames(data),var_all)
# #########################################################
# Learning Bay
# #########################################################
library(openxlsx)
white_list <- read.xlsx("whitelist.xlsx",1)
white_list <- as.matrix(white_list)
black_list <- read.xlsx("whitelist.xlsx",2)

names(black_list)[1] <- c("var_all")
black_list <- merge(black_list,lable,all.x = T,by="var_all")
black_list <- black_list[,-1]
names(black_list) <- c("var_all","from")
black_list <- merge(black_list,lable,all.x = T,by="var_all")
black_list <- black_list[,-1]
names(black_list) <- c("from","to")
black_list <- subset(black_list,from!="a5"&from!="a6"&from!="a15"&from!="a1")


black_list <- as.matrix(black_list)
var_drop <- c("a5","a6","a15","a1")
data <- data[,!colnames(data)%in%var_drop]
#structure learning
structure_bay <- bnlearn::hc(data,blacklist = black_list)
#parameter learning
#fit <- bn.fit(structure_bay, data, cluster = cl);class(structure_bay);class(fit)

# #########################################################
# Bayesian inference
# #########################################################
#Missing data
data_missing <- data_cohort_cox1[,var_all]
colnames(data_missing) <- paste0("a",1:ncol(data_missing))
rownames(data_missing) <- 1:nrow(data_missing)

library(parallel);library(bnlearn)
cl <- makeCluster(2)
data_missing <- data_missing[,!colnames(data_missing)%in%var_drop] 
fit <- bn.fit(structure_bay, data_missing, cluster = cl)

data_evid <- data_missing
vars_evidence <- colnames(data_evid)

colnames(data_event)[1:6] <- c("a61","a51","a31","a27","a67","a63")

vars_event <- var_cox
vars_event <- colnames(data_event)[1:6]

data_csbn <- merge(data_evid, data_event[,c("a61","a51","a31","a27","a67","a63","Ft_select_3","Ft_select_5","Ft_select_1")], by=c("a61","a51","a31","a27","a67","a63"),all.x = T)

data_csbn$check_na <- ifelse(is.na(data_csbn$a61)|is.na(data_csbn$a51)|is.na(data_csbn$a31)|is.na(data_csbn$a27)|is.na(data_csbn$a67)|is.na(data_csbn$a63),"MISSING","COMPLETE")
table(data_csbn$check_na)

mysubset <- function(df, ssubset) {
  subset(df, eval(parse(text=ssubset)))
}

data_bay_summ <- data.frame(matrix(NA,0,0))
data_bay_people <- data.frame(matrix(NA,0,0))

#save(data_csbn, data_event, vars_evidence, vars_event, fit, file="~/wanglijie/Bay/UKB_BAY.RData"  )
#load(file="~/wanglijie/Bay/UKB_BAY.RData")
load("~/wanglijie/Bay/UKB_BAY.RData")
data_csbn <- subset(data_csbn,data_csbn$check_na!="COMPLETE") 
data_csbn <- data_csbn[c(1:2),]
library(parallel);library(bnlearn)
cl <- makeCluster(16)
clusterSetRNGStream(cl, iseed = 1)
mysubset <- function(df, ssubset) {
  subset(df, eval(parse(text=ssubset)))
}

data_bay_summ <- data.frame(matrix(NA,0,0))
data_bay_people <- data.frame(matrix(NA,0,0))


timestart <- Sys.time()
####################################process()要处理一行的
#######################################################
reprocess <- function(x,vars_event,str1){
  str3 <- paste("(", vars_event, "=='",
                sapply(x[vars_event], as.character), "')",
                sep = "", collapse = " & ")
  set.seed(1)
  xxx <- eval(parse(text=paste("cpquery(fit,event =(",str3,"), evidence = (",str1,"),n=100000000,cl)")))
  return(xxx)
}
process <- function(x,vars_evidence,vars_event,data_event){
  if (x[length(x)] == "MISSING"){
    #evidence
    vars_evidence_nomi <- vars_evidence[which(is.na(x[vars_evidence])==F)]
    str1 <- paste("(", vars_evidence_nomi, "=='",
                  sapply(x[vars_evidence_nomi], as.character), "')",
                  sep = "", collapse = " & ")
    #event
    vars_event_nomi <- vars_event[which(is.na(x[vars_event])==F)]
    
    str2 <- paste("(", vars_event_nomi, "=='",
                  sapply(x[vars_event_nomi], as.character), "')",
                  sep = "", collapse = " & ")
    data_event_1 <- mysubset(data_event, str2)
    ####################################################################################
    aabb <- apply(data_event_1,MARGIN = 1,reprocess,vars_event,str1)#aabb是随便写的名字
    ####################################################################################
    #abcd <- length(aabb)
    bbcc <-  sum(aabb*data_event_1[,"Ft_select_1"])                 #bbcc是随便写的名字
    
  }
  else{ bbcc <-  x["Ft_select_1"]
  }
  return(bbcc)
}

aaa16 <- apply(data_csbn, MARGIN=1,process,vars_evidence,vars_event,data_event)
timeend<-Sys.time()
runningtime16<-timeend-timestart
print(runningtime16)   
aaa16
# #########################################################
# Plot the graph #http://106.15.72.70:3838/shinyBN/
# #########################################################



# #########################################################
# Validation
# #########################################################
ares.coxs <- roc(test4$es5,pp5,plot=TRUE,legacy.axes=T,print.thres=T,print.auc=TRUE)##0.745
ci(ares.coxs)#0.7155-0.7744
AIC(res.coxs)#4855.864