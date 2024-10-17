#---------------------------------------------------------------- CSF  --------------------------------------------------------------#
CSF_HTE <- function(data=dat_all[which(dat_all$sim==1),],testdat=mydata$data,time.interest=13,k.folds=10,response_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),tree.unbiased=TRUE,use.grad = TRUE,pre.ctree.control.alpha=0.05,pre.ctree.control.testtype='Bonferroni',rulefit_type="both",ntrees=500,maxdepth=4,normalize=TRUE,nfolds.rulefit.cv=10,removecomplements=FALSE,glmnet.alpha=1, glmnet.relax=F,glmnet.gamma = c(0, 0.25, 0.5, 0.75, 1),glmnet.intercept=FALSE,glmnet.linearterm=F,glmnet.linearterm.X=paste0("V",c(1:10))){
  data<-data[complete.cases(data),]

  cs.forest.prob <- causal_survival_forest(data[,-which(names(data)%in%c("Treatment","Time","Event","sim"))], data$Time, data$Treatment, data$Event,target = "survival.probability",horizon = time.interest,num.trees = 500,mtry = sqrt(15),min.node.size = 15,num.threads=1) # cannot tune tree depth
  csf.prob = predict(cs.forest.prob, data[,-which(names(data)%in%c("Treatment","Time","Event","sim"))])
  data$pseudo<-csf.prob$predictions
  
  tryCatch({
    modeld.pre <- pre::pre(pseudo ~., data = data[,-which(colnames(data)%in%c("Treatment","Time","Event","sim"))],tree.unbiased = tree.unbiased ,maxdepth = maxdepth,nfolds = nfolds.rulefit.cv,type=rulefit_type,normalize = normalize,use.grad = use.grad,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype),nlambda=20,type.measure ="mse",relax = glmnet.relax,gamma =glmnet.gamma,intercept=glmnet.intercept, removecomplements=removecomplements) # ,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype)
  },warning=function(e){
    modeld.pre <<- "error in ctree"
  },
  error=function(e){
    modeld.pre <<- "no candidate rules"
  })
  
  if (length(modeld.pre)==1){
    imps.pre.std <- "no models"
  }else{
    tryCatch({
      imps.pre.std <- pre::importance(modeld.pre, standardize = T, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.std <<- "all shrunk to 0"
    }) 
    
  }
  
  if (length(modeld.pre)==1){
    imps.pre.nonstd <- "no models"
  }else{
    tryCatch({
      imps.pre.nonstd <- pre::importance(modeld.pre, standardize = F, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.nonstd <<- "all shrunk to 0"
    })
    
  }
  
  if (length(modeld.pre)==1){
    pred.pre <- rep("no models",dim(testdat)[1])
  }else if(length(modeld.pre)>1 & length(imps.pre.std)==1){
    pred.pre <- rep("all shrunk to 0",dim(testdat)[1])
  }else{
    pred.pre <- predict(modeld.pre,newdata = testdat[,-which(names(testdat)%in%c("Treatment","Time","Event"))])
  }

  return(list(diff=pred,pre_result=modeld, final_ensemble_result= final.model.ensemble, num.initial.rule=dim(initial_rules_mat)[2],diff.pre=pred.pre,imp.pre.std=imps.pre.std,imp.pre.nonstd=imps.pre.nonstd,pre_result_fit=modeld.pre)) 
  
}

#-------------------------------------------------------- BAFT ----------------------------------------------------------#
baft_HTE_T <- function(data=dat_all[which(dat_all$sim==1),],testdat=mydata$data,time.interest=6.7,k.folds=10,noise.var=T,response_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),tree.unbiased=TRUE,use.grad = TRUE,pre.ctree.control.alpha=0.05,pre.ctree.control.testtype='Bonferroni',rulefit_type="both",ntrees=500,maxdepth=4,normalize=TRUE,nfolds.rulefit.cv=10,removecomplements=FALSE,glmnet.alpha=1, glmnet.relax=F,glmnet.gamma = c(0, 0.25, 0.5, 0.75, 1),glmnet.intercept=FALSE,glmnet.linearterm=F,glmnet.linearterm.X=paste0("V",c(1:10))){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  
  # Fit AFTrees model 
  AFTrees1 = AFTrees(
    x.train = data1[,-c(which(names(data1)%in%c("Time","Event")))],
    y.train =  data1$Time,
    status = data1$Event,
    nskip = 100,
    ndpost = 1000,
    ntree=500,
    x.test = data[,-c(which(names(data)%in%c("Time","Event","sim","Treatment")))],
    nonparametric = T
  )
  
  AFTrees0 = AFTrees(
    x.train = data0[,-c(which(names(data0)%in%c("Time","Event")))],
    y.train =  data0$Time,
    status = data0$Event,
    nskip = 100,
    ndpost = 1000,
    ntree=500,
    x.test = data[,-c(which(names(data)%in%c("Time","Event","sim","Treatment")))],
    nonparametric = T
  )
  
  # Calculate survival probability
  AFTrees1_all_survival_prob <- AFTrees_SurvivalProb(object = AFTrees1, test.only = T,time.points=rep(time.interest,2))
  AFTrees0_all_survival_prob <- AFTrees_SurvivalProb(object = AFTrees0, test.only = T,time.points=rep(time.interest,2)) 
  AFTrees1_survival_prob <- AFTrees1_all_survival_prob$Surv.test[,1]
  AFTrees0_survival_prob <- AFTrees0_all_survival_prob$Surv.test[,1]
  
  data$pseudo <- AFTrees1_survival_prob - AFTrees0_survival_prob
  
  tryCatch({
    modeld.pre <- pre::pre(pseudo ~., data = data[,-which(colnames(data)%in%c("Treatment","Time","Event","sim"))],tree.unbiased = tree.unbiased ,maxdepth = maxdepth,nfolds = nfolds.rulefit.cv,type=rulefit_type,normalize = normalize,use.grad = use.grad,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype),nlambda=20,type.measure ="mse",relax = glmnet.relax,gamma =glmnet.gamma,intercept=glmnet.intercept, removecomplements=removecomplements) # ,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype)
  },warning=function(e){
    modeld.pre <<- "error in ctree"
  },
  error=function(e){
    modeld.pre <<- "no candidate rules"
  })
  
  if (length(modeld.pre)==1){
    imps.pre.std <- "no models"
  }else{
    tryCatch({
      imps.pre.std <- pre::importance(modeld.pre, standardize = T, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.std <<- "all shrunk to 0"
    })
    
  }
  
  if (length(modeld.pre)==1){
    imps.pre.nonstd <- "no models"
  }else{
    tryCatch({
      imps.pre.nonstd <- pre::importance(modeld.pre, standardize = F, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.nonstd <<- "all shrunk to 0"
    })
    
  }
  
  if (length(modeld.pre)==1){
    pred.pre <- rep("no models",dim(testdat)[1])
  }else if(length(modeld.pre)>1 & length(imps.pre.std)==1){
    pred.pre <- rep("all shrunk to 0",dim(testdat)[1])
  }else{
    pred.pre <- predict(modeld.pre,newdata = testdat[,-which(names(testdat)%in%c("Treatment","Time","Event"))])
  }

  return(list(diff=pred,pre_result=modeld, final_ensemble_result= final.model.ensemble, num.initial.rule=dim(initial_rules_mat)[2],diff.pre=pred.pre,imp.pre.std=imps.pre.std,imp.pre.nonstd=imps.pre.nonstd,pre_result_fit=modeld.pre)) 
  
}

#-------------------------------------------------------------- DR-learner ---------------------------------------------------------#
HTE_DR_Kennedy <- function(data=dat_all[which(dat_all$sim==1),],testdat=mydata$data,time.interest=6.7,k.folds=10,noise.var=T,est_pi=1,propensity_method="correct_logistic",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,response_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),tree.unbiased=TRUE,use.grad = TRUE,pre.ctree.control.alpha=0.05,pre.ctree.control.testtype='Bonferroni',rulefit_type="both",ntrees=500,maxdepth=4,normalize=TRUE,nfolds.rulefit.cv=10,removecomplements=FALSE,glmnet.alpha=1, glmnet.relax=F,glmnet.gamma = c(0, 0.25, 0.5, 0.75, 1),glmnet.intercept=FALSE,glmnet.linearterm=TRUE,glmnet.linearterm.X=paste0("V",c(1:10))){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  
  if (noise.var==F){
    rfsrc_data1 <-
      grf::survival_forest(X=as.matrix(data1[,which(names(data1) %in% response_X)]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
    rfsrc_data0 <-
      grf::survival_forest(X=as.matrix(data0[,which(names(data0) %in% response_X)]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
  }else if (noise.var==T){
    rfsrc_data1 <-
      grf::survival_forest(X=as.matrix(data1[,-which(names(data1) %in% c("Time","Event"))]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
    rfsrc_data0 <-
      grf::survival_forest(X=as.matrix(data0[,-which(names(data0) %in% c("Time","Event"))]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
  }
  
  if (est_pi==1){
    if (noise.var==F){
      if (propensity_method=="correct_logistic"){
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
        data$pi.hat <- predict(rf,type='response')
      }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,-which(names(data)%in%c("Time","Event","Treatment","sim"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = 500,num.threads = 1)
        data$pi.hat<-predict(rf)$predictions[,2]
      }else if(propensity_method=="lasso_logistic"){
        rf <- cv.glmnet(y=data$Treatment, x=as.matrix(data[,which(names(data)%in%c(propensity_method_X))]),family = "binomial", type.measure = "class",alpha=1,nlambda = 20)
        data$pi.hat <- predict(rf,as.matrix(data[,which(names(data)%in%c(propensity_method_X))]), type = "response") # This predict() function is for glmnet
        
      }
      
    }else if (noise.var==T){
      if (propensity_method=="correct_logistic"){ 
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,-which(names(data)%in%c("Time","Event","sim"))])
        data$pi.hat <- predict(rf,type='response')
      }else if (propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,-which(names(data)%in%c("Time","Event","Treatment","sim"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = 500,num.threads = 1)
        data$pi.hat<-predict(rf)$predictions[,2]
      }else if(propensity_method=="lasso_logistic"){
        rf <- cv.glmnet(y=data$Treatment, x=as.matrix(data[,which(names(data)%in%c(propensity_method_X))]),family = "binomial", type.measure = "class",alpha=1,nlambda = 20)
        data$pi.hat <- predict(rf,as.matrix(data[,which(names(data)%in%c(propensity_method_X))]), type = "response") # This predict() function is for glmnet
        
      }
      
    }
    
  }else if (est_pi==0){
    data$pi.hat <- rep(propensity,dim(data)[1])
  }
  
  if (noise.var==F){
    data$predict_rfsrc_median0 <- predict(rfsrc_data0, newdata=data[,which(names(data)%in%response_X)], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions 
    data$predict_rfsrc_median1 <- predict(rfsrc_data1, newdata=data[,which(names(data)%in%response_X)], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions 
    
  }else if(noise.var==T){
    data$predict_rfsrc_median0 <- predict(rfsrc_data0, newdata=data[,-which(names(data)%in%c("Treatment","Time","Event","sim","pi.hat"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions 
    data$predict_rfsrc_median1 <- predict(rfsrc_data1, newdata=data[,-which(names(data)%in%c("Treatment","Time","Event","sim","pi.hat","predict_rfsrc_median0"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions 
    
  }
  
  data$I<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
  
  data$pseudo <-as.vector(((data$Treatment-data$pi.hat)/(data$pi.hat*(1-data$pi.hat)))*(data$I-data$Treatment*data$predict_rfsrc_median1-(1-data$Treatment)*data$predict_rfsrc_median0) + data$predict_rfsrc_median1-data$predict_rfsrc_median0 )
  
  U <- pmin(data$Time, time.interest)                      
  fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
  C.hat <- rep(NA, length(fold.id))
  for (z in 1:k.folds) {
    c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
    kmc <- summary(c.fit, times = U[fold.id == z])
    C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
  }
    
  C.hat.cen <- C.hat[!is.na(data$I)]
  
  weight.cen <- 1/ C.hat.cen
  
  data.complete <- data[!is.na(data$pseudo),-which(names(data)%in%c("Treatment","Time","Event","sim","pi.hat","predict_rfsrc_median0","predict_rfsrc_median1","I"))] 
  
  tryCatch({ 
    modeld.pre <- pre::pre(pseudo ~., data = data.complete,weights = weight.cen,tree.unbiased = tree.unbiased ,maxdepth = maxdepth,nfolds = nfolds.rulefit.cv,type=rulefit_type,normalize = normalize,use.grad = use.grad,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype),nlambda=20,type.measure ="mse",relax = glmnet.relax,gamma =glmnet.gamma,intercept=glmnet.intercept, removecomplements=removecomplements) # ,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype)
  },warning=function(e){
    modeld.pre <<- "error in ctree"
  },
  error=function(e){
    modeld.pre <<- "no candidate rules"
  })
  
  if (length(modeld.pre)==1){
    imps.pre.std <- "no models"
  }else{
    tryCatch({
      imps.pre.std <- pre::importance(modeld.pre, standardize = T, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.std <<- "all shrunk to 0"
    })

  }
  
  if (length(modeld.pre)==1){
    imps.pre.nonstd <- "no models"
  }else{
    tryCatch({
      imps.pre.nonstd <- pre::importance(modeld.pre, standardize = F, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.nonstd <<- "all shrunk to 0"
    })
    
  }

  if (length(modeld.pre)==1){
    pred.pre <- rep("no models",dim(testdat)[1])
  }else if(length(modeld.pre)>1 & length(imps.pre.std)==1){
    pred.pre <- rep("all shrunk to 0",dim(testdat)[1])
  }else{
    pred.pre <- predict(modeld.pre,newdata = testdat[,-which(names(testdat)%in%c("Treatment","Time","Event"))])
  }

  return(list(diff=pred,pre_result=modeld, final_ensemble_result= final.model.ensemble, num.initial.rule=dim(initial_rules_mat)[2],diff.pre=pred.pre,imp.pre.std=imps.pre.std,imp.pre.nonstd=imps.pre.nonstd,pre_result_fit=modeld.pre)) 
  
}



#------------------------------------------------- R-learner ------------------------------------------------#
rsf_HTE_R <- function(data=dat_all[which(dat_all$sim==1),],testdat=mydata$data,time.interest=4,k.folds=10,noise.var=F,est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,response_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),tree.unbiased=TRUE,use.grad = TRUE,pre.ctree.control.alpha=0.05,pre.ctree.control.testtype='Bonferroni',rulefit_type="rules",ntrees=500,maxdepth=5,normalize=TRUE,nfolds.rulefit.cv=10,removecomplements=TRUE,glmnet.alpha=1, glmnet.relax=F,glmnet.gamma = c(0, 0.25, 0.5, 0.75, 1),glmnet.intercept=FALSE,glmnet.linearterm=TRUE,glmnet.linearterm.X=paste0("V",c(1:10))){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  
  if (est_pi==1){
    
    if (noise.var==F){
      if (propensity_method=="correct_logistic"){ 
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
        data$pred.pi <- predict(rf,type='response')
      }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,-which(names(data)%in%c("Treatment","Time","Event","sim"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = 500,num.threads = 1)
        data$pred.pi<-predict(rf)$predictions[,2]
      }else if(propensity_method=="lasso_logistic"){
        rf <- cv.glmnet(y=data$Treatment, x=as.matrix(data[,which(names(data)%in%c(propensity_method_X))]),family = "binomial", type.measure = "class",alpha=1,nlambda = 20)
        data$pred.pi <- predict(rf,as.matrix(data[,which(names(data)%in%c(propensity_method_X))]), type = "response") # This predict() function is for glmnet
        
      }
      
    }else if(noise.var==T){
      if (propensity_method=="correct_logistic"){ # not cv version
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,-which(names(data)%in%c("Time","Event","sim"))])
        data$pred.pi <- predict(rf,type='response')
      }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,-which(names(data)%in%c("Treatment","Time","Event","sim"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = 500,num.threads = 1)
        data$pred.pi<-predict(rf)$predictions[,2]
      }else if(propensity_method=="lasso_logistic"){
        rf <- cv.glmnet(y=data$Treatment, x=as.matrix(data[,which(names(data)%in%c(propensity_method_X))]),family = "binomial", type.measure = "class",alpha=1,nlambda = 20)
        data$pred.pi <- predict(rf,as.matrix(data[,which(names(data)%in%c(propensity_method_X))]), type = "response") # This predict() function is for glmnet
        
      }
      
    }
    
  }else if(est_pi==0){
    data$pred.pi <- rep(propensity, dim(data)[1])
  }
  
  if(noise.var==F){
    rfsrc_data1 <-
      grf::survival_forest(X=as.matrix(data1[,which(names(data1) %in% response_X)]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
    rfsrc_data0 <-
      grf::survival_forest(X=as.matrix(data0[,which(names(data0) %in% response_X)]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
  }else if(noise.var==T){
    rfsrc_data1 <-
      grf::survival_forest(X=as.matrix(data1[,-which(names(data1) %in% c("Time","Event"))]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
    rfsrc_data0 <-
      grf::survival_forest(X=as.matrix(data0[,-which(names(data0) %in% c("Time","Event"))]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
  }
  
  if (noise.var==F){
    data$s0_hat <- predict(rfsrc_data0, newdata=data[,which(names(data)%in%response_X)], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions # summary of s0_hat on training data: 0.6 to 0.8
    data$s1_hat <- predict(rfsrc_data1, newdata=data[,which(names(data)%in%response_X)], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions # summary of s0_hat on training data: 0.5 to 0.9
    
  }else if(noise.var==T){
    data$s0_hat <- predict(rfsrc_data0, newdata=data[,-which(names(data)%in%c("Treatment","Time","Event","sim","pred.pi"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions # summary of s0_hat on training data: 0.6 to 0.8
    data$s1_hat <- predict(rfsrc_data1, newdata=data[,-which(names(data)%in%c("Treatment","Time","Event","sim","pred.pi","s0_hat"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions # summary of s0_hat on training data: 0.5 to 0.9
    
  }
  data$m_hat <- data$pred.pi*data$s1_hat+(1-data$pred.pi)*data$s0_hat
  
  data$I<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
  data$pseudo <- as.vector(((data$I) - data$m_hat) / (data$Treatment - data$pred.pi) )

  U <- pmin(data$Time, time.interest)                        
  fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
  C.hat <- rep(NA, length(fold.id))
  for (z in 1:k.folds) {
    c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
    kmc <- summary(c.fit, times = U[fold.id == z])
    C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
  }

  weight <- (1/C.hat)*((data$Treatment - data$pred.pi)^2)
  weight.cen <-  weight[!is.na(data$I)]
  
  data.complete <- data[!is.na(data$pseudo),-which(names(data)%in%c("Treatment","Time","Event","sim","pred.pi","s0_hat","s1_hat","m_hat","I"))] 
  
  tryCatch({ 
    modeld.pre <- pre::pre(pseudo ~., data = data.complete,weights = weight.cen,tree.unbiased = tree.unbiased ,maxdepth = maxdepth,nfolds = nfolds.rulefit.cv,type=rulefit_type,normalize = normalize,use.grad = use.grad,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype),nlambda=20,type.measure ="mse",relax = glmnet.relax,gamma =glmnet.gamma,intercept=glmnet.intercept, removecomplements=removecomplements) # ,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype)
  },warning=function(e){
    modeld.pre <<- "error in ctree"
  },
  error=function(e){
    modeld.pre <<- "no candidate rules"
  })
  
  if (length(modeld.pre)==1){
    imps.pre.std <- "no models"
  }else{
    tryCatch({
      imps.pre.std <- pre::importance(modeld.pre, standardize = T, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.std <<- "all shrunk to 0"
    })

  }
  
  if (length(modeld.pre)==1){
    imps.pre.nonstd <- "no models"
  }else{
    tryCatch({
      imps.pre.nonstd <- pre::importance(modeld.pre, standardize = F, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.nonstd <<- "all shrunk to 0"
    })
    
  }

  if (length(modeld.pre)==1){
    pred.pre <- rep("no models",dim(testdat)[1])
  }else if(length(modeld.pre)>1 & length(imps.pre.std)==1){
    pred.pre <- rep("all shrunk to 0",dim(testdat)[1])
  }else{
    pred.pre <- predict(modeld.pre,newdata = testdat[,-which(names(testdat)%in%c("Treatment","Time","Event"))])
  }

  return(list(diff=pred,pre_result=modeld, final_ensemble_result=final.model.ensemble, num.initial.rule=dim(initial_rules_mat)[2],diff.pre=pred.pre,imp.pre.std=imps.pre.std,imp.pre.nonstd=imps.pre.nonstd,pre_result_fit=modeld.pre)) 
}

#------------------------------------------------------------ DEA-learner --------------------------------------------------------------#
rsf_HTE_D <- function(data=dat_all[which(dat_all$sim==1),],testdat=mydata$data,time.interest=11.8,k.folds=10,EA=F,noise.var=T,est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,response_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),tree.unbiased=TRUE,use.grad = TRUE,pre.ctree.control.alpha=0.05,pre.ctree.control.testtype='Bonferroni',rulefit_type="both",ntrees=500,maxdepth=4,normalize=TRUE,nfolds.rulefit.cv=10,removecomplements=FALSE,glmnet.alpha=1, glmnet.relax=F,glmnet.gamma = c(0, 0.25, 0.5, 0.75, 1),glmnet.intercept=FALSE,glmnet.linearterm=TRUE,glmnet.linearterm.X=paste0("V",c(1:10))){ # EA version is not coded yet, need to figure out mu(X)
  data<-data[complete.cases(data),]

  if (est_pi==1){
    
    if(noise.var==F){
      if (propensity_method=="correct_logistic"){ # not cv version
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
        data$pred.pi <- predict(rf,type='response')
      }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,-which(names(data)%in%c("Treatment","Time","Event","sim"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = 500,num.threads = 1)
        data$pred.pi<-predict(rf)$predictions[,2]
      }else if(propensity_method=="lasso_logistic"){
        rf <- cv.glmnet(y=data$Treatment, x=as.matrix(data[,which(names(data)%in%c(propensity_method_X))]),family = "binomial", type.measure = "class",alpha=1,nlambda = 20)
        data$pred.pi <- predict(rf,as.matrix(data[,which(names(data)%in%c(propensity_method_X))]), type = "response") # This predict() function is for glmnet
        
      }
      
    }else if(noise.var==T){
      if (propensity_method=="correct_logistic"){ # not cv version
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,-which(names(data)%in%c("Time","Event","sim"))])
        data$pred.pi <- predict(rf,type='response')
      }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,-which(names(data)%in%c("Treatment","Time","Event","sim"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = 500,num.threads = 1)
        data$pred.pi<-predict(rf)$predictions[,2]
      }else if(propensity_method=="lasso_logistic"){
        rf <- cv.glmnet(y=data$Treatment, x=as.matrix(data[,which(names(data)%in%c(propensity_method_X))]),family = "binomial", type.measure = "class",alpha=1,nlambda = 20)
        data$pred.pi <- predict(rf,as.matrix(data[,which(names(data)%in%c(propensity_method_X))]), type = "response") # This predict() function is for glmnet
        
      }
      
    }
    
  }else if(est_pi==0){
    data$pred.pi <- rep(propensity,dim(data)[1])
  }
  
  if(noise.var==F){
    rfsrc_data <-
      grf::survival_forest(X=as.matrix(data[,which(names(data) %in% response_X)]),Y=as.vector(data$Time),D=as.vector(data$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
  }else if(noise.var==T){
    rfsrc_data <-
      grf::survival_forest(X=as.matrix(data[,-which(names(data) %in% c("Treatment","Time","Event","sim","pred.pi"))]),Y=as.vector(data$Time),D=as.vector(data$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
  }
  
  data$m_hat <- predict(rfsrc_data,failure.times=rep(time.interest,dim(data)[1]),prediction.times="time")$predictions 
  
  data$I<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
  if (EA==F){
    data$pseudo <- as.vector(2*(2*data$Treatment-1)*data$I )
  }else if (EA==T){
    data$pseudo <- as.vector(2*(2*data$Treatment-1)*(data$I-data$m_hat))
  }
  
  U <- pmin(data$Time, time.interest)                         
  fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
  C.hat <- rep(NA, length(fold.id))
  for (z in 1:k.folds) {
    c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
    kmc <- summary(c.fit, times = U[fold.id == z])
    C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
  }
    
  weight <- (1/C.hat)*(2*data$Treatment-1)*((data$Treatment-data$pred.pi)/(4*data$pred.pi*(1-data$pred.pi)))
  weight.cen <- weight[!is.na(data$I)]
  
  data.complete <- data[!is.na(data$pseudo),-which(names(data)%in%c("Treatment","Time","Event","sim","pred.pi","m_hat","I"))] 

    tryCatch({ 
    modeld.pre <- pre::pre(pseudo ~., data = data.complete,weights = weight.cen,tree.unbiased = tree.unbiased ,maxdepth = maxdepth,nfolds = nfolds.rulefit.cv,type=rulefit_type,normalize = normalize,use.grad = use.grad,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype),nlambda=20,type.measure ="mse",relax = glmnet.relax,gamma =glmnet.gamma,intercept=glmnet.intercept, removecomplements=removecomplements) # ,tree.control=ctree_control(alpha = pre.ctree.control.alpha, testtype=pre.ctree.control.testtype)
  },warning=function(e){
    modeld.pre <<- "error in ctree"
  },
  error=function(e){
    modeld.pre <<- "no candidate rules"
  })
  
  if (length(modeld.pre)==1){
    imps.pre.std <- "no models"
  }else{
    tryCatch({
      imps.pre.std <- pre::importance(modeld.pre, standardize = T, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.std <<- "all shrunk to 0"
    })

  }
  
  if (length(modeld.pre)==1){
    imps.pre.nonstd <- "no models"
  }else{
    tryCatch({
      imps.pre.nonstd <- pre::importance(modeld.pre, standardize = F, round = 4L,plot = FALSE)
    },warning=function(w){
      imps.pre.nonstd <<- "all shrunk to 0"
    })
    
  }

  if (length(modeld.pre)==1){
    pred.pre <- rep("no models",dim(testdat)[1])
  }else if(length(modeld.pre)>1 & length(imps.pre.std)==1){
    pred.pre <- rep("all shrunk to 0",dim(testdat)[1])
  }else{
    pred.pre <- predict(modeld.pre,newdata = testdat[,-which(names(testdat)%in%c("Treatment","Time","Event"))])
  }

  return(list(diff=pred,pre_result=modeld, final_ensemble_result=final.model.ensemble, num.initial.rule=dim(initial_rules_mat)[2],diff.pre=pred.pre,imp.pre.std=imps.pre.std,imp.pre.nonstd=imps.pre.nonstd,pre_result_fit=modeld.pre)) 
}



