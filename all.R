library(MASS)
library(glmnet)
library(lpSolve)
library(ggplot2)




"gds" <-
  function(y,X,lambda=1,tol=0.001,maxit=10,family=binomial,rn=0,int=F){
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    # Add column to X for intercept.
    if (int)
      X <- cbind(1,X)
    weights <- 1
    p <- ncol(X)
    nobs <- length(y)
    # Set initial values for X*beta 
    if (family()$family!="gaussian")
      eval(family()$initialize)
    else mustart =y
    
    muzero <- mustart
    xbeta <- family()$linkfun(muzero)
    # Initialize the beta vector
    beta.old <- rep(2,p)
    beta.new <- rep(1,p)
    num.it <- 1
    # Start iterative weighted Dantzig Algorithm
    while(sum(abs(beta.new-beta.old))/(tol+sum(abs(beta.old)))>tol &
          num.it<=maxit){
      beta.old <- beta.new
      V <- family()$var(muzero)
      Z <- -xbeta-(y-muzero)/V
      # Zstar and Xstar correspond to the weighted new response and design
      # matrix. 
      Zstar <- sqrt(V)*Z
      Xstar <- sqrt(V)*X
      # Rescale the design matrix to have columns of norm one.
      colscale <- apply(X,2,function(x){sqrt(sum(x^2))})
      Xstar <- t(t(Xstar)/colscale)
      
      # linearoptim computes the linear programming solution for the new
      # beta estimate. 
      beta.new<-linearoptim(Zstar,Xstar,lambda,lambdap=sqrt(2*log(p)),rn)/colscale
      
      # Update X*beta and muhat.
      xbeta <- as.vector(X%*%beta.new)
      muzero <- family()$linkinv(xbeta)
      if (sum(abs(beta.new))==0) warning("Beta all zero")
      num.it <- num.it+1
    }
    list(beta=beta.new,iterations=num.it-1,L1=sum(abs(beta.new*colscale)),lambda=lambda)
  }


"linearoptim" <-
  function(Zstar,Xstar,lambda,lambdap,rn){
    p <- ncol(Xstar)
    f.obj <- rep(1,2*p)
    XX1 <- t(Xstar)%*%Xstar-rn
    XX2 <- t(Xstar)%*%Xstar+rn
    b1 <- lambdap*lambda-as.vector(t(Xstar)%*%Zstar)
    b2 <- lambdap*lambda+as.vector(t(Xstar)%*%Zstar)
    f.con1 <- cbind(XX1,-XX2)
    f.con2=cbind(-XX2,XX1)
    f.con <- rbind(f.con1,f.con2)
    f.rhs <- c(b1,b2)
    f.dir <- rep("<=",2*p)
    lpout <- lp("min",f.obj,f.con,f.dir,f.rhs)
    gamma <- lpout$solution[1:p]-lpout$solution[-(1:p)]
    if (lpout$status==2)
      print("Warning: No Solution Found")
    gamma}



#hcs.cv(ytrain,xtrain,lambda=0.4)


"hcs.cv" <-
  function(Y,X,lambda,cv=5,s=NULL,int=F,family=binomial,lambda2=.001,trace=T,rn=0){
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (int)
      X <- cbind(1,X)
    n <- length(Y)
    # Determine the randomization of the data.
    if (is.null(s))
      s <- sample(n,replace=F)
    lambdalen <- length(lambda)
    deviance=matrix(0,lambdalen,cv)
    error_class=matrix(0,lambdalen,cv)
    beta=matrix(0,lambdalen,ncol(X))
    l1=rep(0,lambdalen)
    for (j in 1:lambdalen){
      for (i in 1:cv){
        # Produce the left out portion of the data.
        test.index <- s[(1+(i-1)*floor(n/cv)):(i*floor(n/cv))]
        if (i==cv)
          test.index <- s[(1+(i-1)*floor(n/cv)):n]
        testY <- Y[test.index]
        testX <- X[test.index,]
        # Produce the training data.
        trainY <- Y[-test.index]
        trainX <- X[-test.index,]
        # Compute the DD estimates for beta on the training data.
        fit <-gds(trainY,trainX,lambda[j],int=F,family=family,rn=0,maxit=10)
        trainfit=fit$beta
        l1[j]=fit$L1
        
        muhat <- family()$linkinv(testX%*%trainfit)
        ypred_class=ifelse(muhat>=0.5,1,0)
        error_class[j,i]=mean((abs(ypred_class-testY)))
        # Compute the deviance on the left out portion of the data.
        deviance[j,i]=mean(family()$dev(testY,muhat,1))
        #error[j] <- error[j]+sum(family()$dev(testY,muhat,1))
        beta[j,]=trainfit
      }
      #   if (trace)
      #    print(paste("Finished lambda = ",lambda[j]))
    }
    # Print the optimal value of lambda.
    
    deviance1=apply(deviance,1,mean)
    deviance_sd=apply(deviance,1,sd)
    deviance_sd1=deviance_sd[which.min(deviance1)]
    
    
    error_class1=apply(error_class,1,mean)
    sd_class=apply(error_class,1,sd)
    sd_class1=sd_class[which.min(error_class1)]
    
    if (trace)
      print(paste("Min Deviance lambda = ",lambda[which.min(deviance1)],
                  "Min Classification lambda = ",lambda[which.min(error_class1)]))
    list(lambda=lambda,l1=l1,beta=beta,
         deviance=deviance,deviance1=deviance1,sdd=deviance_sd,sdd1=deviance_sd1,
         lambda1=lambda[which.min(deviance1)],beta1=beta[which.min(deviance1),],
         error=error_class,error1=error_class1,sd=sd_class,sd1=sd_class1,
         optimallambda=lambda[which.min(error_class1)],bestbeta=beta[which.min(error_class1),])
  }





"mhcs.cv" <-
  function(Y,X,lambda=optimallambda,cv=5,s=NULL,int=F,family=binomial,trace=T,rn){
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (int)
      X <- cbind(1,X)
    n <- length(Y)
    # Determine the randomization of the data.
    if (is.null(s))
      s <- sample(n,replace=F)
    lambdalen <- length(rn)
    #error <- rep(0,rnlen)
    
    deviance=matrix(0,lambdalen,cv)
    error_class=matrix(0,lambdalen,cv)
    beta=matrix(0,lambdalen,ncol(X))
    l1=rep(0,lambdalen)
    for (j in 1:lambdalen){
      for (i in 1:cv){
        # Produce the left out portion of the data.
        test.index <- s[(1+(i-1)*floor(n/cv)):(i*floor(n/cv))]
        if (i==cv)
          test.index <- s[(1+(i-1)*floor(n/cv)):n]
        testY <- Y[test.index]
        testX <- X[test.index,]
        # Produce the training data.
        trainY <- Y[-test.index]
        trainX <- X[-test.index,]
        # Compute the DD estimates for beta on the training data.
        fit <-gds(trainY,trainX,lambda,rn=rn[j],int=F,family=family,maxit=10)
        trainfit=fit$beta
        l1[j]=fit$L1
        
        muhat <- family()$linkinv(testX%*%trainfit)
        ypred_class=ifelse(muhat>=0.5,1,0)
        error_class[j,i]=mean((abs(ypred_class-testY)))
        # Compute the deviance on the left out portion of the data.
        deviance[j,i]=mean(family()$dev(testY,muhat,1))
        #error[j] <- error[j]+sum(family()$dev(testY,muhat,1))
        beta[j,]=trainfit
      }
      # if (trace)
      #   print(paste("Finished rn = ",rn[j]))
    }
    # Print the optimal value of lambda.
    
    deviance1=apply(deviance,1,mean)
    deviance_sd=apply(deviance,1,sd)
    deviance_sd1=deviance_sd[which.min(deviance1)]
    
    error_class1=apply(error_class,1,mean)
    sd_class=apply(error_class,1,sd)
    sd_class1=sd_class[which.min(error_class1)]
    
    if (trace)
      print(paste("Min Deviance rn = ",rn[which.min(deviance1)],
                  "Min Classification rn = ",rn[which.min(error_class1)]))
    list(rn=rn,l1=l1,beta=beta,
         deviance=deviance,deviance1=deviance1,sdd=deviance_sd,sdd1=deviance_sd1,
         rn1=rn[which.min(deviance1)],beta1=beta[which.min(deviance1),],
         error=error_class,error1=error_class1,sd=sd_class,sd1=sd_class1,
         optimalrn=rn[which.min(error_class1)],bestbeta=beta[which.min(error_class1),])
  }







family <- get("binomial", mode = "function", envir = parent.frame())

lambda=exp(seq(-6,0.5,0.15))
s0=10
b=1
p=200
rho=0
n=400
ratio=0.5
#ratio_noise=0.3
nrepeat=10
ratio_noise=c(0.1,0.3,0.5,1)
beta=c(rep(b,s0),rep(0,p-s0))


rho=0








"gendata1"=function(s0,p,rho,n,b){
  u=rep(0,p)
  beta=c(rep(b,s0),rep(0,p-s0))
  sigma=(1-rho)*diag(p)+matrix(rho,p,p)
  x=mvrnorm(n,u,sigma)
  pr = 1/(1+exp(-x%*%beta))  
  y = rbinom(n,1,pr) 
  list(x=x,y=y,df=data.frame(x,y))
}


'gendata4'=function(s0,p,rho,n,b){
  
  u=rep(0,p)
  beta=c(rep(b,s0),rep(0,p-s0))
  a=matrix(rep(1:p),p,p)
  sigma=rho^abs(a-t(a))
  x=mvrnorm(n,u,sigma)
  #x12=mvrnorm(n,u1,sigma)
  pr = 1/(1+exp(-x%*%beta))  
  y = rbinom(n,1,pr) 
  list(x=x,y=y,df=data.frame(x,y))
  
}


"genx"=function(df,ratio,ratio_noise,p){
  
  x=df$x
  y=df$y
  df=df$df
  sample= sample.int(n = nrow(df), size = floor(ratio*nrow(df)), replace = F)
  train <- df[sample, ]
  test  <- df[-sample, ]
  xtrain=as.matrix(train[,1:p])
  ytrain=as.vector(train[,p+1])
  xtest=as.matrix(test[,1:p])
  ytest=as.vector(test[,(p+1)])
  
  a=list()
  a$xtrain=xtrain
  a$xtest=xtest
  a$ytrain=ytrain
  a$ytest=ytest
  for (i in 1:length(ratio_noise)){
    sample_noise=sample(nrow(x)*ncol(x), nrow(x)*ncol(x)*ratio_noise[i])
    #xnoise
    xnoise=x
    addnoise=rnorm(nrow(x)*ncol(x)*ratio_noise[i],0,1)
    xnoise[sample_noise]=x[sample_noise]+addnoise
    xtrain_noise=as.matrix(xnoise[sample,])
    xtest_noise=as.matrix(xnoise[-sample,])
    a$xtrain_noise[[i]]=xtrain_noise
    a$xtest_noise[[i]]=xtest_noise
    #xmiss
    xmiss=x
    xmiss[sample_noise]=0
    xtrain_miss=as.matrix(xmiss[sample,])
    xtest_miss=as.matrix(xmiss[-sample,])
    a$xtrain_miss[[i]]=xtrain_miss
    a$xtest_miss[[i]]=xtest_miss
  }
  a
}





"genres"=function(nrepeat,ratio_noise,data){
  
  rlasso=list()
  rridge=list()
  rhcs=list()
  rmhcs=list()
  rlasso$non[[1]]=matrix(0,nrepeat,6)
  rridge$non[[1]]=matrix(0,nrepeat,6)
  rhcs$non[[1]]=matrix(0,nrepeat,6)
  rmhcs$non[[1]]=matrix(0,nrepeat,6)
  
  for (i in 1:(length(ratio_noise)-1)){
    rlasso$miss[[i]]=matrix(0,nrepeat,8)
    rridge$miss[[i]]=matrix(0,nrepeat,8)
    rhcs$miss[[i]]=matrix(0,nrepeat,8)
    rmhcs$miss[[i]]=matrix(0,nrepeat,8)
  }
  
  for (i in 1:length(ratio_noise)){
    rlasso$noise[[i]]=matrix(0,nrepeat,8)
    rridge$noise[[i]]=matrix(0,nrepeat,8)
    rhcs$noise[[i]]=matrix(0,nrepeat,8)
    rmhcs$noise[[i]]=matrix(0,nrepeat,8)
  }
  
  for (j in 1:nrepeat){
    
    print(j)
    Sys.sleep(0.01)
    
    #l=genx(s0,p,rho,n,ratio,ratio_noise)
    l=data[[j]]
    
    rlasso$non[[1]][j,]=result_lasso(xtrain=l$xtrain,xtest=l$xtest,ytrain=l$ytrain,ytest=l$ytest)$result
    
    rridge$non[[1]][j,]=result_ridge(xtrain=l$xtrain,xtest=l$xtest,ytrain=l$ytrain,ytest=l$ytest)$result
    
    f0=result_HCS_MHCS(xtrain=l$xtrain,xtest=l$xtest,ytrain=l$ytrain,ytest=l$ytest)
    rhcs$non[[1]][j,]=f0$HCS
    rmhcs$non[[1]][j,]=f0$MHCS
    
    
    for(i in 1:(length(ratio_noise)-1)){
      rlasso$miss[[i]][j,]=result_lasso3(xtrain=l$xtrain_miss[[i]],xtestz=l$xtest,xtestw=l$xtest_noise[[i]], ytrain=l$ytrain, ytest=l$ytest)
      rridge$miss[[i]][j,]=result_ridge3(xtrain=l$xtrain_miss[[i]],xtestz=l$xtest,xtestw=l$xtest_miss[[i]],ytrain=l$ytrain,ytest=l$ytest)
      
      f1=result_HCS_MHCS3(xtrain=l$xtrain_miss[[i]],xtestz=l$xtest,xtestw=l$xtest_miss[[i]],ytrain=l$ytrain,ytest=l$ytest)
      rhcs$miss[[i]][j,]=f1$HCS
      rmhcs$miss[[i]][j,]=f1$MHCS
    }
    
    for(i in 1:length(ratio_noise)){
      rlasso$noise[[i]][j,]=result_lasso3(xtrain=l$xtrain_noise[[i]],xtestz=l$xtest,xtestw=l$xtest_noise[[i]],ytrain=l$ytrain,ytest=l$ytest)
      rridge$noise[[i]][j,]=result_ridge3(xtrain=l$xtrain_noise[[i]],xtestz=l$xtest,xtestw=l$xtest_noise[[i]],ytrain=l$ytrain,ytest=l$ytest)
      
      
      f2=result_HCS_MHCS3(xtrain=l$xtrain_noise[[i]],xtestz=l$xtest,xtestw=l$xtest_noise[[i]],ytrain=l$ytrain,ytest=l$ytest)
      rhcs$noise[[i]][j,]=f2$HCS
      rmhcs$noise[[i]][j,]=f2$MHCS
    }
  }
  
  
  
  
  aa3=list(lasso=rlasso,ridge=rridge,hcs=rhcs,mhcs=rmhcs)
  a=list()
  a$non=cbind('lasso_mean'=apply(aa3$lasso$non[[1]],2,mean),'lasso_sd'=apply(aa3$lasso$non[[1]],2,sd),
              'ridge_mean'=apply(aa3$ridge$non[[1]],2,mean),'ridge_sd'=apply(aa3$ridge$non[[1]],2,sd),
              'hcs_mean'= apply(aa3$hcs$non[[1]],2,mean),'hcs_sd'=apply(aa3$hcs$non[[1]],2,sd),
              'mhcs_mean'=apply(aa3$mhcs$non[[1]],2,mean),'mhcs_sd'=apply(aa3$mhcs$non[[1]],2,sd))
  for (i in 1:(length(ratio_noise)-1)){
    a$miss[[i]]=cbind('lasso_mean'=apply(aa3$lasso$miss[[i]],2,mean),'lasso_sd'=apply(aa3$lasso$miss[[i]],2,sd),
                      'ridge_mean'=apply(aa3$ridge$miss[[i]],2,mean),'ridge_sd'=apply(aa3$ridge$miss[[i]],2,sd),
                      'hcs_mean'=apply(aa3$hcs$miss[[i]],2,mean),'hcs_sd'=apply(aa3$hcs$miss[[i]],2,sd),
                      'mhcs_mean'=apply(aa3$mhcs$miss[[i]],2,mean),'mhcs_sd'=apply(aa3$mhcs$miss[[i]],2,sd))
    }
  for (i in 1:length(ratio_noise)){
    a$noise[[i]]=cbind('lasso_mean'=apply(aa3$lasso$noise[[i]],2,mean),'lasso_sd'=apply(aa3$lasso$noise[[i]],2,sd),
                       'ridge_mean'=apply(aa3$ridge$noise[[i]],2,mean),'ridge_sd'=apply(aa3$ridge$noise[[i]],2,sd),
                       'hcs_mean'=apply(aa3$hcs$noise[[i]],2,mean),'hcs_sd'=apply(aa3$hcs$noise[[i]],2,sd),
                       'mhcs_mean'= apply(aa3$mhcs$noise[[i]],2,mean),'mhcs_sd'=apply(aa3$mhcs$noise[[i]],2,sd))
  }
  list(result_summary=a,result=aa3)
}





"l1"=function(x,beta=beta){
  beta_unit=beta/sqrt(sum(t(beta)%*%beta))
  l1=sum(abs(x/max(exp(-10),sqrt(sum(t(x)%*%x)))-beta_unit))
  l2=sum((x/max(exp(-10),sqrt(sum(t(x)%*%x)))-beta_unit)^2)
  c(l1,l2)}

"result_lasso"=function(xtrain,ytrain,xtest,ytest){
  #lasso result for non-measurement error data
  fit_lasso=cv.glmnet(xtrain,ytrain,intercept=FALSE,family="binomial")
  
  muhat_lasso=predict(fit_lasso,type="response",newx=xtest,s=fit_lasso$lambda.min)
  deviance_lasso=sum(family()$dev(ytest,muhat_lasso,1))
  deviance_lasso=deviance_lasso/((1-ratio)*n)
  
  ypred_lasso=predict(fit_lasso,type="class",newx=xtest,s=fit_lasso$lambda.min)
  error_lasso=sum(abs(ytest-as.numeric(ypred_lasso)))
  error_lasso=error_lasso/((1-ratio)*n)
  
  
  beta_LASSO=coef(fit_lasso,s=fit_lasso$lambda.min)[2:(p+1)]
  #nnzero(beta_LASSO)
  
  l1_lasso=l1(beta_LASSO,beta)[1]
  l2_lasso=l1(beta_LASSO,beta)[2]
  
  FN_LASSO=1-nnzero(beta_LASSO[1:s0])/s0
  FP_LASSO=nnzero(beta_LASSO[(1+s0):p])/(p-s0)
  list( result=c(deviance_lasso,error_lasso,l1_lasso,l2_lasso,FN_LASSO,FP_LASSO),fit=fit_lasso)
}




"result_lasso3"=function(xtrain,ytrain,xtestw,xtestz,ytest){
  r1=result_lasso(xtrain,ytrain,xtestz,ytest)
  fit_lasso=r1$fit
  r2=result_lasso2(xtestw,ytest,fit_lasso)
  c(r2,r1$result)
}

#result_ridge3(xtrain,ytrain,xtestw=xtestw,xtestz=xtest,ytest)

"result_lasso2"=function(xtest,ytest,fit_lasso){
  muhat_lasso=predict(fit_lasso,type="response",newx=xtest,s=fit_lasso$lambda.min)
  deviance_lasso=sum(family()$dev(ytest,muhat_lasso,1))
  deviance_lasso=deviance_lasso/((1-ratio)*n)
  
  ypred_lasso=predict(fit_lasso,type="class",newx=xtest,s=fit_lasso$lambda.min)
  error_lasso=sum(abs(ytest-as.numeric(ypred_lasso)))
  error_lasso=error_lasso/((1-ratio)*n)
  
  c(deviance_lasso,error_lasso)
}



"result_ridge3"=function(xtrain,ytrain,xtestw,xtestz,ytest){
  r1=result_ridge(xtrain,ytrain,xtestz,ytest)
  fit_lasso=r1$fit
  r2=result_lasso2(xtestw,ytest,fit_lasso)
  c(r2,r1$result)
}

"result_ridge"=function(xtrain,ytrain,xtest,ytest){
  #lasso result for non-measurement error data
  fit_lasso=cv.glmnet(xtrain,ytrain,intercept=FALSE,alpha=0,family="binomial")
  
  muhat_lasso=predict(fit_lasso,type="response",newx=xtest,s=fit_lasso$lambda.min)
  deviance_lasso=sum(family()$dev(ytest,muhat_lasso,1))
  deviance_lasso=deviance_lasso/((1-ratio)*n)
  
  ypred_lasso=predict(fit_lasso,type="class",newx=xtest,s=fit_lasso$lambda.min)
  error_lasso=sum(abs(ytest-as.numeric(ypred_lasso)))
  error_lasso=error_lasso/((1-ratio)*n)
  
  
  beta_LASSO=coef(fit_lasso,s=fit_lasso$lambda.min)[2:(p+1)]
  
  l1_lasso=l1(beta_LASSO,beta)[1]
  l2_lasso=l1(beta_LASSO,beta)[2]
  
  FN_LASSO=1-nnzero(beta_LASSO[1:s0])/s0
  FP_LASSO=nnzero(beta_LASSO[(1+s0):p])/(p-s0)
  list( result=c(deviance_lasso,error_lasso,l1_lasso,l2_lasso,FN_LASSO,FP_LASSO),fit=fit_lasso)
}






"result_HCS_MHCS"=function(xtrain,ytrain,xtest,ytest){
  # a1=seq(-5,0.5, 0.15)
  # lambda=exp(a1)
  
  #  rn=c(0,lambda)
  
  
  # hcs=hcs.cv(ytrain,xtrain,cv=5,rn=0,int=F,lambda=lambda)
  #  mhcs=mhcs.cv(ytrain,xtrain,cv=5,int=F,rn=c(0,lambda),lambda=hcs$optimallambda)
  
  # beta_HCS=hcs$bestbeta
  
  #fit_HCS=gds(ytrain,xtrain,lambda=hcs$optimallambda,rn=0,int=F,maxit=10)
  #beta_HCS=fit_HCS$beta
  
  
  fit_HCS=gds(ytrain,xtrain,lambda=0.02,rn=0,int=F,maxit=10)
  beta_HCS=fit_HCS$beta
  
  
  muhat_HCS=family()$linkinv(xtest%*%beta_HCS)
  ypred_HCS=ifelse(muhat_HCS>=0.5,1,0)
  error_HCS=sum((abs(ypred_HCS-ytest)))
  error_HCS=error_HCS/((1-ratio)*n)
  deviance_HCS=sum(family()$dev(ytest,muhat_HCS,1))
  deviance_HCS=deviance_HCS/((1-ratio)*n)
  
  
  #beta_MHCS=mhcs$bestbeta
  #fit_MHCS=gds(ytrain,xtrain,lambda=hcs$optimallambda,rn=mhcs$optimalrn,int=F, maxit=10)
  #beta_MHCS=fit_MHCS$beta
  
  fit_MHCS=gds(ytrain,xtrain,lambda=0.2,rn=0.004,int=F, maxit=10)
  beta_MHCS=fit_MHCS$beta
  muhat_MHCS=family()$linkinv(xtest%*%beta_MHCS)
  ypred_MHCS=ifelse(muhat_MHCS>=0.5,1,0)
  error_MHCS=sum((abs(ypred_MHCS-ytest)))
  error_MHCS=error_MHCS/((1-ratio)*n)
  deviance_MHCS=sum(family()$dev(ytest,muhat_MHCS,1))
  deviance_MHCS=deviance_MHCS/((1-ratio)*n)
  
  
  l1_hcs=l1(beta_HCS,beta)[1]
  l2_hcs=l1(beta_HCS,beta)[2]
  
  l1_mhcs=l1(beta_MHCS,beta)[1]
  l2_mhcs=l1(beta_MHCS,beta)[2]
  
  FN_HCS=1-nnzero(beta_HCS[1:s0])/s0
  FN_MHCS=1-nnzero(beta_MHCS[1:s0])/s0
  
  
  
  FP_HCS=nnzero(beta_HCS[(1+s0):p])/(p-s0)
  FP_MHCS=nnzero(beta_MHCS[(s0+1):p])/(p-s0)
  
  HCS=c(deviance_HCS,error_HCS,l1_hcs,l2_hcs,FN_HCS,FP_HCS)
  MHCS=c(deviance_MHCS,error_MHCS,l1_mhcs,l2_hcs,FN_MHCS,FP_MHCS)
  list(HCS=HCS,MHCS=MHCS,beta_HCS=beta_HCS,beta_MHCS=beta_MHCS)
  
}


"result_HCS_MHCS2"=function(beta_HCS,beta_MHCS,xtest,ytest){
  
  muhat_HCS=family()$linkinv(xtest%*%beta_HCS)
  ypred_HCS=ifelse(muhat_HCS>=0.5,1,0)
  error_HCS=sum((abs(ypred_HCS-ytest)))
  error_HCS=error_HCS/((1-ratio)*n)
  deviance_HCS=sum(family()$dev(ytest,muhat_HCS,1))
  deviance_HCS=deviance_HCS/((1-ratio)*n)
  
  
  muhat_MHCS=family()$linkinv(xtest%*%beta_MHCS)
  ypred_MHCS=ifelse(muhat_MHCS>=0.5,1,0)
  error_MHCS=sum((abs(ypred_MHCS-ytest)))
  error_MHCS=error_MHCS/((1-ratio)*n)
  deviance_MHCS=sum(family()$dev(ytest,muhat_MHCS,1))
  deviance_MHCS=deviance_MHCS/((1-ratio)*n)
  
  HCS=c(deviance_HCS,error_HCS)
  MHCS=c(deviance_MHCS,error_MHCS)
  list(HCS=HCS,MHCS=MHCS)}

"result_HCS_MHCS3"=function(xtrain,ytrain,xtestz,xtestw,ytest){
  r1=result_HCS_MHCS(xtrain,ytrain,xtestz,ytest)
  r2=result_HCS_MHCS2(r1$beta_HCS,r1$beta_MHCS,xtestw,ytest)
  list(HCS=c(r2$HCS,r1$HCS),MHCS=c(r2$MHCS,r1$MHCS))
}











"route_lambda1"=function(fit,nam='',name1=''){
  lambda=fit$lambda
  # diag(eg0$values)
  tu=fit$beta# %*%sqrt(diag(eg3$values))%*%eg3$vectors
  xrange <- range(lambda) 
  yrange <- range(tu) 
  color=rainbow(200)
  #mypath <- file.path("Desktop","MATHLatex","images",paste("Beta Route by lambda ", nam, ".png", sep = ""))
  #png(filename=mypath,width=7,height=5,units='in',res=300)
  plot(xrange, yrange, type="n", xlab="log(lambda)",
       ylab="beta",log = "x" ,col = "#2E9FDF") 
  for (i in (s0+1):ncol(tu)){
    # for (i in c(1:ncol(tu))){#[-rr[1:5]]
    lines((lambda),tu[,i],col='lightgrey')}
  for (i in 1:s0){
    # for (i in c(rr3[1:5])){
    lines((lambda),tu[,i],col=i)
  }
  #lines(lambda,tu[,i],col=color[i])
  
  #}
  abline(v =c(f[[1]]$hcs$lambda1,fit$lambda1) , col=c("black","blue"), lwd=1, lty=2)
 # abline(v =c(f[[1]]$hcs$lambda1,fit$lambda1,fit$optimallambda) , col=c("black","blue",'red'), lwd=1, lty=2)
  legend("topright", legend=c('f0-lambda*','Deviance-lambda*'),col=c("black", "blue"),box.lty=0,lty=c(2,2),cex=0.8)
  #title(paste("Beta Route by lambda",name1))
  #dev.off()
}








#fit=fset1[[1]]$mhcs

"route_gamma1"=function(fit,nam='',name1=''){
  rn=fit$rn
  rn[1]=0.5*rn[2]
  #tu=fit$beta 
  tu=fit$beta#%*%diag(eg2$values)%*%eg2$vectors
  xrange <- range(rn) 
  yrange <- range(tu) 
  color=rainbow(200)
  
  # mypath <- file.path("Desktop","MATHLatex","images",paste("Beta Route by Gamma ", nam, ".png", sep = ""))
  #png(filename=mypath,width=7,height=5,units='in',res=300)
  
  
  plot(xrange, yrange, type="n", xlab="log(Gamma)",
       ylab="beta",log = "x" ) 
  for (i in 1:ncol(tu)){lines((rn),tu[,i],col='lightgrey')}
  for (i in 1:10){lines((rn),tu[,i],col=i)}
  abline(v =c(max(fit$optimalrn,rn[1]),max(fit$rn1,rn[1])) , col=c("red","black"), lwd=1, lty=2)
  legend("topright", legend=c('CE-gamma*','Deviance-gamma*'),col=c("red", "black"),box.lty=0,lty=c(2,2),cex=0.8)
  #title(paste("Beta Route by Gamma",name1))
  #dev.off()
}












"cvplot_lambda1"=function(fit){
  deviance1=fit$deviance1
  # error2=fit$error1
  sdd=(fit$sdd)
  sdd1=(fit$sdd1)
  # sd=(fit$sd)/5
  # sd1=(fit$sd1)/5
  lambda1=fit$lambda1
  #lambda2=fit$optimallambda
  
  df=data.frame(fit$lambda,deviance1,sdd)
  #mypath <- file.path("Desktop","MATHLatex","images",paste("Cross-Validation of HCS ", nam, ".png", sep = ""))
  
  ggplot(df, aes(x=log(fit$lambda), y=deviance1))+ 
    geom_line(col="lightgrey",linetype=2) +
    # ylim(0.5,5)+
    geom_errorbar(aes(ymin=(deviance1-sdd), ymax=(deviance1+sdd)), col="pink",width=.05,
                  position=position_dodge(1),size=1.2)+
    geom_point(col="grey")+
    labs(x="log(lambda)", y = "Cross Entropy")+
    geom_vline(xintercept = c(log(lambda1),log(lambda1+sdd1)),col=c("black","blue"),linetype=2)+
    ggtitle("Cross-Validation of HCS")
  #ggsave(temp, file=mypath, width = 7, height = 5, units = "in")
}













"cvplot_lambda2"=function(fit){
  error1=fit$error1
  sd=(fit$sd)/5
  sd1=(fit$sd1)/5
  
  lambda2=fit$optimallambda
  
  df2=data.frame(fit$lambda,error1,sd)
  #mypath <- file.path("Desktop","MATHLatex","images",paste("Cross-Validation of HCS ", nam, ".png", sep = ""))
  ggplot(df2, aes(x=log(fit$lambda), y=error1))+ 
    geom_line(col="lightgrey",linetype=2) +
    # ylim(0.5,5)+
    geom_errorbar(aes(ymin=(error1-sd), ymax=(error1+sd)), col="pink",width=.05,
                  position=position_dodge(1),size=1.2)+
    geom_point(col="grey")+
    labs(x="log(lambda)", y = "Classification Error")+
    geom_vline(xintercept = c(log(lambda2),log(lambda2+sd1)),col=c("black","blue"),linetype=2)+
    ggtitle("Cross-Validation of HCS")
  #ggsave(temp, file=mypath, width = 7, height = 5, units = "in")
}












"cvplot_gamma1"=function(fit){
  error1=fit$error1
  sd=(fit$sd)^2
  rn=fit$rn
  rn[1]=rn[2]*0.5
  df=data.frame(rn,error1,sd)
  
  # mypath <- file.path("Desktop","MATHLatex","images",paste("Cross-Validation of MHCS ", nam, ".png", sep = ""))
  
  
  ggplot(df, aes(x=log(rn), y=error1))+ 
    geom_line(col="lightgrey",linetype=2) +
    geom_errorbar(aes(ymin=error1-sd, ymax=error1+sd), col="pink",width=.05,
                  position=position_dodge(1),size=1.2)+
    geom_point(col="grey")+
    #ylim(0.5,5)+
    labs(x="log(gamma)", y = "Classification Error")+
    geom_vline(xintercept = max(log(fit$optimalrn),log(rn[1])),col="black",linetype=2)+
    ggtitle("Cross-Validation of MHCS ")
  #ggsave(temp, file=mypath, width = 7, height = 5, units = "in")
}




't0'=function(a){
  a
  a[,2] <- paste0('&(', a[,2],')&')
  a[,4]=paste0('&(', a[,4],')&')
  a[,6]=paste0('&(', a[,6],')&')
  a[,8]=paste0('&(', a[,8],')\\')
  
  a0=c(#'\textbf{$Deviance(W_{test})$}&',
    #'\textbf{$CE(W_{test})$}&',
    '\textbf{$Deviance$}&',
    '\textbf{$CE$}&',
    '\textbf{$L_1$}&',
    '\textbf{$L2$}&',
    '\textbf{$FP$}&',
    '\textbf{$FN$}&')
  a=cbind(a0,a)
  a
  ax2=c(a[2,],a[1,],a[3,],a[4,],a[5,],a[6,])
  i=1
  ax3=''
  for (i in 1:length(ax2)){
    ax3=paste0(ax3,ax2[i])
    
    i=i+1}
  ax3}


't1'=function(a){
  a
  a[,2] <- paste0('&(', a[,2],')&')
  a[,4]=paste0('&(', a[,4],')&')
  a[,6]=paste0('&(', a[,6],')&')
  a[,8]=paste0('&(', a[,8],')\\')
  
  a0=c('\textbf{$Deviance(W_{test})$}&',
       '\textbf{$CE(W_{test})$}&',
       '\textbf{$Deviance(Z_{test})$}&',
       '\textbf{$CE(Z_{test})$}&',
       '\textbf{$L_1$}&',
       '\textbf{$L2$}&',
       '\textbf{$FP$}&',
       '\textbf{$FN$}&')
  a=cbind(a0,a)
  a
  ax2=c(a[4,],a[3,],a[2,],a[1,],a[5,],a[6,],a[7,],a[8,])
  i=1
  ax3=''
  for (i in 1:length(ax2)){
    ax3=paste0(ax3,ax2[i])
    
    i=i+1}
  
  ax3
  
}

'table'=function(r){
  
  rmiss=list()
  for (i in 1:3){
    a=round(r$result_summary$miss[[i]],2)
    rmiss[[i]]=t1(a)}
  
  rnoise=list()
  for (i in 1:4){
    a=round(r$result_summary$noise[[i]],2)
    rnoise[[i]]=t1(a)}
  
  a=round(r$result_summary$non,2)
  rnon=t0(a)
  
  list(rnon=rnon,rmiss=rmiss,rnoise=rnoise)
}

