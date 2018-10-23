

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
  
  
  fit_HCS=gds(ytrain,xtrain,lambda=0.05,rn=0,int=F,maxit=10)
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
  
  fit_MHCS=gds(ytrain,xtrain,lambda=0.05,rn=0.005,int=F, maxit=10)
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
  
#r1=result_HCS_MHCS(data1$xtrain,data1$ytrain,data1$xtest,data1$ytest)
#result_HCS_MHCS3(xtrain,ytrain,xtestz=xtest,xtestw,ytest)

