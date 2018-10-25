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
    lines(lambda,tu[,i],col='lightgrey')}
  for (i in 1:s0){
    # for (i in c(rr3[1:5])){
    lines(lambda,tu[,i],col=i)
  }
  #lines(lambda,tu[,i],col=color[i])
  
  #}
  abline(v =c(fit$lambda1,fit$optimallambda) , col=c("black",'red'), lwd=1, lty=2)
  legend("topright", legend=c('CE-lambda*','Deviance-lambda*'),col=c("red", "black"),box.lty=0,lty=c(2,2),cex=0.8)
  #title(paste("Beta Route by lambda",name1))
  #dev.off()
}










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
       ylab="eigen vector",log = "x" ) 
  for (i in 1:ncol(tu)){lines(rn,tu[,i],col='lightgrey')}
  for (i in 1:s0){
    lines(rn,tu[,i],col=color[i])}
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




