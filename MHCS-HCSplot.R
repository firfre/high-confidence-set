
#route_lambda(f1)

"route_lambda"=function(fit,nam='',name1='',Type){
  lambda=fit$lambda
  tu=fit$beta 
  xrange <- range(lambda) 
  yrange <- range(tu) 
  color=rainbow(s0)
  mypath <- file.path("Desktop","MATHLatex","images",paste0(Type,"lb", nam, ".png", sep = ""))
  png(filename=mypath,width=7,height=5,units='in',res=300)
  plot(xrange, yrange, type="n", xlab="log(lambda)",
       ylab="beta",col = "#2E9FDF",log = "x")
  for (i in (s0+1):ncol(tu)){
    lines(lambda,tu[,i],col='lightgrey')}
  for (i in 1:s0){
    lines(lambda,tu[,i],col=color[i])}
  abline(v =c(fit$optimallambda,fit$lambda1 ), col=c("blue",'red'), lwd=1, lty=2)
  title(paste("Beta Route by lambda",name1,Type))
  legend('topright', legend=c('lambda*_CE','lambda*_Deviance'),lty=c(2,2),
         col=c("blue", "red"),box.lty=0,cex=0.8)
  dev.off()
}



"route_gamma"=function(fit,nam,name1,Type){
  rn=fit$rn
  rn[1]=0.5*rn[2]
  tu=fit$beta 
  xrange <- range(rn) 
  yrange <- range(tu) 
  color=rainbow(s0)
  
  mypath <- file.path("Desktop","MATHLatex","images",paste0(Type,"gb", nam, ".png", sep = ""))
  png(filename=mypath,width=7,height=5,units='in',res=300)
  
  
  plot(xrange, yrange, type="n", xlab="log(Gamma)",
       ylab="beta",log = "x" ) 
  for (i in (s0+1):ncol(tu)){
    lines(rn,tu[,i],col='lightgrey')}
  for (i in 1:s0){
    lines(rn,tu[,i],col=color[i])}
  abline(v =c(max(fit$optimalrn,rn[1]),max(fit$rn1,rn[1])) , col=c("blue",'red'), lwd=1, lty=2)
  title(paste("Beta Route by Gamma",name1,Type))
  legend('topright', legend=c('gamma*_CE','gamma*_Deviance'),lty=c(2,2),
         col=c("blue", "red"),box.lty=0,cex=0.8)
  dev.off()
}



"cvplot_gamma"=function(fit,nam,name1,Type){
  error1=fit$error1
  sd=(fit$sd)/5
  sd1=(fit$sd1)/5
  rn=fit$rn
  rn[1]=rn[2]*0.5
  df=data.frame(rn,error1,sd)
  
   mypath1 <- file.path("Desktop","MATHLatex","images",paste0(Type,"gc", nam, ".png", sep = ""))
   
  
    temp=ggplot(df, aes(x=log(rn), y=error1))+ 
    geom_line(col="lightgrey",linetype=2) +
    geom_errorbar(aes(ymin=error1-sd, ymax=error1+sd), col="pink",width=.05,
                  position=position_dodge(1),size=1.2)+
    geom_point(col="grey")+
    #ylim(0.5,5)+
    labs(x="log(gamma)", y = "Cross Entropy")+
    geom_vline(xintercept = c(max(log(fit$optimalrn),log(rn[1])),log(fit$optimalrn+sd1)),col=c("black","blue"),linetype=2)+
    ggtitle(paste("Cross-Validation of MHCS ",name1,Type))
    ggsave(temp, file=mypath1, width = 7, height = 5, units = "in")
    
    
    
    dev.off()
    
    
    deviance1=fit$deviance1
    # error2=fit$error1
    sdd=(fit$sdd)/5
    sdd1=(fit$sdd1)/5
    # sd=(fit$sd)/5
    # sd1=(fit$sd1)/5
    rn1=fit$rn1
    #lambda2=fit$optimallambda
    
    df=data.frame(fit$rn,deviance1,sdd)
    mypath2 <- file.path("Desktop","MATHLatex","images",paste0(Type, "gd", nam, ".png", sep = ""))
    
    ggplot(df, aes(x=log(fit$rn), y=deviance1))+ 
      geom_line(col="lightgrey",linetype=2) +
      # ylim(0.5,5)+
      geom_errorbar(aes(ymin=(deviance1-sdd), ymax=(deviance1+sdd)), col="pink",width=.05,
                    position=position_dodge(1),size=1.2)+
      geom_point(col="grey")+
      labs(x="log(rn)", y = "Cross Entropy")+
      geom_vline(xintercept = c(log(rn1),log(rn1+sdd1)),col=c("black","blue"),linetype=2)+
      ggtitle(paste("Cross Entropy Cross-Validation of MHCS",name1,Type))
    ggsave(temp, file=mypath2, width = 7, height = 5, units = "in")
    dev.off()
  }


"cvplot_lambda"=function(fit,nam,name1,Type){
  error1=fit$error1
  sd=(fit$sd)/5
  sd1=fit$sd1/5
  lambda=fit$lambda
  df=data.frame(lambda,error1,sd)
  
  mypath1 <- file.path("Desktop","MATHLatex","images",paste0(Type,"lc", nam, ".png", sep = ""))
  
  temp=ggplot(df, aes(x=log(lambda), y=error1))+ 
    geom_line(col="lightgrey",linetype=2) +
    # ylim(0.5,5)+
    geom_errorbar(aes(ymin=(error1-sd), ymax=(error1+sd)), col="pink",width=.05,
                  position=position_dodge(1),size=1.2)+
    geom_point(col="grey")+
    labs(x="log(lambda)", y = "Classification Error")+
    geom_vline(xintercept = c(log(fit$optimallambda),log(fit$optimallambda+sd1)),col=c("black","blue"),linetype=2)+
    ggtitle(paste("Classification Error Cross-Validation of HCS",name1,Type))
    ggsave(temp, file=mypath1, width = 7, height = 5, units = "in")
    
    dev.off()
    
    
    deviance1=fit$deviance1
    # error2=fit$error1
    sdd=(fit$sdd)/5
    sdd1=(fit$sdd1)/5
    # sd=(fit$sd)/5
    # sd1=(fit$sd1)/5
    lambda1=fit$lambda1
    #lambda2=fit$optimallambda
    
    df=data.frame(fit$lambda,deviance1,sdd)
    mypath2 <- file.path("Desktop","MATHLatex","images",paste0(Type, "ld", nam, ".png", sep = ""))
    
    ggplot(df, aes(x=log(fit$lambda), y=deviance1))+ 
      geom_line(col="lightgrey",linetype=2) +
      # ylim(0.5,5)+
      geom_errorbar(aes(ymin=(deviance1-sdd), ymax=(deviance1+sdd)), col="pink",width=.05,
                    position=position_dodge(1),size=1.2)+
      geom_point(col="grey")+
      labs(x="log(lambda)", y = "Cross Entropy")+
      geom_vline(xintercept = c(log(lambda1),log(lambda1+sdd1)),col=c("black","blue"),linetype=2)+
      ggtitle(paste("Cross Entropy Cross-Validation of HCS",name1,Type))
    ggsave(temp, file=mypath2, width = 7, height = 5, units = "in")
    dev.off()
    }













quartz.save("Desktop/MATHLatex/images/test.png")



name0=c('0','10m','30m','50m',
        '10p','30p','50p')
name=c('without Measurement Error','with 10% Missing Value','with 30% Missing Value','with 50% Missing Value',
       'with 10% Perturbation','with 30% Perturbation','with 50% Perturbation')





  
  
  
  
  


cvplot_lambda(f1)
cvplot_gamma(f2)
x=list()
data=data2[[1]]

x[[1]]=data$xtrain
x[[2]]=data$xtrain_miss[[1]]
x[[3]]=data$xtrain_miss[[2]]
x[[4]]=data$xtrain_miss[[3]]
x[[5]]=data$xtrain_noise[[1]]
x[[6]]=data$xtrain_noise[[2]]
x[[7]]=data$xtrain_noise[[3]]
ytrain=data$ytrain


fset4=list()

lambda=exp(seq(-5.5,0.5,0.15))
for (i in c(1:7)){f1=hcs.cv(Y=ytrain,X=x[[i]],rn=0,lambda=lambda,cv=5)
  f2=mhcs.cv(int=F,Y=ytrain,X=x[[i]],rn=c(0,lambda),lambda=f1$lambda1,cv=5)
  fset4[[i]]=list(hcs=f1,mhcs=f2)
 }

Type="t3"

for (i in 1:7){
  f0=f[[i]]$mhcs
  f00=f[[i]]$hcs
  name1=name[i]
  nam=name0[i]
  route_lambda(f00,nam,name1,Type)
  route_gamma(f0,nam,name1,Type)
  
  cvplot_lambda(f00,nam,name1,Type)
  cvplot_gamma(f0,nam,name1,Type)
  }





