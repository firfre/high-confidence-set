n=400
p=200
lambda=exp(seq(-6,0.5,0.15))


s0=10
b=1

ratio=0.5
#ratio_noise=0.3
nrepeat=10
ratio_noise=c(0.1,0.3,0.5,1)
beta=c(rep(b,s0),rep(0,p-s0))

#generate data set


data1=list()
df=list()
#set.seed(100)
for (j in 1:nrepeat){df[[j]]=gendata1(n=n,p=p,s0=s0,rho=0,b=b)
data1[[j]]=genx(df=df[[j]],p=p,ratio=ratio,ratio_noise=ratio_noise)}


nrepeat=10
data4=list()
df=list()
set.seed(108)
for (j in 1:nrepeat){df[[j]]=gendata4(n=200,p=200,s0=s0,rho=0.5,b=1)
data4[[j]]=genx(df=df[[j]],p=200,ratio=0.5,ratio_noise=ratio_noise)}

nrepeat=10
data2=list()
df=list()
set.seed(222)
for (j in 1:nrepeat){df[[j]]=gendata1(n=200,p=200,s0=s0,rho=0.5,b=1)
data2[[j]]=genx(df=df[[j]],p=200,ratio=0.5,ratio_noise=ratio_noise)}










r1=genres(nrepeat=nrepeat,ratio_noise=ratio_noise,data=data1)
r1$result_summary



round(r1$result_summary$noise[[4]],2)

r1$result_summary




table(r1) 


print.noquote(table(r1) ) 








data=data1[[1]]
x=list()
x[[1]]=data$xtrain

x[[2]]=data$xtrain_miss[[1]]
x[[3]]=data$xtrain_miss[[2]]
x[[4]]=data$xtrain_miss[[3]]
#x[[5]]=data$xtrain_miss[[4]]

x[[5]]=data$xtrain_noise[[1]]
x[[6]]=data$xtrain_noise[[2]]
x[[7]]=data$xtrain_noise[[3]]
x[[8]]=data$xtrain_noise[[4]]
ytrain=data$ytrain



fset1=list()

f0=hcs.cv(Y=ytrain,X=x[[1]],rn=0,lambda=lambda,cv=5)


for (i in 1:8){
  
 # f=mhcs.cv(int=F,Y=ytrain,X=x[[i]],rn=c(0,lambda),lambda=f1$lambda1,cv=5)
  
  f1=hcs.cv(Y=ytrain,X=x[[i]],rn=0,lambda=lambda,cv=5)
  f2=mhcs.cv(int=F,Y=ytrain,X=x[[i]],rn=c(0,lambda),lambda=fset1[[1]]$hcs$lambda1,cv=5)
  fset1[[i]]=list(hcs=f1,mhcs=f2)
}

f=fset1
for (i in 1:8){
  f0=f[[i]]$mhcs
  f00=f[[i]]$hcs
  route_lambda1(f00)
  route_gamma1(f0)
  
  #cvplot_lambda1(f00)
  #cvplot_lambda2(f00)
  #cvplot_gamma1(f0)
}












fit=fset1[[7]]$hcs
cvplot_lambda1(fit)
for (i in 1:7){fit=fset1[[i]]$hcs
cvplot_lambda1(fit)}

for (i in 1:7)
{cvplot_lambda2(fset1[[i]]$hcs)}

for (i in 1:7)
{cvplot_gamma1(fset1[[i]]$mhcs)}