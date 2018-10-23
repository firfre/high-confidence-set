

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
    
   
  
                 
   




















"hcs0.cv" <-
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
    error =matrix(0,lambdalen,cv)
    beta=matrix(0,lambdalen,ncol(X))
    l1=rep(0,lambdalen)
    for (j in 1:lambdalen){
      for (i in 1:cv){
        # Produce the left out portion of the data.
        test.index <-  s[(1+(i-1)*floor(n/cv)):(i*floor(n/cv))]
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
        
        # Compute the deviance on the left out portion of the data.
        muhat <- family()$linkinv(testX%*%trainfit)
        error[j,i]=mean(family()$dev(testY,muhat,1))
        #error[j] <- error[j]+sum(family()$dev(testY,muhat,1))
        beta[j,]=trainfit
      }
      if (trace)
        print(paste("Finished lambda = ",lambda[j]))
    }
    # Print the optimal value of lambda.
    
    error1=apply(error,1,mean)
    sd=apply(error,1,sd)
    sd1=sd[which.min(error1)]
    if (trace)
      print(paste("Min lambda = ",lambda[which.min(error1)]))
    list(lambda=lambda,l1=l1,beta=beta,
         error=error,error1=error1,sd=sd,sd1=sd1,
         optimallambda=lambda[which.min(error1)],bestbeta=beta[which.min(error1),])}

"gds2" <-
  function(y,X,lambda=1,tol=0.001,maxit=50,family=binomial,rn=0,int=F){
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    # Add column to X for intercept.
    if (int)
      X <- cbind(1,X)
    X=xtrain
    y=ytrain
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
    while(
      sum(abs(beta.new-beta.old))/(tol+sum(abs(beta.old)))
      >tol &
      num.it<=maxit){
      beta.old <- beta.new
      V <- family()$var(muzero)
      Z <- xbeta+(y-muzero)/V
      # Zstar and Xstar correspond to the weighted new response and design
      # matrix. 
      Zstar <- sqrt(V)*Z
      Xstar <- sqrt(V)*X
      # Rescale the design matrix to have columns of norm one.
      colscale <- apply(X,2,function(x){sqrt(sum(x^2))})
      Xstar <- t(t(Xstar)/colscale)
      
      # linearoptim computes the linear programming solution for the new
      # beta estimate. 
      
      a=dantzig(Xstar,Zstar,lambda=0.0001,nlambda=100)
      #lambda=0.1
      #lambdap=sqrt(2*log(p))
      beta.new =dantzig.selector(a$lambdalist,a$BETA0,lambda*lambdap)/colscale
      #a$lambdalist[length(a$lambdalist)]
      # Update X*beta and muhat.
      xbeta <- as.vector(X%*%beta.new)
      muzero <- family()$linkinv(xbeta)
      nnzero(beta.new)
      if (sum(abs(beta.new))==0) warning("Beta all zero")
      num.it <- num.it+1
    }
    nnzero(beta.new)
    beta.new[ which(!beta.new == 0)]
    list(beta=beta.new,iterations=num.it-1,L1=sum(abs(beta.new*colscale)),lambda=lambda)
  }


gds(y=d3$ytrain,X=d3$xtrain,lambda=0.4,tol=0.001,maxit=10,family=binomial,rn=0,int=F)
