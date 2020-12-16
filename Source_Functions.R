##########################################################################
###### Generate the design matrix
# simulate A1 of the matrix (rho)^|i-j|
# A1 is p*p precision matrix
##########################################################################
A1gen<-function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}
##########################################################################
# simulate the matrix with the form A2[i,j]=(abs(i-j)==0)+ 0.4*(abs(i-j)==1)+0.2*(abs(i-j)==2)+ 0.2*(abs(i-j)==3)+ 0.1*(abs(i-j)==4)
##########################################################################
A2gen<-function(p){
  A2=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A2[i,j]<-(abs(i-j)==0)+ 0.4*(abs(i-j)==1)+0.2*(abs(i-j)==2)+ 0.2*(abs(i-j)==3)+ 0.1*(abs(i-j)==4)
    }
  }
  A2
}
##########################################################################
###### The following function computes the lasso estimator
# Compute the Lasso estimator:
# - If lambda is given, use glmnet and standard Lasso
# - If lambda is set to the character string "CV", then glmnet with
#   lambda selected by cross-validation is used
# - If lambda is not given or is set to NULL, use square root Lasso
##########################################################################
Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (is.null(lambda)) {
    lambda <- sqrt(qnorm(1 - (0.1 / p)) / n)
    outLas <- slim(X, y, lambda = lambda, method = "lq", q = 2,
                   verbose = FALSE)
    # Objective : sqrt(RSS/n) + lambda * penalty
    c(as.vector(outLas$intercept), as.vector(outLas$beta))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "scalreg") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else {
    outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }
  
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

##########################################################################
###### The following function is to compute weight vector for all dimensions
##########################################################################
opt.weight<-function(Gamma,delta){
  L<-dim(Gamma)[2]
  opt.weight<-rep(NA, L)
  opt.reward<-NA
  # Problem definition
  v<-Variable(L)
  Diag.matrix<-diag(eigen(Gamma)$values)
  for(ind in 1:L){
    Diag.matrix[ind,ind]<-max(Diag.matrix[ind,ind],0.001)
  }
  Gamma.positive<-eigen(Gamma)$vectors%*%Diag.matrix%*%t(eigen(Gamma)$vectors)
  objective <- Minimize(quad_form(v,Gamma.positive+diag(delta,L)))
  constraints <- list(v >= 0, sum(v)== 1)
  prob.weight<- Problem(objective, constraints)
  if(is_dcp(prob.weight)){
    # Problem solution
    result<- solve(prob.weight)
    opt.status<-result$status
    opt.sol<-result$getValue(v)
    #opt.reward<-t(opt.sol)%*%Gamma%*%opt.sol
    for(l in 1:L){
      opt.weight[l]<-opt.sol[l]*(abs(opt.sol[l])>10^{-8})
    }
  }
  v<-Variable(L)
  objective<-Minimize(2*t(v)%*%Gamma.positive%*%opt.weight-t(opt.weight)%*%Gamma.positive%*%opt.weight)
  constraints<-list(v >= 0, sum(v)== 1)
  delta.optim<-Problem(objective, constraints)
  result<- solve(delta.optim)
  opt.reward<-result$value
  returnList <- list("weight" = opt.weight,
                     "reward" = opt.reward
  )
  return(returnList)
}

##########################################################################
###### The following function maps the index of lower triangular matrix to its vectorized version
index.map<-function(L,l,k){
  return((2*L-k)*(k-1)/2+l)
}
##########################################################################



cov.inner<-function(Var.vec,Pred.mat.total,l1,k1,l2,k2){
  n.all<-dim(Pred.mat.total)[1]
  index.set.l1<-seq((l1-1)*n+1,l1*n)
  index.set.k1<-seq((k1-1)*n+1,k1*n)
  index.set.l2<-seq((l2-1)*n+1,l2*n)
  index.set.k2<-seq((k2-1)*n+1,k2*n)
  var1<-Var.vec[l1]*mean(Pred.mat.total[index.set.l1,k1]*Pred.mat.total[index.set.l2,k2])*(l2==l1)
  var1<-var1+Var.vec[l1]*mean(Pred.mat.total[index.set.l1,k1]*Pred.mat.total[index.set.k2,l2])*(k2==l1)
  var1<-var1+Var.vec[k1]*mean(Pred.mat.total[index.set.k1,l1]*Pred.mat.total[index.set.l2,k2])*(l2==k1)
  var1<-var1+Var.vec[k1]*mean(Pred.mat.total[index.set.k1,l1]*Pred.mat.total[index.set.k2,l2])*(k2==k1)
  var2<-mean((diag(Pred.mat.total[,k1]%*%t(Pred.mat.total[,l1]))-mean(Pred.mat.total[,k1]*Pred.mat.total[,l1]))*(diag(Pred.mat.total[,k2]%*%t(Pred.mat.total[,l2]))-mean(Pred.mat.total[,k2]*Pred.mat.total[,l2])))
  var<-var1/n+var2/n.all
  return((var))
}



################## these are the codes for the covariate shift setting
################################################################################################
############### compute the projection direction for estimating the regression covariance matrix
################################################################################################
proj.direction<-function(Xc,loading,maxiter=6,resol=1.25){
  n<-dim(Xc)[1]
  p<-dim(Xc)[2]
  loading.norm<-sqrt(sum(loading^2))
  sigma.hat <- (1/n)*(t(Xc)%*%Xc);
  if ((n>=6*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  if ((n>=6*p)&&(tmp>=1e-4)){
    direction <- solve(sigma.hat)%*%loading
  }else{
    #    if(n>0.5*p){
    step.vec<-rep(NA,3)
    for(t in 1:3){
      index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
      Direction.Est.temp<-Direction_searchtuning_lin(Xc[index.sel,],loading,mu=NULL, resol, maxiter)
      step.vec[t]<-Direction.Est.temp$step
    }
    step<-getmode(step.vec)
    Direction.Est<-Direction_fixedtuning_lin(Xc,loading,mu=sqrt(2.01*log(p)/n)*resol^{-(step-1)})
    while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
      #print(paste("step is", step))
      step<-step-1
      Direction.Est <- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(pp) / n) * resol^{-(step - 1)})
    }
    #    }
    #    else{
    ### for option 2
    #      Direction.Est<-Direction_searchtuning_lin(Xc,loading,mu=NULL, resol, maxiter)
    #      step<-Direction.Est$step
    #print(paste("step is", step))
    #    }
    #print(paste("step is", step))
    direction<-loading.norm*Direction.Est$proj
  }
  return(direction)
}
#############################################################################################
####### compute the estimator of the regression covariance matrix in the covariate shift setting
#############################################################################################
Gamma.shift<-function(plug.in,X.l,X.k,omega.l,omega.k,Y.l,Y.k,Pred.l,Pred.k){
  u.lk<-proj.direction(X.l,omega.k)
  u.kl<-proj.direction(X.k,omega.l)
  n.k<-nrow(X.k)
  n.l<-nrow(X.l)
  prop.est<-plug.in+t(u.kl)%*%t(X.k)%*%(Y.k-Pred.k)/n.k+t(u.lk)%*%t(X.l)%*%(Y.l-Pred.l)/n.l
  returnList <- list("est" = prop.est,
                     "proj.lk" = u.lk,
                     "proj.kl" = u.kl
  )
  return(returnList)
}
##############################################################################################################################################
#################### compute the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
#################### the covariance for the targeted distribution is unknown
##############################################################################################################################################
cov.inner.shift<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Pred.mat.target,Proj.array){
  N<-dim(Pred.mat.target)[1]
  Sigma.est.l1<-(1/dim(X.l1)[1])*(t(X.l1)%*%X.l1)
  Sigma.est.k1<-(1/dim(X.k1)[1])*(t(X.k1)%*%X.k1)
  var1<-0
  if(l2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[l2,k2,]/dim(X.l1)[1]
  }
  if(k2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[k2,l2,]/dim(X.l1)[1]
  }
  if(l2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[l2,k2,]/dim(X.k1)[1]
  }
  if(k2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[k2,l2,]/dim(X.k1)[1]
  }
  var2<-mean((diag(Pred.mat.target[,k1]%*%t(Pred.mat.target[,l1]))-mean(Pred.mat.target[,k1]*Pred.mat.target[,l1]))*(diag(Pred.mat.target[,k2]%*%t(Pred.mat.target[,l2]))-mean(Pred.mat.target[,k2]*Pred.mat.target[,l2])))
  var<-var1+var2/N
  return((var))
}
##############################################################################################################################################
#################### compute the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
#################### the covariance for the targeted distribution is known
##############################################################################################################################################
cov.inner.shift.known<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Proj.array){
  Sigma.est.l1<-(1/dim(X.l1)[1])*(t(X.l1)%*%X.l1)
  Sigma.est.k1<-(1/dim(X.k1)[1])*(t(X.k1)%*%X.k1)
  var1<-0
  if(l2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[l2,k2,]/dim(X.l1)[1]
  }
  if(k2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[k2,l2,]/dim(X.l1)[1]
  }
  if(l2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[l2,k2,]/dim(X.k1)[1]
  }
  if(k2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[k2,l2,]/dim(X.k1)[1]
  }
  var<-var1
  return((var))
}


