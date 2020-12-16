# Maximin-Inference
### Read the packages
```R
library(MASS);
library(Matrix);
library(glmnet);
library(intervals);
library(CVXR);
```
### Load the source functions
```R
source('LF_Inference.R', encoding = 'UTF-8')
source('Source_Functions.R', encoding = 'UTF-8')
```
### Example
```R
### sample size n
n=500
##### the dimension p
p=500
##### L is the number of groups
L<-3
##### shift known is set as true or false
##### shift known=FALSE indicates the the covariance matrix for the targeted distribution is unknown
##### shift known=TRUE indicates the the covariance matrix for the targeted distribution is known
shift.known=FALSE
##### the amount of unlabelled data for the targeted distribution
N=2000
###### the ridge penalty level (the only tuning parameter decided by the user)
###### delta=0 corresponds to the regular non-penalized maximin effect
###### delta>0 is the general ridge type maximin effect
###### a larger delta leads to a more stable inference targets
###### the reward optimality decreases with a larger delta
delta<-0
### tuning parameter for LASSO estimator (we are using CV.min by default)
lam.value<-"CV.min"
### the sampling number (defaulted as 500)
sampling.number=500
########################
### coefficients (setting 3 in the paper)
rho=0.6
Cov<-(A1gen(rho,p))
###### generate the multiple regression vectors 
  data<-read.table('Coef.txt',header=TRUE,sep=" ")
  Coef.matrix<-as.matrix(data)[,1:L]
###### generate the xnew  
  New_Obs<- read.csv(file="New_Obs_Logistic1.csv", header=TRUE, sep=",")
  loading<- New_Obs[,1]/5
###### generate the covariance matrix for the targeted population
Cov.target<-Cov
for(i in 1:6){
  for(j in 1:6){
    if(i!=j){
      Cov.target[i,j]<-0.75
    }
  }
}
diag(Cov.target)<-diag(Cov.target)+0.1
##### compute the true regression covariance matrix
Gamma<-t(Coef.matrix)%*%Cov.target%*%Coef.matrix
##### compute the true optimal aggregation weight with the chosen delta
weight.sol<-opt.weight(Gamma,delta)
##### weight.opt is an L dim vector
weight.opt<-weight.sol$weight
###### assign the weight to compute the betamin 
beta.maxmin<-Coef.matrix%*%weight.opt
###### assign the loading and compute the true value
true.val<-sum(beta.maxmin*loading)
### generate the data
multivnorm=matrix(mvrnorm(n,mu=rep(0,L),Sigma=diag(L)),nrow=n,ncol=L)
### generate the covariate data
X.all<-mvrnorm(n*L, rep(0, p), Cov)
### generate the unlabelled data for the targeted population
X.target<-mvrnorm(N, rep(0, p), Cov.target)
### store the intermediate estimators
  Coef.est<-matrix(0,p+1,L)
  Pred.mat<-matrix(0,n,L)
  Pred.mat.total<-matrix(0,n*L,L)
  Pred.mat.target<-matrix(0,N,L)
  Proj.mat<-matrix(NA,nrow=p+1, ncol=L)
  Point.vec<-rep(NA,L)
  Var.vec<-rep(NA,L)
  SE.vec<-rep(NA,L)
  Y<-matrix(NA,nrow=n,ncol=L)
  for(l in 1:L){
    index.set<-seq((l-1)*n+1,l*n)
    X<-X.all[index.set,]
    col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
    X.norm <- X %*% diag(col.norm);
    Y[,l]<- X%*%Coef.matrix[,l] + multivnorm[,l]
    Coef.est[,l]<- Lasso (X.norm,Y[,l],lambda=lam.value,intercept=TRUE)
    Coef.est[-1,l]<-Coef.est[-1,l]*col.norm;
    Pred.mat[,l]<-X%*%Coef.est[-1,l]+Coef.est[1,l]
    Pred.mat.total[,l]<-X.all%*%Coef.est[-1,l]+Coef.est[1,l]
    Pred.mat.target[,l]<-X.target%*%Coef.est[-1,l]+Coef.est[1,l]
    est<-LF(X,Y[,l],loading,intercept=TRUE,init.Lasso=Coef.est[,l])
    SE.vec[l]<-est$se
    Point.vec[l]<-est$prop.est
    supp.l=which(abs(Coef.est[,l])>0.01)
    n.eff<-max(0.9*n,n-length(supp.l))
    Var.vec[l]<-sum((Y[,l]-Pred.mat[,l])^2)/n.eff
  }
  X.target.b<-cbind(1,X.target)
###### Gamma.plugin is the plug-in estimator of the regression covariance matrix  
  if(shift.known==TRUE){
    Sigma.target.est<-matrix(0,p+1,p+1)
    Sigma.target.est[1,1]=1
    Sigma.target.est[-1,-1]=Cov.target
    Gamma.plugin<-t(Coef.est)%*%Sigma.target.est%*%Coef.est
  }else{
    Sigma.target.est<-(1/dim(X.target.b)[1])*(t(X.target.b)%*%X.target.b)
    Gamma.plugin<-t(Pred.mat.target)%*%Pred.mat.target/dim(Pred.mat.target)[1]
  }
  Omega.est<-Sigma.target.est%*%Coef.est
###### conduct the bias correction
###### Gamma.prop is our proposed bias corrected estimator
Gamma.prop<-Gamma.plugin
  Proj.array<-array(NA, dim=c(L,L,p+1))
  for(l in 1:L){
    for(k in l:L){
      index.set.l<-seq((l-1)*n+1,l*n)
      index.set.k<-seq((k-1)*n+1,k*n)
      X.l<-cbind(1,X.all[index.set.l,])
      X.k<-cbind(1,X.all[index.set.k,])
      output<-Gamma.shift(Gamma.plugin[l,k],X.l,X.k,Omega.est[,l],Omega.est[,k],Y[,l],Y[,k],Pred.mat[,l],Pred.mat[,k])
      Gamma.prop[l,k]<-output$est
      Proj.array[l,k,]<-output$proj.lk
      Proj.array[k,l,]<-output$proj.kl  
    }
  }
  for(l in 2:L){
    for(k in 1:(l-1)){
      Gamma.prop[l,k]<-Gamma.prop[k,l]
    }
  }
###### compuate the mean and covariance matrix for the sampling distribution
  gen.mu<-Gamma.prop[lower.tri(Gamma.prop, diag = TRUE)]
  #Gamma.est.mat[i,]<- gen.mu
  gen.dim<-L*(L+1)/2
  gen.Cov<-matrix(NA,nrow=gen.dim,ncol=gen.dim)
  for(k1 in 1:L){
    for(l1 in k1:L){
      index1<-index.map(L,l1,k1)
      for(k2 in 1:L){
        for(l2 in k2:L){
          index2<-index.map(L,l2,k2)
          index.set.l1<-seq((l1-1)*n+1,l1*n)
          index.set.k1<-seq((k1-1)*n+1,k1*n)
          index.set.l2<-seq((l2-1)*n+1,l2*n)
          index.set.k2<-seq((k2-1)*n+1,k2*n)
          X.l1<-cbind(1,X.all[index.set.l1,])
          X.k1<-cbind(1,X.all[index.set.k1,])
          X.l2<-cbind(1,X.all[index.set.l2,])
          X.k2<-cbind(1,X.all[index.set.k2,])
          if(shift.known==TRUE){
            gen.Cov[index1,index2]<-cov.inner.shift.known(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Proj.array)
          }else{
            gen.Cov[index1,index2]<-cov.inner.shift(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Pred.mat.target,Proj.array)
          }
        }
      }
    }
  }
  gen.Cov
###### we conduct the sample procedure with sampling size 500
  tau<-0.2
  gen.Cov<-gen.Cov+diag(max(tau*diag(gen.Cov),1/n),dim(gen.Cov)[2])
  gen.size=500
  gen.samples=matrix(mvrnorm(gen.size,mu=gen.mu,Sigma=gen.Cov),nrow=gen.size,ncol=gen.dim)
  gen.weight.mat<-matrix(NA,nrow =gen.size, ncol=L)
  gen.est<-matrix(NA,nrow =gen.size, ncol=L)
  solution<-opt.weight(Gamma.prop,delta)
  weight.vector<-solution$weight
###### point.est is the point estimator of the maximin effect  
  point.est<-sum(Point.vec*weight.vector)
###### construct CI for the maximin effect by sampling   
    for(g in 1: gen.size){
      gen.matrix<-matrix(NA,nrow=L,ncol=L)
      gen.matrix[lower.tri(gen.matrix, diag = TRUE)]<-gen.samples[g,]
      for(l in 1:L){
        for(k in 2:L){
          gen.matrix[l,k]<-gen.matrix[k,l]
        }
      }
      gen.solution<-opt.weight(gen.matrix,delta)
      gen.weight.vector<-gen.solution$weight
      gen.weight.mat[g,]<- gen.weight.vector
      gen.point<-sum(Point.vec*gen.weight.vector)
      gen.se<-sqrt(sum(gen.weight.vector^2*SE.vec^2))
      gen.est[g,1]<-gen.point
      gen.est[g,2]<-gen.se
    }
    CI.original<-cbind(gen.est[,1]-1.96*gen.est[,2],gen.est[,1]+1.96*gen.est[,2])
    CI<-na.omit(CI.original)
    CI.lower<-min(CI[,1])
    CI.upper<-max(CI[,2])
    uni<- Intervals(CI)
 ###### construct the confidence interval by taking a union
   CI.union<-as.matrix(interval_union(uni))
 ###### contain indicates whether the constructed CI covers the true.val
 ###### det indicates whether the constructed CI covers zero
 ###### leng.CI represents the CI length
    contain<-FALSE
    det<-FALSE
    leng.CI<-0
    for(g in 1: dim(CI.union)[1]){
      if((CI.union[g,1]<true.val)*(CI.union[g,2]>true.val)){
        contain<-TRUE
      }
      if((CI.union[g,1]<0)*(CI.union[g,2]>0)==FALSE){
        det<-TRUE
      }
      leng.CI<-leng.CI+(CI.union[g,2]-CI.union[g,1])
    }
    
point.est   
CI.lower
CI.upper
CI.union
contain
det
```
