library(hslm)
data(diabetes_x_train)
data(diabetes_y_train)
hs_res <- hslm(diabetes_y_train, diabetes_x_train)
colMeans(hs_res$beta[-(1:1000),]) # Mean of beta parameters with burnin (1000) removed
hslm()


hs_res <- hslm(y, X_data)
y
X_data

diabetes_y_train


dim(diabetes_x_train)

hs_res<-horseshoe(diabetes_y_train,diabetes_x_train,sigma_hat = 1)
colMeans(hs_res$beta[-(1:1000),])
theta
hslm()

dim(diabetes_x_train)
diabetes_x_train<-as.matrix(diabetes_x_train)
diabetes_y_train<- as.matrix(diabetes_y_train)
diabe_lm<-lm(diabetes_y_train~diabetes_x_train+0)
summary(diabe_lm)$coef[,1]




#######generate data
set.seed(1200)
ns<-200
index_0<-sample(1:2,prob=c(0.05,0.95),size=ns,replace=TRUE)
theta_mix<-cbind(as.numeric(rt(ns,df=2))*3,0)
theta<-as.numeric()
for(i in 1:ns){
  theta[i]<-theta_mix[i,index_0[i]]
}
theta
#X_data<-matrix(rnorm(500*ns,mean=0,sd=100),ncol=ns,nrow=500)
X_data<-diag(1,200)

y<-rnorm(200,mean=X_data%*%theta,sd=rep(1,200))
y<-as.matrix(y,nrow=200)
theta


#1.
#test MLE
X_data<-as.matrix(X_data)
y<-as.matrix(y)
dim(X_data)
theta
beta<-solve(t(X_data)%*%X_data)%*%t(X_data)%*%y

sum((beta-theta)^2)

glm<-lm(y~X_data+0)
sum((summary(glm)$coef[,1]-theta)^2)

#2.
#test what I write
var(y)
tt<-horseshoe(y=y,x=X_data,iter=2000,sigma_hat = 1)

###sum of square error of horseshoe estimator
beta_hl<-apply(tt$beta[1000:2000,],2,mean)
sum((beta_hl-theta)^2)

#3.
#test horseshoe
#convert it to data frame
y<-as.data.frame(y)
X_data<-as.data.frame(X_data)
y<-rnorm(200,mean=as.matrix(X_data)%*%theta,sd=1)
y<-as.data.frame(y)
test_hlsm<-hslm(y,X_data)
sum((colMeans(test_hlsm$beta[-(1:1000),])-theta)^2)


tt<-horseshoe(y=y,diabetes_x_train,sigma_hat = 1)



rchisq(1,df=200)













#####Data set1 w=0.05, kuosi=2,  X=I200*200
set.seed(2100)

for(dat_iter in 1:500){
  
  mle<-as.numeric()
  horseshoe_sq<-as.numeric()
  
  ns<-200
  index_0<-sample(1:2,prob=c(0.05,0.95),size=ns,replace=TRUE)
  theta_mix<-cbind(as.numeric(rt(ns,df=2))*3,0)
  theta<-as.numeric()
  for(i in 1:ns){
    theta[i]<-theta_mix[i,index_0[i]]
  }
  theta
  #X_data<-matrix(rnorm(500*ns,mean=0,sd=100),ncol=ns,nrow=500)
  X_data<-diag(1,200)
  
  y<-rnorm(200,mean=X_data%*%theta,sd=rep(1,200))
  y<-as.matrix(y,nrow=200)
  theta
  
  #test MLE
  beta<-solve(t(X_data)%*%X_data)%*%t(X_data)%*%y
  mle[dat_iter]<-sum((beta-theta)^2)
  
  #test what I write
  
  tt<-horseshoe(y=y,x=X_data,iter=2000,sigma_hat = 1)
  
  ###sum of square error of horseshoe estimator
  beta_hl<-apply(tt$beta[1000:2000,],2,mean)
  horseshoe_sq[dat_iter]<-sum((beta_hl-theta)^2)
  
}
  
mean_mle<-mean(mle)
mean_horseshoe<-mean(horseshoe_sq)




#####Data set2 w=0.05, kuosi=2,  X=I200*200
set.seed(2100)

for(dat_iter in 1:500){
  
  mle<-as.numeric()
  horseshoe_sq<-as.numeric()
  
  ns<-200
  index_0<-sample(1:2,prob=c(0.05,0.95),size=ns,replace=TRUE)
  theta_mix<-cbind(as.numeric(rt(ns,df=2))*3,0)
  theta<-as.numeric()
  for(i in 1:ns){
    theta[i]<-theta_mix[i,index_0[i]]
  }
  theta
  X_data<-matrix(rnorm(200*ns,mean=0,sd=1),ncol=ns,nrow=200)
  #X_data<-diag(1,200)
  
  y<-rnorm(200,mean=X_data%*%theta,sd=rep(1,200))
  y<-as.matrix(y,nrow=200)
  theta
  
  #test MLE
  beta<-solve(t(X_data)%*%X_data)%*%t(X_data)%*%y
  mle[dat_iter]<-sum((beta-theta)^2)
  
  #test what I write
  
  tt<-horseshoe(y=y,x=X_data,iter=2000,sigma_hat = 1)
  
  ###sum of square error of horseshoe estimator
  beta_hl<-apply(tt$beta[1000:2000,],2,mean)
  horseshoe_sq[dat_iter]<-sum((beta_hl-theta)^2)
  cat("data number:", dat_iter)
  
  
}

mean_mle<-mean(mle)
mean_horseshoe<-mean(horseshoe_sq)

