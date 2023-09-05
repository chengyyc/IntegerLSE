#Yingyu Cheng, 09/06/2021
##The code consists of three sections: (i) functions of Hassibi's method (with LLL) algorithms and other basic functions we will use; 
##(ii) functions of 9 methods and prediction; (iii) a simple illustrative example using simulation.

##To apply the code, please read README.txt first.


library(MASS)
library(matlib)
library(expm)
library(pracma)
library(rstan)
library(rstanarm)
library(stringr)
library(pROC)

#-------- Section 1: functions of Hassibi's method (with LLL) algorithms --------

#Part 1: LLL algotithm function
LLL<-function(B,delta) {
  n<-ncol(B)-1
  B_star<-GramSchmidt(B,normalize=FALSE)
  mu<-function(i,j) {
    v<-B[,i]
    u<-B_star[,j]
    return(dot(v,u)/dot(u,u))
  }
  k<-1
  while (k<=n) {
    for (j in (k-1):0) {
      if (abs(mu(k+1,j+1))>0.5) {
        B[,k+1]<-B[,k+1]-round(mu(k+1,j+1))*B[,j+1]
        B_star<-GramSchmidt(B,normalize=FALSE)
      }
    }
    if (dot(B_star[,k+1],B_star[,k+1])>=(delta-mu(k+1,k)^2)*dot(B_star[,k],B_star[,k]))
      k<-k+1
    else {
      t<-B[,k]
      B[,k]<-B[,k+1]
      B[,k+1]<-t
      B_star<-GramSchmidt(B,normalize=FALSE)
      k<-max(k-1,1)
    }
  }
  return(B)
}

#Part 2: Computing Upper and Lower Bounds on P_c (IV-C of Hassibi's paper)
##Calculating P_c_up
P_c_up_cal<-function(G,q) {
  a_q<-pi^(q/2)/gamma(q/2+1)
  P_c_up<-pchisq((det(G)/a_q)^(2/q),q)
  return(P_c_up)
}
##Calculating a lower bound of d_min
d_min_lower_cal<-function(G) {
  G_LLL<-LLL(G,0.75)
  G_gramsch<-GramSchmidt(G_LLL,normalize=FALSE)
  d_min_lower<-sqrt(min(apply(G_gramsch^2,2,sum)))
  return(d_min_lower)
}
##Calculating P_c_low
P_c_low_cal<-function(d,q) {
  P_c_low<-pchisq(d^2/4,q)
  return(P_c_low)
}

#Part 3: Suboptimal Polynomial-Time Algorithms (V-B of Hassibi's paper)
##Calculating z_sub
z_sub_cal<-function(G,delta,y_tilde) {
  if (ncol(G)>1) G_bar<-LLL(G,delta) #delta is usually chosen as 0.75
  if (ncol(G)==1) G_bar<-G
  F<-solve(G)%*%G_bar
  z_sub<-F%*%round(solve(G_bar)%*%y_tilde)
  return(z_sub)
}

#Part 4: Searching for Integer Points Inside an Ellipsoid (V-C of Hassibi's paper)
##Calculating r,t,l,u
r_calculate<-function(r,G,y_tilde,res,z,i) {
  if (i==1) return(res)
  else {
    sum<-0
    for (k in 1:(i-1)) sum<-sum+G[i-1,k]*z[k]
    value<-sqrt(r[i-1]^2-(y_tilde[i-1]-sum)^2)
    return(value)
  }
}
t_calculate<-function(G,y_tilde,z,i) {
  if (i==1) return(y_tilde[1])
  else {
    sum<-0
    for (k in 1:(i-1)) sum<-sum+G[i,k]*z[k]
    value<-y_tilde[i]-sum
    return(value)
  }
}
l_calculate<-function(G,t,r,i){
  l_i<-ceiling((t[i]-sign(G[i,i])*r[i])/G[i,i])
  return(l_i)
}
u_calculate<-function(G,t,r,i){
  u_i<-floor((t[i]+sign(G[i,i])*r[i])/G[i,i])
  return(u_i)
}
##Search the element in S closest to the given value x
closest<-function(x,S) {
  S_new<-abs(S-x)
  index<-which(S_new==min(S_new))
  if (length(index)>1) return(S[index[1]])
  else return(S[index])
}
##Searching function
SearchForInt<-function(G,y_tilde,res,z,q) {
  QR_inv_G_t<-qr(solve(t(G))) #Find unitary transformation Theta such that Theta%*%G is lower triangular
  G_lowertri<-t(solve(qr.R(QR_inv_G_t)))
  Theta<-t(qr.Q(QR_inv_G_t))
  y_tilde_new<-Theta%*%y_tilde
  G_new<-G_lowertri
  r<-rep(0,q)
  t<-rep(0,q)
  l<-rep(0,q)
  u<-rep(0,q)
  S<-list(0)
  i<-1
  while (1) {
    r[i]<-r_calculate(r,G_new,y_tilde_new,res,z,i)
    t[i]<-t_calculate(G_new,y_tilde_new,z,i)
    l[i]<-l_calculate(G_new,t,r,i)
    u[i]<-u_calculate(G_new,t,r,i)
    S[[i]]<-c(l[i]:u[i])
    while (l[i]>u[i] || length(S[[i]])==0) {
      if (i==1) return("NoInteger")
      i<-i-1
    }
    z[i]<-closest((max(S[[i]])+min(S[[i]]))/2,S[[i]])
    bridge<-S[[i]]
    bridge<-bridge[bridge!=z[i]]
    S[[i]]<-bridge
    if (i==q) return(z)
    else i<-i+1
  }
}

#Part 5: Global optimization algorithm (V-D of Hassibi's paper)
global_opt<-function(q,y_tilde,G,delta,d_min_lower) {
  z_ML<-rep(0,q)
  z_sub<-z_sub_cal(G,delta,y_tilde)
  if (sum((y_tilde-(G%*%z_sub))^2)<d_min_lower/2) {
    z_ML<-z_sub
    return(round(z_ML))
  }
  else {
    r<-sum((y_tilde-G%*%z_sub)^2)
    z_star<-z_sub
    times<-0
    while (1) {
      z_search<-SearchForInt(G,y_tilde,r,z_star,q)
      if (z_search=="NoInteger") {
        z_ML<-z_star
        return(round(z_ML))
      }
      z_star<-z_search
      r<-sum((y_tilde-G%*%z_star)^2)
      if (r<d_min_lower/2) {
        z_ML<-z_star
        return(round(z_ML))
      }
      times<-times+1
      if (times>100) {
        z_ML<-z_star
        return(round(z_ML))
      }
    }
  }
}


#Part6: other functions we will use
##Function of calculating the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
##Function of calculating the most frequent row
getmode_row<-function(x) {
  count<-table(apply(x, 1, paste, collapse=" "))
  row_mode<-count[which.max(count)]
  index<-names(row_mode)
  our_result<-strsplit(index,split=" ")
  return(as.integer(our_result[[1]]))
}
##Function generating continuous predictors
continuous_predictors<-function(q,T,rho) {
  mu_b<-rep(0,q)
  sigma_b=1 #now the distribution for each x_i is N(0,1); you can change the variance (sigma_b) to other values
  cov_b<-rho*sigma_b^2
  Sigma_b<-matrix(0,ncol=q,nrow=q)
  for (i in 1:q) {
    for (j in 1:q) {
      if (i==j) Sigma_b[i,j]<-sigma_b^2
      if (i!=j) Sigma_b[i,j]<-cov_b
    }
  }
  predictors<-mvrnorm(T,mu_b,Sigma_b)
  return(predictors)
}
##Function generating binary predictors
binary_predictors<-function(q,T,rho) {
  latent_predictors<-continuous_predictors(q,T,rho)
  predictors<-1*(latent_predictors>0)
  return(predictors)
}
##Function generating continuous and binary predictors
con_bin_predictors<-function(q,qb,T,rho) {
  all_predictors<-continuous_predictors(q,T,rho)
  con_predictors<-all_predictors[,1:(q-qb)]
  latent_predictors<-all_predictors[,(q-qb+1):q]
  bin_predictors<-1*(latent_predictors>0)
  predictors<-cbind(con_predictors,bin_predictors)
  return(predictors)
}


#-------- Section 2: functions of 9 methods (3 benchmarks and 6 proposed methods) and prediction funtion --------

#Prediction function, where obj is any output of the functions below
pred_BIL<-function(obj,x_test) {
  para<-obj$coefficients
  ori<-para[1]+as.matrix(x_test)%*%para[-1]
  pro<-exp(ori)/(1+exp(ori))
  y<-as.numeric(pro>0.5)
  return(y)
}

#Method 1: optimizing result (MAP)
BILOpt<-function(x,y,prior) {
  model_Opt<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,algorithm="optimizing",iter=10000)
  
  para<-model_Opt$coefficients
  x_Opt<-para[1]
  z_Opt<-round(para[-1])
  para_Opt<-c(x_Opt,z_Opt)
  
  ori_Opt<-x_Opt+as.matrix(x)%*%z_Opt
  pro_Opt<-exp(ori_Opt)/(1+exp(ori_Opt))
  y_Opt<-as.numeric(pro_Opt>0.5)
  
  result_Opt<-list()
  result_Opt$coefficients<-para_Opt
  result_Opt$estimated_probability<-pro_Opt
  result_Opt$classification_result<-y_Opt
  
  return(result_Opt)
}

#Method 2: directly rounding Bayesian logistic estimator (mean)
BILBayMean<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  
  para_Bay<-model_Bay$coefficients
  x_BayMean<-para_Bay[1]
  para_Bay2<-apply(post_sam,2,mean)
  z_BayMean<-round(para_Bay2[-1])
  para_BayMean<-c(x_BayMean,z_BayMean)
  
  ori_BayMean<-x_BayMean+as.matrix(x)%*%z_BayMean
  pro_BayMean<-exp(ori_BayMean)/(1+exp(ori_BayMean))
  y_BayMean<-as.numeric(pro_BayMean>0.5)
  
  result_BayMean<-list()
  result_BayMean$coefficients<-para_BayMean
  result_BayMean$estimated_probability<-pro_BayMean
  result_BayMean$classification_result<-y_BayMean
  
  return(result_BayMean)
}

#Method 3: directly rounding Bayesian logistic estimator (median)
BILBayMedian<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  
  para_Bay<-model_Bay$coefficients
  x_BayMedian<-para_Bay[1]
  z_BayMedian<-round(para_Bay[-1])
  para_BayMedian<-c(x_BayMedian,z_BayMedian)
  
  ori_BayMedian<-x_BayMedian+as.matrix(x)%*%z_BayMedian
  pro_BayMedian<-exp(ori_BayMedian)/(1+exp(ori_BayMedian))
  y_BayMedian<-as.numeric(pro_BayMedian>0.5)
  
  result_BayMedian<-list()
  result_BayMedian$coefficients<-para_BayMedian
  result_BayMedian$estimated_probability<-pro_BayMedian
  result_BayMedian$classification_result<-y_BayMedian
  
  return(result_BayMedian)
}

#Method 4: project by rounding then obtain median
BILRoundMedian<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  post_sam_trunc<-post_sam[,-1]
  
  para_Bay<-model_Bay$coefficients
  x_RoundMedian<-para_Bay[1]
  z_round<-round(post_sam_trunc)
  z_RoundMedian<-apply(z_round,2,median)
  para_RoundMedian<-c(x_RoundMedian,z_RoundMedian)
  
  ori_RoundMedian<-x_RoundMedian+as.matrix(x)%*%z_RoundMedian
  pro_RoundMedian<-exp(ori_RoundMedian)/(1+exp(ori_RoundMedian))
  y_RoundMedian<-as.numeric(pro_RoundMedian>0.5)
  
  result_RoundMedian<-list()
  result_RoundMedian$coefficients<-para_RoundMedian
  result_RoundMedian$estimated_probability<-pro_RoundMedian
  result_RoundMedian$classification_result<-y_RoundMedian
  
  return(result_RoundMedian)
}

#Method 5: project by rounding then obtain mode
BILRoundMode<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  post_sam_trunc<-post_sam[,-1]
  
  para_Bay<-model_Bay$coefficients
  x_RoundMode<-para_Bay[1]
  z_round<-round(post_sam_trunc)
  z_RoundMode<-apply(z_round,2,getmode)
  para_RoundMode<-c(x_RoundMode,z_RoundMode)
  
  ori_RoundMode<-x_RoundMode+as.matrix(x)%*%z_RoundMode
  pro_RoundMode<-exp(ori_RoundMode)/(1+exp(ori_RoundMode))
  y_RoundMode<-as.numeric(pro_RoundMode>0.5)
  
  result_RoundMode<-list()
  result_RoundMode$coefficients<-para_RoundMode
  result_RoundMode$estimated_probability<-pro_RoundMode
  result_RoundMode$classification_result<-y_RoundMode
  
  return(result_RoundMode)
}

#Method 6: project by rounding then obtain mode of row
BILRoundModeRow<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  post_sam_trunc<-post_sam[,-1]
  
  para_Bay<-model_Bay$coefficients
  x_RoundModeRow<-para_Bay[1]
  z_round<-round(post_sam_trunc)
  z_RoundModeRow<-getmode_row(z_round)
  para_RoundModeRow<-c(x_RoundModeRow,z_RoundModeRow)
  
  ori_RoundModeRow<-x_RoundModeRow+as.matrix(x)%*%z_RoundModeRow
  pro_RoundModeRow<-exp(ori_RoundModeRow)/(1+exp(ori_RoundModeRow))
  y_RoundModeRow<-as.numeric(pro_RoundModeRow>0.5)
  
  result_RoundModeRow<-list()
  result_RoundModeRow$coefficients<-para_RoundModeRow
  result_RoundModeRow$estimated_probability<-pro_RoundModeRow
  result_RoundModeRow$classification_result<-y_RoundModeRow
  
  return(result_RoundModeRow)
}

#Method 7: project by LLL searching algorithm then obtain median
BILLLLMedian<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  post_sam_trunc<-post_sam[,-1]
  
  q<-ncol(x)
  sam_cov<-cov(post_sam_trunc)
  G_Bay<-sqrtm(solve(sam_cov))$B
  y_tilde_Bay<-G_Bay%*%t(post_sam_trunc)
  z_sub_Bay<-z_sub_cal(G_Bay,0.75,y_tilde_Bay) #delta is usually chosen as 0.75
  d_min_lower_Bay<-d_min_lower_cal(G_Bay)
  z_lll<-global_opt(q,y_tilde_Bay,G_Bay,0.75,d_min_lower_Bay)
  z_lll<-t(z_lll)
  z_LLLMedian<-apply(z_lll,2,median)
  para_Bay<-model_Bay$coefficients
  x_LLLMedian<-para_Bay[1]
  para_LLLMedian<-c(x_LLLMedian,z_LLLMedian)
  
  ori_LLLMedian<-x_LLLMedian+as.matrix(x)%*%z_LLLMedian
  pro_LLLMedian<-exp(ori_LLLMedian)/(1+exp(ori_LLLMedian))
  y_LLLMedian<-as.numeric(pro_LLLMedian>0.5)
  
  result_LLLMedian<-list()
  result_LLLMedian$coefficients<-para_LLLMedian
  result_LLLMedian$estimated_probability<-pro_LLLMedian
  result_LLLMedian$classification_result<-y_LLLMedian
  
  return(result_LLLMedian)
}

#Method 8: project by LLL searching algorithm then obtain mode
BILLLLMode<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  post_sam_trunc<-post_sam[,-1]
  
  q<-ncol(x)
  sam_cov<-cov(post_sam_trunc)
  G_Bay<-sqrtm(solve(sam_cov))$B
  y_tilde_Bay<-G_Bay%*%t(post_sam_trunc)
  z_sub_Bay<-z_sub_cal(G_Bay,0.75,y_tilde_Bay) #delta is usually chosen as 0.75
  d_min_lower_Bay<-d_min_lower_cal(G_Bay)
  z_lll<-global_opt(q,y_tilde_Bay,G_Bay,0.75,d_min_lower_Bay)
  z_lll<-t(z_lll)
  z_LLLMode<-apply(z_lll,2,getmode)
  para_Bay<-model_Bay$coefficients
  x_LLLMode<-para_Bay[1]
  para_LLLMode<-c(x_LLLMode,z_LLLMode)
  
  ori_LLLMode<-x_LLLMode+as.matrix(x)%*%z_LLLMode
  pro_LLLMode<-exp(ori_LLLMode)/(1+exp(ori_LLLMode))
  y_LLLMode<-as.numeric(pro_LLLMode>0.5)
  
  result_LLLMode<-list()
  result_LLLMode$coefficients<-para_LLLMode
  result_LLLMode$estimated_probability<-pro_LLLMode
  result_LLLMode$classification_result<-y_LLLMode
  
  return(result_LLLMode)
}

#Method 9: project by LLL searching algorithm then obtain mode of row
BILLLLModeRow<-function(x,y,prior) {
  model_Bay<-stan_glm(y~.,data=cbind(x,y),family=binomial,refresh=0,prior=prior,chains=1,iter=20000)
  post_sam<-as.matrix(model_Bay)
  post_sam_trunc<-post_sam[,-1]
  
  q<-ncol(x)
  sam_cov<-cov(post_sam_trunc)
  G_Bay<-sqrtm(solve(sam_cov))$B
  y_tilde_Bay<-G_Bay%*%t(post_sam_trunc)
  z_sub_Bay<-z_sub_cal(G_Bay,0.75,y_tilde_Bay) #delta is usually chosen as 0.75
  d_min_lower_Bay<-d_min_lower_cal(G_Bay)
  z_lll<-global_opt(q,y_tilde_Bay,G_Bay,0.75,d_min_lower_Bay)
  z_lll<-t(z_lll)
  z_LLLModeRow<-getmode_row(z_lll)
  para_Bay<-model_Bay$coefficients
  x_LLLModeRow<-para_Bay[1]
  para_LLLModeRow<-c(x_LLLModeRow,z_LLLModeRow)
  
  ori_LLLModeRow<-x_LLLModeRow+as.matrix(x)%*%z_LLLModeRow
  pro_LLLModeRow<-exp(ori_LLLModeRow)/(1+exp(ori_LLLModeRow))
  y_LLLModeRow<-as.numeric(pro_LLLModeRow>0.5)
  
  result_LLLModeRow<-list()
  result_LLLModeRow$coefficients<-para_LLLModeRow
  result_LLLModeRow$estimated_probability<-pro_LLLModeRow
  result_LLLModeRow$classification_result<-y_LLLModeRow
  
  return(result_LLLModeRow)
}


#-------- Section 3: A simple illustrative example using simulation--------

#options(digits=3) #you can change to other digits you want to save

#Settings
n<-50 #number of repeatings
q<-10 #number of predictors (not including intercept)
#qb<-5 #number of binary predictors
N<-1000 #sample size
rho<-0 #correlation between predictors
x0<-matrix(1,ncol=1) #intercept
constant_range<-c(-3,-2,-1,0,1,2,3) #if use integer coefficients
my_prior<-normal(location=rep(0,q),scale=rep(100,q)) #diffused normal prior

#Result recording: 9 methods in total
##for diffused normal prior
z_true<-matrix(0,ncol=q,nrow=n)

MSE_Bay_opt<-rep(0,n)
MSE_Bay_median<-rep(0,n)
MSE_Bay_mean<-rep(0,n)
MSE_round_median<-rep(0,n)
MSE_round_mode<-rep(0,n)
MSE_round_mode_row<-rep(0,n)
MSE_lll_median<-rep(0,n)
MSE_lll_mode<-rep(0,n)
MSE_lll_mode_row<-rep(0,n)

auc_true_in<-rep(0,n)
auc_Bay_opt_in<-rep(0,n)
auc_Bay_median_in<-rep(0,n)
auc_Bay_mean_in<-rep(0,n)
auc_round_median_in<-rep(0,n)
auc_round_mode_in<-rep(0,n)
auc_round_mode_row_in<-rep(0,n)
auc_lll_median_in<-rep(0,n)
auc_lll_mode_in<-rep(0,n)
auc_lll_mode_row_in<-rep(0,n)

auc_true<-rep(0,n)
auc_Bay_opt<-rep(0,n)
auc_Bay_median<-rep(0,n)
auc_Bay_mean<-rep(0,n)
auc_round_median<-rep(0,n)
auc_round_mode<-rep(0,n)
auc_round_mode_row<-rep(0,n)
auc_lll_median<-rep(0,n)
auc_lll_mode<-rep(0,n)
auc_lll_mode_row<-rep(0,n)

#Simulation process
for (l in 1:n) {
  
  #-------- 1st Part: Generating data --------
  A<-matrix(1,ncol=1,nrow=N)
  
  ###predictors
  B<-continuous_predictors(q,N,rho)
  #you can also use other functions to generate binary predictors or half continuous half binary predictors
  #B<-con_bin_predictors(q,qb,N,rho)
  
  ###Coefficient
  z<-sample(constant_range,q,replace=TRUE)
  #you can also use other functions to generate coefficients of different distribution
  #z<-runif(q,-3,3)
  #z<-rnorm(q,0,1.5)
  z_true[l,]<-z
  
  ###Response
  ori<-A%*%x0+B%*%z
  pro<-exp(ori)/(1+exp(ori))
  y<-rbinom(N,1,pro)
  
  ###Combining predictors and response
  mydata<-as.data.frame(cbind(y,B))
  
  ###Training set
  mydata_train<-mydata
  y_train<-mydata$y
  A_train<-A
  B_train<-B
  x_train<-mydata[,-1]
  
  ###Testing set
  A_test<-matrix(1,ncol=1,nrow=N*0.3)
  B_test<-continuous_predictors(q,N*0.3,rho)
  ori_test<-A_test%*%x0+B_test%*%z
  pro_test<-exp(ori_test)/(1+exp(ori_test))
  yy_test<-rbinom(N*0.3,1,pro_test)
  mydata_test<-as.data.frame(cbind(yy_test,B_test))
  colnames(mydata_test)[1]<-"y"
  x_test<-mydata_test[,-1]
  y_test<-mydata_test$y
  
  #-------- 2nd Part: In-sample Prediction --------
  
  ##3 benchmarks
  method1<-BILOpt(x_train,y_train,my_prior)
  pred1<-pred_BIL(method1,x_train)
  method2<-BILBayMean(x_train,y_train,my_prior)
  pred2<-pred_BIL(method2,x_train)
  method3<-BILBayMedian(x_train,y_train,my_prior)
  pred3<-pred_BIL(method3,x_train)
  ##6 proposed methods
  method4<-BILRoundMedian(x_train,y_train,my_prior)
  pred4<-pred_BIL(method4,x_train)
  method5<-BILRoundMode(x_train,y_train,my_prior)
  pred5<-pred_BIL(method5,x_train)
  method6<-BILRoundModeRow(x_train,y_train,my_prior)
  pred6<-pred_BIL(method6,x_train)
  method7<-BILLLLMedian(x_train,y_train,my_prior)
  pred7<-pred_BIL(method7,x_train)
  method8<-BILLLLMode(x_train,y_train,my_prior)
  pred8<-pred_BIL(method8,x_train)
  method9<-BILLLLModeRow(x_train,y_train,my_prior)
  pred9<-pred_BIL(method9,x_train)
  
  #MSE
  MSE_Bay_opt[l]<-mean((method1$coefficients[-1]-z)^2)
  MSE_Bay_median[l]<-mean((method2$coefficients[-1]-z)^2)
  MSE_Bay_mean[l]<-mean((method3$coefficients[-1]-z)^2)
  MSE_round_median[l]<-mean((method4$coefficients[-1]-z)^2)
  MSE_round_mode[l]<-mean((method5$coefficients[-1]-z)^2)
  MSE_round_mode_row[l]<-mean((method6$coefficients[-1]-z)^2)
  MSE_lll_median[l]<-mean((method7$coefficients[-1]-z)^2)
  MSE_lll_mode[l]<-mean((method8$coefficients[-1]-z)^2)
  MSE_lll_mode_row[l]<-mean((method9$coefficients[-1]-z)^2)
  
  #AUC
  roc_Bay_opt_in<-roc(y_train,pred1,quiet=TRUE)
  roc_Bay_median_in<-roc(y_train,pred2,quiet=TRUE)
  roc_Bay_mean_in<-roc(y_train,pred3,quiet=TRUE)
  roc_round_median_in<-roc(y_train,pred4,quiet=TRUE)
  roc_round_mode_in<-roc(y_train,pred5,quiet=TRUE)
  roc_round_mode_row_in<-roc(y_train,pred6,quiet=TRUE)
  roc_lll_median_in<-roc(y_train,pred7,quiet=TRUE)
  roc_lll_mode_in<-roc(y_train,pred8,quiet=TRUE)
  roc_lll_mode_row_in<-roc(y_train,pred9,quiet=TRUE)
  
  auc_Bay_opt_in[l]<-as.numeric(auc(roc_Bay_opt_in))
  auc_Bay_median_in[l]<-as.numeric(auc(roc_Bay_median_in))
  auc_Bay_mean_in[l]<-as.numeric(auc(roc_Bay_mean_in))
  auc_round_median_in[l]<-as.numeric(auc(roc_round_median_in))
  auc_round_mode_in[l]<-as.numeric(auc(roc_round_mode_in))
  auc_round_mode_row_in[l]<-as.numeric(auc(roc_round_mode_row_in))
  auc_lll_median_in[l]<-as.numeric(auc(roc_lll_median_in))
  auc_lll_mode_in[l]<-as.numeric(auc(roc_lll_mode_in))
  auc_lll_mode_row_in[l]<-as.numeric(auc(roc_lll_mode_row_in))
  
  
  #-------- 3rd Part: Out-of-sample Prediction --------
  
  ##3 benchmarks
  pred1_out<-pred_BIL(method1,x_test)
  pred2_out<-pred_BIL(method2,x_test)
  pred3_out<-pred_BIL(method3,x_test)
  ##6 proposed methods
  pred4_out<-pred_BIL(method4,x_test)
  pred5_out<-pred_BIL(method5,x_test)
  pred6_out<-pred_BIL(method6,x_test)
  pred7_out<-pred_BIL(method7,x_test)
  pred8_out<-pred_BIL(method8,x_test)
  pred9_out<-pred_BIL(method9,x_test)
  
  #AUC
  roc_Bay_opt<-roc(y_test,pred1_out,quiet=TRUE)
  roc_Bay_median<-roc(y_test,pred2_out,quiet=TRUE)
  roc_Bay_mean<-roc(y_test,pred3_out,quiet=TRUE)
  roc_round_median<-roc(y_test,pred4_out,quiet=TRUE)
  roc_round_mode<-roc(y_test,pred5_out,quiet=TRUE)
  roc_round_mode_row<-roc(y_test,pred6_out,quiet=TRUE)
  roc_lll_median<-roc(y_test,pred7_out,quiet=TRUE)
  roc_lll_mode<-roc(y_test,pred8_out,quiet=TRUE)
  roc_lll_mode_row<-roc(y_test,pred9_out,quiet=TRUE)
  
  auc_Bay_opt[l]<-as.numeric(auc(roc_Bay_opt))
  auc_Bay_median[l]<-as.numeric(auc(roc_Bay_median))
  auc_Bay_mean[l]<-as.numeric(auc(roc_Bay_mean))
  auc_round_median[l]<-as.numeric(auc(roc_round_median))
  auc_round_mode[l]<-as.numeric(auc(roc_round_mode))
  auc_round_mode_row[l]<-as.numeric(auc(roc_round_mode_row))
  auc_lll_median[l]<-as.numeric(auc(roc_lll_median))
  auc_lll_mode[l]<-as.numeric(auc(roc_lll_mode))
  auc_lll_mode_row[l]<-as.numeric(auc(roc_lll_mode_row))
  
  print(l)
}

##Calculate the average result of these n repeatings
AMSE_Bay_opt<-mean(MSE_Bay_opt)
AMSE_Bay_median<-mean(MSE_Bay_median)
AMSE_Bay_mean<-mean(MSE_Bay_mean)
AMSE_round_median<-mean(MSE_round_median)
AMSE_round_mode<-mean(MSE_round_mode)
AMSE_round_mode_row<-mean(MSE_round_mode_row)
AMSE_lll_median<-mean(MSE_lll_median)
AMSE_lll_mode<-mean(MSE_lll_mode)
AMSE_lll_mode_row<-mean(MSE_lll_mode_row)

Aauc_Bay_opt_in<-mean(auc_Bay_opt_in)
Aauc_Bay_median_in<-mean(auc_Bay_median_in)
Aauc_Bay_mean_in<-mean(auc_Bay_mean_in)
Aauc_round_median_in<-mean(auc_round_median_in)
Aauc_round_mode_in<-mean(auc_round_mode_in)
Aauc_round_mode_row_in<-mean(auc_round_mode_row_in)
Aauc_lll_median_in<-mean(auc_lll_median_in)
Aauc_lll_mode_in<-mean(auc_lll_mode_in)
Aauc_lll_mode_row_in<-mean(auc_lll_mode_row_in)

sd_auc_Bay_opt_in<-sqrt(var(auc_Bay_opt_in))/sqrt(n)
sd_auc_Bay_median_in<-sqrt(var(auc_Bay_median_in))/sqrt(n)
sd_auc_Bay_mean_in<-sqrt(var(auc_Bay_mean_in))/sqrt(n)
sd_auc_round_median_in<-sqrt(var(auc_round_median_in))/sqrt(n)
sd_auc_round_mode_in<-sqrt(var(auc_round_mode_in))/sqrt(n)
sd_auc_round_mode_row_in<-sqrt(var(auc_round_mode_row_in))/sqrt(n)
sd_auc_lll_median_in<-sqrt(var(auc_lll_median_in))/sqrt(n)
sd_auc_lll_mode_in<-sqrt(var(auc_lll_mode_in))/sqrt(n)
sd_auc_lll_mode_row_in<-sqrt(var(auc_lll_mode_row_in))/sqrt(n)

Aauc_Bay_opt<-mean(auc_Bay_opt)
Aauc_Bay_median<-mean(auc_Bay_median)
Aauc_Bay_mean<-mean(auc_Bay_mean)
Aauc_round_median<-mean(auc_round_median)
Aauc_round_mode<-mean(auc_round_mode)
Aauc_round_mode_row<-mean(auc_round_mode_row)
Aauc_lll_median<-mean(auc_lll_median)
Aauc_lll_mode<-mean(auc_lll_mode)
Aauc_lll_mode_row<-mean(auc_lll_mode_row)

sd_auc_Bay_opt<-sqrt(var(auc_Bay_opt))/sqrt(n)
sd_auc_Bay_median<-sqrt(var(auc_Bay_median))/sqrt(n)
sd_auc_Bay_mean<-sqrt(var(auc_Bay_mean))/sqrt(n)
sd_auc_round_median<-sqrt(var(auc_round_median))/sqrt(n)
sd_auc_round_mode<-sqrt(var(auc_round_mode))/sqrt(n)
sd_auc_round_mode_row<-sqrt(var(auc_round_mode_row))/sqrt(n)
sd_auc_lll_median<-sqrt(var(auc_lll_median))/sqrt(n)
sd_auc_lll_mode<-sqrt(var(auc_lll_mode))/sqrt(n)
sd_auc_lll_mode_row<-sqrt(var(auc_lll_mode_row))/sqrt(n)

#Output the results
c(AMSE_Bay_opt,AMSE_Bay_median,AMSE_Bay_mean,AMSE_round_median,AMSE_round_mode,AMSE_round_mode_row,AMSE_lll_median,AMSE_lll_mode,AMSE_lll_mode_row)
c(Aauc_Bay_opt_in,Aauc_Bay_median_in,Aauc_Bay_mean_in,Aauc_round_median_in,Aauc_round_mode_in,Aauc_round_mode_row_in,Aauc_lll_median_in,Aauc_lll_mode_in,Aauc_lll_mode_row_in)
c(sd_auc_Bay_opt_in,sd_auc_Bay_median_in,sd_auc_Bay_mean_in,sd_auc_round_median_in,sd_auc_round_mode_in,sd_auc_round_mode_row_in,sd_auc_lll_median_in,sd_auc_lll_mode_in,sd_auc_lll_mode_row_in)
c(Aauc_Bay_opt,Aauc_Bay_median,Aauc_Bay_mean,Aauc_round_median,Aauc_round_mode,Aauc_round_mode_row,Aauc_lll_median,Aauc_lll_mode,Aauc_lll_mode_row)
c(sd_auc_Bay_opt,sd_auc_Bay_median,sd_auc_Bay_mean,sd_auc_round_median,sd_auc_round_mode,sd_auc_round_mode_row,sd_auc_lll_median,sd_auc_lll_mode,sd_auc_lll_mode_row)
