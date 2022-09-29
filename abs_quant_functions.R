abs_quant_cal<-function(m_sim,p_retained,B=1000) {
  
  rep_num<-length(m_sim)
  
  # record the results of multinomial distributed Qs
  lambda_multinom_boots<-matrix(nrow = rep_num,ncol=B)
  ## average m_sim to get mu_est for BinomVar
  mu_est<-mean(m_sim)
  
  for (i in 1:rep_num){
    pi_i<-1-exp(-mu_est/p_retained[i])
    x_multinom_boots=rbinom(B,p_retained[i],pi_i)
    lambda_multinom_boots[i,]<--log(1-x_multinom_boots/p_retained[i])
  }
  
  
  result_list<-list(# bootstrapping method result (conditioned on a fixed m)
    lambda_multinom_boots=lambda_multinom_boots  # multinomial sampling
  )
  return(result_list)
  
}

abs_quant_bootstrap<-function(x_num,partitionnr,B=1000){
  # simulation results
  rep_num<-length(x_num)
  
  glmm_est<-numeric(rep_num)
  
  lambda_sim<--log(1-x_num/partitionnr)
  m_sim<-partitionnr*lambda_sim
  
  result<-abs_quant_cal(m_sim,partitionnr,B)
  
  lambda_multinom<-result$lambda_multinom_boots
  
  #  empirical mean of ms as an estimate of lambdas
  lambda_ave_sim<-mean(lambda_sim)
  
  ##NonPVar
  lambda_ave_cond_total_var<-var(lambda_sim)/rep_num
  
  # bootstrapping total variance
  lambda_ave_multinom_var<-mean(apply(lambda_multinom,1,var))/rep_num
  
  if (!is.na(lambda_ave_cond_total_var)){
    ci_cond<-c(lambda_ave_sim-qt(0.975,rep_num-1)*sqrt(lambda_ave_cond_total_var),lambda_ave_sim+qt(0.975,rep_num-1)*sqrt(lambda_ave_cond_total_var))
    
  } else {ci_cond<-NA}
  
  if (!is.na(lambda_ave_multinom_var)){
    ci_multinom<-c(lambda_ave_sim-qnorm(0.975)*sqrt(lambda_ave_multinom_var),lambda_ave_sim+qnorm(0.975)*sqrt(lambda_ave_multinom_var))
    
  } else {ci_multinom<-NA}
  
  
  ## Delta
  pos_prop<-x_num/partitionnr
  pos_prop_ave<-mean(pos_prop)
  
  pos_prop_var<-pos_prop*(1-pos_prop)/partitionnr
  pos_prop_var_ave<-mean(pos_prop_var)/rep_num
  
  delta_abs_quant_estimate<--log(1-pos_prop_ave)
  
  delta_abs_quant_CIU<--log(1-pos_prop_ave-1.96*sqrt(pos_prop_var_ave))
  delta_abs_quant_CIL<--log(1-pos_prop_ave+1.96*sqrt(pos_prop_var_ave))
  
  delta_abs_quant_CI<-c(delta_abs_quant_estimate,delta_abs_quant_CIL,delta_abs_quant_CIU)
  
  # glmm
  # different wells
  ## GLMM
  data_test<-data.frame(name='same gene', well='same well',pos=x_num,neg=partitionnr-x_num,target=T)
  
  db_glmm<-generate.df(data_test)
  model.rand<-glmer(Y~1+(1|repl),data=db_glmm,family=binomial(link='cloglog'),nAGQ=1,verbose=F)
  
  glmm_est[1]<-(exp(model.rand@beta[1]+0.5*summary(model.rand)$varcor$repl[1]))
  glmm_est[2]<-(exp(model.rand@beta[1]+0.5*summary(model.rand)$varcor$repl[1]-1.96*sqrt(vcov(model.rand)[1,1])))
  glmm_est[3]<-(exp(model.rand@beta[1]+0.5*summary(model.rand)$varcor$repl[1]+1.96*sqrt(vcov(model.rand)[1,1])))
  
  glmm_var<-vcov(model.rand)[1,1]
  
  result_abs_quant<-list(lambda_ave_sim=lambda_ave_sim,
                         ci_cond=ci_cond,
                         ci_multinom=ci_multinom,
                         delta_abs_quant_CI=delta_abs_quant_CI,
                         glmm_est=glmm_est)
  return(result_abs_quant)
}