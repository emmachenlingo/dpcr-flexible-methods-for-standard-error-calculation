cnv_cal<-function(m_A_sim,m_B_sim,partitionnrA,partitionnrB,Nb=2,B=1000) {
  
  rep_num<-length(m_A_sim)
  # record the results of multinomial distributed Qs
  lambda_A_multinom_boots<-matrix(nrow = rep_num,ncol=B)
  lambda_B_multinom_boots<-matrix(nrow = rep_num,ncol=B)
  
  mu_A_est<-mean(m_A_sim)
  mu_B_est<-mean(m_B_sim)
  
  for (i in 1:rep_num){
    pi_A_i<-1-exp(-mu_A_est/partitionnrA[i])
    x_A_multinom_boots=rbinom(B,partitionnrA[i],pi_A_i)
    lambda_A_multinom_boots[i,]<--log(1-x_A_multinom_boots/partitionnrA[i])
    
    pi_B_i<-1-exp(-mu_B_est/partitionnrB[i])
    x_B_multinom_boots=rbinom(B,partitionnrB[i],pi_B_i)
    lambda_B_multinom_boots[i,]<--log(1-x_B_multinom_boots/partitionnrB[i])
  }
  
  result_list<-list(
    lambda_A_multinom_boots=lambda_A_multinom_boots,
    lambda_B_multinom_boots=lambda_B_multinom_boots
    # binomial sampling
  )
  return(result_list)
  
}


cnv_bootstrap_singleplex<-function(x_A_num,x_B_num,partitionnrA,partitionnrB,Nb=2,B=1000){
  # simulation results
  rep_num<-length(x_A_num)
  
  glmm_cnv<-numeric(rep_num)
  
  lambda_A_sim<--log(1-x_A_num/partitionnrA)
  m_A_sim<-partitionnrA*lambda_A_sim
  
  lambda_B_sim<--log(1-x_B_num/partitionnrB)
  m_B_sim<-partitionnrB*lambda_B_sim
  
  result<-cnv_cal(m_A_sim,m_B_sim,partitionnrA,partitionnrB,Nb,B)
  
  lambda_A_multinom<-result$lambda_A_multinom_boots
  lambda_B_multinom<-result$lambda_B_multinom_boots
  
  #  empirical mean of ms as an estimate of lambdas
  lambda_ave_A_sim<-mean(lambda_A_sim)
  lambda_ave_B_sim<-mean(lambda_B_sim)
  
  cnv_ave_sim<-mean(log(lambda_A_sim)-log(lambda_B_sim))
  
  ##NonPVar
  lambda_ave_A_cond_total_var_log<-var(log(lambda_A_sim))/rep_num
  lambda_ave_B_cond_total_var_log<-var(log(lambda_B_sim))/rep_num
  cnv_ave_cond_total_var<-lambda_ave_A_cond_total_var_log+lambda_ave_B_cond_total_var_log
  
  
  # bootstrapping total variance
  lambda_ave_A_multinom_mean_log<-apply(log(lambda_A_multinom),1,mean)
  lambda_ave_A_multinom_var_log<-apply(log(lambda_A_multinom),1,var)
  
  lambda_ave_B_multinom_mean_log<-apply(log(lambda_B_multinom),1,mean)
  lambda_ave_B_multinom_var_log<-apply(log(lambda_B_multinom),1,var)
  
  cnv_ave_multinom_var<-(sum(lambda_ave_A_multinom_var_log)+sum(lambda_ave_B_multinom_var_log))/rep_num^2
  
  
  if (!is.na(cnv_ave_cond_total_var)){
    ci_cond<-exp(c(cnv_ave_sim-qt(0.975,2*rep_num-1)*sqrt(cnv_ave_cond_total_var),cnv_ave_sim+qt(0.975,2*rep_num-1)*sqrt(cnv_ave_cond_total_var)))*Nb
    
  } else {ci_cond<-NA}
  
  if (!is.na(cnv_ave_multinom_var)){
    ci_multinom<-exp(c(cnv_ave_sim-qnorm(0.975)*sqrt(cnv_ave_multinom_var),cnv_ave_sim+qnorm(0.975)*sqrt(cnv_ave_multinom_var)))*Nb
    
  } else {ci_multinom<-NA}
  
  
  # delta
  ## delta method
  pos_partition_aveA<-mean(x_A_num)
  partition_aveA<-mean(partitionnrA)
  
  pos_partition_aveB<-mean(x_B_num)
  partition_aveB<-mean(partitionnrB)
  
  lambda_aveA<--log(1-pos_partition_aveA/partition_aveA)
  lambda_aveB<--log(1-pos_partition_aveB/partition_aveB)
  
  cnv_delta_ave<-(lambda_aveA/lambda_aveB)*Nb
  
  var_log_cnv<-((1-exp(-lambda_aveA))/(partition_aveA*(lambda_aveA^2)*exp(-lambda_aveA))+(1-exp(-lambda_aveB))/(partition_aveB*(lambda_aveB^2)*exp(-lambda_aveB)))/rep_num
  
  delta_cnv_CIU<-exp(log(lambda_aveA/lambda_aveB)+1.96*sqrt(var_log_cnv))*Nb
  delta_cnv_CIL<-exp(log(lambda_aveA/lambda_aveB)-1.96*sqrt(var_log_cnv))*Nb
  delta_cnv_CI<-c(cnv_delta_ave,delta_cnv_CIL,delta_cnv_CIU)
  
  # glmm
  # different wells
  target_data<-data.frame(gene='target gene', well='target well',pos=x_A_num,neg=partitionnrA-x_A_num,target=T)
  reference_data<-data.frame(gene='reference gene', well='reference well',pos=x_B_num,neg=partitionnrB-x_B_num,target=F)
  data_combined<-rbind(target_data,reference_data)
  db_glmm<-generate.df(data_combined)
  model.rand<-glmer(Y~dummy(target,"no")+(1|well),data=db_glmm,family=binomial(link='cloglog'),nAGQ=1,verbose=F)
  
  glmm_cnv[1]<-(exp(-(model.rand@beta[2])+summary(model.rand)$varcor$well[1]))*Nb
  glmm_cnv[2]<-(exp(-(model.rand@beta[2])+summary(model.rand)$varcor$well[1]-1.96*sqrt(vcov(model.rand)[2,2])))*Nb
  glmm_cnv[3]<-(exp(-(model.rand@beta[2])+summary(model.rand)$varcor$well[1]+1.96*sqrt(vcov(model.rand)[2,2])))*Nb
  
  
  result_cnv<-list(cnv_sim=exp(cnv_ave_sim)*Nb,
                   ci_cond=ci_cond,
                   ci_multinom=ci_multinom,
                   delta_cnv_CI=delta_cnv_CI,
                   glmm_cnv=glmm_cnv)
  return(result_cnv)
}


cnv_bootstrap_duplex<-function(x_A_num,x_B_num,partitionnrA,partitionnrB,Nb=2,B=1000){
  # simulation results
  rep_num<-length(x_A_num)
  
  glmm_cnv<-numeric(rep_num)
  
  lambda_A_sim<--log(1-x_A_num/partitionnrA)
  m_A_sim<-partitionnrA*lambda_A_sim
  
  lambda_B_sim<--log(1-x_B_num/partitionnrB)
  m_B_sim<-partitionnrB*lambda_B_sim
  
  result<-cnv_cal(m_A_sim,m_B_sim,partitionnrA,partitionnrB,Nb,B)
  
  lambda_A_multinom<-result$lambda_A_multinom_boots
  lambda_B_multinom<-result$lambda_B_multinom_boots
  cnv_multinom<-(lambda_A_multinom/lambda_B_multinom)*Nb
  
  #  empirical mean of ms as an estimate of lambdas
  lambda_ave_A_sim<-mean(lambda_A_sim)
  lambda_ave_B_sim<-mean(lambda_B_sim)
  
  cnv_sim<-(lambda_A_sim/lambda_B_sim)*Nb
  cnv_ave_sim<-mean(cnv_sim)
  
  
  ##NonPVar
  cnv_ave_cond_total_var<-var(cnv_sim)/rep_num
  
  # bootstrapping total variance
  cnv_ave_multinom_var<-mean(apply(cnv_multinom,1,var))/rep_num
  
  if (!is.na(cnv_ave_cond_total_var)){
    ci_cond<-c(cnv_ave_sim-qt(0.975,rep_num-1)*sqrt(cnv_ave_cond_total_var),cnv_ave_sim+qt(0.975,rep_num-1)*sqrt(cnv_ave_cond_total_var))
    
  } else {ci_cond<-NA}
  
  if (!is.na(cnv_ave_multinom_var)){
    ci_multinom<-c(cnv_ave_sim-qnorm(0.975)*sqrt(cnv_ave_multinom_var),cnv_ave_sim+qnorm(0.975)*sqrt(cnv_ave_multinom_var))
    
  } else {ci_multinom<-NA}
  
  
  # delta
  ## delta method
  pos_partition_aveA<-mean(x_A_num)
  partition_aveA<-mean(partitionnrA)
  
  pos_partition_aveB<-mean(x_B_num)
  partition_aveB<-mean(partitionnrB)
  
  lambda_aveA<--log(1-pos_partition_aveA/partition_aveA)
  lambda_aveB<--log(1-pos_partition_aveB/partition_aveB)
  
  cnv_delta_ave<-(lambda_aveA/lambda_aveB)*Nb
  
  var_log_cnv<-((1-exp(-lambda_aveA))/(partition_aveA*(lambda_aveA^2)*exp(-lambda_aveA))+(1-exp(-lambda_aveB))/(partition_aveB*(lambda_aveB^2)*exp(-lambda_aveB)))/rep_num
  
  delta_cnv_CIU<-exp(log(lambda_aveA/lambda_aveB)+1.96*sqrt(var_log_cnv))*Nb
  delta_cnv_CIL<-exp(log(lambda_aveA/lambda_aveB)-1.96*sqrt(var_log_cnv))*Nb
  delta_cnv_CI<-c(cnv_delta_ave,delta_cnv_CIL,delta_cnv_CIU)
  
  # glmm (please refer to Matthijs Vynck's paper Flexible analysis of digital PCR experiments using generalized linear mixed models)
  # same well
  target_data<-data.frame(gene='target gene', well='same well',pos=x_A_num,neg=partitionnrA-x_A_num,target=T)
  reference_data<-data.frame(gene='reference gene', well='same well',pos=x_B_num,neg=partitionnrB-x_B_num,target=F)
  data_combined<-rbind(target_data,reference_data)
  db_glmm<-generate.df(data_combined)
  model.rand<-glmer(Y~dummy(target,"no")+(1|well),data=db_glmm,family=binomial(link='cloglog'),nAGQ=1,verbose=F)
  
  glmm_cnv[1]<-(exp(-(model.rand@beta[2])))*Nb
  glmm_cnv[2]<-(exp(-(model.rand@beta[2])-1.96*sqrt(vcov(model.rand)[2,2])))*Nb
  glmm_cnv[3]<-(exp(-(model.rand@beta[2])+1.96*sqrt(vcov(model.rand)[2,2])))*Nb
  
  glmm_var<-vcov(model.rand)[2,2]
  
  
  result_cnv<-list(cnv_sim=cnv_ave_sim,
                   ci_cond=ci_cond,
                   ci_multinom=ci_multinom,
                   delta_cnv_CI=delta_cnv_CI,
                   glmm_cnv=glmm_cnv)
  return(result_cnv)
}