MutLoad_cal<-function(m_A_sim,m_B_sim,p_retained,B=1000) {
  
  rep_num<-length(m_A_sim)
  
  # record the results of multinomial distributed Qs
  ml_multinom_boots<-matrix(nrow = rep_num,ncol=B)
  
  mu_A_est<-mean(m_A_sim)
  mu_B_est<-mean(m_B_sim)
  
  for (i in 1:rep_num){
    # random sampling from Qs (assuming a multinomial distribution of Qs,maybe can use multinomial function to sample instead)
    pi_A_i<-1-exp(-mu_A_est/p_retained[i])
    x_A_multinom_boots=rbinom(B,p_retained[i],pi_A_i)
    m_A_multinom_boots<-p_retained[i]*(-log(1-x_A_multinom_boots/p_retained[i]))
    
    pi_B_i<-1-exp(-mu_B_est/p_retained[i])
    x_B_multinom_boots=rbinom(B,p_retained[i],pi_B_i)
    m_B_multinom_boots<-p_retained[i]*(-log(1-x_B_multinom_boots/p_retained[i]))
    
    ml_multinom_boots[i,]<-m_B_multinom_boots/(m_A_multinom_boots+m_B_multinom_boots)
    
  }
  
  
  result_list<-list(
    ml_multinom_boots=ml_multinom_boots
    # binomial sampling
  )
  return(result_list)
  
}


ml_bootstrap<-function(x_A_num,x_B_num,partitionnr,B=1000){
  # simulation results
  rep_num<-length(x_A_num)
  
  ml_sim<-numeric(rep_num)
  
  lambda_A_sim<--log(1-x_A_num/partitionnr)
  m_A_sim<-partitionnr*lambda_A_sim
  
  lambda_B_sim<--log(1-x_B_num/partitionnr)
  m_B_sim<-partitionnr*lambda_B_sim
  
  ml_sim<-m_B_sim/(m_A_sim+m_B_sim)
  
  result<-MutLoad_cal(m_A_sim,m_B_sim,partitionnr,B)
  
  ml_multinom<-result$ml_multinom_boots
  
  ml_ave_cond_total_var<-var(ml_sim)/rep_num
  
  # binomial assumption
  ml_ave_multinom_var<-mean(apply(ml_multinom,1,var))/rep_num
  
  
  # empirical mean of ms as an estimate of lambdas
  ml_ave_sim<-mean(ml_sim)
  
  
  if (!all(x_B_num==0)){
    ci_cond<-ml_ave_sim+c(qt(0.025,rep_num-1),qt(0.975,rep_num-1))*sqrt(ml_ave_cond_total_var)
    ci_binom<-ml_ave_sim+c(-1.96,1.96)*sqrt(ml_ave_multinom_var)
  } else {ci_cond=c(0,0)
  ci_binom=c(0,0)}
  
  
  result_ml<-list(ml_sim=ml_ave_sim,
                  ci_cond=ci_cond,
                  ci_binom=ci_binom)
  
  return(result_ml)
  
  
}