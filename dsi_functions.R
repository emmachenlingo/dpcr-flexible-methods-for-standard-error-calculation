dsi_cal<-function(Q1_num,Q2_num,Q3_num,Q4_num,partitionnr,B=1000) {
  
  rep_num<-length(Q1_num)
  DSI_multinom_boots<-matrix(0,nrow=rep_num,ncol=B)
  
  for (i in 1:rep_num){
    # random sampling from Qs (assuming a multinomial distribution of Qs,maybe can use multinomial function to sample instead)
    random_droplet=rmultinom(B, size = partitionnr[i], prob = c(Q1_num[i]/partitionnr[i],Q2_num[i]/partitionnr[i],Q3_num[i]/partitionnr[i],Q4_num[i]/partitionnr[i]))
    Q1_num_multinom_boots=random_droplet[1,]
    Q2_num_multinom_boots=random_droplet[2,]
    Q3_num_multinom_boots=random_droplet[3,]
    Q4_num_multinom_boots=random_droplet[4,]
    
    
    x_A_num_multinom_boots<-partitionnr[i]*Q1_num_multinom_boots/(Q1_num_multinom_boots+Q3_num_multinom_boots)
    x_B_num_multinom_boots<-partitionnr[i]*Q4_num_multinom_boots/(Q4_num_multinom_boots+Q3_num_multinom_boots)
    x_AB_num_multinom_boots<-Q2_num_multinom_boots-Q1_num_multinom_boots*Q4_num_multinom_boots/Q3_num_multinom_boots
    
    lambda_AB_multinom_boots<--log(1-x_AB_num_multinom_boots/partitionnr[i])
    m_AB_multinom_boots<-partitionnr[i]*lambda_AB_multinom_boots
    m_A_multinom_boots<-partitionnr[i]*(-log(1-x_A_num_multinom_boots/partitionnr[i]))
    m_B_multinom_boots<-partitionnr[i]*(-log(1-x_B_num_multinom_boots/partitionnr[i]))
    
    DSI_multinom_boots[i,]<-((m_A_multinom_boots+m_B_multinom_boots)/2)/((m_A_multinom_boots+m_B_multinom_boots)/2+m_AB_multinom_boots) 
  }
  
  
  
  
  result_list<-list(
    DSI_multinom_boots=DSI_multinom_boots
    # multinomial sampling
  )
  return(result_list)
  
}


dsi_bootstrap<-function(Q1_num,Q2_num,Q3_num,Q4_num,partitionnr,B=1000){
  # simulation results
  rep_num<-length(Q1_num)
  
  # number of retained partitions
  x_A_num<-partitionnr*Q1_num/(Q1_num+Q3_num)
  x_B_num<-partitionnr*Q4_num/(Q4_num+Q3_num)
  x_AB_num<-Q2_num-Q1_num*Q4_num/Q3_num
  
  lambda_AB_sim<--log(1-x_AB_num/partitionnr)
  m_AB_sim<-partitionnr*lambda_AB_sim
  m_A_sim<-partitionnr*(-log(1-x_A_num/partitionnr))
  m_B_sim<-partitionnr*(-log(1-x_B_num/partitionnr))
  
  DSI_sim<-((m_A_sim+m_B_sim)/2)/((m_A_sim+m_B_sim)/2+m_AB_sim)
  
  result<-dsi_cal(Q1_num,Q2_num,Q3_num,Q4_num,partitionnr)
  DSI_multinom<-result$DSI_multinom_boots
  
  DSI_ave_cond_total_var<-var(DSI_sim)/rep_num
  
  # binomial assumption
  DSI_ave_multinom_var<-mean(apply(DSI_multinom,1,var))/rep_num
  
  
  # empirical mean of ms as an estimate of lambdas
  DSI_ave_sim<-mean(DSI_sim)
  
  ci_cond<-DSI_ave_sim+c(qt(0.025,rep_num-1),qt(0.975,rep_num-1))*sqrt(DSI_ave_cond_total_var)
  ci_binom<-DSI_ave_sim+c(-1.96,1.96)*sqrt(DSI_ave_multinom_var)
  
  
  
  result_dsi<-list(DSI_sim=DSI_ave_sim,
                   ci_cond=ci_cond,
                   ci_binom=ci_binom)
  
  return(result_dsi)
}