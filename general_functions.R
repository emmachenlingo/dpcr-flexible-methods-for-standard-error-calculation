## general codes (include the bootstrap procedure when no replicates are available)
est.Bootstrapping <- function(data,Nb=2,estimator='abs quant',replicates=FALSE,B=1000) {
  uniqueID<-unique(data$SampleID)
  sample_num<-length(uniqueID)
  if (!replicates) {
    result_boots<-matrix(0,nrow=sample_num,ncol=B)
    result<-numeric(sample_num)
    
    if (estimator=='abs quant') {
      for (i in 1:length(uniqueID)){
        x1<-data[data$SampleID==uniqueID[i],'Positives']
        p<-data[data$SampleID==uniqueID[i],'Positives']+data[data$SampleID==uniqueID[i],'Negatives']
        if (!(x1%%1==00)){
          return(cat("Error: x1 should be the number of positive partitions. Please input integers"))
        }
        lambda1<-(-log(1-x1/p))
        result[i]<-lambda1
        m1_round<-round(lambda1*p)
        
        for (j in 1:B) {
          x1_boots<-sample(1:p,m1_round,replace=T)
          x1_set_boots<-length(unique(x1_boots))
          lambda1_boots<-(-log(1-x1_set_boots/p))
          result_boots[i,j]<-lambda1_boots
        }
      }
      
    }
    
    if (estimator=='copy number') {
      for (i in 1:length(uniqueID)){
        x1<-data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Positives']
        x2<-data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Positives']
        p1<-data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Positives']+data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Negatives']
        p2<-data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Positives']+data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Negatives']
        
        if (!((x1%%1==00) & (x2%%1==0))){
          return(cat("Error: x1 and x2 should be the number of positive partitions. Please input integers"))
        }
        lambda1<-(-log(1-x1/p1))
        lambda2<-(-log(1-x2/p2))
        result[i]<-lambda1/lambda2*Nb
        m1_round<-round(lambda1*p1)
        m2_round<-round(lambda2*p2)
        
        for (j in 1:B){
          x1_boots<-sample(1:p1,m1_round,replace=T)
          x1_set_boots<-length(unique(x1_boots))
          lambda1_boots<-(-log(1-x1_set_boots/p1))
          x2_boots<-sample(1:p2,m2_round,replace=T)
          x2_set_boots<-length(unique(x2_boots))
          lambda2_boots<-(-log(1-x2_set_boots/p2))
          result_boots[i,j]<-lambda1_boots/lambda2_boots*Nb
        }
      }
      
      
    }
    
    if (estimator=='frac abund') {
      for (i in 1:length(uniqueID)){
        x1<-data[data$SampleID==uniqueID[i],'Mutant_Positives']
        x2<-data[data$SampleID==uniqueID[i],'Wild_Positives']
        p<-data[data$SampleID==uniqueID[i],'Analyzed']
        
        if (!((x1%%1==00) & (x2%%1==0))){
          return(cat("Error: x1 and x2 should be the number of positive partitions. Please input integers"))
        }
        lambda1<-(-log(1-x1/p))
        lambda2<-(-log(1-x2/p))
        result[i]<-lambda1/(lambda1+lambda2)
        m1_round<-round(lambda1*p)
        m2_round<-round(lambda2*p)
        
        for (j in 1:B){
          x1_boots<-sample(1:p,m1_round,replace=T)
          x1_set_boots<-length(unique(x1_boots))
          lambda1_boots<-(-log(1-x1_set_boots/p))
          x2_boots<-sample(1:p,m2_round,replace=T)
          x2_set_boots<-length(unique(x2_boots))
          lambda2_boots<-(-log(1-x2_set_boots/p))
          result_boots[i,j]<-lambda1_boots/(lambda1_boots+lambda2_boots)
        }
      }
      
    }
    
    if (estimator=='DSI') {
      for (i in 1:length(uniqueID)){
        x1<-data[data$SampleID==uniqueID[i],'Single_Positive_A']
        x2<-data[data$SampleID==uniqueID[i],'Double_Positive']
        x3<-data[data$SampleID==uniqueID[i],'Negative']
        x4<-data[data$SampleID==uniqueID[i],'Single_Positive_B']
        p<-x1+x2+x3+x4
        
        if (!((x1%%1==00) & (x2%%1==0) & (x3%%1==0)  & (x4%%1==0))){
          return(cat("Error: x1 and x4 should be the number of single positive partitions.\n x2 should be the number of double positive partitions and x3 should be the number of empty partitions.\n Please input integers"))
        }   
        
        x_A<-p*x1/(x1+x3)
        x_B<-p*x4/(x4+x3)
        x_AB<-x2-x1*x4/x3
        lambda1<-(-log(1-x_A/p))
        lambda2<-(-log(1-x_B/p))
        lambda3<-(-log(1-x_AB/p))
        result[i]<-((lambda1+lambda2)/2)/((lambda1+lambda2)/2+lambda3)
        m1_round<-round(lambda1*p)
        m2_round<-round(lambda2*p)
        m3_round<-round(lambda3*p)
        
        for (j in 1:B) {
          x1_boots<-sample(1:p,m1_round,replace=T)
          x1_set_boots<-length(unique(x1_boots))
          x2_boots<-sample(1:p,m2_round,replace=T)
          x2_set_boots<-length(unique(x2_boots))
          x3_boots<-sample(1:p,m3_round,replace=T)
          x3_set_boots<-length(unique(x3_boots))
          
          ## A, no B and no intact
          Q1_boots<-setdiff(setdiff(x1_boots,x2_boots),x3_boots)
          Q1_num_boots<-length(Q1_boots)
          
          ## B, no A and no intact
          Q4_boots<-setdiff(setdiff(x2_boots,x1_boots),x3_boots)
          Q4_num_boots<-length(Q4_boots)
          
          ## no A, no B and no intact
          Q3_boots<-setdiff(1:p,union(union(x1_boots,x2_boots),x3_boots))
          Q3_num_boots<-length(Q3_boots)
          
          # the rest should be Q2
          Q2_boots<-setdiff(1:p,union(union(Q1_boots,Q3_boots),Q4_boots))
          Q2_num_boots<-length(Q2_boots)
          
          x_A_boots<-p*Q1_num_boots/(Q1_num_boots+Q3_num_boots)
          x_B_boots<-p*Q4_num_boots/(Q4_num_boots+Q3_num_boots)
          x_AB_boots<-Q2_num_boots-Q1_num_boots*Q4_num_boots/Q3_num_boots
          
          lambda1_boots<-(-log(1-x_A_boots/p))
          lambda2_boots<-p*(-log(1-x_B_boots/p))
          lambda3_boots<-p*(-log(1-x_AB_boots/p))
          
          result_boots[i,j]<-((lambda1_boots+lambda2_boots)/2)/((lambda1_boots+lambda2_boots)/2+lambda3_boots)
        }
      }
      
    }
    
    var_boots<-apply(result_boots,1,var)
    CI_lower<-result-1.96*sqrt(var_boots)
    CI_upper<-result+1.96*sqrt(var_boots)
    df_new<-data.frame(SampleID=uniqueID,Estimates=round(result,4),CIL=round(CI_lower,4),CIU=round(CI_upper,4))
    return(df_new)
  } else {
    
    if (estimator=='abs quant') {
      result_sample_list<-vector(mode='list',length(uniqueID))
      for (i in 1:length(uniqueID)){
        result_sample<-matrix(0,nrow=4,ncol=3)
        x_num<-data[data$SampleID==uniqueID[i],'Positives']
        x_num_new<-x_num[!is.na(x_num)]
        partitionnr<-data[data$SampleID==uniqueID[i],'Positives']+data[data$SampleID==uniqueID[i],'Negatives']
        
        results<-abs_quant_bootstrap(x_num_new,partitionnr,B=1000)
        result_sample[1,1]<-results$lambda_ave_sim
        result_sample[1,2]<-results$ci_cond[1]
        result_sample[1,3]<-results$ci_cond[2]
        
        result_sample[2,1]<-results$lambda_ave_sim
        result_sample[2,2]<-results$ci_multinom[1]
        result_sample[2,3]<-results$ci_multinom[2]
        
        result_sample[3,1]<-results$delta_abs_quant_CI[1]
        result_sample[3,2]<-results$delta_abs_quant_CI[2]
        result_sample[3,3]<-results$delta_abs_quant_CI[3]
        
        result_sample[4,1]<-results$glmm_est[1]
        result_sample[4,2]<-results$glmm_est[2]
        result_sample[4,3]<-results$glmm_est[3]
        
        result_sample_list[[i]]<-result_sample
        
      }
      
      
      df_new<-NULL
      for (j in 1:length(result_sample_list)){
        df_new<-rbind(df_new,cbind(uniqueID[j],result_sample_list[[j]],c('NonPVar','BinomVar','Delta','GLMM')))
      }
      df_new<-as.data.frame(df_new)
      colnames(df_new)<-c('SampleID','Estimates','CIL','CIU','Method')
      df_new$Estimates<-as.numeric(df_new$Estimates)
      df_new$Method<-factor(df_new$Method,levels=c('NonPVar','BinomVar','Delta','GLMM'))
      df_new$CIL<-as.numeric(df_new$CIL)
      df_new$CIU<-as.numeric(df_new$CIU)
    }
    
    if (estimator=='copy number singleplex') {
      result_sample_list<-vector(mode='list',length(uniqueID))
      for (i in 1:length(uniqueID)){
        result_sample<-matrix(0,nrow=4,ncol=3)
        x_A_num<-data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Positives']
        x_A_num_new<-x_A_num[!is.na(x_A_num)]
        x_B_num_new<-data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Positives']
        partitionnrA<-data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Positives']+data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Negatives']
        partitionnrB<-data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Positives']+data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Negatives']
        
        results<-cnv_bootstrap_singleplex(x_A_num_new,x_B_num_new,partitionnrA,partitionnrB,Nb=2,B=1000)
        result_sample[1,1]<-results$cnv_sim
        result_sample[1,2]<-results$ci_cond[1]
        result_sample[1,3]<-results$ci_cond[2]
        
        result_sample[2,1]<-results$cnv_sim
        result_sample[2,2]<-results$ci_multinom[1]
        result_sample[2,3]<-results$ci_multinom[2]
        
        result_sample[3,1]<-results$delta_cnv_CI[1]
        result_sample[3,2]<-results$delta_cnv_CI[2]
        result_sample[3,3]<-results$delta_cnv_CI[3]
        
        result_sample[4,1]<-results$glmm_cnv[1]
        result_sample[4,2]<-results$glmm_cnv[2]
        result_sample[4,3]<-results$glmm_cnv[3]
        
        result_sample_list[[i]]<-result_sample
        
      }
      
      
      df_new<-NULL
      for (j in 1:length(result_sample_list)){
        df_new<-rbind(df_new,cbind(uniqueID[j],result_sample_list[[j]],c('NonPVar','BinomVar','Delta','GLMM')))
      }
      df_new<-as.data.frame(df_new)
      colnames(df_new)<-c('SampleID','Estimates','CIL','CIU','Method')
      df_new$Estimates<-round(as.numeric(df_new$Estimates),4)
      df_new$Method<-factor(df_new$Method,levels=c('NonPVar','BinomVar','Delta','GLMM'))
      df_new$CIL<-round(as.numeric(df_new$CIL),4)
      df_new$CIU<-round(as.numeric(df_new$CIU),4)
    }
    
    if (estimator=='copy number duplex') {
      result_sample_list<-vector(mode='list',length(uniqueID))
      for (i in 1:length(uniqueID)){
        result_sample<-matrix(0,nrow=4,ncol=3)
        x_A_num<-data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Positives']
        x_A_num_new<-x_A_num[!is.na(x_A_num)]
        x_B_num<-data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Positives']
        x_B_num_new<-x_B_num[!is.na(x_B_num)]
        partitionnrA<-data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Positives']+data[(data$SampleID==uniqueID[i])&(data$Type=='Target'),'Negatives']
        partitionnrB<-data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Positives']+data[(data$SampleID==uniqueID[i])&(data$Type=='Reference'),'Negatives']
        
        results<-cnv_bootstrap_duplex(x_A_num_new,x_B_num_new,partitionnrA,partitionnrB,Nb=2,B=1000)
        result_sample[1,1]<-results$cnv_sim
        result_sample[1,2]<-results$ci_cond[1]
        result_sample[1,3]<-results$ci_cond[2]
        
        result_sample[2,1]<-results$cnv_sim
        result_sample[2,2]<-results$ci_multinom[1]
        result_sample[2,3]<-results$ci_multinom[2]
        
        result_sample[3,1]<-results$delta_cnv_CI[1]
        result_sample[3,2]<-results$delta_cnv_CI[2]
        result_sample[3,3]<-results$delta_cnv_CI[3]
        
        result_sample[4,1]<-results$glmm_cnv[1]
        result_sample[4,2]<-results$glmm_cnv[2]
        result_sample[4,3]<-results$glmm_cnv[3]
        
        result_sample_list[[i]]<-result_sample
        
      }
      df_new<-NULL
      for (j in 1:length(result_sample_list)){
        df_new<-rbind(df_new,cbind(uniqueID[j],result_sample_list[[j]],c('NonPVar','BinomVar','Delta','GLMM')))
      }
      df_new<-as.data.frame(df_new)
      colnames(df_new)<-c('SampleID','Estimates','CIL','CIU','Method')
      df_new$Estimates<-round(as.numeric(df_new$Estimates),4)
      df_new$Method<-factor(df_new$Method,levels=c('NonPVar','BinomVar','Delta','GLMM'))
      df_new$CIL<-round(as.numeric(df_new$CIL),4)
      df_new$CIU<-round(as.numeric(df_new$CIU),4)
      
    }
    
    if (estimator=='frac abund') {
      result_sample_list<-vector(mode='list',length(uniqueID))
      for (i in 1:length(uniqueID)){
        result_sample<-matrix(0,nrow=2,ncol=3)
        x_A_num<-data[data$SampleID==uniqueID[i],'Mutant_Positives']
        x_A_num_new<-x_A_num[!is.na(x_A_num)]
        x_B_num<-data[data$SampleID==uniqueID[i],'Wild_Positives']
        x_B_num_new<-x_B_num[!is.na(x_B_num)]
        
        
        partitionnr<-data[data$SampleID==uniqueID[i],'Analyzed']
        
        results<-ml_bootstrap(x_A_num_new,x_B_num_new,partitionnr,B=1000)
        result_sample[1,1]<-results$ml_sim
        result_sample[1,2]<-results$ci_cond[1]
        result_sample[1,3]<-results$ci_cond[2]
        
        result_sample[2,1]<-results$ml_sim
        result_sample[2,2]<-results$ci_binom[1]
        result_sample[2,3]<-results$ci_binom[2]
        
        result_sample_list[[i]]<-result_sample
        
      }
      df_new<-NULL
      for (j in 1:length(result_sample_list)){
        df_new<-rbind(df_new,cbind(uniqueID[j],result_sample_list[[j]],c('NonPVar','BinomVar')))
      }
      df_new<-as.data.frame(df_new)
      colnames(df_new)<-c('SampleID','Estimates','CIL','CIU','Method')
      df_new$Estimates<-round(as.numeric(df_new$Estimates),4)
      df_new$Method<-factor(df_new$Method,levels=c('NonPVar','BinomVar'))
      df_new$CIL<-round(as.numeric(df_new$CIL),4)
      df_new$CIU<-round(as.numeric(df_new$CIU),4)
      
    }
    
    if (estimator=='DSI') {
      result_sample_list<-vector(mode='list',length(uniqueID))
      for (i in 1:length(uniqueID)){
        result_sample<-matrix(0,nrow=2,ncol=3)
        Q1_num<-data[data$SampleID==uniqueID[i],'Single_Positive_A']
        Q2_num<-data[data$SampleID==uniqueID[i],'Double_Positive']
        Q3_num<-data[data$SampleID==uniqueID[i],'Negative']
        Q4_num<-data[data$SampleID==uniqueID[i],'Single_Positive_B']
        
        
        partitionnr<-Q1_num+Q2_num+Q3_num+Q4_num
        
        results<-dsi_bootstrap(Q1_num,Q2_num,Q3_num,Q4_num,partitionnr,B=1000)
        result_sample[1,1]<-results$DSI_sim
        result_sample[1,2]<-results$ci_cond[1]
        result_sample[1,3]<-results$ci_cond[2]
        
        result_sample[2,1]<-results$DSI_sim
        result_sample[2,2]<-results$ci_binom[1]
        result_sample[2,3]<-results$ci_binom[2]
        
        result_sample_list[[i]]<-result_sample
        
      }
      df_new<-NULL
      for (j in 1:length(result_sample_list)){
        df_new<-rbind(df_new,cbind(uniqueID[j],result_sample_list[[j]],c('NonPVar','BinomVar')))
      }
      df_new<-as.data.frame(df_new)
      colnames(df_new)<-c('SampleID','Estimates','CIL','CIU','Method')
      df_new$Estimates<-round(as.numeric(df_new$Estimates),4)
      df_new$Method<-factor(df_new$Method,levels=c('NonPVar','BinomVar'))
      df_new$CIL<-round(as.numeric(df_new$CIL),4)
      df_new$CIU<-round(as.numeric(df_new$CIU),4)
      
    }
    
    return(df_new)
    
    
  }
  
}