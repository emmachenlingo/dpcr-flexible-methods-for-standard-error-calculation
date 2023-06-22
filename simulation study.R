library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)

## 2 dimensions
molecules<-c(100,1000,10000)
dsi_s<-c(0.25,0.5,0.75)
dsi_new<-NULL
dsi_bruner_new<-NULL
for (m in molecules){
  for (dsi in dsi_s){
    dsi_2dim<-numeric(100)
    dsi_bruner<-numeric(100)
    for (l in 1:100){
      print(l)
      tryCatch( { 
        p_singles<-c()
        n_targets<-2
        Qs<-c()
        ## total number of partitions
        ##Qs, 2 dimension as an example
        ## Q1: only positive in dimension 1
        ## Q2: only positive in dimension 2
        ## Q3: positive in 1,2
        ## Q4: empty partitions
        dsi<-dsi
        tot_molecules<-m
        p_tot<-20000
        sheared<-dsi*tot_molecules
        ## simulation study, to pair them up (start from the most complicated, the rest will be paired up automatically). not necessary to pair up even, just take half
        ## combinations: 
        combios<-rbind(c(1,0),c(0,1))
        p_combios<-1/nrow(combios)
        sampled_combios<-combios[sample(1:nrow(combios),sheared,replace = TRUE),]
        sheared_combios<-matrix(1,nrow=sheared,ncol=n_targets)-sampled_combios
        
        par_pos<-sample(1:p_tot,nrow(sampled_combios)+nrow(sheared_combios)+tot_molecules*(1-dsi),replace=TRUE)
        par_empty<-setdiff(c(1:p_tot), par_pos)
        non_emp_pars<-rbind(sampled_combios,sheared_combios,matrix(1,nrow=tot_molecules*(1-dsi),ncol=n_targets))
        Qs[2^n_targets]<-length(unique(par_empty))
        
        non_emp_pars_pos<-data.frame(cbind(par_pos,non_emp_pars))
        non_emp_pars_pos_agg<-data.frame(non_emp_pars_pos %>% 
                                           group_by(par_pos) %>% 
                                           summarise(Frequency_Var1 = sum(V2),
                                                     Frequency_Var2 = sum(V3)))
        
        ind_mat<-rbind(c(1,0),c(0,1),c(1,1))
        
        
        for (j in 1:(2^(n_targets)-1)){
          Qs[j]<-length(non_emp_pars_pos_agg[(((1*as.logical(non_emp_pars_pos_agg[,2]))==ind_mat[j,1]) & ((1*as.logical(non_emp_pars_pos_agg[,3]))==ind_mat[j,2])),1])
        }
        
        Q_len<-length(Qs)
        for (i in 1:n_targets){
          p_singles[i]<-Qs[i]/(Qs[i]+Qs[Q_len])
        }
        
        ## for the 2 dimensions
        p_doubles<-((Qs[3]*Qs[4]-Qs[1]*Qs[2])/Qs[4])/p_tot
        

        lambda1<-(-log(1-p_singles[1]))
        lambda2<-(-log(1-p_singles[2]))
        
        lambda12<-(-log(1-p_doubles))
        
        dsi_2dim[l]<-((lambda1+lambda2)/2)/((lambda1+lambda2)/2+lambda12)
        dsi_bruner[l]<-((Qs[1]+Qs[2])/2)/((Qs[1]+Qs[2])/2+Qs[3])
      }, error = function(e) {
        dsi_2dim[l]<-NA
        dsi_bruner[l]<-NA
      })   
      
    }
    dsi_2dim_tmp<-cbind(intactness=1-dsi_2dim,molecules=m,true_intactness=1-dsi)
    dsi_new<-rbind(dsi_new,dsi_2dim_tmp)
    dsi_bruner_tmp<-cbind(intactness=1-dsi_bruner,molecules=m,true_intactness=1-dsi)
    dsi_bruner_new<-rbind(dsi_bruner_new,dsi_bruner_tmp)
    
  }
  
}

dsi_bruner_new<-data.frame(dsi_bruner_new)
dsi_bruner_new<-na.omit(dsi_bruner_new)
dsi_bruner_new<-dsi_bruner_new[dsi_bruner_new$intactness!=0,]
dsi_bruner_new$molecules<-as.factor(dsi_bruner_new$molecules)
dsi_bruner_new$true_dsi<-as.factor(dsi_bruner_new$true_intactness)

true_025<-dsi_bruner_new[dsi_bruner_new$true_intactness==0.25,]
true_05<-dsi_bruner_new[dsi_bruner_new$true_intactness==0.5,]
true_075<-dsi_bruner_new[dsi_bruner_new$true_intactness==0.75,]

plot1 <- ggplot(data = true_025, aes(molecules, intactness)) + geom_boxplot()+geom_hline(yintercept = 0.25, linetype = "dashed", color = "blue")+theme_classic() 
plot2 <- ggplot(data = true_05, aes(molecules, intactness)) + geom_boxplot()+geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")+theme_classic() 
plot3 <- ggplot(data = true_075, aes(molecules, intactness)) + geom_boxplot()+geom_hline(yintercept = 0.75, linetype = "dashed", color = "blue")+theme_classic() 

# Combine the plots horizontally
p<-grid.arrange(plot1, plot2, plot3, nrow = 1)

ggsave("dsi_bruner.png",plot = p,width = 20,
       height = 7, dpi = 350)





## higher order multiplexing
## broken ones 
molecules<-c(100,1000,10000)
dsi_s<-c(0.25,0.5,0.75)
dsi_new<-NULL
for (m in molecules){
  for (dsi in dsi_s){
    dsi_3dim<-numeric(100)
    for (l in 1:100){
      print(l)
      tryCatch( { 
        p_singles<-c()
        p_doubles<-c()
        n_targets<-3
        Qs<-c()
        ## total number of partitions
        ##Qs, 3 dimension as an example
        ## Q1: only positive in dimension 1
        ## Q2: only positive in dimension 2
        ## Q3: only positive in dimension 3
        ## Q4: positive in 1,2
        ## Q5: positive in 1,3
        ## Q6: positive in 2,3
        ## Q7: positive in 1,2,3
        ## Q8: empty partitions
        dsi<-dsi
        tot_molecules<-m
        p_tot<-20000
        sheared<-dsi*tot_molecules
        ## simulation study, to pair them up (start from the most complicated, the rest will be paired up automatically). not necessary to pair up even, just take half
        ## combinations: 
        combios<-rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(0,1,1))
        p_combios<-1/nrow(combios)
        sampled_combios<-combios[sample(1:nrow(combios),sheared,replace = TRUE),]
        sheared_combios<-matrix(1,nrow=sheared,ncol=n_targets)-sampled_combios
        sheared_combios_imp<-((sheared_combios[,1]==1) & (sheared_combios[,3]==1))
        sheared_combios<-sheared_combios[!sheared_combios_imp,]
        sheared_combios<-rbind(sheared_combios,matrix(c(1,0,0),nrow=sum(sheared_combios_imp),ncol=3,byrow=TRUE),matrix(c(0,0,1),nrow=sum(sheared_combios_imp),ncol=3,byrow=TRUE))
        
        par_pos<-sample(1:p_tot,nrow(sampled_combios)+nrow(sheared_combios)+tot_molecules*(1-dsi),replace=TRUE)
        par_empty<-setdiff(c(1:p_tot), par_pos)
        non_emp_pars<-rbind(sampled_combios,sheared_combios,matrix(1,nrow=tot_molecules*(1-dsi),ncol=n_targets))
        Qs[8]<-length(unique(par_empty))
        
        non_emp_pars_pos<-data.frame(cbind(par_pos,non_emp_pars))
        non_emp_pars_pos_agg<-data.frame(non_emp_pars_pos %>% 
                                           group_by(par_pos) %>% 
                                           summarise(Frequency_Var1 = sum(V2),
                                                     Frequency_Var2 = sum(V3),
                                                     Frequency_Var3 = sum(V4)))
        
        ind_mat<-rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
        
        
        for (j in 1:(2^(n_targets)-1)){
          Qs[j]<-length(non_emp_pars_pos_agg[(((1*as.logical(non_emp_pars_pos_agg[,2]))==ind_mat[j,1]) & ((1*as.logical(non_emp_pars_pos_agg[,3]))==ind_mat[j,2]) & ((1*as.logical(non_emp_pars_pos_agg[,4]))==ind_mat[j,3])),1])
        }
        
        Q_len<-length(Qs)
        for (i in 1:n_targets){
          p_singles[i]<-Qs[i]/(Qs[i]+Qs[Q_len])
        }
        
        ## for the 2 dimensions
        doubles<-function(x,ratio,single_targets){
          #ratio-(single_targets[1]*single_targets[2]+
          #         x*(1+single_targets[1]+single_targets[2]+single_targets[1]*single_targets[2]))/((1-single_targets[1])*(1-single_targets[2])*(1-x))
          ratio-(single_targets[1]*single_targets[2]+x-single_targets[1]*single_targets[2]*x)/((1-single_targets[1])*(1-single_targets[2])*(1-x))
        }
        
        tribles_3dim<-function(x,empty_prop,p_singles,p_doubles){
          empty_prop-prod(1-p_singles)*prod(1-p_doubles)*(1-x)
        }
        
        ## for 3 dimensions
        k=n_targets
        p_doubles<-c()
        for (i in 1:(n_targets-1)){
          for (j in (i+1):n_targets){
            k=k+1
            if (i==1 & j==3){
            p_doubles<-c(p_doubles,0)}
            else{
            roots<-uniroot(doubles, interval=c(0,1),ratio=Qs[k]/Qs[Q_len],single_targets=c(p_singles[i],p_singles[j]))
            p_doubles<-c(p_doubles,roots$root)}
          }
          
        }
        
        p_tribles_3dim<-uniroot(tribles_3dim, interval=c(0,1),empty_prop=Qs[Q_len]/p_tot,p_singles=p_singles,p_doubles=p_doubles)$root
        
        lambda1<-(-log(1-p_singles[1]))
        lambda2<-(-log(1-p_singles[2]))
        lambda3<-(-log(1-p_singles[3]))
        lambda12<-(-log(1-p_doubles[1]))
        lambda13<-(-log(1-p_doubles[2]))
        lambda23<-(-log(1-p_doubles[3]))
        
        lambda123<-(-log(1-p_tribles_3dim))
        
        dsi_3dim[l]<-(lambda12+lambda23+(lambda1+lambda3-lambda12-lambda23+lambda2)/3)/(lambda12+lambda23+(lambda1+lambda3-lambda12-lambda23+lambda2)/3+lambda123)
      }, error = function(e) {
        dsi_3dim[l]<-NA
      })   
      
    }
    dsi_3dim_tmp<-cbind(intactness=1-dsi_3dim,molecules=m,true_intactness=1-dsi)
    dsi_new<-rbind(dsi_new,dsi_3dim_tmp)
  }
  
}


## 4 dimensions
molecules<-c(100,1000,10000)
dsi_s<-c(0.25,0.5,0.75)
dsi_new<-NULL
for (m in molecules){
  for (dsi in dsi_s){
    dsi_4dim<-numeric(100)
    p_singles_mat<-NULL
    p_doubles_mat<-NULL
    p_triples_mat<-NULL
    p_quadras_vec<-NULL
    Qs_mat<-NULL
    for (l in 1:100){
      print(l)
      
      tryCatch( { 
        p_singles<-c()
        p_doubles<-c()
        n_targets<-4
        Qs<-c()
        ## total number of partitions
        ##Qs, 3 dimension as an example
        ## Q1: only positive in dimension 1
        ## Q2: only positive in dimension 2
        ## Q3: only positive in dimension 3
        ## Q4: positive in 1,2
        ## Q5: positive in 1,3
        ## Q6: positive in 2,3
        ## Q7: positive in 1,2,3
        ## Q8: empty partitions
        dsi<-dsi
        tot_molecules<-m
        p_tot<-20000
        sheared<-dsi*tot_molecules
        ## simulation study, to pair them up (start from the most complicated, the rest will be paired up automatically). not necessary to pair up even, just take half
        ## combinations: 
        combios<-rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1),c(1,1,0,0),c(0,1,1,0),c(0,0,1,1),c(1,1,1,0),c(0,1,1,1))
        p_combios<-1/nrow(combios)
        sampled_combios<-as.matrix(combios[sample(1:nrow(combios),sheared,replace = TRUE),])
        sheared_combios<-matrix(1,nrow=sheared,ncol=n_targets)-sampled_combios
        types<-apply(sheared_combios,1,function(x) paste0(x[1],x[2],x[3],x[4]))
        sheared_combios_imp<-types %in% c('1011','1001','1010','0101','1101')
        sheared_combios_imp1<-(types == '1011')
        sheared_combios_imp2<-(types == '1001')
        sheared_combios_imp3<-(types == '1010')
        sheared_combios_imp4<-(types == '0101')
        sheared_combios_imp5<-(types == '1101')
        sheared_combios<-sheared_combios[!sheared_combios_imp,]
        sheared_combios<-rbind(sheared_combios,matrix(c(1,0,0,0),nrow=sum(sheared_combios_imp1),ncol=4,byrow=TRUE),matrix(c(0,0,1,1),nrow=sum(sheared_combios_imp1),ncol=4,byrow=TRUE),
                               matrix(c(1,0,0,0),nrow=sum(sheared_combios_imp2),ncol=4,byrow=TRUE),matrix(c(0,0,0,1),nrow=sum(sheared_combios_imp2),ncol=4,byrow=TRUE),
                               matrix(c(1,0,0,0),nrow=sum(sheared_combios_imp3),ncol=4,byrow=TRUE),matrix(c(0,0,1,0),nrow=sum(sheared_combios_imp3),ncol=4,byrow=TRUE),
                               matrix(c(0,1,0,0),nrow=sum(sheared_combios_imp4),ncol=4,byrow=TRUE),matrix(c(0,0,0,1),nrow=sum(sheared_combios_imp4),ncol=4,byrow=TRUE),
                               matrix(c(1,1,0,0),nrow=sum(sheared_combios_imp5),ncol=4,byrow=TRUE),matrix(c(0,0,0,1),nrow=sum(sheared_combios_imp5),ncol=4,byrow=TRUE))
        
        
        par_pos<-sample(1:p_tot,nrow(sampled_combios)+nrow(sheared_combios)+tot_molecules*(1-dsi),replace=TRUE)
        par_empty<-setdiff(c(1:p_tot), par_pos)
        non_emp_pars<-rbind(sampled_combios,sheared_combios,matrix(1,nrow=tot_molecules*(1-dsi),ncol=n_targets))
        Qs[2^(n_targets)]<-length(unique(par_empty))
        
        
        non_emp_pars_pos<-data.frame(cbind(par_pos,non_emp_pars))
        non_emp_pars_pos_agg<-data.frame(non_emp_pars_pos %>% 
          group_by(par_pos) %>% 
          summarise(Frequency_Var1 = sum(V2),
                    Frequency_Var2 = sum(V3),
                    Frequency_Var3 = sum(V4),
                    Frequency_Var4 = sum(V5)))
        ## indicator matrix
        ind_mat<-rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1),c(1,1,0,0),c(1,0,1,0),c(1,0,0,1),c(0,1,1,0),c(0,1,0,1),
                       c(0,0,1,1),c(1,1,1,0),c(1,1,0,1),c(1,0,1,1),c(0,1,1,1),c(1,1,1,1))
        
        
        for (j in 1:(2^(n_targets)-1)){
          Qs[j]<-length(non_emp_pars_pos_agg[(((1*as.logical(non_emp_pars_pos_agg[,2]))==ind_mat[j,1]) & ((1*as.logical(non_emp_pars_pos_agg[,3]))==ind_mat[j,2]) & ((1*as.logical(non_emp_pars_pos_agg[,4]))==ind_mat[j,3]) & ((1*as.logical(non_emp_pars_pos_agg[,5]))==ind_mat[j,4])),1])
        }
        
        # Qs[11:14]<-822
        Qs_mat<-rbind(Qs_mat,Qs)
        
        ## for the doubles 
        Q_len<-length(Qs)
        for (i in 1:n_targets){
          p_singles[i]<-Qs[i]/(Qs[i]+Qs[Q_len])
        }
        
        p_singles_mat<-rbind(p_singles_mat,p_singles)
        
        ## for the 2 dimensions
        doubles<-function(x,ratio,single_targets){
          #ratio-(single_targets[1]*single_targets[2]+
          #         x*(1+single_targets[1]+single_targets[2]+single_targets[1]*single_targets[2]))/((1-single_targets[1])*(1-single_targets[2])*(1-x))
          ratio-(single_targets[1]*single_targets[2]+x-single_targets[1]*single_targets[2]*x)/((1-single_targets[1])*(1-single_targets[2])*(1-x))
        }
        
        
        ## for 3 dimensions
        k=n_targets
        p_doubles<-c()
        for (i in 1:(n_targets-1)){
          for (j in (i+1):n_targets){
            k=k+1
            if ((i==1 & j==3)|(i==1 & j==4) | (i==2 & j==4)){
              p_doubles<-c(p_doubles,0)}
            else{
            roots<-uniroot(doubles, interval=c(0,1),ratio=Qs[k]/Qs[Q_len],single_targets=c(p_singles[i],p_singles[j]))
            p_doubles<-c(p_doubles,roots$root)
            }
          }
          
        }
        
        p_doubles_mat<-rbind(p_doubles_mat,p_doubles)
        
        triples<-function(x,ratio,single_targets,double_targets){
          x1<-prod(single_targets)
          x2<-single_targets[1]*double_targets[3]
          x3<-single_targets[2]*double_targets[2]
          x4<-single_targets[3]*double_targets[1]

          ratio-(x1*10^8+x2*10^8+x3*10^8+x4*10^8+x*10^8-sum(x1*c(x2*10^8,x3*10^8,x4*10^8,x*10^8))-sum(x2*c(x3*10^8,x4*10^8,x*10^8))-sum(x3*c(x4*10^8,x*10^8))-x4*x*10^8+sum(x1*c(x2*x3*10^8,x2*x4*10^8,x2*x*10^8,x3*x4*10^8,x3*x*10^8,x4*x*10^8))+
                   sum(x2*c(x3*x4*10^8,x3*x*10^8,x4*x*10^8))+x3*x4*x*10^8-sum(x1*c(x2*x3*x4*10^8,x2*x3*x*10^8,x2*x4*x*10^8,x3*x4*x*10^8))-x2*x3*x4*x*10^8+x1*x2*x3*x4*x*10^8)/(10^8*prod(1-single_targets)*prod(1-double_targets)*(1-x))
                   
        }
        
        p_triples<-c()
        k=k+1
        roots<-uniroot(triples, interval=c(0,1),ratio=Qs[k]/Qs[Q_len],single_targets=c(p_singles[1],p_singles[2],p_singles[3]),double_targets=c(p_doubles[1],p_doubles[2],p_doubles[4]))
        p_triples<-c(p_triples,roots$root)
        
        k=k+1
        # roots<-uniroot(triples, interval=c(0,1),ratio=Qs[k]/Qs[Q_len],single_targets=c(p_singles[1],p_singles[2],p_singles[4]),double_targets=c(p_doubles[1],p_doubles[3],p_doubles[5]))
        p_triples<-c(p_triples,0)
        
        k=k+1
        # roots<-uniroot(triples, interval=c(0,1),ratio=Qs[k]/Qs[Q_len],single_targets=c(p_singles[1],p_singles[3],p_singles[4]),double_targets=c(p_doubles[2],p_doubles[3],p_doubles[6]))
        p_triples<-c(p_triples,0)
        
        k=k+1
        roots<-uniroot(triples, interval=c(0,1),ratio=Qs[k]/Qs[Q_len],single_targets=c(p_singles[2],p_singles[3],p_singles[4]),double_targets=c(p_doubles[4],p_doubles[5],p_doubles[6]))
        p_triples<-c(p_triples,roots$root)
        
        p_triples_mat<-rbind(p_triples_mat,p_triples)
        
        quadras<-function(x,empty_prop,p_singles,p_doubles,p_triples){
          empty_prop-prod(1-p_singles)*prod(1-p_doubles)*prod(1-p_triples)*(1-x)
        }
        
        p_quadras<-uniroot(quadras, interval=c(0,1),empty_prop=Qs[Q_len]/p_tot,p_singles=p_singles,p_doubles=p_doubles,p_triples=p_triples)$root
        p_quadras_vec<-c(p_quadras_vec,p_quadras)
        
        lambda1<-(-log(1-p_singles[1]))
        lambda2<-(-log(1-p_singles[2]))
        lambda3<-(-log(1-p_singles[3]))
        lambda4<-(-log(1-p_singles[4]))
        
        lambda12<-(-log(1-p_doubles[1]))
        lambda13<-(-log(1-p_doubles[2]))
        lambda14<-(-log(1-p_doubles[3]))
        lambda23<-(-log(1-p_doubles[4]))
        lambda24<-(-log(1-p_doubles[5]))
        lambda34<-(-log(1-p_doubles[6]))
        
        lambda123<-(-log(1-p_triples[1]))
        lambda124<-(-log(1-p_triples[2]))
        lambda134<-(-log(1-p_triples[3]))
        lambda234<-(-log(1-p_triples[4]))
        
        lambda1234<-(-log(1-p_quadras))
        
        
        dsi_4dim[l]<-(lambda123+lambda234+(lambda12+lambda34)/2+lambda23+lambda2)/(lambda123+lambda234+(lambda12+lambda34)/2+lambda23+lambda2+lambda1234)
      }, error = function(e) {
        dsi_4dim[l]<-NA
      })
    }
    dsi_4dim_tmp<-cbind(intactness=1-dsi_4dim,molecules=m,true_intactness=1-dsi)
    dsi_new<-rbind(dsi_new,dsi_4dim_tmp)
  }
  
}




