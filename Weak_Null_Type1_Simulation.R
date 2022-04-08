library(gurobi)

simu_weak_type1<-function(I, p_0, p_1, Uni_low, Uni_upp, alpha, gap, S){
  run_time_vec<-rep(0, S)
  total_time_vec<-rep(0, S)
  WA_vec<-rep(0, S)
  pvalue_pri<-rep(0, S)
  pvalue_opt_vec<-rep(0, S)
  SW_00_vec<-rep(0, S)
  SW_01_vec<-rep(0, S)
  SW_10_vec<-rep(0, S)
  SW_11_vec<-rep(0, S)
  total_sample_vec<-rep(0, S)
  total_treated_vec<-rep(0, S)
  rej_ind_vec<-rep(0, S)
  opt_ind_vec<-rep(0, S)
  result_matrix<-matrix(0, nrow = S, ncol = 13)
  for (d in 1:S){
    start_time <- Sys.time()
    
    n<-rep(0, I)
    m<-rep(0, I)
    
    for (i in 1:I){
      m[i]=sample(c(Uni_low:Uni_upp), size = 1)
      n[i]=sample(c(Uni_low:Uni_upp), size = 1)+m[i]
    }
    
    total_sample_vec[d]=sum(n)
    total_treated_vec[d]=sum(m)
    
    #Input
    Q<-rep(0, sum(n))
    Q_sub<-rep(0, sum(n))
    Z<-rep(0, sum(n))
    Y<-rep(0, sum(n))
    c=qchisq(1-alpha, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE) #=qnorm(0.975)^2
    
    N<-rep(0, I)
    for (i in 1:I){
      N[i]=sum(n[1:i])
    }
    
    Q[1:N[1]]=1 # The index of the stratum 
    Q_sub[1:N[1]]=c(1:N[1]) # The index of the subject in the stratum
    Z[1:m[1]]=1
    Z[(m[1]+1):N[1]]=0
    Z[1:N[1]]=sample(Z[1:N[1]], size = N[1])
    
    if (I>1){
      for (i in 2:I){
        Q[(N[i-1]+1):N[i]]=i
        Q_sub[(N[i-1]+1):N[i]]=c(1:n[i])
        Z[(N[i-1]+1):(N[i-1]+m[i])]=1
        Z[(N[i-1]+m[i]+1):(N[i])]=0
        Z[(N[i-1]+1):N[i]]=sample(Z[(N[i-1]+1):N[i]], size = n[i])
      }
    }
    
    for (s in 1:sum(n)){
      if (Z[s]==1){
        Y[s]=rbinom(1, 1, prob = p_1)
      } else {
        Y[s]=rbinom(1, 1, prob = p_0)
      }
    }
    
    T_Neyman=0
    V_null_1=0
    V_null_2=0
    for (i in 1:I){
      T_Neyman=T_Neyman+(n[i]/sum(n))*((sum(Z[Q==i]*Y[Q==i]))/m[i]-(sum((1-Z[Q==i])*Y[Q==i]))/(n[i]-m[i]))
      V_T=sum(Z[Q==i]*(Y[Q==i]-(1/m[i])*sum(Z[Q==i]*Y[Q==i]))^2)
      V_C=sum((1-Z[Q==i])*(Y[Q==i]-(1/(n[i]-m[i]))*sum((1-Z[Q==i])*Y[Q==i]))^2)      
      V_null_1=V_null_1+(n[i]/sum(n))^2*((1/(m[i]*(m[i]-1)))*V_T+(1/((n[i]-m[i])*(n[i]-m[i]-1)))*V_C)
      V_null_2=V_null_2+(n[i]/sum(n))^2*((1/m[i])*var(Y[Q==i & Z==1])+(1/(n[i]-m[i]))*var(Y[Q==i & Z==0]))
    }    
    T_Neyman_normalized<-((T_Neyman)^2)/V_null_1
    pvalue<-pchisq(T_Neyman_normalized, df=1, lower.tail = FALSE)
    pvalue_pri[d]<-pvalue
    if (pvalue<alpha){
      rej_ind_vec[d]=1
    }
    
    K_1<-rep(0, 4*I)
    K_2<-rep(0, 4*I)
    for (i in 1:I){
      K_1[4*i-3]=i
      K_1[4*i-2]=i
      K_1[4*i-1]=i
      K_1[4*i]=i
      K_2[4*i-3]=1
      K_2[4*i-2]=2
      K_2[4*i-1]=3
      K_2[4*i]=4
    }
    
    q<-rep(0, 4*I)
    for (i in 1:I){
      q[4*i-3]=-1/sum(n)
      q[4*i-2]=1/sum(n)
      q[4*i-1]=-1/sum(n)
      q[4*i]=1/sum(n)
    }
    
    
    Q_1<-matrix(0, nrow=4*I, ncol=4*I)
    for (s in 1:(4*I)){
      for (t in 1:(4*I)){
        if (K_1[s]==K_1[t]){
          if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
            Q_1[s,t]=(n[K_1[s]]^2)/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])^2)+c*(n[K_1[s]]^2)/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])^2*(n[K_1[s]]-m[K_1[s]]-1))
          } else if ((K_2[s]==3 & K_2[t]==3)|(K_2[s]==3 & K_2[t]==4)|(K_2[s]==4 & K_2[t]==3)|(K_2[s]==4 & K_2[t]==4)){
            Q_1[s,t]=(n[K_1[s]]^2)/((sum(n))^2*(m[K_1[s]])^2)+c*(n[K_1[s]]^2)/((sum(n))^2*(m[K_1[s]])^2*(m[K_1[s]]-1))
          } else {
            Q_1[s,t]=-(n[K_1[s]]^2)/((sum(n))^2*m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))
          }
        } else {
          if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
            Q_1[s,t]=(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])*(n[K_1[t]]-m[K_1[t]]))
          } else if ((K_2[s]==1 & K_2[t]==3)|(K_2[s]==1 & K_2[t]==4)|(K_2[s]==2 & K_2[t]==3)|(K_2[s]==2 & K_2[t]==4)){
            Q_1[s,t]=-(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])*m[K_1[t]])
          } else if ((K_2[s]==3 & K_2[t]==1)|(K_2[s]==3 & K_2[t]==2)|(K_2[s]==4 & K_2[t]==1)|(K_2[s]==4 & K_2[t]==2)){
            Q_1[s,t]=-(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*m[K_1[s]]*(n[K_1[t]]-m[K_1[t]]))
          } else {
            Q_1[s,t]=(n[K_1[s]]*n[K_1[t]])/((sum(n))^2*m[K_1[s]]*m[K_1[t]])
          }
        }
      }
    }
    q_1<-rep(0, 4*I)
    for (s in 1:(4*I)){
      if (K_2[s]==1 | K_2[s]==2){
        q_1[s]=-c*(n[K_1[s]]^2)/((sum(n))^2*(n[K_1[s]]-m[K_1[s]])*(n[K_1[s]]-m[K_1[s]]-1))
      } else {
        q_1[s]=-c*(n[K_1[s]]^2)/((sum(n))^2*m[K_1[s]]*(m[K_1[s]]-1))
      }
    }
    
    ##################################################
    # Initiate a model
    model = list()
    model$modelsense = 'max'
    model$vtype = rep('I', 4*I)
    L<-rep(0, 4*I)
    U<-rep(0, 4*I)
    for (i in 1:I){
      U[4*i-3]=sum((1-Z[Q==i])*(1-Y[Q==i]))
      U[4*i-2]=sum((1-Z[Q==i])*Y[Q==i])
      U[4*i-1]=sum(Z[Q==i]*(1-Y[Q==i]))
      U[4*i]=sum(Z[Q==i]*Y[Q==i])
    }
    model$lb<-L
    model$ub<-U
    
    ###################################################
    #Objective function
    model$obj = q  #linear term
    model$objcon = sum(1-Y)/sum(n)     #constant term
    
    ###########################################################
    # Linear constraint
    model$A=matrix(0, nrow=1, ncol=4*I)
    
    ###########################################################
    # Quadratic constraint
    if (pvalue<alpha){
      model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')
    } else {
      model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '>')
    }
    
    #if (verbose) flag = 1
    #else flag = 0
    
    #https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
    res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0, TimeLimit = 100))
    
    opt_ind_vec[d]=as.numeric(res$status=="OPTIMAL")
    if (opt_ind_vec[d]==1){
      
      WA_vec[d]<-res$objval
      
      end_time <- Sys.time()
      
      run_time_vec[d] <- res$runtime
      total_time_vec[d] <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      #https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_results.html
      
      Up_opt=res$x
      T_Neyman_opt=0
      V_null_opt=0
      for (i in 1:I){
        T_Neyman_opt=T_Neyman_opt+(n[i]/sum(n))*((Up_opt[4*i-1]+Up_opt[4*i])/m[i]-(Up_opt[4*i-3]+Up_opt[4*i-2])/(n[i]-m[i]))
        V_null_opt_1=(Up_opt[4*i-1]+Up_opt[4*i])/(m[i]*(m[i]-1))
        V_null_opt_2=((Up_opt[4*i-1]+Up_opt[4*i])^2)/((m[i]^2)*(m[i]-1))
        V_null_opt_3=(Up_opt[4*i-3]+Up_opt[4*i-2])/((n[i]-m[i])*(n[i]-m[i]-1))
        V_null_opt_4=((Up_opt[4*i-3]+Up_opt[4*i-2])^2)/(((n[i]-m[i])^2)*(n[i]-m[i]-1))
        V_null_opt=V_null_opt+(n[i]/sum(n))^2*(V_null_opt_1-V_null_opt_2+V_null_opt_3-V_null_opt_4)
      }    
      T_Neyman_normalized_opt<-((T_Neyman_opt)^2)/V_null_opt
      pvalue_opt_vec[d]<-pchisq(T_Neyman_normalized_opt, df=1, lower.tail = FALSE) #Conservative pvalue
      
      ##############Sensitivity Weights########################
      SW=matrix(0, nrow = 2, ncol = 2)
      colnames(SW)<-c("False Positives", "False Negatives")
      row.names(SW)<-c("Treated", "Control")
      SW_00=0
      SW_01=0
      SW_10=0
      SW_11=0
      
      for (i in 1:I){
        SW_00=SW_00+res$x[4*i-3]
        SW_01=SW_01+U[4*i-2]-res$x[4*i-2]
        SW_10=SW_10+res$x[4*i-1]
        SW_11=SW_11+U[4*i]-res$x[4*i]
      }
      
      SW[2,2]=SW_00/(SW_00+SW_01+SW_10+SW_11)
      SW[2,1]=SW_01/(SW_00+SW_01+SW_10+SW_11)
      SW[1,2]=SW_10/(SW_00+SW_01+SW_10+SW_11)
      SW[1,1]=SW_11/(SW_00+SW_01+SW_10+SW_11)
      
      SW_00_vec[d]<-SW[2,2]
      SW_01_vec[d]<-SW[2,1]
      SW_10_vec[d]<-SW[1,2]
      SW_11_vec[d]<-SW[1,1]
    } else {
      WA_vec[d]=-1
    }
    
    if (d%%100==0){
      print(d)
    }
    
  }
  
  result_matrix[,1]=run_time_vec
  result_matrix[,2]=total_time_vec
  result_matrix[,3]=WA_vec
  result_matrix[,4]=pvalue_pri
  result_matrix[,5]=pvalue_opt_vec
  result_matrix[,6]=SW_00_vec
  result_matrix[,7]=SW_01_vec
  result_matrix[,8]=SW_10_vec
  result_matrix[,9]=SW_11_vec
  result_matrix[,10]=total_sample_vec
  result_matrix[,11]=total_treated_vec
  result_matrix[,12]=rej_ind_vec
  result_matrix[,13]=opt_ind_vec
  
  colnames(result_matrix)<-c("run time", "total time", "WA", "pvalue_pri", "pvalue_opt", "SW_00", "SW_01", "SW_10", "SW_11", "total_sample", "total_treated", "rej or not", "opt or not")
  
  return(result_matrix)
}

res_weak_type1_2000_3_4<-simu_weak_type1(40, 0.3, 0.4, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_3_6<-simu_weak_type1(40, 0.3, 0.6, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_3_8<-simu_weak_type1(40, 0.3, 0.8, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_6_7<-simu_weak_type1(40, 0.6, 0.7, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_6_8<-simu_weak_type1(40, 0.6, 0.8, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_6_9<-simu_weak_type1(40, 0.6, 0.9, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_9_2<-simu_weak_type1(40, 0.9, 0.2, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_9_4<-simu_weak_type1(40, 0.9, 0.4, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_2000_9_6<-simu_weak_type1(40, 0.9, 0.6, 10, 40, 0.05, 0.00001, 1000)

res_weak_type1_10000_3_4<-simu_weak_type1(200, 0.3, 0.4, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_3_6<-simu_weak_type1(200, 0.3, 0.6, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_3_8<-simu_weak_type1(200, 0.3, 0.8, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_6_7<-simu_weak_type1(200, 0.6, 0.7, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_6_8<-simu_weak_type1(200, 0.6, 0.8, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_6_9<-simu_weak_type1(200, 0.6, 0.9, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_9_2<-simu_weak_type1(200, 0.9, 0.2, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_9_4<-simu_weak_type1(200, 0.9, 0.4, 10, 40, 0.05, 0.00001, 1000)
res_weak_type1_10000_9_6<-simu_weak_type1(200, 0.9, 0.6, 10, 40, 0.05, 0.00001, 1000)

