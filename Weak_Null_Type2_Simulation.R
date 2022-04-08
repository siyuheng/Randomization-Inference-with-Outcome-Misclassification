library(gurobi)
library(mgcv)

simu_weak_type2<-function(I, p_0, p_1, Uni_low, Uni_upp, alpha, gap, S){
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
    
    Q[1:N[1]]=1
    Q_sub[1:N[1]]=c(1:N[1])
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
    T_Neyman_normalized<-(T_Neyman^2)/V_null_1
    pvalue<-pchisq(T_Neyman_normalized, df=1, lower.tail = FALSE)
    pvalue_pri[d]<-pvalue
    if (pvalue<alpha){
      rej_ind_vec[d]=1
    }
    
    I_info<-matrix(0, nrow = I, ncol = 4) # 1-00, 2-01, 3-10, 4-11
    for (i in 1:I){
      I_info[i,1]=sum((1-Z[Q==i])*(1-Y[Q==i]))
      I_info[i,2]=sum((1-Z[Q==i])*Y[Q==i])
      I_info[i,3]=sum(Z[Q==i]*(1-Y[Q==i]))
      I_info[i,4]=sum(Z[Q==i]*Y[Q==i])
    }
    S_unique<-uniquecombs(I_info,ordered=FALSE) # Number of unique 2 by 2 tables
    S=nrow(S_unique)
    P_s<-rep(0, S)
    for (s in 1:S){
      count=0
      for (i in 1:I){
        if (sum(as.numeric(I_info[i,]==S_unique[s,]))==4){
          count=count+1
        }
      }
      P_s[s]=count
    }
    N_s<-rep(0, S)
    for (s in 1:S){
      N_s[s]=(S_unique[s,1]+1)*(S_unique[s,2]+1)*(S_unique[s,3]+1)*(S_unique[s,4]+1)
    }
    Delta<-matrix(0, nrow = sum(N_s), ncol = 4)
    for (s in 1:S){
      C_00<-seq(from = 0, to = S_unique[s,1], by = 1)
      C_01<-seq(from = 0, to = S_unique[s,2], by = 1)
      C_10<-seq(from = 0, to = S_unique[s,3], by = 1)
      C_11<-seq(from = 0, to = S_unique[s,4], by = 1)
      C_all<-as.matrix(expand.grid(C_00, C_01, C_10, C_11))
      if (s == 1){
        Delta[c(1:N_s[1]),]=C_all
      } else {
        end_1<-sum(N_s[1:(s-1)])+1
        end_2<-sum(N_s[1:s])
        Delta[c(end_1:end_2), ]=C_all
      }
    }
    N_s_index<-rep(0, sum(N_s))
    for (s in 1:S){
      if (s == 1){
        N_s_index[1:N_s[1]]=1
      } else {
        end_1<-sum(N_s[1:(s-1)])+1
        end_2<-sum(N_s[1:s])
        N_s_index[end_1:end_2]=s
      }
    }
    n_s<-rep(0, S)
    m_s<-rep(0, S)
    for (s in 1:S){
      n_s[s]=sum(S_unique[s, ])
      m_s[s]=S_unique[s,3]+S_unique[s,4]
    }
    
    ##################################################
    # Initiate a model
    model = list()
    model$modelsense = 'max'
    model$vtype = rep('I', sum(N_s))
    L<-rep(0, sum(N_s))
    model$lb<-L
    
    ###################################################
    #Objective function
    q<-rep(0, sum(N_s))
    for (i in 1:sum(N_s)){
      q[i]=(Delta[i,2]+Delta[i,4]-Delta[i,1]-Delta[i,3])/(sum(n))
    }
    model$obj = q  #linear term
    model$objcon = sum(1-Y)/sum(n)     #constant term
    
    ###########################################################
    # Linear constraint
    A_linear=matrix(0, nrow=length(N_s), ncol=sum(N_s))
    for (i in 1:length(N_s)){
      if (i == 1){
        A_linear[i, c(1:N_s[1])]=rep(1, N_s[1])
      } else {
        end_1<-sum(N_s[1:(i-1)])+1
        end_2<-sum(N_s[1:i])
        A_linear[i, c(end_1:end_2)]=rep(1, N_s[i])
      }
    }
    model$A=A_linear
    model$rhs=P_s
    model$sense='='
    
    ###########################################################
    # Quadratic constraint
    Q_1<-matrix(0, nrow = sum(N_s), ncol = sum(N_s))
    Q_1_i<-matrix(0, nrow = sum(N_s), ncol = 1)
    Q_1_j<-matrix(0, nrow = 1, ncol = sum(N_s))
    for (i in 1:sum(N_s)){
      Q_1_i[i,1]=(n_s[N_s_index[i]]*(Delta[i, 3]+Delta[i,4]))/(sum(n)*m_s[N_s_index[i]])-(n_s[N_s_index[i]]*(Delta[i, 1]+Delta[i,2]))/(sum(n)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]))
      Q_1_j[1,i]=(n_s[N_s_index[i]]*(Delta[i, 3]+Delta[i,4]))/(sum(n)*m_s[N_s_index[i]])-(n_s[N_s_index[i]]*(Delta[i, 1]+Delta[i,2]))/(sum(n)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]))
    }  
    Q_1<-Q_1_i%*%Q_1_j
    q_1<-rep(0, sum(N_s))
    for (i in 1:sum(N_s)){
      q_1_1=(Delta[i,3]+Delta[i,4])/(m_s[N_s_index[i]]*(m_s[N_s_index[i]]-1))
      q_1_2=((Delta[i,3]+Delta[i,4])^2)/(m_s[N_s_index[i]]^2*(m_s[N_s_index[i]]-1))
      q_1_3=(Delta[i,1]+Delta[i,2])/((n_s[N_s_index[i]]-m_s[N_s_index[i]])*(n_s[N_s_index[i]]-m_s[N_s_index[i]]-1))
      q_1_4=((Delta[i,1]+Delta[i,2])^2)/(((n_s[N_s_index[i]]-m_s[N_s_index[i]])^2)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]-1))
      q_1[i]=-c*((n_s[N_s_index[i]]/sum(n))^2)*(q_1_1-q_1_2+q_1_3-q_1_4)
    }
    
    if (pvalue<alpha){
      model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')
    } else {
      model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '>')
    }
    
    #if (verbose) flag = 1
    #else flag = 0
    #https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
    res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0, TimeLimit=100))
    
    opt_ind_vec[d]=as.numeric(res$status=="OPTIMAL")
    if (opt_ind_vec[d]==1){
      
      WA_vec[d]<-res$objval
      
      end_time <- Sys.time()
      
      run_time_vec[d] <- res$runtime
      total_time_vec[d] <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      #https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_results.html
      
      d_opt=res$x
      T_Neyman_opt=0
      V_null_opt=0
      # Delta 1-00, 2-01, 3-10, 4-11
      for (i in 1:sum(N_s)){
        T_Neyman_opt=T_Neyman_opt+d_opt[i]*((n_s[N_s_index[i]]/(sum(n)*m_s[N_s_index[i]]))*(Delta[i,3]+Delta[i,4])-(n_s[N_s_index[i]]/(sum(n)*(n_s[N_s_index[i]]-m_s[N_s_index[i]])))*(Delta[i,1]+Delta[i,2]))
        V_null_opt_1=(Delta[i,3]+Delta[i,4])/(m_s[N_s_index[i]]*(m_s[N_s_index[i]]-1))
        V_null_opt_2=((Delta[i,3]+Delta[i,4])^2)/((m_s[N_s_index[i]]^2)*(m_s[N_s_index[i]]-1))
        V_null_opt_3=(Delta[i,1]+Delta[i,2])/((n_s[N_s_index[i]]-m_s[N_s_index[i]])*(n_s[N_s_index[i]]-m_s[N_s_index[i]]-1))
        V_null_opt_4=((Delta[i,1]+Delta[i,2])^2)/(((n_s[N_s_index[i]]-m_s[N_s_index[i]])^2)*(n_s[N_s_index[i]]-m_s[N_s_index[i]]-1))
        V_null_opt=V_null_opt+d_opt[i]*((n_s[N_s_index[i]]/sum(n))^2)*(V_null_opt_1-V_null_opt_2+V_null_opt_3-V_null_opt_4)
      }
      T_Neyman_normalized_opt<-(T_Neyman_opt^2)/V_null_opt
      pvalue_opt_vec[d]<-pchisq(T_Neyman_normalized_opt, df=1, lower.tail = FALSE) #Conservative pvalue
      
      ##############Sensitivity Weights########################
      SW=matrix(0, nrow = 2, ncol = 2)
      colnames(SW)<-c("False Positives", "False Negatives")
      row.names(SW)<-c("Treated", "Control")
      SW_00=0
      SW_01=0
      SW_10=0
      SW_11=0
      
      # S_unique 1-00, 2-01, 3-10, 4-11
      for (i in 1:sum(N_s)){
        SW_00=SW_00+res$x[i]*Delta[i,1]
        SW_01=SW_01+res$x[i]*(S_unique[N_s_index[i], 2]-Delta[i,2])
        SW_10=SW_10+res$x[i]*Delta[i,3]
        SW_11=SW_11+res$x[i]*(S_unique[N_s_index[i], 4]-Delta[i,4])
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

res_weak_type2_2000_3_4<-simu_weak_type2(400, 0.3, 0.4, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_3_6<-simu_weak_type2(400, 0.3, 0.6, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_3_8<-simu_weak_type2(400, 0.3, 0.8, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_6_7<-simu_weak_type2(400, 0.6, 0.7, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_6_8<-simu_weak_type2(400, 0.6, 0.8, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_6_9<-simu_weak_type2(400, 0.6, 0.9, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_9_2<-simu_weak_type2(400, 0.9, 0.2, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_9_4<-simu_weak_type2(400, 0.9, 0.4, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_2000_9_6<-simu_weak_type2(400, 0.9, 0.6, 2, 3, 0.05, 0.00001, 1000)

res_weak_type2_10000_3_4<-simu_weak_type2(2000, 0.3, 0.4, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_3_6<-simu_weak_type2(2000, 0.3, 0.6, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_3_8<-simu_weak_type2(2000, 0.3, 0.8, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_6_7<-simu_weak_type2(2000, 0.6, 0.7, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_6_8<-simu_weak_type2(2000, 0.6, 0.8, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_6_9<-simu_weak_type2(2000, 0.6, 0.9, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_9_2<-simu_weak_type2(2000, 0.9, 0.2, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_9_4<-simu_weak_type2(2000, 0.9, 0.4, 2, 3, 0.05, 0.00001, 1000)
res_weak_type2_10000_9_6<-simu_weak_type2(2000, 0.9, 0.6, 2, 3, 0.05, 0.00001, 1000)

