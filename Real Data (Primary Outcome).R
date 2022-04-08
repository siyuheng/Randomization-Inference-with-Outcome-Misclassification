#install.packages('sas7bdat')

library(sas7bdat)
library(gurobi)

known<-read.sas7bdat("/Users/heng6/Desktop/Measurement Error/PCPT/data/known.sas7bdat")
site<-read.sas7bdat("/Users/heng6/Desktop/Measurement Error/PCPT/data/study_site.sas7bdat")

table(site$ino)
unique(site$ino)
table(site$TRNO)

strata_no<-rep(0, length(unique(site$ino)))
N<-rep(0, length(unique(site$ino)))
N_t<-rep(0, length(unique(site$ino)))
N_c<-rep(0, length(unique(site$ino)))
N_same<-rep(0, length(unique(site$ino)))
Y_ava<-rep(0, length(unique(site$ino)))
for (i in 1:length(strata_no)){
  N[i]<-sum(site$ino==i)
  N_t[i]<-sum(site$ino==i & site$TRNO==1)
  N_c[i]<-sum(site$ino==i & site$TRNO==2)
}
for (i in 1:length(strata_no)){
  for (j in 1:length(strata_no)){
    if (N[i]==N[j] & N_t[i]==N_t[j] & N_c[i]==N_c[j]){
      N_same[i]=N_same[i]+1
    }
  }
}
strata_info<-data.frame(strata_no, N, N_t, N_c, N_same)
table(N_same)
table(known$pca)

#######################################################

data <- merge(x = site, y = known, by = "ALTPATID")
data_no_missing <- data[which(data$known=="Y"), ]

strata_no_ava<-rep(0, length(unique(data_no_missing$ino)))
N_ava<-rep(0, length(unique(data_no_missing$ino)))
N_t_ava<-rep(0, length(unique(data_no_missing$ino)))
N_c_ava<-rep(0, length(unique(data_no_missing$ino)))
N_same_ava<-rep(0, length(unique(data_no_missing$ino)))
Y_ava<-rep(0, length(unique(data_no_missing$ino)))
for (i in 1:length(strata_no_ava)){
  strata_no_ava[i]=unique(data_no_missing$ino)[i]
  N_ava[i]<-sum(data_no_missing$ino==unique(data_no_missing$ino)[i])
  N_t_ava[i]<-sum(data_no_missing$ino==unique(data_no_missing$ino)[i] & data_no_missing$TRNO==1)
  N_c_ava[i]<-sum(data_no_missing$ino==unique(data_no_missing$ino)[i] & data_no_missing$TRNO==2)
}
for (i in 1:length(strata_no_ava)){
  for (j in 1:length(strata_no_ava)){
    if (N_ava[i]==N_ava[j] & N_t_ava[i]==N_t_ava[j] & N_c_ava[i]==N_c_ava[j]){
      N_same_ava[i]=N_same_ava[i]+1
    }
  }
}
strata_info_ava<-data.frame(strata_no_ava, N_ava, N_t_ava, N_c_ava, N_same_ava)
table(N_same_ava)

sum(strata_info_ava$N_ava) #9060 study subjects
sum(strata_info_ava$N_t_ava) #4368 treated
sum(strata_info_ava$N_c_ava) #4692 control

sum(as.numeric(data_no_missing$pca=="Y" & data_no_missing$TRNO==1)) # 803
sum(as.numeric(data_no_missing$TRNO==1)) # 4368
sum(as.numeric(data_no_missing$pca=="Y" & data_no_missing$TRNO==2)) # 1147
sum(as.numeric(data_no_missing$TRNO==2)) # 4692

sum(as.numeric(data_no_missing$pca=="Y" & data_no_missing$TRNO==1))/sum(as.numeric(data_no_missing$TRNO==1)) # 0.1838
sum(as.numeric(data_no_missing$pca=="Y" & data_no_missing$TRNO==2))/sum(as.numeric(data_no_missing$TRNO==2)) # 0.2445

strata_info_ava$Y_t<-rep(0, nrow(strata_info_ava))
strata_info_ava$Y_c<-rep(0, nrow(strata_info_ava))
for (i in 1:nrow(strata_info_ava)){
  strata_info_ava$Y_t[i]=sum(as.numeric(data_no_missing$ino==strata_info_ava$strata_no_ava[i] & data_no_missing$TRNO == 1 & data_no_missing$pca=="Y"))
  strata_info_ava$Y_c[i]=sum(as.numeric(data_no_missing$ino==strata_info_ava$strata_no_ava[i] & data_no_missing$TRNO == 2 & data_no_missing$pca=="Y"))
}
strata_info_ava$N_same_ava<-rep(0, nrow(strata_info_ava))
for (i in 1:nrow(strata_info_ava)){
  for (j in 1:nrow(strata_info_ava)){
    if (strata_info_ava$N_ava[i]==strata_info_ava$N_ava[j] & strata_info_ava$N_t_ava[i]==strata_info_ava$N_t_ava[j] & strata_info_ava$Y_t[i]==strata_info_ava$Y_t[j] & strata_info_ava$Y_c[i]==strata_info_ava$Y_c[j]){
      strata_info_ava$N_same_ava[i]=strata_info_ava$N_same_ava[i]+1
    }
  }
}

data_no_missing$pca[which(data_no_missing$pca=="Y")]=1
data_no_missing$pca[which(data_no_missing$pca=="N")]=0
data_no_missing$TRNO[which(data_no_missing$TRNO==1)]=1
data_no_missing$TRNO[which(data_no_missing$TRNO==2)]=0
data_no_missing$index=0
data_no_missing$pca=as.numeric(data_no_missing$pca)

sum(as.numeric(data_no_missing$pca == 1 & data_no_missing$TRNO==1))/sum(as.numeric(data_no_missing$TRNO==1))
sum(as.numeric(data_no_missing$pca == 1 & data_no_missing$TRNO==0))/sum(as.numeric(data_no_missing$TRNO==0))

########################################

strata_info_ava_final<-strata_info_ava[-which(strata_info_ava$N_t_ava==0 | strata_info_ava$N_c_ava==0), ]
data_no_missing_no_zeroNtNc <- data_no_missing[which(data_no_missing$ino %in% strata_info_ava_final$strata_no_ava),]

strata_info_ava_final$Y_t<-rep(0, nrow(strata_info_ava_final))
strata_info_ava_final$Y_c<-rep(0, nrow(strata_info_ava_final))
for (i in 1:nrow(strata_info_ava_final)){
  strata_info_ava_final$Y_t[i]=sum(as.numeric(data_no_missing_no_zeroNtNc$ino==strata_info_ava_final$strata_no_ava[i] & data_no_missing_no_zeroNtNc$TRNO == 1 & data_no_missing_no_zeroNtNc$pca=="Y"))
  strata_info_ava_final$Y_c[i]=sum(as.numeric(data_no_missing_no_zeroNtNc$ino==strata_info_ava_final$strata_no_ava[i] & data_no_missing_no_zeroNtNc$TRNO == 2 & data_no_missing_no_zeroNtNc$pca=="Y"))
}
strata_info_ava_final$N_same_ava<-rep(0, nrow(strata_info_ava_final))
for (i in 1:nrow(strata_info_ava_final)){
  for (j in 1:nrow(strata_info_ava_final)){
    if (strata_info_ava_final$N_ava[i]==strata_info_ava_final$N_ava[j] & strata_info_ava_final$N_t_ava[i]==strata_info_ava_final$N_t_ava[j] & strata_info_ava_final$Y_t[i]==strata_info_ava_final$Y_t[j] & strata_info_ava_final$Y_c[i]==strata_info_ava_final$Y_c[j]){
      strata_info_ava_final$N_same_ava[i]=strata_info_ava_final$N_same_ava[i]+1
    }
  }
}

data_no_missing_no_zeroNtNc$pca[which(data_no_missing_no_zeroNtNc$pca=="Y")]=1
data_no_missing_no_zeroNtNc$pca[which(data_no_missing_no_zeroNtNc$pca=="N")]=0
data_no_missing_no_zeroNtNc$TRNO[which(data_no_missing_no_zeroNtNc$TRNO==1)]=1
data_no_missing_no_zeroNtNc$TRNO[which(data_no_missing_no_zeroNtNc$TRNO==2)]=0
data_no_missing_no_zeroNtNc$index=0
data_no_missing_no_zeroNtNc$pca=as.numeric(data_no_missing_no_zeroNtNc$pca)

sum(as.numeric(data_no_missing_no_zeroNtNc$pca == 1 & data_no_missing_no_zeroNtNc$TRNO==1))/sum(as.numeric(data_no_missing_no_zeroNtNc$TRNO==1))
sum(as.numeric(data_no_missing_no_zeroNtNc$pca == 1 & data_no_missing_no_zeroNtNc$TRNO==0))/sum(as.numeric(data_no_missing_no_zeroNtNc$TRNO==0))

########################################
count=0
for (i in 1:nrow(strata_info_ava_final)){
  count=count+1
  data_no_missing_no_zeroNtNc$index[which(data_no_missing_no_zeroNtNc$ino==unique(data_no_missing_no_zeroNtNc$ino)[i])]=count
}

table(data_no_missing_no_zeroNtNc$index)

gap=0.00001
I=max(data_no_missing_no_zeroNtNc$index)
n<-rep(0, I)
m<-rep(0, I)

for (i in 1:I){
  n[i]=sum(as.numeric(data_no_missing_no_zeroNtNc$index==i))
  m[i]=sum(as.numeric(data_no_missing_no_zeroNtNc$index==i & data_no_missing_no_zeroNtNc$TRNO==1))
}

#Input 

Q<-rep(0, sum(n))
Q_sub<-rep(0, sum(n))
Z<-data_no_missing_no_zeroNtNc$TRNO
Y<-data_no_missing_no_zeroNtNc$pca
sum(Z*Y)/sum(Z)
sum((1-Z)*Y)/sum(1-Z)
c=qchisq(0.95, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE) #=qnorm(0.975)^2

N<-rep(0, I)
for (i in 1:I){
  N[i]=sum(n[1:i])
}

Q[1:N[1]]=1
Q_sub[1:N[1]]=c(1:N[1])

if (I>1){
  for (i in 2:I){
    Q[(N[i-1]+1):N[i]]=i
    Q_sub[(N[i-1]+1):N[i]]=c(1:n[i])
  }
}


T_MH=sum(Z*Y)
E_null=0
V_null=0
for (i in 1:I){
  E_null=E_null+(m[i]/n[i])*sum(Y[Q==i])
  V_null=V_null+(m[i]*sum(Y[Q==i])*(n[i]-sum(Y[Q==i]))*(n[i]-m[i]))/(n[i]^2*((n[i])-1))
}
T_MH_normalized<-((T_MH-E_null)^2)/V_null
pvalue<-pchisq(T_MH_normalized, df=1, lower.tail = FALSE)

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
        Q_1[s,t]=(m[K_1[s]]/n[K_1[s]])^2+c*(m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/(n[K_1[s]]^2*(n[K_1[s]]-1))
      } else if ((K_2[s]==3 & K_2[t]==3)|(K_2[s]==3 & K_2[t]==4)|(K_2[s]==4 & K_2[t]==3)|(K_2[s]==4 & K_2[t]==4)){
        Q_1[s,t]=(1-m[K_1[s]]/n[K_1[s]])^2+c*(m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/(n[K_1[s]]^2*(n[K_1[s]]-1))
      } else {
        Q_1[s,t]=-(m[K_1[s]]/n[K_1[s]])*(1-m[K_1[s]]/n[K_1[s]])+c*(m[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/(n[K_1[s]]^2*(n[K_1[s]]-1))
      }
    } else {
      if ((K_2[s]==1 & K_2[t]==1)|(K_2[s]==1 & K_2[t]==2)|(K_2[s]==2 & K_2[t]==1)|(K_2[s]==2 & K_2[t]==2)){
        Q_1[s,t]=(m[K_1[s]]/n[K_1[s]])*(m[K_1[t]]/n[K_1[t]])
      } else if ((K_2[s]==1 & K_2[t]==3)|(K_2[s]==1 & K_2[t]==4)|(K_2[s]==2 & K_2[t]==3)|(K_2[s]==2 & K_2[t]==4)){
        Q_1[s,t]=-(m[K_1[s]]/n[K_1[s]])*(1-m[K_1[t]]/n[K_1[t]])
      } else if ((K_2[s]==3 & K_2[t]==1)|(K_2[s]==3 & K_2[t]==2)|(K_2[s]==4 & K_2[t]==1)|(K_2[s]==4 & K_2[t]==2)){
        Q_1[s,t]=-(1-m[K_1[s]]/n[K_1[s]])*(m[K_1[t]]/n[K_1[t]])
      } else {
        Q_1[s,t]=(1-m[K_1[s]]/n[K_1[s]])*(1-m[K_1[t]]/n[K_1[t]])
      }
    }
  }
}
q_1<-rep(0, 4*I)
for (s in 1:(4*I)){
  q_1[s]=-c*((m[K_1[s]]*n[K_1[s]]*(n[K_1[s]]-m[K_1[s]]))/((n[K_1[s]])^2*(n[K_1[s]]-1)))
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
if (pvalue<0.05){
  model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '<')
} else {
  model$quadcon[[1]]=list(Qc = Q_1, q=q_1, rhs = 0, sense = '>')
}

#if (verbose) flag = 1
#else flag = 0

#https://www.gurobi.com/documentation/9.1/refman/parameters.html#sec:Parameters
res = gurobi(model, params = list(NonConvex = 2, MIPGapAbs = gap, OutputFlag = 0, TimeLimit = 100))
res$objval
res$objval*sum(n)
sum(n)-res$objval*sum(n)

#https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_results.html

Up_opt=res$x
T_MH_opt=0
E_null_opt=0
V_null_opt=0
for (i in 1:I){
  T_MH_opt=T_MH_opt+Up_opt[4*i-1]+Up_opt[4*i]
  E_null_opt=E_null_opt+(m[i]/n[i])*(Up_opt[4*i-3]+Up_opt[4*i-2]+Up_opt[4*i-1]+Up_opt[4*i])
  V_null_opt=V_null_opt+(m[i]*(Up_opt[4*i-3]+Up_opt[4*i-2]+Up_opt[4*i-1]+Up_opt[4*i])*(n[i]-(Up_opt[4*i-3]+Up_opt[4*i-2]+Up_opt[4*i-1]+Up_opt[4*i]))*(n[i]-m[i]))/(n[i]^2*((n[i])-1))
}
T_MH_normalized_opt<-((T_MH_opt-E_null_opt)^2)/V_null_opt
pvalue_opt<-pchisq(T_MH_normalized_opt, df=1, lower.tail = FALSE)

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

