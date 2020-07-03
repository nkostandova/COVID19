max_time = 200
sd.dw = 0

x[1,] <- round(ICs)

S = x[, 1:Ncomp] 
St = x[, (Ncomp+1):(2*Ncomp)]
E = x[, (2*Ncomp+1):(3*Ncomp)] 
Et = x[, (3*Ncomp+1):(4*Ncomp)]
A = x[, (4*Ncomp+1):(5*Ncomp)]
At = x[, (5*Ncomp+1):(6*Ncomp)]
Im = x[, (6*Ncomp+1):(7*Ncomp)]
Imt = x[, (7*Ncomp+1):(8*Ncomp)]
Is = x[, (8*Ncomp+1):(9*Ncomp)]
Ist = x[, (9*Ncomp+1):(10*Ncomp)]
R =  x[, (10*Ncomp+1):(11*Ncomp)]

# incidence
incid_A = x[, (11*Ncomp+1):(12*Ncomp)]
incid_Im = x[, (12*Ncomp+1):(13*Ncomp)]
incid_Is = x[, (13*Ncomp+1):(14*Ncomp)]
incid_deaths = x[, (14*Ncomp+1):(15*Ncomp)]
incid_plasma = x[, (15*Ncomp+1):(16*Ncomp)]

P = matrix(NA, nrow = length(time), ncol = 1)
P[1, 1] = p0


out = matrix(NA, nrow = max_time, ncol = Ncomp * dim(nt)[1] + 5*Ncomp + 1 + 1)   # 5 = number of variables in incid matrix; 1 = plasma; 1 = time
temp = (cbind(time[1], S[1,], St[1,], E[1,], Et[1,], A[1,], At[1,], Im[1,], Imt[1,], Is[1,], Ist[1,], R[1,], incid_A[1,], incid_Im[1,], incid_Is[1,], incid_deaths[1,], incid_plasma[1,], P[1,]))
temp2 =  as.vector(temp)
temp2 = temp2[-c(2,36)]
out[1,] = temp2
names(out) = c("time", names(ICs), "plasma_units")



for (it in 1:(max_time-1)){
  
  
  
  WI = W %*% (A[it,]+At[it,]+Im[it, ]+Imt[it,]+Is[it,]+Ist[it,])
  births = rep(0, Ncomp)
  births[1] = mu
  deaths = rep(v, Ncomp)
  deaths_I = rep(v_I, Ncomp)
  deaths_It = rep(v_It, Ncomp)
  
  dw = rtruncnorm(Ncomp, a =0, mean = 1, sd = sd.dw)
  
  # transitions
  foi_prob <- 1 - exp( - beta* WI / N * dw * delta.t)
  exposed_prob <- 1 - exp( - sigma * delta.t)    
  inf_prob <- 1 - exp( - gamma * delta.t)   
  death_prob <- 1 - exp( - deaths * delta.t)
  death_prob_I <- 1 - exp( - deaths_I * delta.t)
  inf_prob_t <- 1 - exp(- gamma_t * delta.t * (1-prop_eti))
  inf_prob_et <- 1 - exp( - sigma * delta.t * prop_eti)
  death_prob_It <- 1 - exp( - deaths_It * delta.t)
  foi_prob_St <- 1 - exp( - ef_est * beta* WI / N * dw * delta.t)
  
  new_St = rep(NA, Ncomp)
  for(i in (1:Ncomp)){
    new_St[i] = ifelse(P[it] * p_use * p_s[i] + SIA_s[it, i]>S[it,i],  S[it, i] * prop_ef, P[it] * p_use * p_s[i] * prop_ef + SIA_s[it, i] * prop_ef)
  }
  p_St = ifelse(new_St == 0, 0, new_St / S[it, ])
  
  new_Et = rep(NA, Ncomp)
  for(i in (1:Ncomp)){
    new_Et[i] = ifelse(P[it] * p_use * p_e[i] + SIA_e[it, i]>E[it,i],  E[it, i] * prop_ef, P[it] * p_use * p_e[i] * prop_ef + SIA_e[it, i] * prop_ef)
  }      
  p_Et = ifelse(new_Et ==0, 0, new_Et / E[it, ])
  
  new_Iat = rep(NA, Ncomp)
  for(i in (1:Ncomp)){
    new_Iat[i] = ifelse(P[it] * p_use * p_ia[i] + SIA_ia[it, i]>A[it,i],  A[it, i] * prop_ef, P[it] * p_use * p_ia[i] * prop_ef + SIA_ia[it, i] * prop_ef)
  } 
  p_Iat = ifelse(new_Iat ==0, 0, new_Iat / A[it, ])
  
  new_Imt = rep(NA, Ncomp)
  for(i in (1:Ncomp)){
    new_Imt[i] = ifelse(P[it] * p_use * p_im[i] + SIA_im[it, i]>Im[it,i],  Im[it, i] * prop_ef, P[it] * p_use * p_im[i] * prop_ef + SIA_im[it, i] * prop_ef)
  }  
  p_Imt = ifelse(new_Imt ==0, 0, new_Imt / Im[it, ])
  
  new_Ist = rep(NA, Ncomp)
  for(i in (1:Ncomp)){
    new_Ist[i] = ifelse(P[it] * p_use * p_is[i] + SIA_is[it, i]>Is[it,i],  Is[it, i] * prop_ef, P[it] * p_use * p_is[i] * prop_ef + SIA_is[it, i] * prop_ef)
  }  
  p_Ist = ifelse(new_Ist ==0, 0, new_Ist / Is[it, ])
  
  
## set up a transition matrix A
nt = rbind(S[it, ], St[it,], E[it, ], Et[it, ], A[it, ], At[it,], Im[it,], Imt[it,], Is[it,], Ist[it,], R[it,])
Tm = matrix(NA, nrow = dim(nt)[1], ncol = Ncomp * dim(nt)[1])
z = rep(0, Ncomp)

# transitions out of S compartment
Tm[1,] = c(1-p_St-t(foi_prob), p_St, t(foi_prob), z, z, z, z, z, z, z, z)
# transitions out of St compartment
Tm[2,] = c(z, 1-t(foi_prob_St), t(foi_prob_St), z, z, z, z, z, z, z, z)
# transitions out of E compartment
Tm[3,] = c(z, z, 1-p_Et-exposed_prob, p_Et, exposed_prob*p_asym, z, exposed_prob*p_mild, z, exposed_prob*p_severe, z, z)
# transitions out of Et compartment
Tm[4,] = c(z, z, z, 1-inf_prob_et-inf_prob_t, inf_prob_et*p_asym, z, inf_prob_et*p_mild, z, inf_prob_et*p_severe, z, inf_prob_t)
# transitions out of Ia compartment
Tm[5,] = c(z, z, z, z, 1-p_Iat-inf_prob, p_Iat, z, z, z, z, inf_prob)
# transitions out of Iat compartment
Tm[6,] = c(z, z, z, z, z, 1-inf_prob_t, z, z, z, z, inf_prob_t)
# transitions out of Im compartment
Tm[7,] = c(z, z, z, z, z, z, 1-p_Imt-inf_prob, p_Imt, z, z, inf_prob)
# transitions out of Imt compartment
Tm[8,] = c(z, z, z, z, z, z, z, 1-inf_prob_t, z, z, inf_prob_t)
# transitions out of Is compartment
Tm[9,] = c(z, z, z, z, z, z, z, z, 1-p_Ist-inf_prob-death_prob_I, p_Ist, inf_prob)
# transitions out of Ist compartment
Tm[10,] = c(z, z, z, z, z, z, z, z, z, 1-inf_prob_t-death_prob_It, inf_prob_t)
# transitions out of R compartment
Tm[11,] = c(z, z, z, z, z, z, z, z, z, z, 1-z)


# do multinomial draws
nstates = dim(nt)[1]
Tm = data.frame(Tm)
TmN = matrix(NA, nrow=nstates, ncol=Ncomp*nstates)
ntp1 = matrix(NA, nrow = nstates, ncol = Ncomp)
incid = matrix(NA, nrow = 5, ncol=Ncomp)
ntp = matrix(NA, nrow = nrow(ntp1)+nrow(incid), ncol = Ncomp)
for (i in 1:Ncomp){
  # separate into matrices for each compartment
  Tm_i = Tm[, seq(i, ncol(Tm), Ncomp)]
  # draw multinomial draws
  TmN_i = matrix(NA, nrow = nstates, ncol = nstates)
  for (j in 1:nrow(TmN_i)){
    size = nt[j, i]
    TmN_i[j,] = rmultinom(1, size, prob = Tm_i[j,])

  }
  ntp1[,i] = t(colSums(TmN_i))
  incid[1,i] = sum(TmN_i[,5])-TmN_i[5,5] # incid_A
  incid[2,i] = sum(TmN_i[,7])-TmN_i[7,7] # incid_Im
  incid[3,i] = sum(TmN_i[,9])-TmN_i[9,9] # incid_Is
  incid[4,i] = nt[9,i]-TmN_i[9,10]-TmN_i[9,11]-TmN_i[9,9]+
    nt[10,i]-TmN_i[10,11]-TmN_i[10,10]  #incid_deaths
  incid[5,i] = sum(TmN_i[,2])-TmN_i[2,2]+sum(TmN_i[,4])-TmN_i[4,4]+sum(TmN_i[,6])-TmN_i[6,6]+sum(TmN_i[,8])-TmN_i[8,8]+sum(TmN_i[,10])-TmN_i[10,10] #incid_plasma
  #ntp[1:nrow(ntp1),i] = ntp1[,i]
  #ntp[(nrow(ntp1)+1):nrow(ntp),i] = incid[,i]
}
ntp1_t = as.vector(t(ntp1))
incid_t = as.vector(t(incid))

P[it+1,] = P[it,] - sum(incid[5,]) + rbinom(n = 1, size = sum(ntp1[11,]), prob = p_give*rate_p)
  
out[it+1,] = cbind(ntp1_t, incid_t, P[it+1,])

S[it+1,] = out[it+1, 1+1:Ncomp] 
St[it+1,] = out[it+1, 1+(Ncomp+1):(2*Ncomp)]
E[it+1,] = out[it+1, 1+(2*Ncomp+1):(3*Ncomp)] 
Et[it+1,] = out[it+1, 1+(3*Ncomp+1):(4*Ncomp)]
A[it+1,] = out[it+1, 1+(4*Ncomp+1):(5*Ncomp)]
At[it+1,] = out[it+1, 1+(5*Ncomp+1):(6*Ncomp)]
Im[it+1,] = out[it+1, 1+(6*Ncomp+1):(7*Ncomp)]
Imt[it+1,] = out[it+1, 1+(7*Ncomp+1):(8*Ncomp)]
Is[it+1,] = out[it+1, 1+(8*Ncomp+1):(9*Ncomp)]
Ist[it+1,] = out[it+1, 1+(9*Ncomp+1):(10*Ncomp)]
R[it+1,] =  out[it+1, 1+(10*Ncomp+1):(11*Ncomp)]

# incidence
incid_A[it+1,] = out[it+1, 1+(11*Ncomp+1):(12*Ncomp)]
incid_Im[it+1,] = out[it+1, 1+(12*Ncomp+1):(13*Ncomp)]
incid_Is[it+1,] = out[it+1, 1+(13*Ncomp+1):(14*Ncomp)]
incid_deaths[it+1,] = out[it+1, 1+(14*Ncomp+1):(15*Ncomp)]
incid_plasma[it+1,] = out[it+1, 1+(15*Ncomp+1):(16*Ncomp)]

P[t+1,] = out[it+1, 1+(16*Ncomp+1)]
}
