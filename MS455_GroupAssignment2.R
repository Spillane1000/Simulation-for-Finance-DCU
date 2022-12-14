# (b) (ii)
## Milstein scheme for Heston stochastic volatility model absolute value for volatity process taken to
## avoid negative volatilty outputs 

HestonMil = function(M,N){
  
  # given parameters
  T = 5 
  h = T/N
  S0 = 100
  V0 = 0.16
  k = 1
  theta = 0.09
  rho = 0.3
  r = 0.03
  xi = 0.2
  
  # initialize matrices for results 
  norm_sample_B = matrix(0L,nrow = M,ncol = N)  
  norm_sample_W = matrix(0L,nrow = M,ncol = N)
  result_Vn = matrix(0L,nrow = M,ncol = N)
  result_Sn = matrix(0L,nrow = M,ncol = N)
  
  
  for (j in 1:M) {
    
    
    
    
    norm_sample_W[j,] = rnorm(N,0,1) # Normal sample for Brownian Motion W
    norm_sample_B[j,] = rnorm(N,0,1) # Normal sample for Brownian Motion B
    
    
    
    
    sample_path_Vn = numeric(N)
    sample_path_Sn = numeric(N)
    
    sample_path_Vn[1] = V0
    sample_path_Sn[1] = S0
    
    # milstein scheme for volatility process, absolute value taken to remove chance of negative volatility
    for (i in 2:N) {
      sample_path_Vn[i] = abs(sample_path_Vn[i-1] + k*(theta - sample_path_Vn[i-1])*h +
                                xi*rho*(sqrt(sample_path_Vn[i-1]))*(sqrt(h)*norm_sample_W[j,i-1]) +
                                xi*(sqrt(1-rho^2))*sqrt(sample_path_Vn[i-1])*(sqrt(h)*norm_sample_B[j,i-1]) +
                                (1/4)*((xi^2)*(rho^2))*(h*(norm_sample_W[j,i-1]^2) - h) +
                                (1/2)*((xi^2)*rho*(sqrt(1-rho^2)))*((norm_sample_W[j,i-1])*norm_sample_B[j,i-1]) - sample_path_Vn[i-1] +
                                (1/4)*((xi^2)*(1-rho^2))*((h*norm_sample_B[j,i-1]^2) - h))
      
      # milstein scheme for stock price process
      sample_path_Sn[i] = sample_path_Sn[i-1] + r*sample_path_Sn[i-1]*h +
        sqrt(sample_path_Vn[i-1])*sample_path_Sn[i-1]*(sqrt(h)*norm_sample_W[j,i-1]) +
        (1/2)*((sample_path_Vn[i-1]*sample_path_Sn[i-1]) + (sample_path_Sn[i-1]/2)*xi*rho)*(h*(norm_sample_W[j,i-1]^2)-h) + 
        (1/4)*(sample_path_Sn[i-1]*xi*(sqrt(1-rho^2)))*((norm_sample_W[j,i-1]*norm_sample_B[j,i-1]) - sample_path_Vn[i-1])
      
    }
    result_Vn[j,] = sample_path_Vn
    result_Sn[j,] = sample_path_Sn
  }
  return(result_Sn)
}


#####################################################################################################################################################################
#####################################################################################################################################################################

# (b) (i)
## Euler-Maruyama scheme for Heston stochastic volatility model absolute value for volatity process taken to
## avoid negative volatilty outputs 

HestonEM = function(M,N){
  
  # given parameters
  T = 5 
  h = T/N
  S0 = 100
  V0 = 0.16
  k = 1
  theta = 0.09
  rho = 0.3
  r = 0.03
  xi = 0.2
  
  # initialize matrices for results 
  norm_sample_B = matrix(0L,nrow = M,ncol = N)  
  norm_sample_W = matrix(0L,nrow = M,ncol = N)
  result_Vn = matrix(0L,nrow = M,ncol = N)
  result_Sn = matrix(0L,nrow = M,ncol = N)
  
  
  for (j in 1:M) {
    
    
    
    
    norm_sample_W[j,] = rnorm(N,0,1) # Normal sample for Brownian Motion W
    norm_sample_B[j,] = rnorm(N,0,1) # Normal sample for Brownian Motion B
    
    
    
    
    sample_path_Vn = numeric(N)
    sample_path_Sn = numeric(N)
    
    sample_path_Vn[1] = V0
    sample_path_Sn[1] = S0
    
    # E.M scheme for volatility process, absolute value used to avoid negative volatilities 
    for (i in 2:N) {
      sample_path_Vn[i] = abs(sample_path_Vn[i-1] + k*(theta - sample_path_Vn[i-1])*h +
                                xi*rho*(sqrt(sample_path_Vn[i-1]))*(sqrt(h)*norm_sample_W[j,i-1]) +
                                xi*(sqrt(1-rho^2))*sqrt(sample_path_Vn[i-1])*(sqrt(h)*norm_sample_B[j,i-1]))
      
      
      # E.M scheme for stock price process
      sample_path_Sn[i] = sample_path_Sn[i-1] + r*sample_path_Sn[i-1]*h +
        sqrt(sample_path_Vn[i-1])*sample_path_Sn[i-1]*(sqrt(h)*norm_sample_W[j,i-1]) 
      
    }
    result_Vn[j,] = sample_path_Vn
    result_Sn[j,] = sample_path_Sn
  }
  return(result_Sn)
}



########################################################################################################################################
########################################################################################################################################

# (c) part (i)
## simulations to estimate stock price at time e

set.seed(4)  # to ensure results are reproducable 


M = c(10,100,250,1000)  # number of paths to simulate

N = c(50,100,1000)   # number of steps for discretisation

estimate = list()   # to store samples of stock price at time e (Euler) 
#each element contains matrix containing samples of stock price with 
#columns representing different values of N


attach(mtcars)          # for plotting empirical distributions not sure which ones he wants plotted
par(mfrow = c(2,2))     


MCE = matrix(0L,ncol = length(M),nrow = length(N))
sample = matrix(0L,nrow = M[1],ncol = length(N)) 


for (i in 1:length(M)) {
  
  
  
  for (k in 1:length(N)) {
    
    
    
    euler = HestonEM(M[i],N[k]) # creates sample of stock price paths change to euler = HestonMil(M[i],N[k]) for Milstein
    
    h = 5/N[k]   # step size
    
    upp = ceiling(exp(1)/h) 
    low = floor(exp(1)/h)
    
    S_e_upp = euler[,upp] 
    S_e_low = euler[,low]
    
    
    
    sample[,k] = S_e_low + ((exp(1)-low*h)/h)*(S_e_upp - S_e_low)  # linear interpolation for stock price at e 
    
    MCE[k,i] = mean(S_e_low + ((exp(1)-low*h)/h)*(S_e_upp - S_e_low))
    
    
    
  }
  estimate[[i]] = sample
  if(i < length(M))sample = matrix(0L,nrow = M[i+1],ncol = length(N))  # if statement just to avoid a warning
}
colnames(MCE) = c("M = 10","M = 100","M = 250","M = 1000")
rownames(MCE) = c("N = 50","N = 100","N = 1000")
View(MCE)


#########################################################################################################################################
#########################################################################################################################################
#(d) (i)
## pricing European call option
set.seed(5)

strike_price = 115
r = 0.03
T = 5


M = c(10,100,250,1000)  # number of paths to simulate
N = c(50,100,1000)   # number of steps for discretisation

option_price_table_european = matrix(0L,nrow = length(N),ncol = length(M))

for (k in 1:length(M)) {
  
  for (i in 1:length(N)) {
    
    
    
    discretisation = HestonEM(M[k],N[i])   # change to HestonMil for milstein scheme
    
    
    S_T = discretisation[,ncol(discretisation)]  #stock price from sample path at time T = 5
    payoff = S_T - strike_price    # payoff for each sample path 
    payoff[payoff < 0] = 0    #removes any negative payoffs ie. where K > S_T (payoff is 0 in these cases)
    MCE_payoff = mean(payoff)   # monte carlo estimate of payoff
    option_price_table_european[i,k] = exp(-1*r*T)*MCE_payoff #discounted option price at time t = 0 
    
  }
}    
colnames(option_price_table_european) = c(paste("M =",M[1]),paste("M = ",M[2]),paste("M = ",M[3]),paste("M = ",M[4]))
rownames(option_price_table_european) = c(paste("N = ",N[1]),paste("N = ",N[2]),paste("N = ",N[3]))
View(option_price_table_european)

#################################################################################################################################
#################################################################################################################################
# (e) (i)
## pricing asian type 1 call option, code is the same as the european except average price of stock over the time
# period is used instead of price at end of period for payoff calculation 



set.seed(8)

strike_price = 115
r = 0.03
T = 5


M = c(10,100,250,1000)  # number of paths to simulate
N = c(50,100,1000)   # number of steps for discretisation

option_price_table_asian = matrix(0L,nrow = length(N),ncol = length(M))

for (k in 1:length(M)) {
  
  for (i in 1:length(N)) {
    
    
    
    discretisation = HestonEM(M[k],N[i])   # change to HestonMil for milstein scheme
    
    
    
    A_T = rowMeans(discretisation)   # arithmetic mean how to estimate the integral?
    payoff = A_T - strike_price
    payoff[payoff < 0] = 0    #removes any negative payoffs ie. where K > S_T (payoff is 0 in these cases)
    MCE_payoff = mean(payoff)   # monte carlo estimate of payoff
    option_price_table_asian[i,k] = exp(-1*r*T)*MCE_payoff #discounted option price at time t = 0 
    
  }
}    
colnames(option_price_table_asian) = c(paste("M =",M[1]),paste("M = ",M[2]),paste("M = ",M[3]),paste("M = ",M[4]))
rownames(option_price_table_asian) = c(paste("N = ",N[1]),paste("N = ",N[2]),paste("N = ",N[3]))
View(option_price_table_asian)