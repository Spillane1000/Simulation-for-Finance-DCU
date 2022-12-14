###################### Question 1 ##########################


# Distribution function of one-sided normal
norm_one_side_dist <- function (Sigma, x){
  F_fist = 2*(pnorm(x/Sigma)) - 1
  return(F_fist)
}


# Density function of one-sided normal
norm_one_side_density <- function (Sigma, x){
  f = sqrt(1/(2*pi))*(2/Sigma)*exp(-(x^2)/(2*(Sigma^2)))  
  return(f)
}




################################################################################################################


## 1 (b)

group_number = 6
sigma = 1+(group_number %% 5)


x = seq(0,6,0.01)
one_sided = (sqrt(1/(2*pi))*(2/sigma)*(exp(-(x^2)/(2*(sigma^2)))))
plot(x,one_sided,type = "l",col = "blue",ylim = c(0,1),ylab = "f(x)",main = "Distribution functions Q1 (b)")
lambda = c(0.25,0.5,1)     
c = (2/sqrt(2*pi))*exp(1/2)
for (i in 1:length(lambda)) {
  lines(x,2*dexp(x,lambda[i]),col = "red")
  
}
lines(x,c*dexp(x,1/sigma),col = "black")
legend(3,0.8,legend = c("One-Sided Normal","Exponential Various Lambda Values","Exponential with optimal C"),
       fill = c("blue","red","black"))



#############################################################################################################


## 1 (e) (f)
group_number = 6
sigma = 1+(group_number %% 5)

sample_list = list()
chi_results_table_list = list()

N_vec <- c(10,50,100,1000,10000)
xaxis_limit <- 11

chi_results_table <- as.data.frame(matrix(rep(0, 4*length(N_vec)),nrow = length(N_vec), ncol = 4))
colnames(chi_results_table) = c("N", "Test Stat", "df", "Crit Value")
attach(mtcars)
par(mfrow = c(2,2))

for (k in 1:length(N_vec)){
  
  N = N_vec[k]
  
  one_sided_norm_sim <- function(N){
    f = numeric(N)
    for( i in 1:N ){
      x=rexp(1, rate = 1/sigma)
      unif = runif(1)
      while(unif > exp((-(x^2)/(2*(sigma^2))) - (1/2) + (x/sigma))){     # wheres formula coming from?
        x = rexp(1, rate = 1/sigma)
        unif = runif(1)
      }
      f[i] = x
    }
    return(f)
  }
  one_sided_norm_sample = one_sided_norm_sim(N)
  sample_list[[k]] = one_sided_norm_sample
  hist(one_sided_norm_sample, breaks=seq(0,max(one_sided_norm_sample)+1, by = 0.4), 
       xlim=c(0,xaxis_limit),freq=FALSE, 
       ylim = c(0,1), ylab = "Density", xlab = " ", 
       main = paste("Empirical + Actual Density Function N =", N))
  x_values=seq(0,xaxis_limit,0.001)
  lines(x_values, norm_one_side_density(sigma, x_values), col = "red")
  lines(x_values, ((2/sqrt(2*pi))*exp(1/2))*dexp(x_values, rate = 1/sigma), col = "blue")
  
  # Take a sample from the set of values generated
  chi_sample = sample(one_sided_norm_sample, N, replace = FALSE)
  
  
  # Loop for different partitions
  par_length = c(0.25,0.5,0.75)
  for (j in 1:length(par_length)) {
    
    
    
    # Define the partition to be used
    part_len = par_length[j]
    part_end = 6
    partition_vec = seq(0,part_end, by = part_len)
    df = length(partition_vec) - 1
    
    #  c(seq(0,part_end/2, by = part_len/2),
    #                 seq(part_end/2+(part_len/2),part_end, by = part_len))
    
    # Count number of elements of sample in each interval of the partition
    count_vec = numeric(length(partition_vec))
    
    for (i in 2:length(partition_vec)){
      count_vec[i-1] = length(chi_sample[chi_sample<=partition_vec[i] & chi_sample>partition_vec[i-1]])
    }
    
    count_vec[length(partition_vec)] = length(chi_sample[chi_sample > partition_vec[length(partition_vec)]])
    
    # Find the expected number of realisations in each of the intervals in the partition
    expect_vec = numeric(length(partition_vec))
    
    for (i in 2:length(partition_vec)){
      expect_vec[i-1] = 
        norm_one_side_dist(sigma,partition_vec[i])-norm_one_side_dist(sigma, partition_vec[i-1])
    }
    expect_vec[length(partition_vec)] = 1-norm_one_side_dist(sigma,partition_vec[length(partition_vec)])
    expect_vec = expect_vec*length(chi_sample)
    
    # Determine the Chi-Squared Statistic
    chi_vec = ((count_vec - expect_vec)^2)/expect_vec
    chisq_stat = sum(chi_vec)
    chisq_stat
    print(paste("Test statistic ==" , chisq_stat, "df ==", df , "crit value ==", 
                qchisq(p=.05, df=df, lower.tail = FALSE)))
    chi_results_table[k,] = c(N, chisq_stat, df, qchisq(p=.05, df=df, lower.tail = FALSE))
    
    chi_results_table_list[[j]] = chi_results_table
  }
  
}
View(chi_results_table_list[[1]])

################################### Question 2 ################################

## 2 (a)

group_number = 6
a = 15 - 2*((group_number)%%7)


N = c(10,25,100,1000)
I = function(x){(x^(3/2))*(1-x)^2}
I_sq = function(x){(x^(3))*(1-x)^4}

mu = integrate(I,lower = 0 , upper = 1)$value
sigma = sqrt(integrate(I_sq,lower = 0 , upper = 1)$value - mu^2)


intv_emp_up =  numeric(length(N))
intv_true_up = numeric(length(N))
intv_emp_lo = numeric(length(N))
intv_true_lo = numeric(length(N))
for (i in 1:length(N)) {
  
  
  sample = runif(N[i])
  MCE = mean((sample^(3/2))*(1-sample)^2)   #((2*sample^(5/2))*(35*sample^2 - 90*sample + 63))/315]
  est_var = sd((sample^(3/2))*(1-sample)^2)
  
  intv_true_up[i] = mu +1.96*sigma/N[i]
  intv_true_lo[i] = mu -1.96*sigma/N[i]
  intv_emp_up[i] = MCE +1.96*est_var/N[i]
  intv_emp_lo[i] = MCE - 1.96*est_var/N[i]
}
plot(N,intv_emp_up, type = "o",col = "blue",ylim = c(0,0.1),ylab = "Integral value estimate")
lines(N,intv_emp_lo,type = "o" , col = "blue")
lines(N,intv_true_lo,type = "o" , col = "red")
lines(N,intv_true_up,type = "o" , col = "red")
legend(400,0.02,legend = c("Theoretical Confidence Interval","Empirical Confidence Interval"),fill = c("red","blue"))



####################################################################################################################################



## 2 (b)



N = c(10,25,100,1000)
mu_beta = 1/3
sigma_beta = sqrt(4/99)

intv_emp_up =  numeric(length(N))
intv_true_up = numeric(length(N))
intv_emp_lo = numeric(length(N))
intv_true_lo = numeric(length(N))


for (i in 1:length(N)) {
  sample_beta = rbeta(N[i],3/2,3)
  MCE_beta = mean(sample_beta)
  est_var = sd(sample_beta)
  
  intv_true_up[i] = mu_beta +1.96*sigma_beta/N[i]
  intv_true_lo[i] = mu_beta -1.96*sigma_beta/N[i]
  intv_emp_up[i] = MCE_beta +1.96*est_var/N[i]
  intv_emp_lo[i] = MCE_beta - 1.96*est_var/N[i]
}
plot(N,intv_emp_up, type = "o",col = "blue",ylim = c(0,0.6),ylab = "Integral value estimate")
lines(N,intv_emp_lo,type = "o" , col = "blue")
lines(N,intv_true_lo,type = "o" , col = "red")
lines(N,intv_true_up,type = "o" , col = "red")
legend(400,0.2,legend = c("Theoretical Confidence Interval","Empirical Confidence Interval"),fill = c("red","blue"))
