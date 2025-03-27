library(ggplot2)
library(ggsci)
library(data.table)

# General case: calculate power of SV test ONLY under the assumption of independent variants, 
# with all variants in mask
calc_power_sv = function(n, h, p, b, asv)
{
  ncp_sv = 2*n*(b^2)*p*(1-p)/(1-h)
  sv_pow = 1 - prod(pchisq(qchisq(asv, df = 1, lower.tail = F), ncp = ncp_sv, df = 1))
  return(sv_pow)
}

# General case: calculate power of burden test ONLY under the assumption of independent variants, 
# with all variants in mask
calc_power_burden = function(n, h, p, b, w, aa)
{
  ncp_b = 2*n*((sum(w*b*p*(1-p)))^2)/((1-h)*(sum((w^2)*p*(1-p))))
  b_pow = pchisq(qchisq(aa, df = 1, lower.tail = F), ncp = ncp_b, df = 1, lower.tail = F)
  return(b_pow)
}

# General case: calculate power of SKAT ONLY under the assumption of independent variants,
# with all variants in mask
calc_power_skat = function(n, h, p, b, w, aa, nsim = 2000000)
{
  ncp_sv = 2*n*(b^2)*p*(1-p)/(1-h)
  scale_chi = 2*n*(1-h)*p*(1-p)*w
  dfs = rep(1, length(ncp_sv))
  rand_chisq = matrix(nrow = nsim, ncol = length(dfs))
  for(i in 1:ncol(rand_chisq))
  {
    rand_chisq[, i] = rchisq(nsim, df = dfs[i], ncp = ncp_sv[i])
  }
  rand_chisq_sum = apply(rand_chisq, 1, function(x) sum(scale_chi * x))
  rand_chisq_null = matrix(rchisq(nsim*length(dfs), df = 1, ncp = 0), nrow = nsim)
  rand_chisq_null_sum = sort(apply(rand_chisq_null, 1, function(x) sum(scale_chi * x)), decreasing = T)
  chisq_crit = rand_chisq_null_sum[nsim*aa]
  s_pow = length(which(rand_chisq_sum > chisq_crit))/nsim
  return(s_pow)
}


# General case: obtain power for all three tests as number of causal variants increases from 1 to v
power_compute = function(mafs, n, h, alpha_sv, alpha_a, iter, seed, eff_model, n_means = c(2, 3, 4, 5), n_sds = c(0.5, 0.5, 0.5, 0.5))
{
  mafs = mafs[mafs < 0.01 & mafs > 0] # keep only rare variants
  v = length(mafs) # number of variants
  sv_power_c = matrix(rep(0, v*iter), nrow = iter)
  bd_power_c = matrix(rep(0, v*iter), nrow = iter)
  sk_power_c = matrix(rep(0, v*iter), nrow = iter)
  causal = matrix(rep(0, v*iter), ncol = v)
  causal_mafs = matrix(rep(0, v*iter), ncol = v)
  
  set.seed(seed)
  #start_time = Sys.time()
  for(i in 1:iter) 
  {
    causal[i, ] = sample(c(1:v), size = v, replace = F) # choose random order variants are added to causal set
    for(c in 1:v)
    {
      #start_time = Sys.time()
      causal_set = rep(0, v)
      causal_set[causal[i, c(1:c)]] = 1 # only first c variants of the random order are in causal set
      causal_mafs[i, c] = sum(mafs[which(causal_set == 1)]) # sum of MAFs of causal variants
      beta = rep(0, v)
      if(eff_model == 1) { # -log10(MAF) model for effect sizes
        beta[causal_set == 1] = -log10(mafs[causal_set == 1]) 
      }
      if(eff_model == 2) { # draw effect sizes from normal distributions depending on MAF bins such that rarer bins have larger normal means
        ind_bin = which(mafs < 0.01 & mafs >= 0.001)
        if(length(ind_bin) > 0) {
          beta[intersect(ind_bin, causal[i, 1:c])] = rnorm(length(intersect(ind_bin, causal[i, 1:c])), mean = n_means[1], sd = n_sds[1])
        }
        ind_bin = which(mafs < 0.001 & mafs >= 0.0001)
        if(length(ind_bin) > 0) {
          beta[intersect(ind_bin, causal[i, 1:c])] = rnorm(length(intersect(ind_bin, causal[i, 1:c])), mean = n_means[2], sd = n_sds[2])
        }
        ind_bin = which(mafs < 0.0001 & mafs >= 0.00001)
        if(length(ind_bin) > 0) {
          beta[intersect(ind_bin, causal[i, 1:c])] = rnorm(length(intersect(ind_bin, causal[i, 1:c])), mean = n_means[3], sd = n_sds[3])
        }
        ind_bin = which(mafs < 0.00001)
        if(length(ind_bin) > 0) {
          beta[intersect(ind_bin, causal[i, 1:c])] = rnorm(length(intersect(ind_bin, causal[i, 1:c])), mean = n_means[4], sd = n_sds[4])
        }
      } 
      const = sqrt(0.5*h/sum((beta^2)*gene_mafs*(1 - gene_mafs))) # scale betas so that total heritability sums to h
      beta = const*beta
      
      # power of SV test
      sv_power_c[i, c] = calc_power_sv(n, h, gene_mafs, beta, alpha_sv) 
      
      # fix weights for variants in mask of aggregation tests; default here are weights recommended in SKAT paper (Wu et al, 2011)
      agg_wts = dbeta(gene_mafs, 1, 25)
      
      # power of burden test
      bd_power_c[i, c] = calc_power_burden(n, h, gene_mafs, beta, agg_wts, alpha_a)
      
      # power of SKAT 
      sk_power_c[i, c] = calc_power_skat(n, h, gene_mafs, beta, agg_wts, alpha_a, nsim = 2*(10^6))
      #Sys.time() - start_time
    }
  }
  #Sys.time() - start_time
  
  # return power for all tests, the random orders for causal sets, and cumulative sums of MAFs of causal variants 
  return(list(SV = sv_power_c,
              Burden = bd_power_c,
              SKAT = sk_power_c,
              Causal = causal,
              Causal_MAFs = causal_mafs))
}

## compute power for parameter choices 
n = 50000 # sample size
h = 0.001 # heritability
# MAFs of 20 rare variants in LINC01305 gene in chr2 of UKB exome sequence data; 
# you can read in another list of MAFs of your choice
gene_mafs = c(1.85e-6, 6.61e-6, 1.19e-5, 3.97e-6, 1.32e-5, 1.45e-5, 6.61e-6, 
          5.29e-6, 3.17e-5, 1.06e-5, 6.61e-6, 1.19e-5, 2.12e-5, 2.64e-6,
          5.29e-6, 2.64e-5, 5.28e-6, 5.29e-6, 2.12e-5, 2.64e-6)
alpha_sv = 5E-8 # level of significance for SV test
alpha_a = 2.5E-6 # level of significance for burden and SKAT tests
iter = 1 # number of random orders for causal sets
seed = 1234 
eff = c(1,2) # MAF-effect size model: 1 for -log10(MAF), 2 for step function normal distributions

# run with -log10(MAF) effect size model
pow_general = power_compute(gene_mafs, n, h, alpha_sv, alpha_a, iter, seed, eff[1])

# run with step function effect size model with n_means and n_sds being the mean and sds of the normal
# distributions for the four MAF bins - [0.001, 0.01), [0.0001, 0.001), [0.00001, 0.0001), and (0, 0.00001)
pow_general = power_compute(gene_mafs, n, h, alpha_sv, alpha_a, iter, seed, eff[2],
                    n_means = c(1, 2, 3, 4), n_sds = c(0.25, 0.25, 0.25, 0.25))

# create scatterplots like those in paper for average power over random causal sets for all three tests
pow = data.frame(Power = c(apply(pow_general$SV, 2, mean),
                           apply(pow_general$Burden, 2, mean),
                           apply(pow_general$SKAT, 2, mean)),
                 Test = rep(c("SV", "Burden", "SKAT"), each = length(gene_mafs)),
                 c = rep(c(1:length(gene_mafs)), times = 3))

# choose breaks for x-axis
br = 2
if(v > 40)
  br = 5
if(v > 100)
  br = 10
if(v > 200)
  br = 20
if(v > 400)
  br = 50

plot1 = ggplot(data = pow, aes(x = c, y = Power, group = Test, color = Test))+geom_line(size = 1.25)+
  geom_point(size=3)+ylim(0,1)+
  theme(plot.title = element_text(hjust = 0.5, size = 28), axis.title=element_text(size=20),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=26), #change legend title font size
        legend.text = element_text(size=20))+
  scale_x_continuous(breaks=seq(br, v, br)) + scale_color_aaas() + xlab("Number of causal variants")
plot1


# Simple case: calculate power of SV, burden, and SKAT tests under the assumption of independent variants,
# with all variants in mask, and also:
# * all variants have equal MAF
# * all causal variants have equal effect sizes
# * all variants have equal weights in the mask
# The function returns power for all three tests as number of causal variants increases from 1 to v.
calc_power = function(n, h, v, alpha1, alpha2)
{
  pow = data.frame(c(1:v), rep(0, v), rep(0, v), rep(0, v))
  colnames(pow) = c("m", "SV", "Burden", "SKAT")
  for(m in 1:v)
  {
    ncp1 = n*h/(m*(1-h))
    ncp2 = n*h*m/(v*(1-h))
    ncp3 = n*h/(1-h)
    pow[m,2] = 1 - (pchisq(qchisq(alpha1, df = 1, lower.tail = F), ncp = ncp1, df = 1)^m)*
      (pchisq(qchisq(alpha1, df = 1, lower.tail = F), ncp = 0, df = 1)^(v-m))
    pow[m,3] = pchisq(qchisq(alpha2, df = 1, lower.tail = F), ncp = ncp2, df = 1, lower.tail = F)
    pow[m,4] = pchisq(qchisq(alpha2, df = v, lower.tail = F), ncp = ncp3, df = v, lower.tail = F)
  }
  pow1 = data.frame(rep(c(1:v), 3),
                    rep(c("SV", "Burden", "SKAT"), each = v),
                    c(pow$SV, pow$Burden, pow$SKAT))
  colnames(pow1) = c("c", "Test", "Power")
  return(pow1)
}

# compute power for parameter choices
n = 50000 # sample size
v = 50 # number of variants
h = 0.001 # heritability

alpha1 = 5*10^(-8) # level of significance for SV test
alpha2 = 2.5*10^(-6) # level of significance for aggregation test

# run for simple case
pow_simple = calc_power(n, h, v, alpha1, alpha2)

# choose breaks for x-axis
br = 2
if(v > 40)
  br = 5
if(v > 100)
  br = 10
if(v > 200)
  br = 20
if(v > 400)
  br = 50

plot2 = ggplot(data = pow_simple, aes(x = c, y = Power, group = Test, color = Test))+geom_line(size = 1.25)+
  geom_point(size=3)+ylim(0,1)+
  theme(plot.title = element_text(hjust = 0.5, size = 28), axis.title=element_text(size=20),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=26), #change legend title font size
        legend.text = element_text(size=20))+
  scale_x_continuous(breaks=seq(br, v, br)) + scale_color_aaas() + xlab("Number of causal variants")
plot2

