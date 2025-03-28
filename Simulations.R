library(data.table)
library(stringr)
library(SKAT)
library(dplyr)

#initialize variables
n = 50000 # sample size
h = 0.002 # heritability
p_caus = c(0.8, 0.3, 0.01) # causal probs for PTV, deleterious missense, and other missense
random_sets = 100 # number of random causal sets
iter = 1 # number of simulation replicates per causal set
seed = 1234 # random seed

#score statistics
score_stat = function(g, y)
{
  return(sum(y*g, na.rm=T))
}

#calculate maf
maf_calc = function(x)
{
  return(min(x, 1-x))
}

#recode genotypes based on minor allele
recode_mac = function(x)
{
  return(ifelse(x == 0, 2, ifelse(x == 2, 0, 1)))
}

# read genotype file
data_full = fread("dummy_gene.txt")

# read annotation file - 3 columns indicating whether variant is a PTV, deleterious missense or other missense variant
annos = fread("annotations.txt")

#choose only rare variants and store annotations for mask
maf_full = colMeans(data_full)/2
inds = which(maf_full > 0.99) # implies variant is rare but minor and major alleles need recoding
if(length(inds) > 0)
{
  for(i in 1:length(inds))
  {
    data_full[,inds[i]] = unlist(sapply(data_full[, inds[i]], recode_mac))
  }
}
maf_full = unlist(sapply(maf_full, maf_calc)) # recalculate MAF 
gc()

inc_inds = as.numeric(which(maf_full < 0.01 & maf_full > 0)) # store indices of rare variants
data_full = data_full %>% select(inc_inds)
annos = annos[inc_inds, ]
maf_rare = maf_full[inc_inds]

# create indicator variables for masks
annos$PTV_Del_Mis = ifelse((annos$PTV + annos$Del_Mis) == 0, 0, 1)
annos$PTV_All_Mis = ifelse((annos$PTV + annos$Del_Mis + annos$Oth_Mis) == 0, 0, 1)

# create empty matrices that will store pvalues and number of variants in each simulation replicate
num_tests = 4
pv_ptv = matrix(rep(0, random_sets * iter * num_tests), nrow = random_sets * iter) # PTV
pv_ptv_del_mis = matrix(rep(0, random_sets * iter * num_tests), nrow = random_sets * iter) # PTV+Missense(Deleterious)
pv_ptv_all_mis = matrix(rep(0, random_sets * iter * num_tests), nrow = random_sets * iter) # PTV+Missense(All)

causal1 = which(annos$PTV == 1) # store indices of PTVs
causal2 = which(annos$Del_Mis == 1) # store indices of deleterious missense variants
causal3 = which(annos$Oth_Mis == 1) # store indices of other missense variants

for(j in 1:random_sets)
{
  #choose causal variants
  set.seed(seed + (j-1))
   
  # sample causal PTVs according to their causal probabilities (p_caus[1])
  causal = causal1[which(rbinom(length(causal1), 1, p_caus[1]) == 1)]
  
  # sample causal deleterious missense variants according to their causal probabilities (p_caus[2])
  causal = c(causal, causal2[which(rbinom(length(causal2), 1, p_caus[2]) == 1)])
  
  # sample causal other missense variants according to their causal probabilities (p_caus[3])
  causal = c(causal, causal3[which(rbinom(length(causal3), 1, p_caus[3]) == 1)])
  
  causal = sort(causal)
  
  # compute effect sizes for causal variants proportional to -log10(MAF)
  beta_true = rep(0, dim(annos)[1])
  if(h > 0)
  {
    const = sqrt(0.5*h/sum(((log10(maf_rare[causal]))^2)*maf_rare[causal]*(1-maf_rare[causal]), na.rm = T))
    beta_true[causal] = -const * log10(maf_rare[causal])
  }
  
  # perform all four tests and return chi-square statistics/pvalues 
  simulate_mask_pheno = function(y, g, n = 100000, h = 0.0005, mask_include)
  {
    # calculate score statistics
    y_centered = y - mean(y)
    sj = apply(g, 2, score_stat, y = y_centered)
    
    # calculate variance of score statistics
    var_sj = suppressWarnings((apply(g, 2, stats::var, na.rm = T)*(n-1))*stats::var(y, na.rm=T)*(n-1)/n)
    var_sj[var_sj == 0] = 1E10 # set variance of "non-variants" arbitrarily high (signifies these are not included in the test in this iteration)
    
    # calculate test statistic and pvalue for SV test
    tj = (sj^2)/var_sj
    sv_pv = pchisq(max(tj), df = 1, lower.tail = F)
    
    # set aggregation test weights to Beta(1,25) as recommended in SKAT paper (Wu et al, 2011)
    agg_wts = c(1, 25)
    null_model = SKAT_Null_Model(y ~ 1)
    
    # find burden test pvalue
    burden_model = SKAT(as.matrix(g[, mask_include]), null_model, r.corr = 1, weights.beta = agg_wts)
    burden_pv = burden_model$p.value
    
    # find SKAT pvalue
    skat_model = SKAT(as.matrix(g[, mask_include]), null_model, weights.beta = agg_wts)
    skat_pv = skat_model$p.value
    
    # find SKAT-O pvalue
    skato_model = SKAT(as.matrix(g[, mask_include]), null_model, weights.beta = agg_wts, method = "optimal.adj")
    skato_pv = skato_model$p.value
    
    # return SV, burden, SKAT, and SKAT-O test pvalues
    return(c(sv_pv, burden_pv, skat_pv, skato_pv))
  }
  
  set.seed(seed + (j-1))
  
  for(i in 1:iter)
  {
    # sample n individuals from full genotype matrix
    test_samples = sample(c(1:dim(data_full)[1]), size = n)
    geno = as.matrix(data_full[test_samples,])
    gc()
    
    maf = colMeans(geno)/2 # compute MAFs of variants in our sample
    
    # create masks
    mask_include1 = which(maf > 0 & maf < 0.01 & annos$PTV == 1) # PTV
    mask_include2 = which(maf > 0 & maf < 0.01 & annos$PTV_Del_Mis == 1) # PTV+Missense(Deleterious)
    mask_include3 = which(maf > 0 & maf < 0.01 & annos$PTV_All_Mis == 1) # PTV+Missense(All)
    
    # simulate normal trait
    y = geno %*% beta_true + rnorm(n, sd = sqrt(1-h))
    
    # run all tests with PTV mask
    if(mask_vars[((j-1)*iter) + i, 1] > 0)
    {
      sim_iter = simulate_mask_pheno(y, geno, n, h, mask_include1)
      pv_ptv[((j-1)*iter) + i,] = sim_iter
    }
    gc()
    
    # run all tests with PTV+Missense(Deleterious) mask
    if(mask_vars[((j-1)*iter) + i, 2] == mask_vars[((j-1)*iter) + i, 1])
    {
      pv_ptv_del_mis[((j-1)*iter) + i,] = pv_ptv[((j-1)*iter) + i,]
    } else {
      sim_iter = simulate_mask_pheno(y, geno, n, h, mask_include2)
      pv_ptv_del_mis[((j-1)*iter) + i,] = sim_iter
    }
    gc()
    
    # run all tests with PTV+Missense(All) mask
    if(mask_vars[((j-1)*iter) + i, 3] == mask_vars[((j-1)*iter) + i, 2])
    {
      pv_ptv_all_mis[((j-1)*iter) + i,] = pv_ptv_del_mis[((j-1)*iter) + i,]
    } else {
      sim_iter = simulate_mask_pheno(y, geno, n, h, mask_include3)
      pv_ptv_all_mis[((j-1)*iter) + i,] = sim_iter
    }
    gc()
  }
}

# store results for all replicates in .txt files
pv_ptv = as.data.frame(pv_ptv)
colnames(pv_ptv) = c("SV", "Burden", "SKAT", "SKATO")
fwrite(pv_ptv, file = "pvalues_ptv.txt", col.names = T, sep = "\t")

pv_ptv_del_mis = as.data.frame(pv_ptv_del_mis)
colnames(pv_ptv_del_mis) = c("SV", "Burden", "SKAT", "SKATO")
fwrite(pv_ptv_del_mis, file = "pvalues_ptv_del_mis.txt", col.names = T, sep = "\t")

pv_ptv_all_mis = as.data.frame(pv_ptv_all_mis)
colnames(pv_ptv_all_mis) = c("SV", "Burden", "SKAT", "SKATO")
fwrite(pv_ptv_all_mis, file = "pvalues_ptv_all_mis.txt", col.names = T, sep = "\t")

# compute power for tests as proportion of significant pvalues
alpha_sv = 5E-8 # level of significance for SV tests
alpha_a = 2.5E-6 # level of significance for aggregation tests

pow_sv = length(which(pv_ptv$SV < alpha_sv))/nrow(pv_ptv)
# power for aggregation tests with PTV mask
pow_burden_ptv = length(which(pv_ptv$Burden < alpha_a))/nrow(pv_ptv)
pow_skat_ptv = length(which(pv_ptv$SKAT < alpha_a))/nrow(pv_ptv)
pow_skato_ptv = length(which(pv_ptv$SKATO < alpha_a))/nrow(pv_ptv)
