
library(tidyverse)
library(Matrix)
library(sparseMVN)
library(spdep)
library(qlcMatrix)
data = readMM('~/sparse._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx')
samples = read.table('~/sparse._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt')
n_samples = dim(data)[1]
n_phenos_to_simulate = 1000
decomp = Cholesky(data)

orig = data.frame(x=data@i + 1, y=data@j + 1, orig=data@x) %>% filter(x != y)
set.seed(663)
orig_uncorr = data.frame(x=round(runif(nrow(orig), max=n_samples)),
                         y=round(runif(nrow(orig), max=n_samples)),
                         orig=0)
orig = orig %>%
  union_all(orig_uncorr %>%
              anti_join(orig, by=c('x', 'y')))

set.seed(663)
phenos = rmvn.sparse(n_phenos_to_simulate, rep(0, n_samples), decomp, prec=FALSE)

output = orig %>%
  rowwise() %>%
  mutate(pheno_corr=cor(phenos[,x], phenos[,y]))

output %>%
  ggplot + aes(x = orig, y = pheno_corr) +
  geom_point() +
  labs(x = 'GRM', y = 'Random pheno correlation')

n_draws = 10
output_df = samples %>%
  transmute(userId=V1)

heritabilities = c(1, 0.5, 0.2, 0.1)
n_heritabilities = length(heritabilities)

set.seed(663)
for (j in 1:n_heritabilities) {
  heritability = heritabilities[j]
  if (heritability == 1) {
    prevalences = c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5)
  } else {
    prevalences = c(1e-3, 1e-2, 0.1, 0.2, 0.5)
  }
  n_prevalences = length(prevalences)
  qnorms = qnorm(prevalences)
  n_phenos = ((n_prevalences + 1)*n_draws)
  start = (j - 1) * n_phenos + 1
  end = j * n_phenos
  use_phenos = t(phenos[start:end,])
  noise_pheno = matrix(rnorm(length(use_phenos)), nrow=nrow(use_phenos))
  use_phenos = use_phenos*sqrt(heritability) + noise_pheno*sqrt(1 - heritability)
  
  test = orig %>%
    rowwise() %>%
    mutate(pheno_corr=cor(use_phenos[x,], use_phenos[y,]))
  
  print(paste('Heritability:', heritability, 
              '; Slope:', round(lm(pheno_corr ~ orig, test)$coefficients[[2]], 4),
              '; Correlation:', round(cor(test$pheno_corr, test$orig), 3)))
  
  for (i in 1:n_prevalences) {
    start = (i - 1) * n_draws + 1
    end = i * n_draws
    out = data.frame(use_phenos[,start:end] < qnorms[i])
    output_df = output_df %>%
      bind_cols(out %>%
                  rename_with(function(x) {paste0('p_', heritability, '_', prevalences[i], '_', gsub('X', '', x))}) %>%
                  mutate_all(as.numeric))
  }
  
  continuous_df = data.frame(use_phenos[,((n_prevalences - 1)*n_draws):(n_prevalences*n_draws)])
  output_df = output_df %>%
    bind_cols(continuous_df %>%
                rename_with(function(x) {paste0('p_', heritability, '_continuous_', gsub('X', '', x))}))
  
}
output_df %>% summary
output_df %>%
  write_tsv('random_phenos2.tsv')
