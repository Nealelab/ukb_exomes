source('~/ukb_exomes/R/constants.R')
gene_sig  = load_ukb_file('gene_sig300k.txt.bgz')
var_sig = load_ukb_file('var_sig300k.txt.bgz')

# Variant Matching
detach(package:plyr)
var_match = get_matched_data(var_sig, freq_col = 'AF')
var_match = var_match %>% mutate(annotation = factor(annotation, levels = annotation_types) )
sig_cnt_summary(var_match, sig_cnt_col = 'sig_pheno_cnt')
print_annotation_wilcoxon_test(var_match, 'sig_pheno_cnt', 'annotation')

p1 = ggplot(var_match, aes(x = sig_pheno_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 5, position = 'dodge') +
  labs(y = 'Variant Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes

p2 = var_match%>%
  na.omit() %>%
  ggplot(., aes(x = sig_pheno_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 1) +
  labs(y = 'Variant Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

# Gene Matching (SKATO)
library(plyr)
gene_match = get_matched_data(gene_sig[gene_sig$result_type=='SKATO', ], freq_col = 'caf')
gene_match = gene_match %>% mutate(annotation = factor(annotation, levels = annotation_types) )
detach(package:plyr)
sig_cnt_summary(gene_match, sig_cnt_col = 'sig_cnt')
print_annotation_wilcoxon_test(gene_match, 'sig_cnt', 'annotation')

p1 = ggplot(gene_match, aes(x = sig_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 5, position = 'dodge') +
  labs(y = 'Gene Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes
p2 = gene_match%>%
  na.omit() %>%
  ggplot(., aes(x = sig_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 1) +
  labs(y = 'Gene Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

# 470 Gene subset - Match by CAF
# gene_sub = read_delim('forKonrad_sig31kDNM_consensus_genes_2021_01_12.txt', delim='\t', col_types = all_cols)
sub = gene_sig %>%
  filter(result_type=='SKATO') %>%
  merge(gene_sub, ., by.x = c('gene_id', 'symbol'), by.y =c('gene_id', 'gene_symbol'), all.x =T)

full_sum = gene_sig %>%
  filter(result_type=='SKATO') %>%
  sig_cnt_summary()

sub_sum = sig_cnt_summary(sub)

sub_lof = gene_sig %>%
  filter(annotation == 'pLoF' & result_type == 'SKATO')%>%
  match_gene_subset_by_caf(gene_subset = gene_sub, sim_times = 1000)

sub_mis = gene_sig %>%
  filter(annotation == 'missense|LC'& result_type == 'SKATO')%>%
  match_gene_subset_by_caf(gene_subset = gene_sub, sim_times = 1000)

sub_syn = gene_sig %>%
  filter(annotation == 'synonymous'& result_type == 'SKATO')%>%
  match_gene_subset_by_caf(gene_subset = gene_sub, sim_times = 1000)

sim_sum1000 = rbind(sub_mis[[2]], sub_lof[[2]], sub_syn[[2]])
sim_sum1000$annotation = factor(rep(c('missense|LC', 'pLoF', 'synonymous'), each=1000), levels=annotation_types)

plt1 = ggplot(sim_sum1000, aes(x = means, color = annotation)) +
  geom_density(alpha = 0.5) + theme_bw() + labs(y = 'Density', x = 'Mean') +
  geom_vline(data = sub_sum, aes(xintercept = means, color = annotation), lty = 2) +
  annotation_color_scale + annotation_fill_scale + themes

plt2 = ggplot(sim_sum1000, aes(x = props, color = annotation)) +
  geom_density(alpha = 0.5 ) + theme_bw() +
  labs(y = 'Density', x = 'Proportion of Genes with 1 + Hit', title = 'Summary Statistics for Significant Count \n 1000 Simulations') +
  geom_vline(data = sub_sum, aes(xintercept = props, color = annotation), lty = 2) +
  scale_x_continuous(breaks = seq(0, 0.2, 0.01)) +
  annotation_color_scale + annotation_fill_scale + themes

# ggpubr::ggarrange(plt2, plt1, nrow=2)
