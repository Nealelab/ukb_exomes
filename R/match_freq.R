source('~/ukb_exomes/R/constants.R')
gene_sig  <- load_ukb_file('gene_sig300k.txt.bgz')
var_sig <- load_ukb_file('var_sig300k.txt.bgz')

# Variant Matching
var_match <- get_matched_data(var_sig, freq_col = 'AF')
detach(package:plyr)
summary <- sig_cnt_summary(var_match, sig_cnt_col = 'sig_pheno_cnt')
summary
print_annotation_wilcoxon_test(var_match,'sig_pheno_cnt','annotation')

p1 <- ggplot(var_match, aes(x = sig_pheno_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 5, position = 'dodge') +
  labs(y = 'Variant Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes
var_match$annotation <- factor(var_match$annotation, levels = annotation_types)
p2 <- ggplot(var_match[complete.cases(var_match),], aes(x = sig_pheno_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 1) +
  labs(y = 'Variant Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)


# Gene Matching (SKATO)
library(plyr)
gene_match <- get_matched_data(gene_sig[gene_sig$result_type=='SKATO',], freq_col = 'caf')
detach(package:plyr)
summary <- sig_cnt_summary(gene_match, sig_cnt_col = 'sig_cnt')
summary
print_annotation_wilcoxon_test(gene_match,'sig_cnt','annotation')

p1 <- ggplot(gene_match, aes(x = sig_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 5, position = 'dodge') +
  labs(y = 'Gene Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes
var_match$annotation <- factor(var_match$annotation, levels = annotation_types)
p2 <- ggplot(gene_match[complete.cases(gene_match),], aes(x = sig_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 1) +
  labs(y = 'Gene Count', x = 'Association Count') +
  theme_bw() +scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

# 470 Gene subset - Match by CAF
# gene_sub <- read_delim('forKonrad_sig31kDNM_consensus_genes_2021_01_12.txt', delim='\t', col_types = all_cols)
sub_mis <- match_gene_subset_by_caf(full0[full0$annotation=='missense|LC',],  gene_subset = gene_sub, sim_times = 1000)
sub_lof <- match_gene_subset_by_caf(full0[full0$annotation=='pLoF',],  gene_subset = gene_sub, sim_times = 1000)
sub_syn <- match_gene_subset_by_caf(full0[full0$annotation=='synonymous',],  gene_subset = gene_sub, sim_times = 1000)

sim_sum1000 <- rbind(sub_mis[[2]],sub_lof[[2]], sub_syn[[2]])
sim_sum1000$annotation <- factor(rep(c('missense|LC', 'pLoF', 'synonymous'),each=1000),levels=annotation_types)


plt1 <- ggplot(sim_sum1000, aes(x = means, color = annotation)) +
  geom_density(alpha = 0.5) + theme_bw() + labs(y = 'Density', x = 'Mean') +
  geom_vline(data = sub_sum, aes(xintercept = means, color = annotation), lty = 2) +
  annotation_color_scale + annotation_fill_scale + themes

plt2 <- ggplot(sim_sum1000, aes(x = props, color = annotation)) +
  geom_density(alpha = 0.5 ) + theme_bw() +
  labs(y = 'Density', x = 'Proportion of Genes with 1 + Hit', title = 'Summary Statistics for Significant Count \n 1000 Simulations') +
  geom_vline(data = sub_sum, aes(xintercept = props, color = annotation), lty = 2) +
  scale_x_continuous(breaks = seq(0, 0.2, 0.01)) +
  annotation_color_scale + annotation_fill_scale + themes

# ggpubr::ggarrange(plt2, plt1, nrow=2)
