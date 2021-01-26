source('~/ukb_exomes/R/constants.R')
gene_sig  <- load_ukb_file('gene_sig300k.txt.bgz')
var_sig <- load_ukb_file('var_sig300k.txt.bgz')

# Variant Matching
match1 <- match_freq(data = var_sig,  annt1 = 'mis', annt2 = 'lof')
match2 <- match_freq(data = var_sig, annt1 = 'mis', annt2 = 'syn')
mis <- match1[[1]]
lof <- match1[[2]]
syn <- match2[[2]]
var_match <- rbind(mis, lof, syn)
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
match1 <- match_freq(data = gene_sig[gene_sig$result_type=='SKATO',],  annt1 = 'mis', annt2 = 'lof', freq_name = 'caf')
match2 <- match_freq(data = gene_sig[gene_sig$result_type=='SKATO',], annt1 = 'mis', annt2 = 'syn', freq_name = 'caf')
mis <- match1[[1]]
lof <- match1[[2]]
syn <- match2[[2]]
gene_match <- rbind(mis, lof, syn)
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