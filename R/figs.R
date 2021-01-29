source('~/ukb_exomes/R/constants.R')

gene_sig  <- load_ukb_file('gene_sig300k.txt.bgz')
var_sig <- load_ukb_file('var_sig300k.txt.bgz')
pheno_sig  <- load_ukb_file('pheno_sig300k.txt.bgz')
pheno_var_sig  <- load_ukb_file('pheno_var_sig300k.txt.bgz')
pheno_annt_gene <- load_ukb_file('pheno_annt_gene300k.txt.bgz')
pheno_annt_var <- load_ukb_file('pheno_annt_var300k.txt.bgz')
var_gene <- load_ukb_file('var_gene_comparison300k.txt.bgz')

gene_sig$interval <- get_freq_interval(gene_sig$caf)
var_sig$interval <- get_freq_interval(var_sig$AF)

pheno_annt_gene <- pheno_annt_gene %>%
  pivot_longer(cols = contains('_sig_cnt_'), names_to = 'labels', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(annotation = str_split(labels, "_sig_cnt_") %>% map_chr(., 1)) %>%
  mutate(result_type = str_split(labels, "_sig_cnt_") %>% map_chr(., 2))

pheno_annt_var <- pheno_annt_var[,-8] %>%
  pivot_longer(cols = contains('_sig_cnt'), names_to = 'labels', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(annotation = stringr::str_split(labels, "_sig_cnt") %>% map_chr(., 1))
pheno_annt_var$annotation[pheno_annt_var$annotation=='missense'] = 'missense|LC'

gene_sig$annotation <- factor(gene_sig$annotation,levels=annotation_types)
var_sig$annotation <- factor(var_sig$annotation,levels=annotation_types)
pheno_annt_gene$annotation <- factor(pheno_annt_gene$annotation,levels=annotation_types)
pheno_annt_var$annotation <- factor(pheno_annt_var$annotation,levels=annotation_types)

pheno_sig$trait_type2 <- factor(pheno_sig$trait_type2, levels=trait_types)
pheno_var_sig$trait_type2 <- factor(pheno_var_sig$trait_type2, levels=trait_types)
pheno_annt_gene$trait_type2 <- factor(pheno_annt_gene$trait_type2, levels=trait_types)
pheno_annt_var$trait_type2 <- factor(pheno_annt_var$trait_type2, levels=trait_types)
var_gene$trait_type2 <- factor(var_gene$trait_type2, levels=trait_types)

# Gene-Level Association Count
plt <- ggplot(gene_sig, aes(x = sig_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) +
  labs(y = 'Gene Count', x = 'Association Count') +
  theme_bw() + xlim(0, 50) + scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_grid(result_type~annotation, labeller = label_type)

# ggplot(gene_sig, aes(x = caf, y = sig_cnt, color = annotation, label = gene_symbol)) +
plt <- ggplot(gene_sig[gene_sig$result_type == 'SKATO', ], aes(x = caf, y = sig_cnt, color = annotation, label = gene_symbol)) +
  geom_point( alpha = 0.5) + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n ( Gene-Level )') +
  geom_hline(yintercept = 50, lty = 2) +
  geom_text_repel(data = gene_sig[gene_sig$sig_cnt>50, ]) +
  scale_x_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)
# facet_wrap(~result_type, nrow = 3)

plt <- ggplot(gene_sig[gene_sig$result_type == 'SKATO', ], aes(x = interval, y = sig_cnt, color = annotation, label = gene_symbol)) +
# ggplot(gene_sig[gene_sig$result_type == 'SKATO', ], aes(x = interval, y = sig_cnt, color = annotation, label = gene_symbol)) +
  geom_boxplot() + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n (SKATO)') +
  # labs(x = 'Cumulative Allele Frequency', y = 'Association Count') +
  geom_text_repel(data = gene_sig[gene_sig$sig_cnt>50 & gene_sig$result_type == 'SKATO', ]) +
  scale_x_discrete(limits = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')) +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  # facet_grid(result_type~annotation, scales = 'free', labeller = label_type)
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)



# Variant-Level Association Count
plt <- ggplot(var_sig[!is.na(var_sig$annotation), ], aes(x = sig_pheno_cnt, color = annotation, fill = annotation)) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) +
  labs(y = 'Variant Count', x = 'Association Count') +
  theme_bw() + xlim(0, 120) +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

# pdf('variant_association_dist.pdf', height=12, width=10)
# print(p)
# dev.off()

var_sig0 <- var_sig[!is.na(var_sig$annotation),]
plt1 <- var_sig0 %>%
  ggplot(aes(x = AF, y = sig_pheno_cnt, color = annotation, label = paste(gene, ':', locus))) +
  geom_point( alpha = 0.5) + theme_bw() +
  labs(x = 'Allele Frequency', y = 'Association Count \n ( Variant-Level )') +
  geom_hline(yintercept = 100, lty = 2) +
  geom_text_repel(data = var_sig[var_sig$sig_pheno_cnt>100, ]) +
  scale_x_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

plt2 <- var_sig0 %>%
  ggplot(aes(x = interval, y = sig_pheno_cnt, color = annotation, label = paste(gene, ':', locus))) +
  geom_boxplot() + theme_bw() +
  labs(x = 'Allele Frequency', y = 'Association Count \n ( Variant-Level )') +
  geom_text_repel(data = var_sig[var_sig$sig_pheno_cnt>100, ]) +
  # scale_x_discrete(limits = c('[0, 0.0001)', '[0.0001, )')) +
  scale_x_discrete(limits = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')) +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

# plt <- ggarrange(plt1, plt2, ncol=2)

# Phenotype-Level Association Count (By Gene)
plt <- ggplot(pheno_sig, aes(x = sig_cnt, color = trait_type2, fill = trait_type2)) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) +
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Gene-Level )') +
  theme_bw() + xlim(0, 60) +
  scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(result_type~trait_type2, labeller = label_type)

plt <- ggplot(pheno_sig, aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode)) +
  geom_point( alpha = 0.5) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Gene-Level )') +
  geom_hline(yintercept = 200, lty = 2) +
  geom_text_repel(data = pheno_sig[pheno_sig$sig_cnt>200, ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

plt <- ggplot(pheno_annt_gene, aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode)) +
  geom_point( alpha = 0.4) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Gene-Level )') +
  geom_hline(yintercept = 150, lty = 2) +
  geom_text_repel(data = pheno_annt_gene[pheno_annt_gene$sig_cnt>150 , ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(annotation~result_type, labeller = label_type)

plt <- ggplot(pheno_sig, aes(x = sig_cnt, color = trait_type2, fill = trait_type2)) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) +
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Gene-Level )') +
  theme_bw() + xlim(0, 60) +
  scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(result_type~trait_type2, labeller = label_type)

# pdf('pheno_gene_association_dist.pdf', height=6, width=10)
# print(plt)
# dev.off()


# Phenotype-Level Association Count (By Variant)
plt <- ggplot(pheno_var_sig, aes(x = sig_var_cnt, color = trait_type2, fill = trait_type2)) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) +
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Variant-Level )') +
  theme_bw() + xlim(0,100) +
  scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(~trait_type2, ncol = 3, labeller = label_type)

plt <- ggplot(pheno_var_sig, aes(x = n_cases, y = sig_var_cnt, color = trait_type2, label = phenocode)) +
  geom_point( alpha = 0.5) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Variant-Level )') +
  geom_hline(yintercept = 40000, lty = 2) +
  geom_text_repel(data = pheno_var_sig[pheno_var_sig$sig_var_cnt>40000, ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes

plt <- ggplot(pheno_annt_var, aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode)) +
  geom_point( alpha = 0.4) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Variant-Level )') +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(annotation~trait_type2, scales='free', labeller = label_type)

# Gene Variant Association Count Comparison
plt <- ggplot(var_gene, aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode)) + geom_point( alpha = 0.5) +
  labs(x = 'Number of Cases', y = 'Association Count') + theme_bw() +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(trait_type2~cnt_type, scales = 'free', nrow = 3, labeller = label_type)

# Cumulative Allele Frequency Comparison
plt1 <- ggplot(gene_sig[gene_sig$sig_cnt > 0, ], aes(x = caf, color = annotation)) +
  stat_ecdf() +xlim(0, 0.001) + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Cumulative Distribution', title = 'Genes with Significant Hits') +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

plt2 <- ggplot(gene_sig[gene_sig$sig_cnt == 0, ], aes(x = caf, color = annotation)) +
  stat_ecdf() + xlim(0, 0.001) + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Cumulative Distribution', title = 'Genes without Significant Hits') +
  annotation_color_scale + annotation_fill_scale +themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

# plt <- ggarrange(plt1, plt2, ncol = 2)

# Correlation: CAF vs. Association Count
detach(package:plyr)
print_freq_sig_cor()
print_freq_sig_cor(test='SKAT')
print_freq_sig_cor(test='Burden Test')
print_freq_sig_cor(data = var_sig, test = 'variant', freq_col='AF', sig_col = 'sig_pheno_cnt')