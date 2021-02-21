source('~/ukb_exomes/R/constants.R')

# --------Preliminary Lambda GC Check----------
lambda_gene_full  = load_ukb_file('lambda_gene_full300k.txt.bgz')
lambda_gene_annt  = load_ukb_file('lambda_gene_annt300k.txt.bgz')
lambda_gene_af  = load_ukb_file('lambda_gene_af300k.txt.bgz')
lambda_gene_af_annt  = load_ukb_file('lambda_gene_af_annt300k.txt.bgz')

lambda_var_full  = load_ukb_file('lambda_var_full300k.txt.bgz')
lambda_var_annt  = load_ukb_file('lambda_var_annt300k.txt.bgz')
lambda_var_af  = load_ukb_file('lambda_var_af300k.txt.bgz')
lambda_var_af_annt  = load_ukb_file('lambda_var_af_annt300k.txt.bgz')

lambda_by_gene  = load_ukb_file('lambda_by_gene300k.txt.bgz')

# --------QCed Lambda GC----------

# Data Cleaning
lambda_gene_full = lambda_gene_full %>%
  pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(result_type = str_split(labels, 'lambda_gc_') %>% map_chr(., 2),)%>%
  mutate(result_type = factor(result_type,levels = result_types),
         trait_type2 = factor(trait_type2, levels=trait_types),)

lambda_gene_annt = lambda_gene_annt %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
         result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2),) %>%
  mutate(annotation = factor(annotation,levels = annotation_types),
         result_type = factor(result_type,levels = result_types),
         trait_type2 = factor(trait_type2, levels=trait_types),)

lambda_gene_af = lambda_gene_af %>%
  pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(result_type = str_split(labels, 'lambda_gc_') %>% map_chr(., 2),) %>%
  mutate(result_type = factor(result_type,levels = result_types),
         trait_type2 = factor(trait_type2, levels=trait_types),)

lambda_gene_af_annt = lambda_gene_af_annt %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
         result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2),) %>%
  mutate(annotation = factor(annotation,levels = annotation_types),
         result_type = factor(result_type,levels = result_types),
         trait_type2 = factor(trait_type2, levels=trait_types),)

lambda_by_gene = lambda_by_gene %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
         result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2),) %>%
  mutate(annotation = factor(annotation,levels = annotation_types),
         result_type = factor(result_type,levels = result_types),)

lambda_var_full = lambda_var_full %>% mutate(trait_type2 = factor(trait_type2, levels=trait_types),)
lambda_var_af = lambda_var_af %>% mutate(trait_type2 = factor(trait_type2, levels=trait_types),)

lambda_var_annt = lambda_var_annt %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),) %>%
  mutate(annotation = factor(annotation,levels = annotation_types),
         trait_type2 = factor(trait_type2, levels=trait_types),)

lambda_var_af_annt = lambda_var_af_annt %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),) %>%
  mutate(annotation = factor(annotation,levels = annotation_types),
         trait_type2 = factor(trait_type2, levels=trait_types),)

# Plotting
lambda_gene_full %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~result_type, nrow=3, labeller = label_type) + themes

lambda_gene_annt %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_grid(result_type~annotation, labeller = label_type) + themes

lambda_gene_af %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(result_type~CAF_range, nrow = 3, labeller = label_type) + themes

lambda_gene_af_annt %>%
  filter(result_type == 'skato') %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_grid(CAF_range~annotation, labeller = label_type) + themes

lambda_var_full %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale + themes

lambda_var_annt %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~annotation, nrow = 3, labeller = label_type) + themes

lambda_var_af %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~AF_range, nrow = 3, labeller = label_type) + themes

lambda_var_af_annt %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(annotation~AF_range, nrow = 3, labeller = label_type) + themes

lambda_by_gene %>%
  mutate(lambda_gc = replace(lambda_gc, lambda_gc>2, 2)) %>%
  filter(trait_type == 'categorical') %>%
  ggplot + aes(x = lambda_gc, color = annotation, fill = annotation) +
  geom_density(alpha = 0.5) + theme_bw() +
  # scale_y_log10(label=comma) +
  labs(x = 'Lambda GC (Categorical)', y = 'Density') +
  geom_vline(xintercept = 1, lty = 2) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(result_type~annotation, scale='free', nrow=3, labeller = label_type)

# --------Significant Hits----------
gene_sig  = load_ukb_file('gene_sig300k.txt.bgz')
var_sig = load_ukb_file('var_sig300k.txt.bgz')
pheno_sig  = load_ukb_file('pheno_sig300k.txt.bgz')
pheno_var_sig  = load_ukb_file('pheno_var_sig300k.txt.bgz')
pheno_annt_gene = load_ukb_file('pheno_annt_gene300k.txt.bgz')
pheno_annt_var = load_ukb_file('pheno_annt_var300k.txt.bgz')
var_gene = load_ukb_file('var_gene_comparison300k.txt.bgz')

gene_sig = gene_sig %>%
  mutate(interval = get_freq_interval(caf), 
         annotation = factor(annotation, levels=annotation_types),
         result_type = factor(result_type, levels=result_types))
var_sig = var_sig %>%
  mutate(interval = get_freq_interval(AF), 
         annotation = factor(annotation, levels=annotation_types))

pheno_annt_gene = pheno_annt_gene %>%
  pivot_longer(cols = contains('_sig_cnt_'), names_to = 'labels', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(annotation = str_split(labels, '_sig_cnt_') %>% map_chr(., 1), 
         result_type = str_split(labels, '_sig_cnt_') %>% map_chr(., 2)) %>%
  mutate(annotation = factor(annotation, levels=annotation_types), 
         trait_type2 = factor(trait_type2, levels=trait_types),
         result_type = factor(result_type, levels=result_types))

pheno_annt_var = pheno_annt_var %>%
  select(., -sig_cnt) %>%
  pivot_longer(cols = contains('_sig_cnt'), names_to = 'labels', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(annotation = stringr::str_split(labels, '_sig_cnt') %>% map_chr(., 1)) %>%
  mutate(annotation = factor(annotation, levels=annotation_types), 
         trait_type2 = factor(trait_type2, levels=trait_types))

pheno_sig = pheno_sig %>%
  mutate(trait_type2 = factor(trait_type2, levels=trait_types),
         result_type = factor(result_type, levels=result_types))

pheno_var_sig = pheno_var_sig %>% mutate(trait_type2 = factor(trait_type2, levels=trait_types))
var_gene = var_gene %>% mutate(trait_type2 = factor(trait_type2, levels=trait_types))

# Gene-Level Association Count
plt = gene_sig %>% filter() %>%
  ggplot + aes(x = sig_cnt, color = annotation, fill = annotation) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_bw() + 
  labs(y = 'Gene Count', x = 'Association Count') +
  xlim(0, 50) + scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_grid(result_type~annotation, labeller = label_type)

plt = gene_sig %>% filter(result_type == 'skato') %>%
  ggplot + aes(x = caf, y = sig_cnt, color = annotation, label = gene_symbol) +
  geom_point(alpha = 0.5) + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n ( Gene-Level )') +
  geom_hline(yintercept = 50, lty = 2) +
  geom_text_repel(data = gene_sig[gene_sig$sig_cnt>50, ]) +
  scale_x_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

plt = gene_sig %>% filter(result_type == 'skato') %>%
  ggplot + aes(x = interval, y = sig_cnt, color = annotation, label = gene_symbol) +
  geom_boxplot() + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n (SKAT-O)') +
  geom_text_repel(data = gene_sig[gene_sig$sig_cnt>50 & gene_sig$result_type == 'skato', ]) +
  scale_x_discrete(limits = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')) +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

# Variant-Level Association Count
plt = var_sig %>% filter(!is.na(annotation)) %>%
  ggplot + aes(x = sig_pheno_cnt, color = annotation, fill = annotation) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_bw() + 
  labs(y = 'Variant Count', x = 'Association Count') +
  xlim(0, 120) + scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

plt1 = var_sig %>% filter(!is.na(annotation)) %>%
  ggplot + aes(x = AF, y = sig_pheno_cnt, color = annotation, label = paste(gene, ':', locus)) +
  geom_point(alpha = 0.5) + theme_bw() +
  labs(x = 'Allele Frequency', y = 'Association Count \n ( Variant-Level )') +
  geom_hline(yintercept = 100, lty = 2) +
  geom_text_repel(data = var_sig[var_sig$sig_pheno_cnt>100 & !is.na(var_sig$annotation),]) +
  scale_x_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

plt2 = var_sig %>% filter(!is.na(annotation)) %>%
  ggplot + aes(x = interval, y = sig_pheno_cnt, color = annotation, label = paste(gene, ':', locus)) +
  geom_boxplot() + theme_bw() +
  labs(x = 'Allele Frequency', y = 'Association Count \n ( Variant-Level )') +
  geom_text_repel(data = var_sig[var_sig$sig_pheno_cnt>100  & !is.na(var_sig$annotation), ]) +
  # scale_x_discrete(limits = c('[0, 0.0001)', '[0.0001, )')) +
  scale_x_discrete(limits = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')) +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

# Phenotype-Level Association Count (By Gene)
plt = pheno_sig %>% filter() %>%
  ggplot + aes(x = sig_cnt, color = trait_type2, fill = trait_type2) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_bw() + 
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Gene-Level )') +
  xlim(0, 60) + scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(result_type~trait_type2, labeller = label_type)

plt = pheno_sig %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Gene-Level )') +
  geom_hline(yintercept = 200, lty = 2) +
  geom_text_repel(data = pheno_sig[pheno_sig$sig_cnt>200, ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

plt = pheno_annt_gene %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.4) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Gene-Level )') +
  geom_hline(yintercept = 150, lty = 2) +
  geom_text_repel(data = pheno_annt_gene[pheno_annt_gene$sig_cnt>150 , ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(annotation~result_type, labeller = label_type)

# Phenotype-Level Association Count (By Variant)
plt = pheno_var_sig %>% filter() %>%
  ggplot + aes(x = sig_var_cnt, color = trait_type2, fill = trait_type2) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_bw() + 
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Variant-Level )') +
  xlim(0, 100) + scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(~trait_type2, ncol = 3, labeller = label_type)

plt = pheno_var_sig %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_var_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Variant-Level )') +
  geom_hline(yintercept = 40000, lty = 2) +
  geom_text_repel(data = pheno_var_sig[pheno_var_sig$sig_var_cnt>40000, ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes

plt = pheno_annt_var %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.4) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Variant-Level )') +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(annotation~trait_type2, scales='free', labeller = label_type)

# Gene Variant Association Count Comparison
plt = var_gene %>% filter() %>% 
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) + 
  geom_point( alpha = 0.5) + theme_bw() +
  labs(x = 'Number of Cases', y = 'Association Count') + 
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(trait_type2~cnt_type, scales = 'free', nrow = 3, labeller = label_type)

# Cumulative Allele Frequency Comparison
plt1 = gene_sig %>% filter(sig_cnt > 0) %>%
  ggplot + aes(x = caf, color = annotation) +
  stat_ecdf() +xlim(0, 0.001) + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Cumulative Distribution', title = 'Genes with Significant Hits') +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

plt2 = gene_sig %>% filter(sig_cnt == 0) %>%
  ggplot + aes(x = caf, color = annotation) +
  stat_ecdf() + xlim(0, 0.001) + theme_bw() +
  labs(x = 'Cumulative Allele Frequency', y = 'Cumulative Distribution', title = 'Genes without Significant Hits') +
  annotation_color_scale + annotation_fill_scale +themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

# Correlation: CAF vs. Association Count
detach(package:plyr)
print_freq_sig_cor()
print_freq_sig_cor(test='skat')
print_freq_sig_cor(test='burden')
print_freq_sig_cor(data = var_sig, test = 'variant', freq_col='AF', sig_col = 'sig_pheno_cnt')


