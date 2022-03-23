source('~/ukb_exomes/R/constants.R')

pheno_corr <- read_delim('/Users/wlu/OneDrive/2020Work/Pleiotropy/subset_data_pheno_correlation_bgz.txt', delim = '\t',
                         col_types = cols(i_pheno = col_character(), j_pheno = col_character()))
corr = load_ukb_file(paste0('correlation_table_phenos_', tranche,'.txt.bgz'), subfolder = 'qc/')
corr2 <- corr %>%
  mutate(i_pheno = str_sub(str_replace_all(i_data, '[^[:alnum:]]',''), 29,33),
         j_pheno = str_sub(str_replace_all(j_data, '[^[:alnum:]]',''), 29,35),) %>%
  mutate(j_pheno = str_replace_all(j_pheno, '[a-z]',''))



# --------Missense Check----------
test = 'skato'
gene_info = load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', subfolder = 'analysis/')
gene_info <- gene_info %>%
  mutate(mis_bin = ntile(gene_info$oe_mis, 10))
gene_sig_after = load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
gene_sig_after = gene_sig_after %>%
  filter(annotation != 'pLoF|missense|LC')

mis_constrained = gene_sig_after %>%
    merge(gene_info[, c('gene_id', 'gene', 'mis_bin')], ., by.x = c('gene_id', 'gene'), by.y =c('gene_id', 'gene_symbol')) %>%
    filter(mis_bin %in% c(1,2)) %>% distinct(gene_id)
matched_constrained_sum = get_subset_matched_data_summary(gene_sig_after, subset = mis_constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = 'all_sig_pheno_cnt', oversample = 1000,
                                                          ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Constrained')

matched_constrained_sum %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot +
    geom_pointrange(aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, group = group, color = annotation, fill=annotation, pch = group),
                    stat = "identity", position = position_dodge(width = 1)) +
    labs(y = 'Proportion', x = NULL, alpha = NULL)  +
    scale_y_continuous(label = label_percent(accuracy = 1), breaks = c(0,0.05,0.10,0.15)) +
    scale_x_discrete(labels = annotation_names, limits = rev(levels(matched_constrained_sum$annotation))) +
    annotation_color_scale + annotation_fill_scale  +
    scale_shape_manual(name = NULL, values = c(1, 16)) +
    coord_flip(ylim = c(0,0.2)) +
    theme(panel.spacing = unit(1, "lines"),
        axis.text= element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face='bold', size=14),
        strip.text.y = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.title = element_text(size = 13))

# --------Preliminary Phenotype Check----------
# pre_qc_pheno = load_ukb_file(paste0('pheno_results.txt.bgz'), subfolder = 'results/')
pre_qc_pheno = read_delim('/Users/wlu/Desktop/ukb_exomes_450k/pre_filter/pheno_text_file/pheno_results_without_pharma_tag_new.txt.bgz', delim = '\t',)
pre_qc_pheno = pre_qc_pheno %>%
  select('n_cases_defined', 'n_controls', 'trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier', 'description', 'description_more', 'coding_description', 'category')
write.table(pre_qc_pheno, paste0('/Users/wlu/Desktop/ukb_exomes_450k/pre_filter/pheno_text_file/ukb_exome_450k_phenotypes_pharma_tag_removed_new.txt'), quote = FALSE, row.names = FALSE, sep = '\t')

# --------Preliminary CAF*n_cases Check----------
lambda_gene_full_before  = load_ukb_file(paste0('lambda_by_pheno_full_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
lambda_gene_full_after  = load_ukb_file(paste0('lambda_by_pheno_full_filtered_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')


# --------Preliminary Lambda GC Check----------
lambda_gene_full_before  = load_ukb_file(paste0('lambda_by_pheno_full_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
lambda_gene_annt_before  = load_ukb_file(paste0('lambda_by_pheno_annt_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
lambda_gene_caf_before  = load_ukb_file(paste0('lambda_by_pheno_freq_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
lambda_by_gene_before  = load_ukb_file(paste0('lambda_by_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
lambda_by_gene_before_coverage  = load_ukb_file(paste0('lambda_by_gene_before_coverage_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')

# Data Cleaning
lambda_gene_full_before = lambda_gene_full_before %>%
  pivot_longer_lambda_data() %>%
  mutate(trait_type2 = factor(trait_type2, levels=trait_types))

lambda_gene_annt_before = lambda_gene_annt_before %>%
  select(-c(7:9)) %>%
  pivot_longer_lambda_data() %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1) ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types),
         trait_type2 = factor(trait_type2, levels=trait_types), )

lambda_gene_caf_before = lambda_gene_caf_before %>%
  pivot_longer_lambda_data() %>%
  mutate(trait_type2 = factor(trait_type2, levels=trait_types),
         CAF_range = factor(CAF_range, levels = caf_types), )

lambda_by_gene_before = lambda_by_gene_before %>%
  select(-c(6:8)) %>%
  pivot_longer_lambda_data() %>%
  mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
         annotation = factor(annotation, levels = annotation_types))

lambda_by_gene_before_coverage = lambda_by_gene_before_coverage %>%
  select(-c(6:8)) %>%
  pivot_longer_lambda_data( ) %>%
  mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
         coverage_int = get_coverage_interval(mean_coverage),
         annotation = factor(annotation, levels = annotation_types),
         coverage_int = factor(coverage_int, levels = c('[0, 10]', '(10, 20]', '(20, 30]', '(30, 40]', '(40, 50]', '(50, )'))) %>%
  filter(trait_type != 'icd10')

## Plotting
p_lambda_gene_full_before = lambda_gene_full_before %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 0.5) + ylim(0, 2) + theme_classic() +
  labs(y = 'Lambda GC', x = '') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~result_type, ncol=3, labeller = label_type) + themes

png(paste0(path_to, 'lambda_gene_full_pre_filter.png'), height = 3, width = 7, units = 'in', res = 300)
print(p_lambda_gene_full_before)
# print(ggarrange(p_lambda_gene_full_before, p_lambda_gene_full_after, labels = c('Before', 'After'), nrow = 2,
#                   label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2)))
dev.off()

p_lambda_gene_caf_bin_before = lambda_gene_caf_before %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 0.5) + ylim(0, 2) + theme_classic() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_grid(result_type~CAF_range, labeller = label_type) + themes

png(paste0(path_to, 'lambda_gene_caf_bin_pre_filter.png'), height = 5, width = 10, units = 'in', res = 300)
print(p_lambda_gene_caf_bin_before)
dev.off()

detach(package:plyr)
gene_cnt_before = lambda_by_gene_before %>%
  filter(trait_type == 'all') %>%
  group_by(result_type, annotation) %>%
  summarise(cnt = sum(lambda_gc>2, na.rm = TRUE))
p_lambda_dist_by_gene_before = lambda_by_gene_before %>%
  mutate(lambda_gc = replace(lambda_gc, lambda_gc>2, 2)) %>%
  filter(trait_type == 'all') %>%
  ggplot + aes(x = lambda_gc, color = annotation, fill = annotation) +
  geom_density(alpha = 0.5) + theme_classic() +
  # scale_y_log10(label=comma) +
  labs(x = 'Gene Lambda GC', y = 'Density') +
  geom_vline(xintercept = 1, lty = 2) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_grid(result_type~annotation, scale='free', labeller = label_type) +
  geom_text(data = gene_cnt_before, aes(label = paste('Genes with lambda>2: ', as.character(cnt)), x = 1.5, y = 2.5, face = 'bold', color = annotation, group = NULL), size = 2)

png(paste0(path_to, 'lambda_dist_by_gene_without_filter.png'), height = 5, width = 9, units = 'in', res = 300)
print(p_lambda_dist_by_gene_before)
dev.off()

p_lambda_gene_annt_before = lambda_gene_annt_before %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size=0.5) + ylim(0, 2) + theme_classic() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_grid(result_type~annotation, labeller = label_type) + themes

png(paste0(path_to, 'lambda_gene_annt_pre_filter.png'), height = 5, width = 7.5, units = 'in', res = 300)
print(p_lambda_gene_annt_before)
dev.off()


# --------QCed Lambda GC----------
lambda_gene_full  = load_ukb_file('lambda_full_filtered_gene_300k.txt.bgz')
lambda_gene_annt  = load_ukb_file('lambda_annt_filtered_gene_300k.txt.bgz')
lambda_gene_caf  = load_ukb_file('lambda_freq_filtered_gene_300k.txt.bgz')
lambda_gene_caf_annt  = load_ukb_file('lambda_freq_annt_filtered_gene_300k.txt.bgz')

lambda_var_full  = load_ukb_file('lambda_full_filtered_SErm_var_300k.txt.bgz')
lambda_var_annt  = load_ukb_file('lambda_annt_filtered_SErm_var_300k.txt.bgz')
lambda_var_af  = load_ukb_file('lambda_freq_filtered_SErm_var_300k.txt.bgz')
lambda_var_af_annt  = load_ukb_file('lambda_freq_annt_SErm_filtered_var_300k.txt.bgz')

lambda_by_gene  = load_ukb_file('lambda_by_gene_filtered_300k.txt.bgz')
lambda_var_ac  = load_ukb_file('lambda_expectedAC_filtered_var_300k.txt.bgz')

# Data Cleaning
lambda_gene_full = lambda_gene_full %>%
  pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(result_type = str_split(labels, 'lambda_gc_') %>% map_chr(., 2), )%>%
  mutate(result_type = factor(result_type, levels = result_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), )

lambda_gene_annt = lambda_gene_annt %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1), 
         result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2), ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types), 
         result_type = factor(result_type, levels = result_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), )

lambda_gene_caf = lambda_gene_caf %>%
  pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(result_type = str_split(labels, 'lambda_gc_') %>% map_chr(., 2), ) %>%
  mutate(result_type = factor(result_type, levels = result_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), 
         CAF_range = factor(CAF_range, levels = caf_types), )

lambda_gene_caf_annt = lambda_gene_caf_annt %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1), 
         result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2), ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types), 
         result_type = factor(result_type, levels = result_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), 
         CAF_range = factor(CAF_range, levels = caf_types), )

lambda_by_gene = lambda_by_gene %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1), 
         result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2), ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types), 
         result_type = factor(result_type, levels = result_types), )

lambda_var_full = lambda_var_full %>% mutate(trait_type2 = factor(trait_type2, levels=trait_types), )
lambda_var_af = lambda_var_af %>%
  mutate(trait_type2 = factor(trait_type2, levels=trait_types), 
         AF_range = factor(AF_range, levels = af_types))

lambda_var_annt = lambda_var_annt %>%
  select(-lambda_gc) %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1), ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), )

lambda_var_af_annt = lambda_var_af_annt %>%
  select(-lambda_gc) %>%
  pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1), ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), 
         AF_range = factor(AF_range, levels = af_types))

lambda_var_ac = lambda_var_ac %>%
  pivot_longer(cols = contains('lambda_gc'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
  mutate(ac_type = str_split(labels, 'lambda_gc') %>% map_chr(., 2), 
         trait_type2 = if_else(trait_type == 'icd_first_occurrence', 'icd10', trait_type))%>%
  mutate(ac_type = factor(ac_type, levels = ac_types), 
         trait_type2 = factor(trait_type2, levels=trait_types), )
lambda_var_ac_cnt = lambda_var_ac %>% group_by(ac_type) %>% summarise(cnt = sum(!is.na(lambda_gc)))

# Plotting


p_lambda_gene_full_after = lambda_gene_full %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + ylim(0, 2) + theme_classic() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~result_type, ncol=3, labeller = label_type) + themes +
  theme(legend.position = 'none')

png('~/Figures/supp_figures/figureS14_lambda_gene_before_after.png', height = 4, width = 6, units = 'in', res = 300)
print(ggarrange(p_lambda_gene_full_before, p_lambda_gene_full_after, labels = c('Before', 'After'), nrow = 2, 
                  label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2)))
dev.off()

p_lambda_gene_annt_after = lambda_gene_annt %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + ylim(0, 2) + theme_classic() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_grid(result_type~annotation, labeller = label_type) + themes

png('~/Figures/lambdaGC/filtered/lambda_gene_annt_filtered.png', height = 5, width = 7.5, units = 'in', res = 300)
print(p_lambda_gene_annt_after)
dev.off()



p_lambda_gene_caf_bin_after = lambda_gene_caf %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + ylim(0, 2) + theme_classic() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_grid(result_type~CAF_range, labeller = label_type) + themes

# Filtered n_var<2 & coverage<20
png('~/Figures/lambdaGC/filtered/lambda_gene_caf_bin_filtered.png', height = 5, width = 9, units = 'in', res = 300)
print(p_lambda_gene_caf_bin_after)
dev.off()

png('~/Figures/supp_figures/figureS13_lambda_gene_caf_bin.png', height = 5, width = 9, units = 'in', res = 300)
print(p_lambda_gene_caf_bin_after)
dev.off()

p_lambda_dist_by_pheno_after = lambda_gene_annt %>%
  ggplot + aes(x = lambda_gc, color = trait_type2, fill = trait_type2) +
  geom_density(alpha = 0.5) + theme_classic() +
  # scale_y_log10(label=comma) +
  labs(x = 'Lambda GC', y = 'Density') +
  geom_vline(xintercept = 1, lty = 2) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(result_type~annotation, scale='free', labeller = label_type)

p_lambda_dist_by_pheno_before = lambda_gene_annt_before %>%
  ggplot + aes(x = lambda_gc, color = trait_type2, fill = trait_type2) +
  geom_density(alpha = 0.5) + theme_classic() +
  # scale_y_log10(label=comma) +
  labs(x = 'Lambda GC', y = 'Density') +
  geom_vline(xintercept = 1, lty = 2) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(result_type~annotation, scale='free', labeller = label_type)

png('~/Figures/lambdaGC/filtered/lambda_dist_by_pheno_filtered.png', height = 5, width = 7.5, units = 'in', res = 300)
print(p_lambda_dist_by_pheno_after)
dev.off()

png('~/Figures/supp_figures/figureS16_lambda_by_pheno_filtered.png', height = 5, width = 7.5, units = 'in', res = 300)
print(p_lambda_dist_by_pheno_after)
dev.off()

png('~/Figures/lambdaGC/without_filter/lambda_dist_by_pheno_without_filter.png', height = 5, width = 9, units = 'in', res = 300)
print(p_lambda_dist_by_pheno_before)
dev.off()

detach(package:plyr)
gene_cnt_after = lambda_by_gene %>%
  filter(trait_type == 'all') %>%
  group_by(result_type, annotation) %>%
  summarise(cnt = sum(lambda_gc>2, na.rm = TRUE))
p_lambda_dist_by_gene_after = lambda_by_gene %>%
  mutate(lambda_gc = replace(lambda_gc, lambda_gc>2, 2)) %>%
  filter(trait_type == 'all') %>%
  ggplot + aes(x = lambda_gc, color = annotation, fill = annotation) +
  geom_density(alpha = 0.5) + theme_classic() +
  # scale_y_log10(label=comma) +
  labs(x = 'Gene Lambda GC', y = 'Density') +
  geom_vline(xintercept = 1, lty = 2) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_grid(result_type~annotation, scale='free', labeller = label_type) +
  geom_text(data = gene_cnt_after, aes(label = paste('Genes with \n lambda>2: ', as.character(cnt)), x = 1.5, y = 2.5, face = 'bold', color = annotation, group = NULL), size = 3)

png('~/Figures/lambdaGC/filtered/lambda_dist_by_gene_filtered.png', height = 5, width = 7.5, units = 'in', res = 300)
print(p_lambda_dist_by_gene_after)
dev.off()

png('~/Figures/supp_figures/figureS17_lambda_by_gene_filtered.png', height = 5, width = 7.5, units = 'in', res = 300)
print(p_lambda_dist_by_gene_after)
dev.off()



p_lambda_var_ac_after = lambda_var_ac %>%
  ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5, size = 2) + ylim(0, 2) + theme_bw() +
  labs(y = 'Lambda GC', x = 'Number of Cases') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_log10(label = comma, limits = c(2, NA)) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~ac_type, nrow = 1, labeller = label_type) + themes +
  geom_text(data = lambda_var_ac_cnt, aes(label = paste('Defined Phenotype: ', as.character(cnt)), x = 100, y = 1.75, face = 'bold', color = NULL, group = NULL), size = 3)

png('~/Figures/lambdaGC/filtered/lambda_var_expectedAC_filtered.png', height = 3, width = 15, units = 'in', res = 300)
print(p_lambda_var_ac_after)
dev.off()

p_lambda_var_ac_dist_after =lambda_var_ac %>%
  ggplot + aes(x = lambda_gc, color = trait_type2, fill = trait_type2) +
  geom_density(alpha = 0.5)  + theme_bw() +
  labs(y = 'Density', x = 'Lambda GC') +
  geom_vline(xintercept = 1, lty = 2) +
  trait_color_scale + trait_fill_scale +
  facet_wrap(~ac_type, nrow = 1, scale = 'free', labeller = label_type) + themes

png('~/Figures/lambdaGC/filtered/lambda_var_expectedAC_dist_filtered.png', height = 3, width = 15, units = 'in', res = 300)
print(p_lambda_var_ac_dist_after)
dev.off()

# --------Preliminary Significant Hits----------
test = 'skato'
gene_sig_before = load_ukb_file('sig_cnt_before_gene_300k.txt.bgz')
var_sig_before = load_ukb_file('sig_cnt_before_var_300k.txt.bgz', force_cols = var_cols)
pheno_sig_before = load_ukb_file('pheno_sig_cnt_before_gene_300k.txt.bgz')
pheno_var_sig_before = load_ukb_file('pheno_sig_cnt_before_var_300k.txt.bgz')

gene_sig_after = load_ukb_file(paste0('sig_cnt_after_gene_300k_', test, '.txt.bgz'))
var_sig_after = load_ukb_file(paste0('sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'), force_cols = var_cols)
pheno_sig_after = load_ukb_file(paste0('pheno_sig_cnt_after_gene_300k_', test, '.txt.bgz'))
pheno_var_sig_after = load_ukb_file(paste0('pheno_sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'))
var_gene_by_pheno = load_ukb_file(paste0('var_gene_comparison_by_pheno_after_300k_', test, '.txt.bgz'))

gene_sig = gene_sig_after %>%
  mutate(interval = get_freq_interval(CAF), 
         annotation = factor(annotation, levels=annotation_types)) %>%
  mutate(interval = factor(interval, levels = c('[2e-05, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')))
var_sig = var_sig_after %>%
  mutate(interval = get_freq_interval(AF), 
         annotation = factor(annotation, levels=annotation_types))%>%
  mutate(interval = factor(interval, levels = c('[2e-05, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')))

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
var_gene_by_pheno = var_gene_by_pheno %>% mutate(trait_type2 = factor(trait_type2, levels=trait_types))
var_gene = var_gene_by_pheno %>%
  pivot_longer(cols = contains('_sig_cnt'), names_to = 'cnt_type', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(cnt_type = factor(cnt_type, levels = c('pheno_gene_var_sig_cnt', 'pheno_gene_sig_cnt', 'pheno_var_sig_cnt', 'pheno_none_sig_cnt'), 
                           labels = c('Gene and Variant', 'Gene Only', 'Variant Only', 'None') ))

# Gene-Level Association Count
p_gene_sig_dist = gene_sig %>% filter() %>%
  ggplot + aes(x = sig_cnt, color = annotation, fill = annotation) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_classic() + 
  labs(y = 'Gene Count', x = 'Association Count') +
  xlim(0, 50) + scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_grid(result_type~annotation, labeller = label_type)

png('~/Figures/sig_cnt_filtered/gene_sig_dist.png', height = 5, width = 9, units = 'in', res = 300)
print(p_gene_sig_dist)
dev.off()

p_gene_sig_scatter = gene_sig %>% filter(result_type == 'skato') %>%
  ggplot + aes(x = CAF, y = sig_cnt, color = annotation, label = gene_symbol) +
  geom_point(alpha = 0.5) + theme_classic() +
  labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n ( Burden Test )') +
  geom_text_repel(data = gene_sig[gene_sig$sig_cnt>50, ]) +
  scale_x_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

png('~/Figures/sig_cnt_filtered/gene_sig_scatter_skato.png', height = 6, width = 5, units = 'in', res = 300)
print(p_gene_sig_scatter)
dev.off()

p_gene_sig_box = gene_sig %>%
  ggplot + aes(x = interval, y = all_sig_pheno_cnt_skato, color = annotation, label = gene_symbol) +
  geom_boxplot() + theme_classic() +
  labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n ( SKAT-O)') +
  # geom_text_repel(data = gene_sig[gene_sig$sig_cnt>50 & gene_sig$result_type == 'skato', ]) +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

png('~/Figures/sig_cnt_figures/gene_sig_boxplot_skato.png', height = 6, width = 5, units = 'in', res = 300)
print(p_gene_sig_box)
dev.off()

var_sum = var_sig_after %>%
  group_by(annotation) %>%
  sig_cnt_summary('all_sig_pheno_cnt')
plt = var_sum %>%
  mutate(annotation = factor(annotation, levels = annotation_types) )%>%
  ggplot + aes(x=annotation, y=prop, ymin=prop-sd, ymax=prop+sd, color=annotation, fill=annotation) +
  geom_pointrange(stat="identity", position=position_dodge(width = 0.8)) +
  labs(y = 'Proportion', x = '')  +
  scale_y_continuous(label=label_percent(accuracy=0.1)) +
  scale_x_discrete(labels = annotation_names) +
  annotation_color_scale + annotation_fill_scale  +
  theme_classic() + themes+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(face='bold'), 
        axis.text= element_text(size = 10), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.95))

png('~/Figures/sig_cnt_figures/var_sig_prop.png', height = 4, width = 6, units = 'in', res = 300)
print(plt)
dev.off()

# Variant-Level Association Count
p_var_sig_dist = var_sig %>% filter(!is.na(annotation)) %>%
  ggplot + aes(x = sig_pheno_cnt, color = annotation, fill = annotation) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_classic() + 
  labs(y = 'Variant Count', x = 'Association Count') +
  xlim(0, 120) + scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

png('~/Figures/sig_cnt_filtered/var_sig_dist.png', height = 6, width = 5, units = 'in', res = 300)
print(p_var_sig_dist)
dev.off()

p_var_sig_scatter = var_sig %>% filter(!is.na(annotation)) %>%
  ggplot + aes(x = AF, y = sig_pheno_cnt, color = annotation, label = paste(gene, ':', locus)) +
  geom_point(alpha = 0.5) + theme_classic() +
  labs(x = 'Allele Frequency', y = 'Association Count \n ( Variant-Level )') +
  # geom_text_repel(data = var_sig[var_sig$sig_pheno_cnt>100 & !is.na(var_sig$annotation), ]) +
  scale_x_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, labeller = label_type)

png('~/Figures/sig_cnt_filtered/var_sig_scatter.png', height = 6, width = 5, units = 'in', res = 300)
print(p_var_sig_scatter)
dev.off()

p_var_sig_boxplot = var_sig %>%
  ggplot + aes(x = interval, y = all_sig_pheno_cnt, color = annotation, label = paste(gene, ':', locus)) +
  geom_boxplot() + theme_classic() +
  labs(x = 'Allele Frequency', y = 'Association Count \n ( Variant-Level )') +
  scale_y_log10(label = comma) +
  annotation_color_scale + annotation_fill_scale + themes +
  facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

png('~/Figures/sig_cnt_figures/var_sig_boxplot.png', height = 6, width = 5, units = 'in', res = 300)
print(p_var_sig_boxplot)
dev.off()

# Phenotype-Level Association Count (By Gene)
p_pheno_gene_sig_dist = pheno_sig %>% filter() %>%
  ggplot + aes(x = sig_cnt, color = trait_type2, fill = trait_type2) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_classic() + 
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Gene-Level )') +
  xlim(0, 60) + scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(result_type~trait_type2, labeller = label_type)

png('~/Figures/sig_cnt_filtered/pheno_gene_sig_dist.png', height = 5, width = 9, units = 'in', res = 300)
print(p_pheno_gene_sig_dist)
dev.off()

p_pheno_gene_sig_scatter = pheno_sig %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + theme_classic() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Gene-Level )') +
 # geom_text_repel(data = pheno_sig[pheno_sig$sig_cnt>200, ]) +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(~result_type, nrow = 3, labeller = label_type)

png('~/Figures/sig_cnt_filtered/pheno_gene_sig_scatter.png', height = 6, width = 5, units = 'in', res = 300)
print(p_pheno_gene_sig_scatter)
dev.off()

p_pheno_gene_annt = pheno_annt_gene %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.4) + theme_classic() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Gene-Level )') +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(annotation~result_type, labeller = label_type)

png('~/Figures/sig_cnt_filtered/pheno_gene_annt.png', height = 5, width = 9, units = 'in', res = 300)
print(p_pheno_gene_annt)
dev.off()

# Phenotype-Level Association Count (By Variant)
p_pheno_var_sig_dist = pheno_var_sig %>% filter() %>%
  ggplot + aes(x = sig_var_cnt, color = trait_type2, fill = trait_type2) +
  geom_histogram(alpha = 0.5, binwidth = 1 ) + theme_classic() + 
  labs(y = 'Phenotype Count', x = 'Association Count \n ( Variant-Level )') +
  xlim(0, 100) + scale_y_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_wrap(~trait_type2, ncol = 3, labeller = label_type)

png('~/Figures/sig_cnt_filtered/pheno_var_sig_dist.png', height = 3, width = 6, units = 'in', res = 300)
print(p_pheno_var_sig_dist)
dev.off()

p_pheno_var_sig_scatter = pheno_var_sig %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_var_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.5) + theme_classic() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Variant-Level )') +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes

png('~/Figures/sig_cnt_filtered/pheno_var_sig_scatter.png', height = 4, width = 6, units = 'in', res = 300)
print(p_pheno_var_sig_scatter)
dev.off()

p_pheno_var_sig_annt = pheno_annt_var %>% filter() %>%
  ggplot + aes(x = n_cases, y = sig_cnt, color = trait_type2, label = phenocode) +
  geom_point(alpha = 0.4) + theme_classic() +
  labs(x = 'Number of Cases', y = 'Association Count \n ( Variant-Level )') +
  scale_x_log10(label = comma) +
  trait_color_scale + trait_fill_scale + themes +
  facet_grid(annotation~trait_type2, scales='free', labeller = label_type)

png('~/Figures/sig_cnt_filtered/pheno_var_sig_annt.png', height = 5, width = 9, units = 'in', res = 300)
print(p_pheno_var_sig_annt)
dev.off()

# --------Summary Statistics----------
gene = load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/') %>% mutate(annotation = factor(annotation, levels = annotation_types))
pheno_var = load_ukb_file(paste0('pheno_sig_cnt_filtered_', test, '_var_', tranche,'.txt.bgz'), subfolder = 'analysis/')
var = load_ukb_file(paste0('variant_qc_metrics_ukb_exomes_', tranche,'.txt.bgz'), subfolder = 'qc/') %>%
  mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))%>%
  mutate(annotation = factor(annotation, levels = annotation_types))
pheno = load_ukb_file(paste0('pheno_qc_metrics_ukb_exomes_', tranche,'.txt.bgz'), subfolder = 'qc/') %>% mutate(trait_type = factor(trait_type, levels = trait_types))

gene %>% group_by(annotation) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt))
gene %>% group_by(annotation) %>% summarise(cnt = n())
gene %>% filter(keep_gene_n_var) %>% group_by(annotation) %>% summarise(cnt = n())
gene %>% filter(keep_gene_n_var & keep_gene_coverage) %>% group_by(annotation) %>% summarise(cnt = n())
gene %>% filter(keep_gene_expected_ac & keep_gene_coverage & keep_gene_n_var) %>% group_by(annotation) %>% summarise(cnt = n())
gene %>% filter(keep_gene_expected_ac & keep_gene_coverage & keep_gene_n_var & keep_gene_skato) %>% group_by(annotation) %>% summarise(cnt = n())

gene %>% filter(keep_gene_n_var) %>% summarise(cnt = n())
gene %>% filter(keep_gene_n_var & keep_gene_coverage) %>% summarise(cnt = n())
gene %>% filter(keep_gene_expected_ac & keep_gene_coverage & keep_gene_n_var) %>% summarise(cnt = n())
gene %>% filter(keep_gene_expected_ac & keep_gene_coverage & keep_gene_n_var & keep_gene_skato) %>% summarise(cnt = n())

gene %>% group_by(annotation) %>% summarise(cnt = n())
as.data.frame(gene %>% filter(keep_gene_n_var) %>% group_by(annotation) %>% summarise(cnt = n()))[,2]/
  as.data.frame(gene %>% group_by(annotation) %>% summarise(cnt = n()))[,2]
as.data.frame(gene %>% filter(keep_gene_n_var & keep_gene_coverage) %>% group_by(annotation) %>% summarise(cnt = n()))[,2]/
  as.data.frame(gene %>% group_by(annotation) %>% summarise(cnt = n()))[,2]
as.data.frame(gene %>% filter(keep_gene_expected_ac & keep_gene_coverage & keep_gene_n_var) %>% group_by(annotation) %>% summarise(cnt = n()))[,2]/
  as.data.frame(gene %>% group_by(annotation) %>% summarise(cnt = n()))[,2]
as.data.frame(gene %>% filter(keep_gene_expected_ac & keep_gene_coverage & keep_gene_n_var & keep_gene_skato) %>% group_by(annotation) %>% summarise(cnt = n()))[,2]/
  as.data.frame(gene %>% group_by(annotation) %>% summarise(cnt = n()))[,2]

var = var %>% mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation))
var %>% group_by(annotation) %>% summarise(cnt = n())
nrow(var %>% filter(keep_var_annt))
var %>% filter(keep_var_annt & keep_var_expected_ac) %>% group_by(annotation) %>% summarise(cnt = n())

pheno %>% group_by(trait_type) %>% summarise(cnt = n())
pheno %>% filter(keep_pheno_skato) %>% group_by(trait_type) %>% summarise(cnt = n())
pheno %>% filter(keep_pheno_skato & keep_pheno_unrelated) %>% group_by(trait_type) %>% summarise(cnt = n())

pheno_var %>% group_by(trait_type) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_var_cnt))


gene_sig_before %>% group_by(annotation) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt_skato))
var_sig_before %>% group_by(annotation) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt))
pheno_sig_before %>% group_by(trait_type) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_gene_cnt_skato))
pheno_var_sig_before %>% group_by(trait_type) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_var_cnt))

gene_sig_before %>% group_by(annotation) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt_burden))
pheno_sig_before %>% group_by(trait_type) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_gene_cnt_burden))

gene_sig_after  %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt))
var_sig_after %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt))
gene_sig_after %>% group_by(annotation) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt))
var_sig_after %>% group_by(annotation) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_pheno_cnt))
pheno_sig_after %>% group_by(trait_type2) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_gene_cnt))
pheno_var_sig_after %>% group_by(trait_type2) %>% summarise(cnt = n(), sig_cnt = sum(all_sig_var_cnt))

gene_sig_long = gene_sig_after %>%
  pivot_longer(cols = contains('_sig_pheno_cnt_'), names_to = 'labels', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(trait_type = str_split(labels, '_sig_pheno_cnt_') %>% map_chr(., 1), 
         result_type = str_split(labels, '_sig_pheno_cnt_') %>% map_chr(., 2), ) %>%
  group_by(annotation, trait_type) %>%
  summarise(cnt = n(), sig_cnt = sum(sig_cnt, na.rm = T))

var_sig_long = var_sig_after %>%
  pivot_longer(cols = contains('_sig_pheno_cnt'), names_to = 'labels', names_repair = 'unique', values_to = 'sig_cnt') %>%
  mutate(trait_type = str_split(labels, '_sig_pheno_cnt') %>% map_chr(., 1)) %>%
  group_by(annotation, trait_type) %>%
  summarise(cnt = n(), sig_cnt = sum(sig_cnt, na.rm = T))

get_two_col_correlation_table(var_sig_after, group_col = 'annotation', col1 = 'AF', col2 = 'all_sig_pheno_cnt')
get_two_col_correlation_table(gene_sig_after, group_col = 'annotation', col1 = 'CAF', col2 = 'all_sig_pheno_cnt')

save_var_gene_comparison_table(filter = T, normalize = T, save_plot = T, output_path = paste0(output, 'var_gene_comparison_after_norm_', test, '_3annt', tranche, '.png'))
save_var_gene_comparison_table(filter = T, normalize = F, save_plot = T, output_path = paste0(output, 'var_gene_comparison_after_count_', test, '_3annt', tranche, '.png'))
save_var_gene_comparison_table(filter = F, normalize = T, save_plot = T, output_path = paste0(output, 'var_gene_comparison_before_norm_', test, '_', tranche, '.png'))
save_var_gene_comparison_table(filter = F, normalize = F, save_plot = T, output_path = paste0(output, 'var_gene_comparison_before_count_', test, '_', tranche, '.png'))