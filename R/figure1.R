source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/Desktop/final_figures/'

gene = load_ukb_file('gene_qc_metrics_ukb_exomes_300k.txt.bgz', subfolder = 'qc/')
var = load_ukb_file('variant_qc_metrics_ukb_exomes_300k.txt.bgz', subfolder = 'qc/')
pheno = load_ukb_file('pheno_qc_metrics_ukb_exomes_300k.txt.bgz', subfolder = 'qc/')

figure1 = function(save_plot = F, output_path){
  figure1_cnt = data.frame(
    type = rep(c('Variants', 'Groups', 'Phenotypes'), 2),
    cnt = as.integer(c(nrow(var), nrow(gene), nrow(pheno),
                       sum(var$keep_var_af & var$keep_var_annt),
                       sum(gene$keep_gene_caf & gene$keep_gene_coverage & gene$keep_gene_n_var & gene[paste0('keep_gene_',test)], na.rm = TRUE),
                       sum(pheno$keep_pheno_unrelated & pheno[paste0('keep_pheno_', test)]))),
    filter = rep(c('Before Filtering', 'After Filtering'), each = 3)
  ) %>% mutate(filter = factor(filter, levels = c('Before Filtering', 'After Filtering')))

  figure1a = save_count_barplot_figure(cnt_data = figure1_cnt, cnt_type = 'Phenotypes',  save_plot = T, output_path = paste0(output, 'figure1a_pheno_cnt.png'))
  figure1b = save_count_barplot_figure(cnt_data = figure1_cnt, cnt_type = 'Variants',  save_plot = T, output_path = paste0(output, 'figure1b_var_cnt.png'))
  figure1c = save_count_barplot_figure(cnt_data = figure1_cnt, cnt_type = 'Groups',  save_plot = T, output_path = paste0(output, 'figure1c_gene_cnt.png'))

  gene_cnt_by_freq = gene %>%
    filter(keep_gene_caf & keep_gene_coverage & keep_gene_n_var & get(paste0('keep_gene_',test))) %>%
    format_count_by_freq_data(freq_col = 'CAF')
  var_cnt_by_freq = var %>%
    filter(keep_var_af & keep_var_annt) %>%
    mutate(annotation = if_else(annotation %in% c('missense', 'LC'), 'missense|LC', annotation)) %>%
    format_count_by_freq_data(freq_col = 'AF')

  figure1d = save_count_by_freq_figure(cnt_data = var_cnt_by_freq, type = 'Variants', save_plot = T, output_path = paste0(output, 'figure1d_var_cnt_by_AF.png'))
  figure1e = save_count_by_freq_figure(cnt_data = gene_cnt_by_freq, type = 'Groups', save_plot = T, output_path = paste0(output, 'figure1e_gene_cnt_by_CAF.png'))

  figure1_top =  egg::ggarrange(figure1a, figure1b, figure1c, labels = c('(A) Phenotypes', '(B) Variants', '(C) Groups'), ncol = 3,
                                label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 1.5, hjust = 0))
  figure1_bottom = egg::ggarrange(figure1d, figure1e, labels = c('(D) Variants by AF Interval', '(E) Groups by CAF Interval'), ncol = 2,
                                  label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 1.5))
  figure = ggpubr::ggarrange(figure1_top, figure1_bottom, nrow = 2)
  if(save_plot){
    png(output_path, height = 6, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure1(save_plot = T, output_path = paste0(output, 'figure1_', test, '.png'))