source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/skato/'

gene_sig_before = load_ukb_file('sig_cnt_before_gene_300k.txt.bgz')
var_sig_before = load_ukb_file('sig_cnt_before_var_300k.txt.bgz',
                        force_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(),
                                          most_severe_consequence = col_character(), mean_proportion = col_double()))
pheno_sig_before = load_ukb_file('pheno_sig_cnt_before_gene_300k.txt.bgz')
pheno_var_sig_before = load_ukb_file('pheno_sig_cnt_before_var_300k.txt.bgz')

gene_sig_after = load_ukb_file(paste0('sig_cnt_after_gene_300k_', test, '.txt.bgz'))
var_sig_after = load_ukb_file(paste0('sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'),
                        force_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(),
                                          most_severe_consequence = col_character(), mean_proportion = col_double()))
pheno_sig_after = load_ukb_file(paste0('pheno_sig_cnt_after_gene_300k_', test, '.txt.bgz'))
pheno_var_sig_after = load_ukb_file(paste0('pheno_sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'))

# genes_to_keep = load_ukb_file(paste0('genes_to_keep_300k_', test, '.txt.bgz'))
# phenos_to_keep = load_ukb_file(paste0('phenos_to_keep_300k_', test, '.txt.bgz'))

figure1 = function(save_plot = F, output_path){
  figure1_cnt = data.frame(
    type = rep(c('Variants', 'Groups', 'Phenotypes'), 2),
    cnt = as.integer(c(nrow(var_sig_before), nrow(gene_sig_before), nrow(pheno_sig_before), nrow(var_sig_after), nrow(gene_sig_after), nrow(pheno_sig_after))),
    filter = rep(c('Before Filtering', 'After Filtering'), each = 3)
  ) %>% mutate(filter = factor(filter, levels = c('Before Filtering', 'After Filtering')))

  figure1a = save_count_barplot_figure(cnt_data = figure1_cnt, cnt_type = 'Variants',  save_plot = T, output_path = paste0(output, 'figure1a_var_cnt.png'))
  figure1b = save_count_barplot_figure(cnt_data = figure1_cnt, cnt_type = 'Groups',  save_plot = T, output_path = paste0(output, 'figure1b_gene_cnt.png'))
  figure1c = save_count_barplot_figure(cnt_data = figure1_cnt, cnt_type = 'Phenotypes',  save_plot = T, output_path = paste0(output, 'figure1c_pheno_cnt.png'))

  gene_cnt_by_freq = format_count_by_freq_data(gene_sig_after, freq_col = 'CAF')
  var_cnt_by_freq = format_count_by_freq_data(var_sig_after, freq_col = 'AF')

  figure1d = save_count_by_freq_figure(cnt_data = var_cnt_by_freq, type = 'Variants', save_plot = T, output_path = paste0(output, 'figure1d_var_cnt_by_AF.png'))
  figure1e = save_count_by_freq_figure(cnt_data = gene_cnt_by_freq, type = 'Groups', save_plot = T, output_path = paste0(output, 'figure1e_gene_cnt_by_CAF.png'))

  figure1_top =  egg::ggarrange(figure1a, figure1b, figure1c, labels = c('(A) Variants', '(B) Groups', '(C) Phenotypes'), ncol = 3,
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