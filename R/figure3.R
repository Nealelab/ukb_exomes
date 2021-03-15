source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

gene_sig_after = load_ukb_file('sig_cnt_after_gene_300k.txt.bgz')
var_sig_after = load_ukb_file('sig_cnt_after_SErm_var_300k.txt.bgz',
                        force_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(),
                                          most_severe_consequence = col_character(), mean_proportion = col_double()))
output = '~/Desktop/'

figure3 = function(save_plot = F, output_path){
  gene_sig_sum = format_sig_cnt_summary_data(gene_sig_after, freq_col = 'CAF', label = '_sig_pheno_cnt_skato')
  var_sig_sum = format_sig_cnt_summary_data(var_sig_after, freq_col = 'AF', label = '_sig_pheno_cnt')

  figure3a = save_prop_by_annt_freq_figure(var_sig_sum, output_path = paste0(output, 'figure3a.png'), save_plot = save_plot)
  figure3b = save_prop_by_annt_freq_figure(gene_sig_sum, output_path = paste0(output, 'figure3b.png'), save_plot = save_plot)
  figure3 = ggarrange(figure3a, figure3b, labels = c('(A) Single-Variant', '(B) SKAT-O'), nrow=2, label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
      png(output_path, height = 6, width = 7.5, units = 'in', res = 300)
      print(figure3)
      dev.off()
  }
  return(figure3)
}

figure3(save_plot = T, paste0(output, 'figure3.png'))
