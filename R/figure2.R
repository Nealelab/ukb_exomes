source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figure2 = function(save_plot = F, output_path) {
  # gene_sig_without_filter  = load_ukb_file('gene_sig300k.txt.bgz')
  gene_sig_full_filtered  = load_ukb_file('sig_cnt_filtered_gene_300k.txt.bgz')
  gene_sig_full_filtered  = format_gene_sig_data(gene_sig_full_filtered)

  plt2 = gene_sig_full_filtered %>%
    filter(result_type == 'skato') %>%
    ggplot + aes(x = interval, y = sig_cnt, color = annotation, label = gene_symbol) +
    geom_boxplot() + theme_classic() +
    labs(x = 'Cumulative Allele Frequency', y = 'Association Count \n (SKAT-O)') +
    scale_x_discrete(limits = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')) +
    scale_y_log10(label = comma) +
    annotation_color_scale + annotation_fill_scale + themes +
    facet_wrap(~annotation, nrow = 3, scales = 'free', labeller = label_type)

  if(save_plot){
    png(output_path, height=6, width=5, units = 'in', res=300)
    print(plt2)
    dev.off()
  }
  return(plt2)
}