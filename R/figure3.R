source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/skato/'

gene_sig_after = load_ukb_file(paste0('sig_cnt_after_gene_300k_', test, '.txt.bgz'))
var_sig_after = load_ukb_file(paste0('sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'),
                        force_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(),
                                          most_severe_consequence = col_character(), mean_proportion = col_double()))

figure3c = function(save_plot=F, output_path){
  pext = var_sig_after %>%
    filter(tx_annotation_csq == most_severe_consequence &
             most_severe_consequence %in% c('splice_acceptor_variant', 'splice_donor_variant') &
             lof == 'HC') %>%
    mutate(mean_prop_bin = factor(get_mean_prop_interval(mean_proportion_expressed), levels = c('[0, 0.2]', '(0.2, 0.8]', '(0.8, 1]'))) %>%
    group_by(mean_prop_bin, most_severe_consequence) %>%
    sig_cnt_summary(sig_col = 'all_sig_pheno_cnt')

  figure = pext %>%
    filter(!is.na(mean_prop_bin)) %>%
    ggplot + aes(x = most_severe_consequence, y = prop, ymin = prop-sd, ymax = prop+sd, color = mean_prop_bin, fill = mean_prop_bin) +
    geom_pointrange(stat = "identity", position = position_dodge(width = 0.6), size = 0.4) +
    labs(y = 'Proportion', x = 'Splice Variant')  +
    scale_y_continuous(label = label_percent(accuracy = 0.1)) +
    scale_x_discrete(labels = c('Acceptor', 'Donor')) +
    scale_color_manual(name = 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
    scale_fill_manual(name = 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
    theme_classic() + themes + theme(panel.spacing = unit(1, "lines"))
  if(save_plot){
    png(output_path, height = 4, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure3 = function(save_plot = F, output_path){
  gene_sig_sum = format_sig_cnt_summary_data(gene_sig_after, freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt')
  var_sig_sum = format_sig_cnt_summary_data(var_sig_after, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt')
  clinvar_sum = var_sig_after %>% filter(!is.na(pathogenicity)) %>%
    mutate(interval = get_freq_interval(AF)) %>%
    mutate(interval = factor(interval, levels = c('[2e-05, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )'))) %>%
    group_by(interval, pathogenicity) %>%
    sig_cnt_summary('all_sig_pheno_cnt')

  figure3a = save_prop_by_annt_freq_figure(var_sig_sum, output_path = paste0(output, 'figure3a.png'), save_plot = save_plot)
  figure3b = save_prop_by_annt_freq_figure(gene_sig_sum, output_path = paste0(output, 'figure3b.png'), save_plot = save_plot)
  figure3c = figure3c(save_plot = T, paste0(output, 'figure3c_splice_variant.png'))
  figure3d = save_group_matched_figure2(clinvar_sum, levels = c('P/LP', 'VUS', 'B/LB'),
                              labels = c( 'Pathogenic/Likely Pathogenic', 'Uncertain significance', 'Benign/Likely Benign'),
                              group_col = 'pathogenicity', x_lab = 'ClinVar Pathogenicity', output_path = paste0(output, 'figure3d_clinvar_pathogenicity.png'), save_plot = F) +
    facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) +
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = 'bold'),
          axis.text= element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95) )
  figure3 = ggarrange(figure3a, figure3b, figure3c, figure3d, labels = c('(A) Single-Variant', '(B) SKAT-O', '(C) Splice Variants', '(D) ClinVar Pathogenicity'), nrow=4, label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
      png(output_path, height = 12, width = 7.5, units = 'in', res = 300)
      print(figure3)
      dev.off()
  }
  return(figure3)
}

figure3(save_plot = T, paste0(output, 'figure3_', test, '.png'))

