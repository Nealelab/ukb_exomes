source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

gene_sig_after = load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
var_sig_after = load_ukb_file(paste0('var_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/', force_cols = var_cols) %>%
  format_variant_call_stats(.)

if(tranche=='500k'){
  gene_sig_after = gene_sig_after %>%
    filter(annotation != 'pLoF|missense|LC')
}

save_prop_by_annt_freq_figure_500k_supp = function(matched_summary, type, output_path, save_plot = F){
  matched_summary = matched_summary %>%
      filter(!(interval %in% c('(0.0001, 0.001]', '(0.001, 0.005]','(0.005, 0.01]'))) %>%
      mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names))
  plt = matched_summary %>%
      ggplot + aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, color = annotation, fill = annotation) +
      geom_pointrange(stat = "identity", position = position_dodge(width = 2)) +
      labs(y = 'Proportion', x = NULL)  +
      scale_y_continuous(label = label_percent(accuracy = 1)) +
      annotation_color_scale + annotation_fill_scale  +
      facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) + theme_classic() + themes+
      theme(panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(face = 'bold'),
            axis.text= element_text(size = 8),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95),
            legend.key.size = unit(0.5, 'cm'),
            axis.title.y = element_text(size = 9, face = 'bold'),
            legend.title = element_text(size = 8, face = 'bold'),
            legend.text = element_text(size = 7),) +
    # scale_x_discrete(labels = c(matched_summary$labels)) +
    geom_text(aes(y = 0.002, label=paste(cnt, type), vjust = -1), size = 1.9)
  if(save_plot){
    png(output_path, height = 4, width = 7.5, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

figureS19 <- function(save_plot=F, save_pdf = F, output_path){
  if(tranche=='500k'){
    gene_sig_sum = format_sig_cnt_summary_data_500k(gene_sig_after, freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt')
    var_sig_sum = format_sig_cnt_summary_data_500k(var_sig_after, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt')
  }else{
    gene_sig_sum = format_sig_cnt_summary_data(gene_sig_after, freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt')
    var_sig_sum = format_sig_cnt_summary_data(var_sig_after, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt')
  }
  figure1 = save_prop_by_annt_freq_figure_500k_supp(var_sig_sum, 'variants' ,  output_path = paste0(output, 'figure3b_supp_var_',tranche, '_', test,'.png'), save_plot = TRUE)
  figure2 = save_prop_by_annt_freq_figure_500k_supp(gene_sig_sum, 'genes', output_path = paste0(output, 'figure3b_supp_gene_',tranche, '_', test,'.png'), save_plot = TRUE)

  figure = ggarrange(figure1, figure2, labels = c('(A) Single-Variant', '(B) SKAT-O'), nrow=2, label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
      png(output_path, height = 6, width = 5, units = 'in', res = 300)
      print(figure)
      dev.off()
  }
  if(save_pdf){
    ggsave(filename=output_path, plot=figure, device=cairo_pdf,
       width = 114/in2mm, height = 140/in2mm, units = "in", dpi = 300)
    # pdf(output_path, height = 140/in2mm, width = 114/in2mm)
    # print(figure)
    # dev.off()
  }
  return(figure)
}

figureS19(save_plot = T, save_pdf = T, paste0(output_path, 'final_pdf_figures/figureS19_common_variant_gene_',tranche, '_', test,'.pdf'))
# figureS19(save_plot = T, paste0(output_path, 'figureS19_common_variant_gene_',tranche, '_', test,'.png'))

