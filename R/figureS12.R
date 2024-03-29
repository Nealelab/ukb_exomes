source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figureS12 = function(save_plot = F, save_pdf = F, output_path){
  lambda_gene_full_before = load_ukb_file(paste0('lambda_by_pheno_full_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
  lambda_gene_full_after = load_ukb_file(paste0('lambda_by_pheno_full_filtered_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')

  lambda_gene_full_before = format_full_lambda_data(lambda_gene_full_before)
  lambda_gene_full_after = format_full_lambda_data(lambda_gene_full_after)

  figureS12a = lambda_gene_full_before %>%
    ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
    geom_point(alpha = 0.5, size = 0.8) + ylim(0, 2) + theme_classic() +
    labs(y = 'Lambda GC', x = '') +
    geom_hline(yintercept = 1, lty = 2) +
    scale_x_log10(label = comma, limits = c(2, NA)) +
    trait_color_scale + trait_fill_scale +
    facet_wrap(~result_type, ncol = 3, labeller = label_type) + themes +
    theme(plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"))+
    theme(axis.text.x = element_text(size = 6, vjust = 0.8))

  figureS12b = lambda_gene_full_after %>%
    ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
    geom_point(alpha = 0.5, size = 0.8) + ylim(0, 2) + theme_classic() +
    labs(y = 'Lambda GC', x = 'Number of Cases') +
    geom_hline(yintercept = 1, lty = 2) +
    scale_x_log10(label = comma, limits = c(2, NA)) +
    trait_color_scale + trait_fill_scale +
    facet_wrap(~result_type, ncol = 3, labeller = label_type) + themes +
    theme(legend.position = 'none',
          plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"),
          axis.text.x = element_text(size = 6, vjust = 0.8))

  figure = ggarrange(figureS12a, figureS12b,  labels = c('(A) Before', '(B) After'), nrow = 2,
                  label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  if(save_pdf){
    pdf(output_path, height = 116/in2mm, width = 174/in2mm)
    print(figure)
    dev.off()
  }

  return(figure)
}

figureS12(save_plot = T, save_pdf = T, output_path = paste0(output_path, 'final_pdf_figures/figureS12_lambda_gene_before_after_filter.pdf'))
# figureS12(save_plot = T, output_path = paste0(output_path, 'figureS12_lambda_gene_before_after_filter.png'))

