source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
output_path = '~/Desktop/'

figureS11 = function(save_plot = F, output_path){
  # Filtered n_var<2 & coverage<20
  lambda_gene_caf = load_ukb_file('lambda_freq_filtered_gene_300k.txt.bgz')
  lambda_gene_caf = lambda_gene_caf %>%
    pivot_longer_lambda_data() %>%
    mutate(trait_type2 = factor(trait_type2, levels = trait_types),
           CAF_range = factor(CAF_range, levels = caf_types))

  figure = lambda_gene_caf %>%
    ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
    geom_point(alpha = 0.5) + ylim(0, 2) + theme_classic() +
    labs(y = 'Lambda GC', x = 'Number of Cases') +
    geom_hline(yintercept = 1, lty = 2) +
    scale_x_log10(label = comma, limits = c(2, NA)) +
    trait_color_scale + trait_fill_scale +
    facet_grid(result_type~CAF_range, labeller = label_type) + themes

  if(save_plot){
    png(output_path, height = 5, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
figureS11(save_plot = T, output_path = paste0(output_path, 'figureS11_lambda_gene_by_caf.png'))


