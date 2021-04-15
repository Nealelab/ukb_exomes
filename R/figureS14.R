source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
output_path = '~/Desktop/'

figureS14 = function(save_plot = F, output_path){
  lambda_gene_annt = load_ukb_file('lambda_annt_filtered_gene_300k.txt.bgz')
  lambda_gene_annt = lambda_gene_annt %>%
    pivot_longer_lambda_data() %>%
    mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
           annotation = factor(annotation,levels = annotation_types),
           trait_type2 = factor(trait_type2, levels = trait_types),)

  figure = lambda_gene_annt %>%
    ggplot + aes(x = lambda_gc, color = trait_type2, fill = trait_type2) +
    geom_density(alpha = 0.5) + theme_classic() +
    # scale_y_log10(label=comma) +
    labs(x = 'Lambda GC', y = 'Density') +
    geom_vline(xintercept = 1, lty = 2) +
    trait_color_scale + trait_fill_scale + themes +
    facet_grid(result_type~annotation, scale = 'free', labeller = label_type)

  if(save_plot){
    png(output_path, height = 5, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
figureS14(save_plot = T, output_path = paste0(output_path, 'figureS14_lambda_dist_by_pheno_filtered.png'))

