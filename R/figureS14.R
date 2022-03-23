source('~/ukb_exomes/R/constants.R')
detach(package:plyr)


figureS14 = function(save_plot = F, output_path){
  lambda_gene_annt = load_ukb_file(paste0('lambda_by_pheno_annt_filtered_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')

  lambda_gene_annt = lambda_gene_annt %>%
    select(-c(7:9)) %>%
    pivot_longer_lambda_data() %>%
    mutate(annotation = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
           annotation = factor(annotation,levels = annotation_types),
           trait_type2 = factor(trait_type2, levels = trait_types),)
  if(tranche=='500k'){
    lambda_gene_annt = lambda_gene_annt %>%
      filter(annotation != 'pLoF|missense|LC')
  }

  figure = lambda_gene_annt %>%
    #filter(trait_type2 != 'categorical' | annotation != 'pLoF'| result_type != 'skat') %>%
    ggplot + aes(x = lambda_gc, color = trait_type2, fill = trait_type2) +
    geom_density(alpha = 0.5) + theme_classic() +
    # geom_histogram(alpha = 0.5, position = 'dodge', bins = 10) + theme_classic() +
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
figureS14(save_plot = T, output_path = paste0(output_path, 'figureS14_lambda_dist_by_pheno_filtered_density.png'))

