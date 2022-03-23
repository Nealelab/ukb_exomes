source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figureS15 = function(save_plot = F, output_path){
  lambda_by_gene = load_ukb_file(paste0('lambda_by_gene_filtered_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
  if(tranche == '500k'){
    lambda_by_gene = lambda_by_gene %>%
      filter(annotation != 'pLoF|missense|LC')
  }
  lambda_by_gene = lambda_by_gene %>%
    select(-c(6:8)) %>%
    pivot_longer_lambda_data()%>%
    mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
           annotation = factor(annotation,levels = annotation_types))

  gene_cnt = lambda_by_gene %>%
    filter(trait_type == 'all') %>%
    group_by(result_type, annotation) %>%
    summarise(cnt = sum(lambda_gc > 2,na.rm = TRUE))

  figure = lambda_by_gene %>%
    mutate(lambda_gc = replace(lambda_gc, lambda_gc > 2, 2)) %>%
    filter(trait_type == 'all') %>%
    ggplot + aes(x = lambda_gc, color = annotation, fill = annotation) +
    geom_density(alpha = 0.5) + theme_classic() +
    # geom_bar(alpha = 0.5, position = position_dodge()) + theme_classic() +
    labs(x = 'Gene Lambda GC', y = 'Density') +
    geom_vline(xintercept = 1, lty = 2) +
    annotation_color_scale + annotation_fill_scale + themes +
    facet_grid(result_type~annotation, scale = 'free', labeller = label_type) +
    geom_text(data = gene_cnt, aes(label = paste('Genes with \n lambda>2: ',as.character(cnt)),
                                         x = 1.5, y = 2.5,  face = 'bold', color = annotation, group = NULL), size = 3)
  if(save_plot){
    png(output_path, height = 5, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
figureS15(save_plot = T, output_path = paste0(output_path, 'figureS15_lambda_dist_by_gene_filtered.png'))

