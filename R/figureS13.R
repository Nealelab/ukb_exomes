source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
output_path = '~/Desktop/'

figureS13 = function(save_plot = T, output_path){
  lambda_by_gene_before_coverage = load_ukb_file('lambda_by_gene_before_coverage_filtered_300k.txt.bgz')
  lambda_by_gene_before_coverage = lambda_by_gene_before_coverage %>%
    pivot_longer(cols = contains('_lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
    mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
           result_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 2),
           coverage_int = get_coverage_interval(mean_coverage)) %>%
    mutate(annotation = factor(annotation, levels = annotation_types),
           result_type = factor(result_type, levels = result_types),
           coverage_int = factor(coverage_int, levels = c('[0, 10]', '(10, 20]', '(20, 30]', '(30, 40]', '(40, 50]', '(50, )'))) %>%
    filter(trait_type != 'icd10')
  figureS13a = lambda_by_gene_before_coverage %>%
    filter(!is.na(coverage_int)) %>%
    filter(result_type == 'skato') %>%
    ggplot + aes(x = coverage_int, y = lambda_gc, color = annotation) +
    geom_boxplot() + geom_hline(aes(yintercept = 1), lty = 2) +
    labs(y = 'Gene Lambda GC', x = '') +
    theme_classic() + ylim(0, 2) +
    annotation_color_scale + annotation_fill_scale + themes +
    facet_grid(trait_type~annotation,scale = 'free', labeller = label_type)
  figureS13b = lambda_by_gene_before_coverage %>%
    filter(!is.na(coverage_int)) %>%
    filter(result_type == 'burden') %>%
    ggplot + aes(x = coverage_int, y = lambda_gc, color = annotation) +
    geom_boxplot() + geom_hline(aes(yintercept = 1), lty = 2) +
    labs(y = 'Gene Lambda GC', x = 'Mean Coverage') +
    theme_classic() + ylim(0, 2) +
    annotation_color_scale + annotation_fill_scale +  themes +
    facet_grid(trait_type~annotation,scale = 'free', labeller = label_type)
  figure = ggarrange(figureS13a, figureS13b,  labels = c('(A) SKAT-O', '(B) Burden Test'), nrow = 2,
                  label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
    png(output_path, height = 10, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}
figureS13(save_plot = T, output_path = paste0(output_path, 'figureS13_lambda_by_gene_by_coverage.png'))