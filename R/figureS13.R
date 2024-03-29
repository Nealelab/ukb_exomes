source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figureS13 = function(save_plot = T, save_pdf = F, output_path){
   lambda_by_gene_before_coverage  = load_ukb_file(paste0('lambda_by_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
  if(tranche == '500k'){
    lambda_by_gene_before_coverage = lambda_by_gene_before_coverage %>%
      filter(annotation != 'pLoF|missense|LC')
  }

  lambda_by_gene_before_coverage = lambda_by_gene_before_coverage %>%
    select(-c(6:8)) %>%
    pivot_longer_lambda_data( ) %>%
    mutate(trait_type = str_split(labels, '_lambda_gc_') %>% map_chr(., 1),
           coverage_int = get_coverage_interval(mean_coverage),
           annotation = factor(annotation, levels = annotation_types),
           coverage_int = factor(coverage_int, levels = c('[0, 10]', '(10, 20]', '(20, 30]', '(30, 40]', '(40, 50]', '(50, )'),
                                 labels = c('[0, 10]', '(10, 20]', '(20, 30]', '(30, 40]', '(40, 50]', paste0('(50,', bquote("\U221E"), ' )')))) %>%
    filter(trait_type != 'icd10')

  figureS13a = lambda_by_gene_before_coverage %>%
    filter(!is.na(coverage_int)) %>%
    filter(result_type == 'skato') %>%
    ggplot + aes(x = coverage_int, y = lambda_gc, color = annotation) +
    geom_boxplot() + geom_hline(aes(yintercept = 1), lty = 2) +
    labs(y = 'Gene Lambda GC', x = '') +
    theme_classic() + ylim(0, 2) +
    annotation_color_scale + annotation_fill_scale + themes +
    theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.8))+
    facet_grid(trait_type~annotation,scale = 'free', labeller = label_type)

  figureS13b = lambda_by_gene_before_coverage %>%
    filter(!is.na(coverage_int)) %>%
    filter(result_type == 'burden') %>%
    ggplot + aes(x = coverage_int, y = lambda_gc, color = annotation) +
    geom_boxplot() + geom_hline(aes(yintercept = 1), lty = 2) +
    labs(y = 'Gene Lambda GC', x = 'Mean Coverage') +
    theme_classic() + ylim(0, 2) +
    annotation_color_scale + annotation_fill_scale +  themes +
    theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.8))+
    facet_grid(trait_type~annotation,scale = 'free', labeller = label_type)
  figure = ggarrange(figureS13a, figureS13b,  labels = c('(A) SKAT-O', '(B) Burden Test'), nrow = 2,
                  label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
    png(output_path, height = 10, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  if(save_pdf){
    ggsave(filename=output_path, plot=figure, device=cairo_pdf,
       width = 174/in2mm, height = 232/in2mm, units = "in", dpi = 300)
    # pdf(output_path, height = 232/in2mm, width = 174/in2mm)
    # print(figure)
    # dev.off()
  }
  return(figure)
}

figureS13(save_plot = T, save_pdf = T, output_path = paste0(path_to, 'final_pdf_figures/figureS13_lambda_by_gene_by_coverage.pdf'))
# figureS13(save_plot = T, output_path = paste0(path_to, 'figureS13_lambda_by_gene_by_coverage.png'))

