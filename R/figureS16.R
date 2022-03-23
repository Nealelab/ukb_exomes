source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
output_path = '~/Desktop/final_figures/'

figureS16 = function(save_plot = F, filter = T, output_path){
  corr = load_ukb_file(paste0('pheno_correlation_before_filter_', tranche,'.txt.bgz'), subfolder = 'qc/')
  count = data.frame(
    # pickle.load(open('pheno_correlated_cnt_list_500k', 'rb'))
    Count = c(1391, 1118, 917, 774, 640, 503, 349, 219, 95), # Generated from https://github.com/Nealelab/ukb_exomes/hail/analysis.py#L30
    Corr = seq(0.1, 0.9, 0.1))

  corr = corr %>% filter(i != j)

  figureS16a = corr %>%
    ggplot + aes(x = entry^2) +
    geom_histogram(binwidth = 0.01, alpha = 0.7) + theme_classic() +
    scale_y_log10(label = comma) + scale_x_continuous(breaks = seq(0, 1, 0.1))+
    labs(y = 'Phenotype pairs', x = expression(bold(paste('Correlation (r'^2, ')')))) +
    theme(axis.title.x = element_text(family = "Arial", face = 'bold'))+
    trait_color_scale + trait_fill_scale + themes
  figureS16b = count %>%
    ggplot + aes(x = factor(Corr), y = Count) +
    geom_bar(stat = "identity", alpha = 0.7) + theme_classic() +
    labs(y = 'Phenotypes removed', x = 'Correlation threshold') +
    geom_text(aes(label = Count), vjust = -0.3, size = 2.5)  + themes
  figure = ggarrange(figureS16a, figureS16b,  labels = c('A', 'B'), nrow = 2,
                     label.args = list(gp = gpar(font = 2, cex = 0.75), vjust = 2))
  if(save_plot){
    png(output_path, height = 5, width = 5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figureS16(save_plot = T, filter = F, output_path = paste0(output_path, 'figureS16_pheno_corr_count.png'))

