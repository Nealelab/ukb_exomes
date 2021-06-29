source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
output_path = '~/Desktop/final_figures/'

rp_gene_p_value_skato  = load_ukb_file('rp_gene_p_value_skato_300k.txt.bgz', subfolder = 'analysis/')
rp_gene_p_value_burden = load_ukb_file('rp_gene_p_value_burden_300k.txt.bgz', subfolder = 'analysis/')
rp_var_p_value  = load_ukb_file('rp_var_p_value300k.txt.bgz', subfolder = 'analysis/')

figureS8 = function(save_plot = F, direction = 'long', output_path){
  var_p_value = format_random_pheno_p_data(rp_var_p_value, test_type = 'Single-Variant')
  gene_p_value_skato = format_random_pheno_p_data(rp_gene_p_value_skato, test_type = 'SKAT-O')
  gene_p_value_burden = format_random_pheno_p_data(rp_gene_p_value_burden, test_type = 'Burden Test')
  var_p_value_filtered = var_p_value %>%
    filter(Pvalue <= 0.01 | (Pvalue > 0.01 & n < 10000) |
             (Pvalue > 0.01 & n >= 10000 & n < 100000 & rank %% 50 == 1) |
              (Pvalue > 0.01 & n >= 100000 & rank %% 1000 == 1))

  freq_names = c('( , 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )')
  names(freq_names) = c(caf_types)
  gene_p_value_skato = gene_p_value_skato %>% mutate(af = as.character(factor(CAF_range, labels = freq_names, levels = caf_types)))
  gene_p_value_burden = gene_p_value_burden %>% mutate(af = as.character(factor(CAF_range, labels = freq_names, levels = caf_types)))
  colnames(gene_p_value_burden)[9] = 'Pvalue'
  names(freq_names) = c(af_types)
  var_p_value_filtered2 = var_p_value_filtered %>% mutate(af = as.character(factor(af, labels = freq_names, levels = af_types)))
  p_value = rbind(gene_p_value_skato %>% select(observed, expected, modifier, af, phenocode, test),
                  gene_p_value_burden %>% select(observed, expected, modifier, af, phenocode, test),
                  var_p_value_filtered2%>% select(observed, expected, modifier, af, phenocode, test)) %>%
    mutate(test = factor(test, levels = c('Single-Variant', 'SKAT-O', 'Burden Test')))

  if(direction == 'long'){
    facet_setting = facet_grid(phenocode~test, labeller = label_type, scales = 'free')
    height = 20
    width = 8
  }else{
    facet_setting = facet_grid(test~phenocode, labeller = label_type, scales = 'free')
    height = 8
    width = 20
  }

  figure = p_value %>%
    filter(is.na(modifier)) %>%
    ggplot + aes(y = observed, x = expected, color = af) +
    geom_point(alpha = 0.5, position = 'jitter')+
    geom_abline(intercept = 0, slope = 1)+
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)")+
    scale_color_brewer(name = 'Freq Interval', palette = 'RdYlBu') +
    theme_classic()+ theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 20, face = 'bold'),
                           axis.title = element_text(color = 'Black', size = 18, face = 'bold'),
                           legend.title = element_text(color = 'Black', size = 17, face = 'bold'),
                           legend.text = element_text(color = 'Black', size = 17),
                           legend.position = 'top', legend.box = 'vertical',
                           strip.text = element_text(color = 'Black', size = 17)) +
    coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +
    guides(color = guide_legend(override.aes = list(size = 5, pch = 19, alpha = 1) ) ) +
    facet_setting

  if(save_plot){
    png(output_path, height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figureS8(save_plot = T, output_path = paste0(output_path, 'figureS8_random_pheno_qq_plot.png'))




