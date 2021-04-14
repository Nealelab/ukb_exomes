source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/Desktop/Figures/skato_collection/'
var_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(), most_severe_consequence = col_character(), mean_proportion = col_double())



var_sig_after = load_ukb_file(paste0('sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'), force_cols = var_cols)

figureS18 = function(save_plot = F, output_path){
  polyphen_sum = var_sig_after %>%
    filter(AF <= 0.01 & !is.na(polyphen2) & polyphen2 != 'unknown') %>%
    mutate(interval = factor(get_freq_interval(AF), levels = af_int)) %>%
    group_by(interval, polyphen2) %>%
    sig_cnt_summary('all_sig_pheno_cnt') %>%
      mutate(no_sig_cnt = cnt - sig_cnt)

  prop_test_pvalue = vector()
  for(i in 1:3){
    prop_test_table = as.data.frame(polyphen_sum) %>%
      filter(interval == af_int[i])
    rownames(prop_test_table) = prop_test_table[, 'polyphen2']
    prop_test_table = prop_test_table %>%
    select(sig_cnt, no_sig_cnt) %>%
    as.matrix()
    test = pairwise.prop.test(prop_test_table)
    temp_p = c(test$p.value[1], test$p.value[2], test$p.value[4])
    prop_test_pvalue = c(prop_test_pvalue, temp_p)
  }

  annt_sig = data.frame(interval = polyphen_sum[, 'interval'],
                        group1 = rep(c('benign', 'benign', 'possibly_damaging'), 3),
                        group2 = rep(c('possibly_damaging', 'probably_damaging', 'probably_damaging'), 3),
                        pvalue = prop_test_pvalue,
                        prop = rep(as.data.frame(polyphen_sum %>% group_by(interval) %>% summarise(max(prop)))[,2], each=3)) %>%
    filter(pvalue < 0.05) %>%
    mutate(sig_label = if_else(pvalue < 0.001, '**', '*'),
           x  = match(group1, polyphen2_levels),
           xend = match(group2, polyphen2_levels),
           y = prop + (x - xend) * 0.002)

  figure = polyphen_sum %>%
    mutate(group = factor(polyphen2, levels = polyphen2_levels, labels = polyphen2_labels)) %>%
    ggplot +
    geom_pointrange(aes(x = group, y = prop, ymin = prop-sd, ymax = prop+sd, color = group),
                    stat = "identity", position = position_dodge(width = 2)) +
    facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) +
    scale_color_manual(name = 'Polyphen2', values = c("#D95F02", "#7570B3", "#1B9E77")) +
    scale_y_continuous(name = 'Proportion', label = label_percent(accuracy = 1)) +
    scale_x_discrete(name = NULL, labels = c('Probably Damaging', 'Possibly Damaging', 'Benign')) +
    theme_classic() + themes +
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = 'bold'),
          axis.text= element_text(size = 9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95) ) +
    geom_segment(data = annt_sig, size = .6, aes(x = x, y = y, xend = xend, yend = y)) +
    geom_text(data = annt_sig, aes(x = (x+xend)/2, y = y, label = sig_label), size = 6)

  if(save_plot){
    png(output_path, height = 4, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figureS18(save_plot = T, output_path = paste0(output, 'figureS18_polyphen.png'))