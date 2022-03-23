source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

var_sig_after = load_ukb_file(paste0('var_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/' ,force_cols = var_cols)

figureS18 = function(save_plot = F, output_path){
  polyphen_sum = var_sig_after %>%
    filter(AF <= 0.1 & AF > 0.0001) %>%
    filter(!is.na(polyphen2) & polyphen2 != 'unknown') %>%
    mutate(interval = factor(get_freq_interval(AF), levels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                             labels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'))) %>%
    group_by(interval, polyphen2) %>%
    sig_cnt_summary('all_sig_pheno_cnt') %>%
      mutate(no_sig_cnt = cnt - sig_cnt)

 annt_sig = polyphen_sum %>%
    select(sig_cnt, no_sig_cnt) %>%
    group_by(interval) %>%
    group_modify(~ broom::tidy(pairwise.prop.test(as.matrix(.)))) %>%
    merge(.,polyphen_sum %>% group_by(interval) %>% summarise(prop = max(prop)), by = 'interval') %>%
    mutate(sig_label = if_else(p.value < 0.001, '**', '*'),
           x  = as.numeric(group1),
           xend = as.numeric(group2),
           y = prop + abs(x + xend-3.15) * 0.015,
           group1 = polyphen2_levels[x],
           group2 = polyphen2_levels[xend]) %>%
    filter(p.value<0.05)

  figure = polyphen_sum %>%
    mutate(group = factor(polyphen2, levels = polyphen2_levels, labels = polyphen2_labels),
           interval = factor(interval, levels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                             labels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'))) %>%
    ggplot +
    geom_pointrange(aes(x = group, y = prop, ymin = prop-sd, ymax = prop+sd, color = group),
                    stat = "identity", position = position_dodge(width = 2)) +
    facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) +
    scale_color_manual(name = 'Polyphen2', values = c("#1B9E77", "#7570B3", "#D95F02")) +
    scale_y_continuous(name = 'Proportion', label = label_percent(accuracy = 1)) +
    scale_x_discrete(name = NULL, labels = polyphen2_labels) +
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

figureS18(save_plot = T, output_path = paste0(output, 'figureS18_polyphen_',tranche, '_', test,'.png'))