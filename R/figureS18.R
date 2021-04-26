source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/Desktop/Figures/skato_collection/'

var_sig_after = load_ukb_file(paste0('sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'), force_cols = var_cols)

figureS18 = function(save_plot = F, output_path){
  polyphen_sum = var_sig_after %>%
    filter(AF <= 0.01 & !is.na(polyphen2) & polyphen2 != 'unknown') %>%
    mutate(interval = factor(get_freq_interval(AF), levels = af_int)) %>%
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
           y = prop + (x - xend) * 0.002,
           group1 = polyphen2_levels[x],
           group2 = polyphen2_levels[xend]) %>%
    filter(p.value<0.05)

  figure = polyphen_sum %>%
    mutate(group = factor(polyphen2, levels = polyphen2_levels, labels = polyphen2_labels)) %>%
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

figureS18(save_plot = T, output_path = paste0(output, 'figureS18_polyphen.png'))