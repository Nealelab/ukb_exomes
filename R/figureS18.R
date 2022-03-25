source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

gene_sig_after = load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
var_sig_after = load_ukb_file(paste0('var_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/', force_cols = var_cols)

if(tranche=='500k'){
  gene_sig_after = gene_sig_after %>%
    filter(annotation != 'pLoF|missense|LC')
}

figureS18 = function(save_plot=F, output_path){
  sum = rbind(
    gene_sig_after %>%
    group_by(annotation) %>% sig_cnt_summary('all_sig_pheno_cnt') %>% mutate(type = "SKAT-O"),
    var_sig_after %>%
    group_by(annotation) %>% sig_cnt_summary('all_sig_pheno_cnt') %>% mutate(type = "Single-Variant")
  ) %>%
    mutate(annotation = factor(annotation, levels = annotation_types))

  figure = sum %>%
    ggplot + aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, color = annotation, fill = annotation) +
    geom_pointrange(stat = "identity", position = position_dodge(width = 2)) +
    labs(y = 'Proportion', x = NULL)  +
    scale_y_continuous(label = label_percent(accuracy = 1)) +
    scale_x_discrete(labels = annotation_names) +
    annotation_color_scale + annotation_fill_scale  +
    theme_classic() + themes+
    theme(panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(face = 'bold'),
            axis.text= element_text(size = 10),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95) ) +
    facet_grid(~type)

  if(save_plot){
      png(output_path, height = 3, width = 5, units = 'in', res = 300)
      print(figure)
      dev.off()
  }
  return(figure)
}
figureS18(save_plot = T, paste0(output_path, 'figureS18_proportion_no_freq_matching_',tranche, '_', test,'.png'))