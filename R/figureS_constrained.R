source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

gene_sig_after = load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')
gene_info = load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', subfolder = 'analysis/')
gene_pheno_group = load_ukb_file(paste0('constrained_gene_con_pheno_', test, '_', tranche,'.txt.bgz'), subfolder = 'analysis/')

if(tranche=='500k'){
  gene_sig_after = gene_sig_after %>%
    filter(annotation != 'pLoF|missense|LC')
}

matched = get_constrained_matched_pheno_group_table(gene_sig_after, gene_info, test, write = TRUE)

figure1 = function(save_plot = F, output_path){
  matched = read_csv('~/ukb_exomes/data/matched_constrained_by_pheno_groups.csv')
  matched= matched %>%
    filter(pheno_group != 'binary') %>%
    mutate(annotation = factor(annotation, levels = annotation_types),
           gene_set_name = factor(gene_set_name, levels = names(gene_list_names)),
           group = factor(group, levels = c('Background', 'Test Set')))

  matched_test = matched %>%
    mutate(no_sig_cnt = cnt - sig_cnt) %>%
    select(sig_cnt, no_sig_cnt, gene_set_name, annotation, pheno_group) %>%
    group_by(annotation, gene_set_name, pheno_group)%>%
    group_modify(~ broom::tidy(fisher.test(as.matrix(.)))) %>%
    merge(.,matched %>% group_by(gene_set_name, pheno_group) %>% summarise(pos = max(prop+sd)), by = c('gene_set_name', 'pheno_group')) %>%
    mutate(sig_label = if_else(p.value < 0.001, '**', if_else(p.value < 0.05, '*', '')),
           annotation = factor(annotation, levels = annotation_types), ) %>%
    filter(p.value< 0.05)
  figure = matched %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot +
    geom_pointrange(aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, group = group, color = annotation, fill=annotation, pch = group),
                    stat = "identity", position = position_dodge(width = 1)) +
    labs(y = 'Proportion', x = NULL, alpha = NULL)  +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = annotation_names, limits = rev(levels(matched$annotation))) +
    annotation_color_scale + annotation_fill_scale  +
    scale_shape_manual(name = NULL, values = c(1, 16)) +
    coord_flip() +
    facet_grid(~pheno_group, labeller = label_type, scale = 'free') + theme_classic() + themes+
    theme(panel.spacing = unit(1, "lines"),
          axis.text.y= element_text(size = 14),
          axis.text.x= element_text(size = 14, angle = 45, vjust = 1, hjust = 0.95),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face='bold', size=14),
          strip.text.y = element_blank(),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.title = element_text(size = 13)) +
    geom_text(data = matched_test, aes(x = 4-as.numeric(annotation), y = pos, label = sig_label), size = 6)
  if(save_plot){
    png(output_path, height = 4, width = 8, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure2 = function(save_plot = F, output_path){
  constrained = gene_pheno_group %>%
    merge(gene_info[, c('gene_id', 'gene', 'oe_lof_upper_bin')], ., by.x = c('gene_id', 'gene'), by.y =c('gene_id', 'gene_symbol')) %>%
    filter(oe_lof_upper_bin == 0) %>% distinct(gene_id)
  groups <- unique(gene_pheno_group$pheno_group)
  gene_pheno_group <- gene_pheno_group %>%
    pivot_wider(names_from = 'pheno_group', values_from = 'pheno_group_sig_cnt')

  matched = data.frame()
  for(i in groups){
    matched_constrained_sum = get_subset_matched_data_summary(gene_pheno_group, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = i, oversample = 1000,
                                                                  ref_label = 'Background', sub_label = 'Test Set') %>% mutate(pheno_group = i)
    matched = matched %>% rbind(matched_constrained_sum)
  }
  matched= matched %>%
    mutate(annotation = factor(annotation, levels = annotation_types),
           group = factor(group, levels = c('Background', 'Test Set')))

  matched_test = matched %>%
    mutate(no_sig_cnt = cnt - sig_cnt) %>%
    select(sig_cnt, no_sig_cnt, annotation, pheno_group) %>%
    group_by(annotation, pheno_group)%>%
    group_modify(~ broom::tidy(fisher.test(as.matrix(.)))) %>%
    merge(.,matched %>% group_by(annotation, pheno_group) %>% summarise(pos = max(prop+sd)), by = c('annotation', 'pheno_group')) %>%
    mutate(sig_label = if_else(p.value < 0.001, '**', if_else(p.value < 0.05, '*', '')),
           annotation = factor(annotation, levels = annotation_types), ) %>%
    filter(p.value< 0.05)
  figure = matched %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot +
    geom_pointrange(aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, group = group, color = annotation, fill=annotation, pch = group),
                    stat = "identity", position = position_dodge(width = 1)) +
    labs(y = 'Proportion', x = NULL, alpha = NULL)  +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = annotation_names,
                     limits = rev(levels(matched$annotation))
                     ) +
    annotation_color_scale + annotation_fill_scale  +
    scale_shape_manual(name = NULL, values = c(1, 16)) +
    coord_flip() +
    # facet_grid(~pheno_group, labeller = label_type, scale = 'free') + theme_classic() + themes+
    facet_grid(~pheno_group, labeller = label_type, scale = 'free') + theme_classic() + themes+
    theme(panel.spacing = unit(1, "lines"),
          axis.text.y= element_text(size = 14),
          axis.text.x= element_text(size = 14, angle = 45, vjust = 1, hjust = 0.95),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face='bold', size=12),
          strip.text.y = element_blank(),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.title = element_text(size = 13))+
    geom_text(data = matched_test, aes(x = 4-as.numeric(annotation), y = pos, label = sig_label), size = 6)
  if(save_plot){
    png(output_path, height = 4, width = 8, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure1 = figure1(save_plot = T, paste0(output_path, 'figureS_constrained_1',tranche, '_', test,'.png'))

figure2 = figure2(save_plot = T, paste0(output_path, 'figureS_constrained_2',tranche, '_', test,'.png'))