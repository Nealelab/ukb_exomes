source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figure3a = function(save_plot=F){
  var_info  = load_ukb_file('var_info300k.txt.bgz')

  var_mis = var_info %>%
    filter(annotation == 'missense|LC' & polyphen2 != 'unknown' & !is.na(polyphen2)) %>%
    group_by(annotation, polyphen2) %>%
    sig_cnt_summary(annt_col = NULL, sig_col = 'sig_pheno_cnt') %>%
    mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt))

  var_non_mis = var_info %>%
    filter(annotation != 'missense|LC') %>%
    sig_cnt_summary(annt_col = 'annotation', sig_col = 'sig_pheno_cnt') %>%
    mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt),
           polyphen2 = '' ) %>%
    select(colnames(var_mis))

  var_sum = rbind(var_mis, var_non_mis) %>%
    mutate(annotation=factor(annotation, levels = annotation_types),
           polyphen2=factor(polyphen2, levels = c('probably_damaging', 'possibly_damaging', 'benign',''),
                            labels = c('Probably Damaging', 'Possibly Damaging', 'Benign', '')))

  plt3a = var_sum %>%
    ggplot + aes(x=annotation, y=prop, ymin=prop-sd, ymax=prop+sd, group=polyphen2, fill=annotation, color = annotation) +
    geom_pointrange(stat="identity", position = position_dodge(1), size=0.4) +
    labs(y = 'Proportion', x = '')  +
    scale_y_continuous(label = label_percent(accuracy=1)) +
    scale_x_discrete(labels = annotation_names) +
    annotation_color_scale + annotation_fill_scale +
    geom_text(data = var_sum, aes(x = annotation, y = prop, label = polyphen2), size = 1.5, vjust = -3, hjust = -0.01, position = position_dodge(1) ) +
    theme_classic() +themes

  if(save_plot){
    pdf('~/Desktop/3a_proportion_var_annotation.pdf', height=2, width=3)
    print(plt3a)
    dev.off()
  }
  return(plt3a)
}

figure3b = function(save_plot=F){
  var_info  = load_ukb_file('var_info300k.txt.bgz')

  blosum_sum = var_info %>%
    filter(annotation == 'missense|LC' & nchar(amino_acids)==3) %>%
    mutate(amino_acid1 = str_split(amino_acids, '/') %>% map_chr(., 1),
           amino_acid2 = str_split(amino_acids, '/') %>% map_chr(., 2),)%>%
    group_by(amino_acid1, amino_acid2) %>%
    sig_cnt_summary(annt_col = NULL, sig_col = 'sig_pheno_cnt')

  data(BLOSUM62)
  blosum_score = as.data.frame(BLOSUM62) %>%
    mutate(amino_acid1=rownames(.))%>%
    pivot_longer(!amino_acid1, names_to = 'amino_acid2', values_to = 'score')

  blosum_info = merge(blosum_sum,blosum_score, by=colnames(blosum_score)[1:2])

  plt3b = blosum_info %>%
    na.omit() %>%
    filter(cnt>=100) %>%
    ggplot + aes(x = score, y = prop) +
    geom_point(position = 'jitter', color = color_mis, alpha=0.5) +
    labs(y = 'Proportion', x = 'BLOSUM62 Score') +
    scale_y_continuous(label=label_percent(accuracy=1))+
    theme_classic() + themes + theme(legend.position = 'none')

  if(save_plot){
    pdf('3b_proportion_blosum.pdf', height=2, width=3)
    print(plt3b)
    dev.off()
  }
  return(plt3b)
}

figure3c = function(save_plot=F){
  var_info  = load_ukb_file('var_info300k.txt.bgz')

  pext = var_info %>%
    filter(tx_annotation_csq == most_severe_consequence &
             most_severe_consequence %in% c('splice_acceptor_variant', 'splice_donor_variant') &
             lof == 'HC') %>%
    mutate(mean_prop_bin=factor(get_mean_prop_interval(mean_proportion), levels = c('[0, 0.2]','(0.2, 0.8]','(0.8, 1]')))%>%
    group_by(mean_prop_bin, most_severe_consequence) %>%
    sig_cnt_summary(annt_col = NULL, sig_col = 'sig_pheno_cnt') %>%
    mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt))

  plt3c = pext %>%
    filter(!is.na(mean_prop_bin)) %>%
    ggplot + aes(x=most_severe_consequence, y=prop, ymin=prop-sd, ymax=prop+sd, color=mean_prop_bin, fill=mean_prop_bin) +
    geom_pointrange(stat="identity", position=position_dodge(width=0.4), size=0.4) +
    labs(y = 'Proportion', x = 'Splice Variant')  +
    scale_y_continuous(label=label_percent(accuracy=1)) +
    scale_x_discrete(labels= c('Acceptor', 'Donor')) +
    scale_color_manual(name= 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
    scale_fill_manual(name= 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
    theme_classic() + themes

  if(save_plot){
    pdf('3c_proportion_pext_splice_var.pdf', height=2, width=3)
    print(plt3c)
    dev.off()
  }
  return(plt3c)
}

# detach(package:plyr)
figure3d = function(save_plot=F){
  gene_info = load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
  gene_sig = load_ukb_file('gene_sig300k.txt.bgz')
  gene_sum = gene_sig %>%
    filter(result_type=='skato') %>%
    merge(gene_info[,c('gene_id', 'gene', 'oe_lof_upper_bin')],.,by.x = c('gene_id','gene'),by.y =c('gene_id','gene_symbol'),all.y = T) %>%
    mutate(oe_lof_upper_decile=if_else(oe_lof_upper_bin==0, 'Constrained', 'Unconstrained'),
           annotation = factor(annotation, levels=annotation_types)) %>%
    group_by(annotation, oe_lof_upper_decile) %>%
    sig_cnt_summary(annt_col = NULL, sig_col = 'sig_cnt') %>%
    mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt))

  plt3d = gene_sum %>%
    filter(!is.na(oe_lof_upper_decile)) %>%
    ggplot + aes(x=oe_lof_upper_decile, y=prop, ymin=prop-sd, ymax=prop+sd,group = oe_lof_upper_decile, color=annotation, fill=annotation) +
    geom_pointrange(stat="identity", position = position_dodge(1),  size=0.4) +
    labs(y = 'Proportion', x = ' ')  +
    scale_y_continuous(label=label_percent(accuracy=1)) +
    annotation_color_scale + annotation_fill_scale  +
    facet_grid(~annotation, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) + theme_classic() + themes+
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_text(size=4))

  if(save_plot){
    pdf('3d_proportion_constrained_var.pdf', height=2, width=3)
    print(plt3d)
    dev.off()
  }
  return(plt3d)
}

figure3 = function(){
  png('figure3.png', height=4, width=6, units = 'in', res=300)
  print(ggarrange(p3a, p3b, p3c, p3d, ncol = 2, nrow=2, labels = c('a', 'b', 'c', 'd'),
                  label.args = list(gp=gpar(font=2, cex=1), vjust=1)))
  dev.off()

  pdf('figure3.pdf', height=4, width=6)
  print(ggarrange(p3a, p3b, p3c, p3d, ncol=2, nrow=2, labels = c('a', 'b', 'c', 'd'),
                  label.args = list(gp=gpar(font=2, cex=1), vjust=1)))
  dev.off()
}
