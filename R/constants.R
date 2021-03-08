packages = c('dplyr', 'GGally', 'reshape2',  'data.table', 'scales', 'Matrix', 'ggplot2', 'extrafont', 'gridExtra', 'grDevices', 'grid',
              'RCurl', 'trelliscopejs', 'tidyverse', 'Hmisc', 'devtools', 'broom', 'plotly', 'slackr', 'magrittr', 'gapminder', 'readr', 
              'purrr', 'skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar', 'BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges', 
              'ggpubr', 'meta', 'tidygraph', 'pbapply', 'RMySQL', 'egg', 'ggwordcloud', 'patchwork', 'ggrastr', 'ggthemes', 'STRINGdb', 'ggrepel', 'LowMACA')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
}

# BiocManager::install('STRINGdb')
# BiocManager::install("LowMACA")
# devtools::install_github('VPetukhov/ggrastr')
# devtools::install_github('hafen/trelliscopejs')
# devtools::install_github('thomasp85/patchwork')
source('~/ukbb_pan_ancestry/constants.R')

annotation_types = c('pLoF', 'missense|LC', 'synonymous')
annotation_names = c('pLoF', 'Missense', 'Synonymous')
result_types = c('skato', 'skat', 'burden')
result_names = c('SKAT-O', 'SKAT', 'Burden Test')
af_types = c('AF:(None, 0.0001]', 'AF:(None, 2e-05]','AF:(2e-05, 0.001]','AF:(2e-05, 0.0001]', 'AF:(0.0001, 0.001]', 'AF:(0.001, 0.01]', 'AF:(0.01, 0.1]', 'AF:(0.1, None]')
af_names = c('AF:( , 0.0001]', 'AF:(, 2e-05]','AF:(2e-05, 0.001]','AF:(2e-5, 0.0001]', 'AF:(0.0001, 0.001]', 'AF:(0.001, 0.01]', 'AF:(0.01, 0.1]', 'AF:(0.1, )')
caf_types = c('CAF:(None, 0.0001]', 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', 'CAF:(0.1, None]')
caf_names = c('CAF:( , 0.0001]', 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', 'CAF:(0.1, )')
ac_types = c('0', '1', '10', '100', '1000', '10000', '100000')
ac_names = c('Expected AC:( , 1]', 'Expected AC:(1, 10]', 'Expected AC:(10, 100]',
             'Expected AC:(100 , 1000]', 'Expected AC:(1000, 10000]', 'Expected AC:(10000, 100000]', 'Expected AC:(100000, )')

names(ac_names) = ac_types
names(af_names) = af_types
names(caf_names) = caf_types
names(result_names) = result_types
names(annotation_names) = annotation_types

annotation_color_scale = scale_colour_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
annotation_fill_scale = scale_fill_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
themes = theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 10, face = 'bold'),
               axis.text = element_text(color = 'Black', size = 6),
              # axis.line = element_line(size=0.5),
               axis.title = element_text(color = 'Black', size = 10, face = 'bold'),
               legend.title = element_text(color = 'Black', size = 10, face = 'bold'),
               legend.text = element_text(color = 'Black', size = 10),
               legend.position = 'top', legend.box = 'vertical',
               strip.text = element_text(color = 'Black', size = 10),
               strip.background = element_rect( color="black", size=0.5, linetype="solid") )
# themes  =  theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 20, face = 'bold'),
#                   axis.title = element_text(color = 'Black', size = 18, face = 'bold'),
#                   legend.title = element_text(color = 'Black', size = 17, face = 'bold'),
#                   legend.text = element_text(color = 'Black', size = 17),
#                   legend.position = 'top', legend.box = 'vertical',
#                   strip.text = element_text(color = 'Black', size = 17))
label_type = labeller(trait_type2=trait_type_names, annotation=annotation_names, result_type=result_names, CAF_range=caf_names, AF_range=af_names, ac_type=ac_names)

get_ukb_data_url = function() {
  return(paste0('https://storage.googleapis.com/ukbb-exome-public/summary_statistics_analysis/'))
}

get_freq_interval = function(freq){
  interval = case_when(
    freq <= 1e-4  ~ '[0, 0.0001]', 
    freq <= 1e-3  ~ '(0.0001, 0.001]', 
    freq <= 1e-2  ~ '(0.001, 0.01]', 
    freq <= 1e-1  ~ '(0.01, 0.1]', 
    freq > 1e-1  ~ '(0.1, )', 
  )
  return(interval)
}

get_mean_prop_interval = function(mean_proportion){
  bin = case_when(
    mean_proportion <= 0.2  ~ '[0, 0.2]',
    mean_proportion <= 0.8  ~ '(0.2, 0.8]',
    mean_proportion <= 1  ~ '(0.8, 1]',
  )
  return(bin)
}

get_coverage_interval = function(coverage){
  interval = case_when(
    coverage <= 10  ~ '[0, 10]',
    coverage <= 20  ~ '(10, 20]',
    coverage <= 30  ~ '(20, 30]',
    coverage <= 40  ~ '(30, 40]',
    coverage <= 50  ~ '(40, 50]',
    coverage > 50  ~ '(50, )',
  )
  return(interval)
}

set_freq_bins = function(data, freq_col='AF', interval=0.001){
  br = seq(0, max(data[, freq_col], na.rm = T), by=interval)
  data = data %>%
    mutate(bins = cut(unlist(get(freq_col)), breaks=br, labels=paste(head(br, -1), br[-1], sep=' - ')), 
           bin_labels = cut(unlist(get(freq_col)), breaks = br, labels=1:(length(br)-1)))
  return(data)
}

check_annotation = function(annotation){
  if(tolower(annotation) %in% c('plof', 'lof')){
    annotation='pLoF'
  }else if(tolower(annotation) %in% c('synonymous', 'syn')){
    annotation='synonymous'
  }else if(tolower(annotation) %in% c('missense', 'mis', 'missense|lc')){
    annotation='missense|LC'
  }
  return(annotation)
}

get_group_matched_data = function(data, ref_group, group1, group2, group_col, freq_col, interval=0.01, seed=1024){
  set.seed(seed)
  match1 = match_group_by_freq_1set(data, group1 = ref_group, group2 = group1, group_col = group_col, freq_col=freq_col, interval = interval)
  match2 = match_group_by_freq_1set(data, group1 = ref_group, group2 = group2, group_col = group_col, freq_col=freq_col, interval = interval)
  matched_data = rbind( data %>% filter(get(group_col) == ref_group),
                        match1 %>% select(1:ncol(data)),
                        match2 %>% select(1:ncol(data)))
  return(matched_data)
}

get_group_matched_data_summary = function(data, ref_group, group1, group2, group_col, freq_col, interval=0.01, seed=1024){
  matched_data = get_group_matched_data(data, ref_group, group1, group2, group_col, freq_col, interval=interval, seed=seed)
  matched_sum = matched_data %>%
    group_by(get(group_col)) %>%
    sig_cnt_summary()
  colnames(matched_sum)[1] = group_col
  return(matched_sum)
}

get_subset_matched_data_summary = function(ref_data, subset, freq_col = 'CAF', id_col, ref_label, sub_label, interval = 0.01,  seed = seed){
  ref_data = as.data.frame(ref_data)
  subset = as.data.frame(subset)
  set.seed(seed)
  match_lof = ref_data %>% filter(annotation == 'pLoF') %>% match_subset_by_freq_1set(subset, freq_col = freq_col, id_col = id_col, interval = interval)
  match_mis = ref_data %>% filter(annotation == 'missense|LC') %>% match_subset_by_freq_1set(subset, freq_col = freq_col, id_col = id_col, interval = interval)
  match_syn = ref_data %>% filter(annotation == 'synonymous') %>% match_subset_by_freq_1set(subset, freq_col = freq_col, id_col = id_col, interval = interval)
  ref_matched_data = rbind(match_lof, match_mis, match_syn) %>% mutate(group=ref_label) %>% select(group, annotation, sig_cnt)
  sub_data = ref_data %>% filter(get(id_col) %in% subset[, id_col])%>% mutate(group=sub_label) %>% select(group, annotation, sig_cnt)
  matched_data = rbind(ref_matched_data, sub_data)
  matched_sum = matched_data %>% group_by(annotation, group) %>% sig_cnt_summary()
  return(matched_sum)
}

save_group_matched_figure = function(matched_summary, output_path){
  plt = matched_summary %>%
    mutate(annotation = factor(annotation, levels=annotation_types)) %>%
    ggplot + aes(x=annotation, y=prop, ymin=prop-sd, ymax=prop+sd, color=annotation) +
    geom_pointrange(stat="identity", position=position_dodge(width=0.4)) +
    labs(y = 'Proportion', x = 'Annotation')  +
    scale_y_continuous(label=label_percent(accuracy=1)) +
    scale_x_discrete(labels = annotation_names) +
    annotation_color_scale + annotation_fill_scale +
    theme_classic() + themes

  png(output_path, height=4, width=6, units = 'in', res=300)
  print(plt)
  dev.off()
}

save_group_matched_figure2 = function(matched_summary, levels, labels, group_col, x_lab,  output_path){
  plt = matched_summary %>%
    mutate(group = factor(get(group_col), levels=levels, labels = labels)) %>%
    ggplot + aes(x=group, y=prop, ymin=prop-sd, ymax=prop+sd, color=group) +
    geom_pointrange(stat="identity", position=position_dodge(width=0.4)) +
    labs(y = 'Proportion', x = x_lab)  +
    scale_color_manual(name = x_lab, values = c("#D95F02", "#7570B3", "#1B9E77")) +
    scale_y_continuous(label=label_percent(accuracy=1)) +
    theme_classic() + themes +
    guides(color = guide_legend(nrow=3) )

  png(output_path, height=4, width=6, units = 'in', res=300)
  print(plt)
  dev.off()
}

save_subset_matched_figure = function(matched_summary, output_path){
  plt = matched_summary %>%
    mutate(annotation = factor(annotation,levels = annotation_types)) %>%
    ggplot + aes(x=group, y=prop, ymin=prop-sd, ymax=prop+sd,group = group, color=annotation, fill=annotation) +
    geom_pointrange(stat="identity", position=position_dodge(width=0.4)) +
    labs(y = 'Proportion', x = ' ')  +
    scale_y_continuous(label=label_percent(accuracy=1)) +
    annotation_color_scale + annotation_fill_scale  +
    facet_grid(~annotation, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) + theme_classic() + themes+
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face='bold'))

  png(output_path, height=4, width=6, units = 'in', res=300)
  print(plt)
  dev.off()
}

sig_cnt_summary = function(data, sig_col='sig_cnt'){
    summary = data %>%
      summarise(mean = mean(get(sig_col), na.rm = TRUE),
                prop = sum(get(sig_col) > 0, na.rm = T) / n(),
                sig_cnt = sum(get(sig_col) > 0, na.rm = T),
                cnt=n()) %>%
      mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt))
  return(summary)
}

print_annotation_wilcoxon_test = function(data, test_col, annt_col){
  lof = data[data[[annt_col]]=='pLoF', ]
  mis = data[data[[annt_col]]=='missense|LC', ]
  syn = data[data[[annt_col]]=='synonymous', ]
  test = wilcox.test(lof[[test_col]], mis[[test_col]], alternative='two.sided')
  print(paste0('pLoF vs. Missense: p-value = ', test$p.value, '; test-statistic = ', test$statistic ))
  test = wilcox.test(lof[[test_col]], syn[[test_col]], alternative='two.sided')
  print(paste0('pLoF vs. Synonymous: p-value = ', test$p.value, '; test-statistic = ', test$statistic))
  test = wilcox.test(mis[[test_col]], syn[[test_col]], alternative='two.sided')
  print(paste0('Missense vs. Synonymous: p-value = ', test$p.value, '; test-statistic = ', test$statistic))
}

print_freq_sig_cor = function(data=gene_sig, test='skato', freq_col='CAF', sig_col='sig_cnt'){
  if( freq_col=='caf'){data = data[data$result_type == test, ]}
  cat('*** Correlation - CAF vs. Association Count (', test, ') *** \n')
  data %>%
    na.omit() %>%
    group_by(annotation) %>%
    summarise(tidy(cor.test(get(freq_col), get(sig_col))))
}

match_subset_by_freq = function(ref_data, subset, freq_col='CAF', id_col='gene_id',  interval=0.01, sim_times=5, seed=12345){
  set.seed(seed)
  sim_data = replicate(sim_times, match_subset_by_freq_1set(ref_data, subset, freq_col, id_col, interval), simplify=FALSE)
  sim_sum = map_dfr(sim_data, sig_cnt_summary, .id = 'simulation')

  return(list(sim_data, sim_sum))
}

match_subset_by_freq_1set = function(ref_data, subset, freq_col='CAF', id_col='gene_id', interval=0.01){
  ref_data = as.data.frame(ref_data)
  subset = as.data.frame(subset)
  ref_data = set_freq_bins(ref_data, freq_col, interval)
  sub_data = ref_data %>% filter(get(id_col) %in% subset[, id_col])
  bin_sum = ref_data %>% filter(get(id_col) %in% subset[, id_col]) %>% count(bin_labels)
  remain_data = ref_data %>% filter(!(get(id_col) %in% subset[, id_col]))
  modifier = if_else(nrow(sub_data)<10000, 10000/nrow(sub_data), 1)

  n = nrow(bin_sum)
  ids = vector()
  pend = 0
  for(i in n:1){
    temp = remain_data %>% filter(bin_labels == bin_sum[i, 1])
    size = bin_sum[i, 2] * modifier
    if(nrow(temp)==0){
      pend = pend + size
      next
    }else{
      id = sample(row.names(temp), size+pend, replace = T)
      ids = c(ids, id)
      pend = 0
    }
  }
  # sims = remain_data %>%
  #   filter(row.names(.) %in% ids)
  sims = remain_data[ids,]
 return(sims)
}

match_group_by_freq = function(data, group1, group2, group_col, freq_col='AF', interval=0.01, sim_times=5, seed=12345){
  set.seed(seed)
  sim_data = replicate(sim_times, match_group_by_freq_1set(data, group1, group2, group_col ,freq_col, interval), simplify=FALSE)
  sim_sum = map_dfr(sim_data, sig_cnt_summary, .id = 'simulation')

  return(list(sim_data, sim_sum))
}

match_group_by_freq_1set = function(data, group1, group2, group_col = 'annotation',freq_col='AF', interval=0.01){
  data = set_freq_bins(data = data, freq_col = freq_col, interval = interval)
  bin1 = data %>% filter(get(group_col) == group1) %>% count(bin_labels)
  data2 = data %>% filter(get(group_col) == group2)

  n = nrow(bin1)
  ids = vector()
  pend = 0
  for(i in n:1){
    temp = data2 %>% filter(bin_labels == unlist(bin1[i, 1]))
    size = unlist(bin1[i, 2])
    if(nrow(temp)==0){
      pend = pend + size
      next
    }else{
      id = sample(row.names(temp), size+pend, replace = T)
      ids = c(ids, id)
      pend = 0
    }
  }
  sims = data2[ids, ]
  return(sims)
}

get_var_gene_overlap_count = function(data = var_gene_by_pheno, pheno_group = 'all', normalize = TRUE){
  print(paste('Significant associations for', pheno_group, 'Phenotypes'))
  if(pheno_group == 'all'){pheno_group = c('categorical', 'continuous')}
  data = data %>%
    mutate(trait_type = if_else(trait_type %in% c('icd10', 'icd_first_occurrence'), 'categorical', trait_type)) %>%
    filter(trait_type %in% pheno_group)
  n = if_else(normalize, nrow(data), as.integer(1))
  print(paste('Both gene and variant: ', sum(data %>% select(pheno_gene_var_sig_cnt))/n))
  print(paste('Gene only: ', sum(data %>% select(pheno_gene_sig_cnt))/n))
  print(paste('Variant only: ', sum(data %>% select(pheno_var_sig_cnt))/n))
  print(paste('Neither gene nor variant: ', sum(data %>% select(pheno_none_sig_cnt))/n))
}

format_gene_sig_data = function(gene_sig){
  gene_sig = gene_sig %>%
    mutate(interval = get_freq_interval(CAF),
           annotation = factor(annotation, levels=annotation_types),
           result_type = factor(result_type, levels=result_types))
  return(gene_sig)
}

