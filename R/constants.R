packages = c('dplyr', 'GGally', 'reshape2', 'data.table', 'scales', 'Matrix', 'ggplot2', 'extrafont', 'gridExtra', 'grDevices', 'grid',
              'RCurl', 'trelliscopejs', 'tidyverse', 'Hmisc', 'devtools', 'broom', 'plotly', 'slackr', 'magrittr', 'gapminder', 'readr', 
              'purrr', 'skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar', 'BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges', 'RColorBrewer', 
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

tranche = '500k'
test = 'skato'
path_to = '~/Desktop/ukb_exomes_450k/'
output = '~/Desktop/ukb_exomes_450k/'
output_path = '~/Desktop/ukb_exomes_450k/'



trait_types = c("continuous", "categorical", "icd10")
trait_type_colors = c('#d56f3e', '#43aa8b', '#b594b6')
names(trait_type_colors) = trait_types
trait_type_names = c('Continuous', 'Binary', 'Disease (ICD)')
names(trait_type_names) = trait_types
trait_type_colors['hide'] = 'gray'
trait_color_scale = scale_color_manual(breaks = trait_types, values=trait_type_colors, name='Trait type', labels=trait_type_names)
trait_fill_scale = scale_fill_manual(breaks = trait_types, values=trait_type_colors, name='Trait type', labels=trait_type_names)

annotation_types = c('pLoF', 'missense|LC', 'synonymous')
annotation_names = c('pLoF', 'Missense', 'Synonymous')
result_types = c('skato', 'skat', 'burden')
result_names = c('SKAT-O', 'SKAT', 'Burden Test')
af_types = c('AF:(None, 0.0001]', 'AF:(0.0001, 0.001]', 'AF:(0.001, 0.01]', 'AF:(0.01, 0.1]', 'AF:(0.1, None]')
af_names = c('AF:(', bquote("-\U221E"), ' , 0.0001]', 'AF:(, 2e-05]', 'AF:(2e-05, 0.001]', 'AF:(2e-5, 0.0001]', 'AF:(0.0001, 0.001]', 'AF:(0.001, 0.01]', 'AF:(0.01, 0.1]', paste0('AF:(0.1, ', bquote("\U221E"), ' )'))
caf_types = c('CAF:(None, 0.0001]', 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', 'CAF:(0.1, None]')
caf_names = c(paste0('CAF:[0 , 0.0001]'), 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', paste0('CAF:(0.1, ', bquote("\U221E"), ' )'))
ac_types = c('(0, 5]', '(5, 50]', '(50, 500]', '(500, 5000]', '(5000, 50000]', '(50000, )')
ac_names = c('(0, 5]', '(5, 50]', '(50, 500]', '(500, 5000]', '(5000, 50000]', paste0('(50000, ', bquote("\U221E"), ' )'))
gene_list_names = c('Constrained' = 'Constrained', 'Developmental Delay' = 'Developmental\nDelay', 'all_ad' = 'Autosomal dominant', 'all_ar' = 'Autosomal recessive', 
                    'fda_approved_drug_targets' = 'FDA approved\ndrug targets', 'gwascatalog' = 'GWAS catalog', 'CEGv2_subset_universe' = 'Cell essential', 'NEGv1_subset_universe' = 'Cell non-essential')
icd_names = c('A' = 'Infectious', 'B' = 'Infectious','C' = 'Neoplasms', 'D' = 'Blood/immune', 'E' = 'Endocrine/metabolic', 'F' = 'Mental/behavioral', 'G' = 'Nervous', 'H1' = 'Eye',
              'H2' =  'Ear', 'I' = 'Circulatory', 'J' = 'Respiratory', 'K' = 'Digestive', 'L' = 'Skin/subcutaneous', 'M' = 'Musculoskeletal', 'N' = 'Genitourinary', 'O' = 'Pregnancy',
              'P' = 'Perinatal', 'Q' = 'Congenital', 'R' = 'Symptoms', 'S' = 'Injury/poison', 'T' = 'Injury/poison', 'V' = 'External causes', 'Y' = 'External causes', 'Z' = 'Health Factors')
random_pheno_subset = c('random_1e-04', 'random_0.001', 'random_0.01', 'random_0.05', 'random_0.1', 'random_0.2', 'random_0.5', 'random_continuous')
polyphen2_levels = c('benign', 'possibly_damaging', 'probably_damaging')
polyphen2_labels = c('Benign', 'Possibly Damaging', 'Probably Damaging')

names(ac_names) = ac_types
names(af_names) = af_types
names(caf_names) = caf_types
names(result_names) = result_types
names(trait_type_names) = trait_types
names(annotation_names) = annotation_types

qual_col_pals = brewer.pal.info %>% filter(category == 'qual')
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(2)
manhattan_color = sample(col_vector[14:29], 16, replace = F)
icd_colors = c('A' = '#B2341E', 'B' = '#B2341E','C' = '#DC632C', 'D' = '#AFA148', 'E' = '#808080', 'F' = '#90BEBE', 'G' = '#B07F80', 'H1' = '#841E62',
               'H2' =  '#8EA870', 'I' = '#364C63', 'J' = '#504F1F', 'K' = '#CA7D05', 'L' = '#709DE2', 'M' = '#002FA7', 'N' = '#FBC9BE', 'O' = '#EDB235',
               'P' = '#4C0009', 'Q' = '#C4AE94', 'R' = '#CCE1F4', 'S' = '#97A1CD', 'T' = '#97A1CD', 'V' = '#A96BA6', 'Y' = '#A96BA6', 'Z' = '#4E8717')
annotation_color_scale = scale_colour_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
annotation_fill_scale = scale_fill_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
themes = theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 10, face = 'bold'),
               axis.text = element_text(color = 'Black', size = 8), 
               axis.title = element_text(color = 'Black', size = 10, face = 'bold'), 
               legend.title = element_text(color = 'Black', size = 10, face = 'bold'), 
               legend.text = element_text(color = 'Black', size = 10), 
               legend.position = 'top', legend.box = 'vertical', 
               strip.text = element_text(color = 'Black', size = 10), 
               strip.background = element_rect( color = "black", size=0.5, linetype="solid") )
label_type = labeller(trait_type2 = trait_type_names, annotation = annotation_names, result_type = result_names, 
                      CAF_range = caf_names, AF_range = af_names, ac_type = ac_names, gene_set_name = gene_list_names)
var_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(), most_severe_consequence = col_character(), mean_proportion = col_double())


check_annotation = function(annotation){
  if(tolower(annotation) %in% c('plof', 'lof')){
    annotation = 'pLoF'
  }else if(tolower(annotation) %in% c('synonymous', 'syn')){
    annotation = 'synonymous'
  }else if(tolower(annotation) %in% c('missense', 'mis', 'missense|lc')){
    annotation = 'missense|LC'
  }
  return(annotation)
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

get_freq_interval = function(freq){
  interval = case_when(
    freq <= 1e-4  ~ '[0, 0.0001]',
    freq <= 1e-3  ~ '(0.0001, 0.001]', 
    freq <= 1e-2  ~ '(0.001, 0.01]', 
    freq <= 1e-1  ~ '(0.01, 0.1]', 
    freq > 1e-1  ~ '(0.1, 1]',
  )
  return(interval)
}

get_freq_interval_500k = function(freq){
  interval = case_when(
    freq <= 1e-4  ~ '[0, 0.0001]',
    freq <= 1e-3  ~ '(0.0001, 0.001]',
    freq <= 5e-3  ~ '(0.001, 0.005]',
    freq <= 1e-2  ~ '(0.005, 0.01]',
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


get_ukb_data_url = function() {
  return(paste0('https://storage.googleapis.com/ukbb-exome-public/',tranche,'/'))
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
    mutate(bins = cut(unlist(get(freq_col)), breaks=br, labels = paste(head(br, -1), br[-1], sep=' - ')), 
           bin_labels = cut(unlist(get(freq_col)), breaks = br, labels = 1:(length(br)-1)))
  return(data)
}

get_group_matched_data = function(data, ref_group, group1, group2, group_col, freq_col, interval = 0.01, seed = 1024){
  set.seed(seed)
  match1 = match_group_by_freq_1set(data, group1 = ref_group, group2 = group1, group_col = group_col, freq_col = freq_col, interval = interval)
  match2 = match_group_by_freq_1set(data, group1 = ref_group, group2 = group2, group_col = group_col, freq_col = freq_col, interval = interval)
  matched_data = rbind( data %>% filter(get(group_col) == ref_group), 
                        match1 %>% select(1:ncol(data)), 
                        match2 %>% select(1:ncol(data)))
  return(matched_data)
}

get_group_matched_data_summary = function(data, ref_group, group1, group2, group_col, freq_col, sig_col, interval = 0.01, by_freq_interval = FALSE, seed = 1024){
  matched_data = get_group_matched_data(data, ref_group, group1, group2, group_col, freq_col, interval = interval, seed = seed)
  matched_data = matched_data %>% mutate(interval = get_freq_interval(get(freq_col)))
  matched_data$sig_cnt = matched_data[, sig_col]

  if(by_freq_interval){
    matched_sum = matched_data %>%
      group_by(get(group_col), interval) %>%
      sig_cnt_summary(sig_col)
  }else{
    print(colnames(matched_data))
    print(group_col)
    matched_sum = matched_data %>%
      group_by(get(group_col)) %>%
      sig_cnt_summary(sig_col)
  }
  colnames(matched_sum)[1] = group_col
  return(matched_sum)
}

get_subset_matched_data_summary = function(ref_data, subset, freq_col = 'CAF', id_col, sig_col, oversample, ref_label, sub_label, interval = 0.01, seed = 1024){
  ref_data = as.data.frame(ref_data)
  subset = as.data.frame(subset)
  ref_data = set_freq_bins(ref_data, freq_col, interval)
  ref_data$sig_cnt = ref_data[ , sig_col]
  sub_data = ref_data %>% filter(get(id_col) %in% subset[, id_col])
  bin_sum = ref_data %>% filter(get(id_col) %in% subset[, id_col]) %>% count(bin_labels)
  remain_data = ref_data %>% filter(!(get(id_col) %in% subset[, id_col]))
  modifier = 1
  if(oversample){modifier = if_else(oversample>nrow(sub_data), oversample/nrow(sub_data), 1)}
  set.seed(seed)
  match_lof = remain_data %>% filter(annotation == 'pLoF') %>% match_subset_by_freq_1set(. , bin_sum, modifier = modifier)
  match_mis = remain_data %>% filter(annotation == 'missense|LC') %>% match_subset_by_freq_1set(., bin_sum, modifier = modifier)
  match_syn = remain_data %>% filter(annotation == 'synonymous') %>% match_subset_by_freq_1set(., bin_sum, modifier = modifier)
  ref_matched_data = rbind(match_lof, match_mis, match_syn) %>% mutate(group=ref_label) %>% select(group, annotation, sig_cnt)
  sub_data = ref_data %>% filter(get(id_col) %in% subset[, id_col]) %>% mutate(group=sub_label) %>% select(group, annotation, sig_cnt)
  matched_data = rbind(ref_matched_data, sub_data)
  matched_sum = matched_data %>% group_by(annotation, group) %>% sig_cnt_summary()
  return(matched_sum)
}

get_two_col_correlation_table = function(data, group_col, col1, col2){
  cat('*** Correlation between', col1, ' and', col2, '*** \n')
  cor_table = data %>% group_by(get(group_col)) %>%
          summarise(tidy(cor.test(get(col1), get(col2))))
  colnames(cor_table)[1] = group_col
  return(cor_table)
}

save_group_matched_figure = function(matched_summary, save_plot = F, output_path){
  plt = matched_summary %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot + aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, color = annotation) +
    geom_pointrange(stat = "identity", position = position_dodge(width = 0.4)) +
    labs(y = 'Proportion', x = 'Annotation')  +
    scale_y_continuous(label = label_percent(accuracy = 0.1)) +
    scale_x_discrete(labels = annotation_names) +
    annotation_color_scale + annotation_fill_scale +
    theme_classic() + themes
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_group_matched_figure2 = function(matched_summary, levels, labels, group_col, x_lab, save_plot = F, output_path){
  plt = matched_summary %>%
    mutate(group = factor(get(group_col), levels = levels, labels = labels)) %>%
    ggplot + aes(x = group, y = prop, ymin = prop-sd, ymax = prop+sd, color = group) +
    geom_pointrange(stat = "identity", position = position_dodge(width = 2)) +
    labs(y = 'Proportion', x = x_lab)  +
    scale_color_manual(name = x_lab, values = c("#D95F02", "#7570B3", "#1B9E77")) +
    scale_y_continuous(label = label_percent(accuracy = 1)) +
    scale_x_discrete(labels = levels) +
    theme_classic() + themes +
    facet_grid(~interval) +
    guides(color = guide_legend(nrow = 3) )
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_subset_matched_figure = function(matched_summary, save_plot = F, output_path){
  plt = matched_summary %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot + aes(x = group, y = prop, ymin = prop-sd, ymax = prop+sd, group = group, color = annotation, fill=annotation) +
    geom_pointrange(stat = "identity", position = position_dodge(width = 0.4)) +
    labs(y = 'Proportion', x = NULL)  +
    scale_y_continuous(label = label_percent(accuracy = 0.1)) +
    annotation_color_scale + annotation_fill_scale  +
    facet_grid(~annotation, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) + theme_classic() + themes+
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(), 
          strip.placement = "outside", 
          strip.text = element_text(face='bold'))
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_subset_matched_figure2 = function(matched_summary, matched_test, save_plot = F, output_path){
  plt = matched_summary %>%
    mutate(annotation = factor(annotation, levels = annotation_types)) %>%
    ggplot +
    geom_pointrange(aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, group = group, color = annotation, fill=annotation, pch = group),
                    stat = "identity", position = position_dodge(width = 1)) +
    labs(y = 'Proportion', x = NULL, alpha = NULL)  +
    scale_y_continuous(label = label_percent(accuracy = 1), breaks = c(0,0.05,0.10,0.15)) +
    scale_x_discrete(labels = annotation_names, limits = rev(levels(matched_summary$annotation))) +
    annotation_color_scale + annotation_fill_scale  +
    scale_shape_manual(name = NULL, values = c(1, 16)) +
    facet_wrap(gene_set_name~., nrow = 4, labeller = label_type) + theme_classic() + themes+
    coord_flip(ylim = c(0,0.2)) +
    theme(panel.spacing = unit(1, "lines"),
        axis.text= element_text(size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face='bold', size=14),
        strip.text.y = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.title = element_text(size = 13)) +
    geom_text(data = matched_test, aes(x = 4-as.numeric(annotation), y = pos, label = sig_label), size = 6)
  if(save_plot){
    png(output_path, height = 10, width = 5, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

sig_cnt_summary = function(data, sig_col = 'sig_cnt'){
    summary = data %>%
      summarise(mean = mean(get(sig_col), na.rm = TRUE), 
                prop = sum(get(sig_col) > 0, na.rm = T) / n(), 
                sem = sd(get(sig_col), na.rm = T)/sqrt(n()), 
                sig_cnt = sum(get(sig_col) > 0, na.rm = T), 
                cnt = n()) %>%
      mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt))
  return(summary)
}

print_annotation_wilcoxon_test = function(data, test_col, annt_col){
  lof = data[data[[annt_col]]=='pLoF', ]
  mis = data[data[[annt_col]]=='missense|LC', ]
  syn = data[data[[annt_col]]=='synonymous', ]
  test = wilcox.test(lof[[test_col]], mis[[test_col]], alternative = 'two.sided')
  print(paste0('pLoF vs. Missense: p-value = ', test$p.value, '; test-statistic = ', test$statistic ))
  test = wilcox.test(lof[[test_col]], syn[[test_col]], alternative = 'two.sided')
  print(paste0('pLoF vs. Synonymous: p-value = ', test$p.value, '; test-statistic = ', test$statistic))
  test = wilcox.test(mis[[test_col]], syn[[test_col]], alternative = 'two.sided')
  print(paste0('Missense vs. Synonymous: p-value = ', test$p.value, '; test-statistic = ', test$statistic))
}

match_subset_by_freq = function(ref_data, subset, freq_col = 'CAF', id_col = 'gene_id', sig_col, oversample=1000, interval = 0.01, sim_times = 5, seed = 12345){
  set.seed(seed)
  ref_data = as.data.frame(ref_data)
  ref_data$sig_cnt = ref_data[, sig_col]
  subset = as.data.frame(subset)
  ref_data = set_freq_bins(ref_data, freq_col, interval)
  sub_data = ref_data %>% filter(get(id_col) %in% subset[, id_col])
  bin_sum = ref_data %>% filter(get(id_col) %in% subset[, id_col]) %>% count(bin_labels)
  remain_data = ref_data %>% filter(!(get(id_col) %in% subset[, id_col]))
  modifier = 1
  if(oversample){modifier = if_else(oversample > nrow(sub_data), oversample/nrow(sub_data), 1)}

  sim_data = replicate(sim_times, match_subset_by_freq_1set(remain_data = remain_data, bin_sum = bin_sum, modifier = modifier), simplify = FALSE)
  sim_sum = map_dfr(sim_data, sig_cnt_summary, .id = 'simulation')

  return(list(sim_data, sim_sum))
}

match_subset_by_freq_1set = function(remain_data, bin_sum, modifier = modifier){
  n = nrow(bin_sum)
  ids = vector()
  pend = 0
  for(i in n:1){
    temp = remain_data %>% filter(bin_labels == bin_sum[i, 1])
    size = bin_sum[i, 2] * modifier
    if(nrow(temp) == 0){
      pend = pend + size
      next
    }else{
      id = sample(row.names(temp), size+pend, replace = T)
      ids = c(ids, id)
      pend = 0
    }
  }
  sims = remain_data[ids, ]
 return(sims)
}

match_group_by_freq = function(data, group1, group2, group_col, freq_col = 'AF', sig_col, interval = 0.01, sim_times = 5, seed = 12345){
  set.seed(seed)
  sim_data = replicate(sim_times, match_group_by_freq_1set(data, group1, group2, group_col , freq_col, interval), simplify = FALSE)
  sim_data$sig_cnt = sim_data[ , sig_col]
  sim_sum = map_dfr(sim_data, sig_cnt_summary, .id = 'simulation')

  return(list(sim_data, sim_sum))
}

match_group_by_freq_1set = function(data, group1, group2, group_col = 'annotation', freq_col = 'AF', interval = 0.01){
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

match_gene_set_multi_sample = function(gene_list, id_col = 'gene_symbol', sig_col, oversample = F, sim_times = 1000, seed = 1024, save_plot = F, output_path, filename){
  sub_sum = gene_universe %>%
    filter(gene_symbol %in% gene_list[ , id_col]) %>%
    group_by(annotation) %>%
    sig_cnt_summary(sig_col)
  sub_lof = gene_universe %>%
    filter(annotation == 'pLoF') %>%
    match_subset_by_freq(subset = gene_list, freq_col = 'CAF', id_col = id_col, sig_col = sig_col, oversample = oversample, interval = 0.01, sim_times = sim_times, seed = seed)

  sub_mis = gene_universe %>%
    filter(annotation == 'missense|LC') %>%
    match_subset_by_freq(subset = gene_list, freq_col = 'CAF', id_col = id_col, sig_col = sig_col, oversample = oversample, interval = 0.01, sim_times = sim_times, seed = seed)

  sub_syn = gene_universe %>%
    filter(annotation == 'synonymous') %>%
    match_subset_by_freq(subset = gene_list, freq_col = 'CAF', id_col = id_col, sig_col = sig_col, oversample = oversample, interval = 0.01, sim_times = sim_times, seed = seed)

  matched_sum = rbind(sub_mis[[2]], sub_lof[[2]], sub_syn[[2]])
  matched_sum$annotation = factor(rep(c('missense|LC', 'pLoF', 'synonymous'), each = sim_times), levels = annotation_types)

  if(save_plot == T){
    match_mean = ggplot(matched_sum, aes(x = mean, color = annotation)) +
      geom_density(alpha = 0.5) + theme_classic() + labs(y = 'Density', x = 'Mean') +
      geom_vline(data = sub_sum, aes(xintercept = mean, color = annotation), lty = 2) +
      annotation_color_scale + annotation_fill_scale + themes

    png(paste0(output_path, 'match_mean_', sim_times, filename, '.png'), height = 2.5, width = 4, units = 'in', res = 300)
    print(match_mean)
    dev.off()

    match_prop = ggplot(matched_sum, aes(x = prop, color = annotation)) +
      geom_density(alpha = 0.5 ) + theme_classic() +
      labs(y = 'Density', x = 'Proportion of Genes with 1 + Hit') +
      geom_vline(data = sub_sum, aes(xintercept = prop, color = annotation), lty = 2) +
      scale_x_continuous(label = label_percent(accuracy = 0.1)) +
      annotation_color_scale + annotation_fill_scale + themes

    png(paste0(output_path, 'match_prop_', sim_times, filename, '.png'), height = 2.5, width = 4, units = 'in', res = 300)
    print(match_prop)
    dev.off()
  }
  return(matched_sum)
}

get_var_gene_overlap_count = function(data = var_gene_by_pheno, type = 'pheno', group = 'all', normalize = T, print = F){
  full_set = group
  label = if_else(type== 'pheno', 'pheno_', '')
  if(type == 'pheno'){
    if(full_set == 'all'){full_set = c('categorical', 'continuous', 'icd10')}
    group = if_else(group == 'icd_first_occurrence', 'icd10', group)
    data = data %>%
      mutate(trait_type = if_else(trait_type == 'icd_first_occurrence', 'icd10', trait_type)) %>%
      filter(trait_type %in% full_set)
  }else{
    if(full_set == 'all'){full_set = c('missense|LC', 'pLoF', 'synonymous')}
    data = data %>% filter(annotation %in% full_set)
  }
  n = if_else(normalize, nrow(data), as.integer(1))
  if(print){
    print(paste('Significant associations for', group, if_else(type=='pheno', 'phenotypes', 'groups')))
    print(paste('Both gene and variant: ', sum(data %>% select(get(paste0(label, 'gene_var_sig_cnt'))))/n))
    print(paste('Gene only: ', sum(data %>% select(get(paste0(label, 'gen_sig_cnt'))))/n))
    print(paste('Variant only: ', sum(data %>% select(get(paste0(label, 'var_sig_cnt'))))/n))
    print(paste('Neither gene nor variant: ', sum(data %>% select(get(paste0(label, 'none_sig_cnt'))))/n))
  }
  return(data.frame(
    significant_by_gene = c(T, F, T, F),
    significant_by_variant = c(T, T, F, F),
    value = c(sum(data %>% select(paste0(label, 'gene_var_sig_cnt')))/n,
              sum(data %>% select(paste0(label, 'var_sig_cnt')))/n,
              sum(data %>% select(paste0(label, 'gene_sig_cnt')))/n,
              sum(data %>% select(paste0(label, 'none_sig_cnt')))/n ),
    category = rep(group, 4),
    label = c('Significant by both', 'Significant by variant', 'Significant by gene', 'None')))
}

format_full_lambda_data = function(data){
  data = data  %>%
    pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
    mutate(result_type = str_split(labels, 'lambda_gc_') %>% map_chr(., 2), )%>%
    mutate(result_type = factor(result_type, levels = result_types), 
           trait_type2 = factor(trait_type2, levels=trait_types), )
  return(data)
}

format_sig_cnt_summary_data = function(data, pheno_group = 'all', freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt'){
  data = data %>%
    mutate(interval = get_freq_interval(get(freq_col))) %>%
    group_by(interval, annotation) %>% sig_cnt_summary(sig_col) %>%
    mutate(annotation = factor(annotation, levels = annotation_types), 
           interval = factor(interval,
                             levels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )'),
                             labels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', if_else(freq_col =='CAF', paste0('(0.1, ', bquote("\U221E"), ' )'),'(0.1, 1]'))))
  return(data)
}

format_sig_cnt_summary_data_500k = function(data, pheno_group = 'all', freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt'){
  data = data %>%
    mutate(interval = get_freq_interval_500k(get(freq_col))) %>%
    group_by(interval, annotation) %>% sig_cnt_summary(sig_col) %>%
    mutate(annotation = factor(annotation, levels = annotation_types),
           interval = factor(interval, levels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.005]', '(0.005, 0.01]', '(0.01, 0.1]', '(0.1, )'),
                             labels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.005]', '(0.005, 0.01]', '(0.01, 0.1]', if_else(freq_col =='CAF', paste0('(0.1, ', bquote("\U221E"), ' )'),'(0.1, 1]'))))
  return(data)
}

format_count_by_freq_data = function(data, freq_col = 'CAF'){
  data = data %>%
    mutate(interval = get_freq_interval(get(freq_col))) %>%
    group_by(interval, annotation) %>% summarise(cnt = n()) %>%
    mutate(interval = factor(interval, levels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, 1]'),
                             labels = c('[0, 0.0001]', '(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', if_else(freq_col =='CAF', paste0('(0.1, ', bquote("\U221E"), ' )'),'(0.1, 1]'))),
           annotation = factor(annotation, levels = annotation_types))
    # mutate(interval = factor(interval, levels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', '(0.1, )'),
    #                          labels = c('(0.0001, 0.001]', '(0.001, 0.01]', '(0.01, 0.1]', paste0('(0.1, ', bquote("\U221E"), ' )'))),
    #        annotation = factor(annotation, levels = annotation_types))
  return(data)
}


format_random_pheno_p_data = function(data, test_type = 'SKAT-O'){
  freq_col = if_else(test_type == 'Single-Variant', 'af', 'CAF_range')
  p_col = if_else(test_type == 'Burden Test', 'Pvalue_Burden', 'Pvalue')
  data = data %>%
    filter(is.na(modifier) & phenocode %in% random_pheno_subset) %>%
    group_by(phenocode, get(freq_col)) %>%
    arrange(get(p_col)) %>%
    add_count(phenocode) %>%
    mutate(observed = -log10(get(p_col)), 
           rank = order(get(p_col)), 
           expected = -(log10(rank / (n+1))), 
           phenocode = factor(phenocode, levels = random_pheno_subset), 
           test = test_type)
  return(data)
}

format_variant_call_stats = function(data){
  data = data %>%
  mutate(AF = as.numeric(strsplit(call_stats, '[:,]') %>% map_chr(., 4)),
         AC = as.numeric(strsplit(call_stats, '[:,]') %>% map_chr(., 2)))
  return(data)
}


save_icd_manhattan_figure = function(data, p_filter = 1e-2, width = 10, spacing = 3, sig_level, save_plot = F, output_path){
   icd_labels = sort(deframe(data %>% distinct(icd10)))
   ctr = ((1:length(icd_labels))-1) * (width +spacing) + width/2
   data = data %>%
     mutate(icd_ind = match(icd10, icd_labels)) %>%
     group_by(icd10) %>%
     arrange(chr, as.integer(position)) %>%
     mutate(pos = seq(0, width, length.out = n()) + (icd_ind-1)*(width+spacing),
            icd10 = factor(icd10, levels = names(icd_names)))
  icd_mh_plt = data %>%
    filter(min_p < p_filter) %>%
    ggplot + aes(x = pos, y = -log10(min_p), color = icd10) +
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_level), lty = 2)+
    scale_x_continuous(breaks = ctr, labels = icd_names[icd_labels]) +
    scale_color_manual(values = icd_colors) +
    scale_y_continuous(trans = gwas_loglog_trans(), breaks = loglog_breaks, name = expression(bold(paste('-log'[10], '(', italic(p), ')'))))+
    labs(x = NULL, y = expression(bold(paste('-log'[10], '(', italic(p), ')')))) +
    theme_classic()   + themes +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95, size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(0.8, 0.1, 0.1, 0.1, "cm"), 
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 12, face = 'bold'),
          legend.text = element_text(size = 12), )
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(icd_mh_plt)
    dev.off()
  }
  return(icd_mh_plt)
}

save_prop_by_annt_freq_figure = function(matched_summary, output_path, save_plot = F){
  plt = matched_summary %>%
      ggplot + aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, color = annotation, fill = annotation) +
      geom_pointrange(stat = "identity", position = position_dodge(width = 2)) +
      labs(y = 'Proportion', x = NULL)  +
      scale_y_continuous(label = label_percent(accuracy = 1)) +
      scale_x_discrete(labels = annotation_names) +
      annotation_color_scale + annotation_fill_scale  +
      facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) + theme_classic() + themes+
      theme(panel.spacing = unit(0, "lines"), 
            strip.background = element_blank(), 
            strip.placement = "outside", 
            strip.text = element_text(face = 'bold'), 
            axis.text= element_text(size = 10), 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95) )
  if(save_plot){
    png(output_path, height = 4, width = 7.5, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_prop_by_annt_freq_figure_500k = function(matched_summary, output_path, save_plot = F){
  plt = matched_summary %>%
      filter(interval %in% c('[0, 0.0001]','(0.0001, 0.001]', '(0.001, 0.005]','(0.005, 0.01]')) %>%
      ggplot + aes(x = annotation, y = prop, ymin = prop-sd, ymax = prop+sd, color = annotation, fill = annotation) +
      geom_pointrange(stat = "identity", position = position_dodge(width = 2)) +
      labs(y = 'Proportion', x = NULL)  +
      scale_y_continuous(label = label_percent(accuracy = 1)) +
      scale_x_discrete(labels = annotation_names) +
      annotation_color_scale + annotation_fill_scale  +
      facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type) + theme_classic() + themes+
      theme(panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(face = 'bold'),
            axis.text= element_text(size = 10),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95) )
  if(save_plot){
    png(output_path, height = 3, width = 5, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_count_barplot_figure = function(cnt_data, cnt_type, save_plot = F, output_path){
  plt = cnt_data %>%
    filter(type == cnt_type) %>%
    ggplot + aes(x = filter, y = cnt, alpha = filter) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.5, fill = '#187bcd') +
    scale_y_continuous(label = comma) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    labs(y = paste('Number of', cnt_type), x = NULL, alpha = NULL) +
    annotation_color_scale + annotation_fill_scale + themes +
    geom_text(aes(label = cnt), vjust = -0.3, size = 3, position = position_dodge(width = 1), color = '#187bcd') +
    theme_classic() + themes + theme(legend.position = 'none')+
    theme(plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"), 
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 10), 
          legend.title = element_text(size = 8, face = 'bold'), 
          legend.text = element_text(size = 8), )
  if(save_plot){
    png(output_path, height = 4, width = 4, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_count_by_freq_figure = function(cnt_data, type, save_plot = F, output_path){
  freq = if_else(type == 'Groups', 'CAF', 'AF')
  if(type == 'Groups'){
    cnt_data = cnt_data %>%
      filter(annotation != 'pLoF|Missense')
  }
  plt = cnt_data %>%
  ggplot + aes(x = interval, y = cnt, color = annotation, fill = annotation) +
  geom_bar(stat ='identity', position = 'dodge') +
  scale_y_continuous(label = comma) +
  labs(y = paste('Number of', type), x = paste(freq, 'Interval'))+
  annotation_color_scale + annotation_fill_scale + themes +
  geom_text(aes(label = cnt, color = annotation), vjust = -0.3, size = 2, position = position_dodge(width = 1)) +
  theme_classic() + themes +
  theme(plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"), 
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 30,  vjust = 1, hjust = 0.95),
        legend.title = element_text(size = 8, face = 'bold'), 
        legend.text = element_text(size = 8))
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(plt)
}

save_var_gene_comparison_table = function(filter = T, normalize = T, save_plot = F, output_path){
  if(filter){
    var_gene = load_ukb_file(paste0('var_gene_comparison_by_pheno_filtered_', test, '_3annt_', tranche,'.txt.bgz'), subfolder = 'analysis/')
  }else{
    var_gene = load_ukb_file(paste0('var_gene_comparison_by_pheno_unfiltered_', test,  '_', tranche, '.txt.bgz'), subfolder = 'analysis/')
  }
  var_gene_summary = rbind(get_var_gene_overlap_count(data = var_gene, group = 'all', normalize = normalize, print = F),
                           get_var_gene_overlap_count(data = var_gene, group = 'icd10', normalize = normalize, print = F),
                           get_var_gene_overlap_count(data = var_gene, group = 'categorical', normalize = normalize, print = F),
                           get_var_gene_overlap_count(data = var_gene, group = 'continuous', normalize = normalize, print = F)) %>%
    mutate(x_offset = rep(c(-0.25, 0.25, 0.25, -0.25), each = 4), 
           y_offset = rep(c(0.25, -0.25, 0.25, -0.25), each = 4), 
           value = round(value, digits = 2))

  figure = var_gene_summary %>%
    ggplot + aes(x = (1 - significant_by_gene) + x_offset, 
                 y = significant_by_variant + y_offset, label = value, fill = category) +
    geom_tile(color = 'darkgray') + theme_void() +
    trait_fill_scale + geom_text() +
    geom_hline(yintercept = 0.5, size = 3, color = 'white') +
    geom_vline(xintercept = 0.5, size = 3, color = 'white') +
    annotate('text', x = -0.55, y = 0, hjust = 0.5, vjust = 0, angle = 90, label = 'Not significant by variant') +
    annotate('text', x = -0.55, y = 1, hjust = 0.5, vjust = 0, angle = 90, label = 'Significant by variant') +
    annotate('text', x = 0, y = 1.55, hjust = 0.5, vjust = 0, label = 'Significant by gene') +
    annotate('text', x = 1, y = 1.55, hjust = 0.5, vjust = 0, label = 'Not significant by gene')+
    theme(legend.title = element_text(size = 8, face = 'bold'), 
          legend.text = element_text(size = 8), )
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

pivot_longer_lambda_data = function(data){
  data = data %>% pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
    mutate(result_type = str_split(labels, 'lambda_gc_') %>% map_chr(., 2),
           result_type = factor(result_type,levels = result_types))
  return(data)
}


get_constrained_matched_pheno_group_table <- function(gene_sig_after, gene_info, test, write = FALSE){
  constrained = gene_sig_after %>%
    merge(gene_info[, c('gene_id', 'gene', 'oe_lof_upper_bin')], ., by.x = c('gene_id', 'gene'), by.y =c('gene_id', 'gene_symbol')) %>%
    filter(oe_lof_upper_bin == 0) %>% distinct(gene_id)
  matched_constrained_sum_all = get_subset_matched_data_summary(gene_sig_after, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = 'all_sig_pheno_cnt', oversample = 1000,
                                                            ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Constrained')
  matched_constrained_sum_con = get_subset_matched_data_summary(gene_sig_after, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = paste0('continuous_sig_pheno_cnt_', test), oversample = 1000,
                                                              ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Constrained')
  matched_constrained_sum_cat = get_subset_matched_data_summary(gene_sig_after, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = paste0('categorical_sig_pheno_cnt_', test), oversample = 1000,
                                                              ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Constrained')
  matched_constrained_sum_icd = get_subset_matched_data_summary(gene_sig_after, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = paste0('icd10_sig_pheno_cnt_', test), oversample = 1000,
                                                              ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Constrained')
  matched_constrained_sum <- matched_constrained_sum_all %>%
    mutate(pheno_group = 'all') %>%
    rbind(matched_constrained_sum_con %>% mutate(pheno_group = 'continuous'))%>%
    rbind(matched_constrained_sum_cat %>% mutate(pheno_group = 'categorical'))%>%
    rbind(matched_constrained_sum_icd %>% mutate(pheno_group = 'icd10'))
  if(write){
    write_csv(matched_constrained_sum, '~/ukb_exomes/data/matched_constrained_by_pheno_groups.csv')
  }
  return(matched_constrained_sum)
}