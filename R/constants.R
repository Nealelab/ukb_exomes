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
af_types = c('AF:(None, 0.0001]', 'AF:(2e-05, 0.0001]', 'AF:(0.0001, 0.001]', 'AF:(0.001, 0.01]', 'AF:(0.01, 0.1]', 'AF:(0.1, None]')
af_names = c('AF:( , 0.0001]', 'AF:(2e-5, 0.0001]', 'AF:(0.0001, 0.001]', 'AF:(0.001, 0.01]', 'AF:(0.01, 0.1]', 'AF:(0.1, )')
caf_types = c('CAF:(None, 0.0001]', 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', 'CAF:(0.1, None]')
caf_names = c('CAF:( , 0.0001]', 'CAF:(0.0001, 0.001]', 'CAF:(0.001, 0.01]', 'CAF:(0.01, 0.1]', 'CAF:(0.1, )')

names(af_names) = af_types
names(caf_names) = caf_types
names(result_names) = result_types
names(annotation_names) = annotation_types

annotation_color_scale = scale_colour_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
annotation_fill_scale = scale_fill_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
themes = theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 10, face = 'bold'),
               axis.text = element_text(color = 'Black', size = 5),
               axis.line = element_line(size=0.5),
               axis.title = element_text(color = 'Black', size = 7, face = 'bold'),
               legend.title = element_text(color = 'Black', size = 5, face = 'bold'),
               legend.text = element_text(color = 'Black', size = 5),
               legend.position = 'top', legend.box = 'vertical',
               strip.text = element_text(color = 'Black', size = 7),
               strip.background = element_rect( color="black", size=0.5, linetype="solid") )
label_type = labeller(trait_type2=trait_type_names, annotation=annotation_names, result_type=result_names, CAF_range=caf_names, AF_range=af_names)

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

get_matched_data = function(data, annt1='lof', freq_col, interval=0.01){
  data = as.data.frame(data)
  match1 = match_annotation_by_freq_1set(data, annt1 = annt1, annt2 = 'mis', freq_col = freq_col, interval=interval)
  match2 = match_annotation_by_freq_1set(data, annt1 = annt1, annt2 = 'syn', freq_col = freq_col, interval=interval)
  ref_data = data %>% filter(annotation == check_annotation(annt1)) %>% set_freq_bins(freq_col = freq_col, interval = interval)
  matched = rbind(ref_data, match1, match2)
  return(matched)
}

sig_cnt_summary = function(data, annt_col='annotation', sig_col='sig_cnt'){
  if(is.null(annt_col)){
    sums = data %>%
      summarise(mean=mean(get(sig_col), na.rm=TRUE),
                median=median(get(sig_col), na.rm=TRUE),
                prop=sum(get(sig_col)>0, na.rm = T)/sum(get(sig_col)>=0, na.rm=T),
                sig_cnt=sum(get(sig_col)>0, na.rm = T),
                cnt=n())
  }else{
    sums = data %>%
      group_by(get(annt_col)) %>%
      summarise(mean=mean(get(sig_col), na.rm=TRUE),
                median=median(get(sig_col), na.rm=TRUE),
                prop=sum(get(sig_col)>0, na.rm = T)/sum(get(sig_col)>=0, na.rm=T),
                sig_cnt=sum(get(sig_col)>0, na.rm = T),
                cnt=n())
    colnames(sums)[1] = annt_col}
  return(sums)
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

print_freq_sig_cor = function(data=gene_sig, test='skato', freq_col='caf', sig_col='sig_cnt'){
  if( freq_col=='caf'){data = data[data$result_type == test, ]}
  cat('*** Correlation - CAF vs. Association Count (', test, ') *** \n')
  data %>%
    na.omit() %>%
    group_by(annotation) %>%
    summarise(tidy(cor.test(get(freq_col), get(sig_col))))
}

match_subset_by_freq = function(ref_data, subset, freq_col='caf', id_col='gene_id', sig_col='sig_cnt',
                                    interval=0.01, sim_times=5, seed=12345){
  ref_data = as.data.frame(ref_data)
  subset = as.data.frame(subset)

  ref_data = set_freq_bins(ref_data, freq_col, interval)
  sub_data = ref_data %>% filter(get(id_col) %in% subset[, id_col])
  remain_data = ref_data %>% filter(!(get(id_col) %in% subset[, id_col]))
  bin_sum = sub_data %>% count(bin_labels)

  set.seed(seed)
  sim_data = replicate(sim_times, match_subset_by_freq_1set(remain_data, bin_sum, sig_col = sig_col), simplify=FALSE)
  sim_sum = map_dfr(sim_data, sig_cnt_summary, .id = 'simulation')

  return(list(sim_data, sim_sum))
}


match_subset_by_freq_1set = function(remain_data, bin_sum, sig_col='sig_cnt'){
  n = nrow(bin_sum)
  ids = vector()
  pend = 0
  for(i in n:1){
    temp = remain_data %>% filter(bin_labels == bin_sum[i, 1])
    size = bin_sum[i, 2]
    if(nrow(temp)==0){
      pend = pend + size
      next
    }else{
      id = sample(row.names(temp), size+pend, replace = T)
      ids = c(ids, id)
      pend = 0
    }
  }
  sims = remain_data %>%
    filter(row.names(.) %in% ids)
 return(sims)
}

match_annotation_by_freq = function(data, annt1, annt2, annt_col = 'annotation', freq_col='AF', interval=0.0001, sim_times=5, seed=12345){
  set.seed(seed)
  sim_data = replicate(sim_times, match_annotation_by_freq_1set(data, annt1, annt2, annt_col ,freq_col, interval), simplify=FALSE)
  sim_sum = map_dfr(sim_data, sig_cnt_summary, .id = 'simulation', annt_col=annt_col)

  return(list(sim_data, sim_sum))
}

match_annotation_by_freq_1set = function(data, annt1, annt2, annt_col = 'annotation',freq_col='AF', interval=0.001){
  data = set_freq_bins(data = data, freq_col = freq_col, interval = interval)
  bin1 = data %>% filter(get(annt_col) == check_annotation(annt1)) %>% count(bin_labels)
  data2 = data %>% filter(get(annt_col) == check_annotation(annt2))

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