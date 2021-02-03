packages <- c('dplyr', 'GGally', 'reshape2',  'data.table', 'scales', 'Matrix', 'ggplot2', 'extrafont', 'gridExtra', 'grDevices',
              'RCurl', 'trelliscopejs', 'tidyverse', 'Hmisc', 'devtools', 'broom', 'plotly', 'slackr', 'magrittr', 'gapminder', 'readr',
              'purrr', 'skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar', 'BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges',
              'ggpubr', 'meta', 'tidygraph', 'pbapply', 'RMySQL', 'egg', 'ggwordcloud', 'patchwork', 'ggrastr', 'ggthemes', 'STRINGdb', 'ggrepel')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
}

# BiocManager::install('STRINGdb')
# devtools::install_github('VPetukhov/ggrastr')
# devtools::install_github('hafen/trelliscopejs')
# devtools::install_github('thomasp85/patchwork')
source('~/ukbb_pan_ancestry/constants.R')

annotation_types <- c('pLoF', 'missense|LC', 'synonymous')
annotation_names <- c('pLoF', 'Missense', 'Synonymous')
result_types <- c('skato', 'skat', 'burden')
result_names <- c('SKAT-O', 'SKAT', 'Burden Test')

names(result_names) <- result_types
names(annotation_names) <- annotation_types
annotation_color_scale <- scale_colour_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
annotation_fill_scale <- scale_fill_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
themes <-   theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 20, face = 'bold'),
                  axis.title = element_text(color = 'Black', size = 18, face = 'bold'),
                  legend.title = element_text(color = 'Black', size = 17, face = 'bold'),
                  legend.text = element_text(color = 'Black', size = 17),
                  legend.position = 'top', legend.box = 'vertical',
                  strip.text = element_text(color = 'Black', size = 17))
label_type <- labeller(trait_type2=trait_type_names, annotation=annotation_names, result_type=result_names)

get_ukb_data_url <- function() {
  return(paste0('https://storage.googleapis.com/ukbb-exome-public/summary_statistics_analysis/'))
}

get_freq_interval <- function(freq){
  interval = case_when(
    freq <= 1e-4  ~ '[0, 0.0001]',
    freq <= 1e-3  ~ '(0.0001, 0.001]',
    freq <= 1e-2  ~ '(0.001, 0.01]',
    freq <= 1e-1  ~ '(0.01, 0.1]',
    freq > 1e-1  ~ '(0.1, )',
  )
  return(interval)
}

set_freq_bins <- function(data, freq_col='AF', interval=0.001){
  br = seq(0,max(data[,freq_col], na.rm = T),by=interval)
  data = data %>%
    mutate(bins = cut(unlist(get(freq_col)), breaks=br, labels=paste(head(br,-1), br[-1], sep=" - ")),
           bin_labels = cut(unlist(get(freq_col)), breaks = br, labels=1:(length(br)-1)))
  return(data)
}

check_annotation <- function(annotation){
  if(tolower(annotation) %in% c('plof','lof')){
    annotation='pLoF'
  }else if(tolower(annotation) %in% c('synonymous','syn')){
    annotation='synonymous'
  }else if(tolower(annotation) %in% c('missense','mis','missense|lc')){
    annotation='missense|LC'
  }else{
    stop('Invalid Annotation Type')
  }
  return(annotation)
}

match_annotation_by_freq <- function(data, annt1, annt2, freq_col='AF', interval=0.001, seed=12345){
  data = set_freq_bins(data = data, freq_col = freq_col, interval = interval)
  data1 = data %>% filter(annotation == check_annotation(annt1))
  data2 = data %>% filter(annotation == check_annotation(annt2))
  bin1 = data1 %>% count(bin_labels)

  n <- nrow(bin1)
  ids <- vector()
  pend <- 0
  set.seed(seed)
  for(i in n:1){
    temp <- data2 %>% filter(bin_labels == unlist(bin1[i,1]))
    size <- unlist(bin1[i,2])
    if(nrow(temp)==0){
      pend = pend + size
      next
    }else{
      id <- sample(row.names(temp),size+pend,replace = T)
      ids <- c(ids,id)
      pend <- 0
    }
  }
  sims <- data2[ids,]
  return(list(data1, sims))
}

get_matched_data <- function(data, freq_col){
  match1 <- match_annotation_by_freq(data,  annt1 = 'mis', annt2 = 'lof', freq_col)
  match2 <- match_annotation_by_freq(data, annt1 = 'mis', annt2 = 'syn', freq_col)
  mis <- match1[[1]]
  lof <- match1[[2]]
  syn <- match2[[2]]
  matched <- rbind(mis, lof, syn)
  return(matched)
}

sig_cnt_summary <- function(data,sig_cnt_col='sig_cnt'){
  sums = data %>%
    na.omit() %>%
    group_by(annotation) %>%
    summarise(means=mean(get(sig_cnt_col),na.rm=TRUE),
              medians=median(get(sig_cnt_col),na.rm=TRUE),
              prop=sum(get(sig_cnt_col)>0,na.rm = T)/sum(get(sig_cnt_col)>=0,na.rm=T))
  return(sums)
}

print_annotation_wilcoxon_test <- function(data, test_col, annt_col){
  lof <- data[data[[annt_col]]=='pLoF',]
  mis <- data[data[[annt_col]]=='missense|LC',]
  syn <- data[data[[annt_col]]=='synonymous',]
  test <- wilcox.test(lof[[test_col]], mis[[test_col]], alternative='two.sided')
  print(paste0('pLoF vs. Missense: p-value = ',test$p.value, '; test-statistic = ', test$statistic ))
  test <- wilcox.test(lof[[test_col]], syn[[test_col]], alternative='two.sided')
  print(paste0('pLoF vs. Synonymous: p-value = ',test$p.value, '; test-statistic = ', test$statistic))
  test <- wilcox.test(mis[[test_col]], syn[[test_col]], alternative='two.sided')
  print(paste0('Missense vs. Synonymous: p-value = ',test$p.value, '; test-statistic = ', test$statistic))
}

print_freq_sig_cor <- function(data=gene_sig, test='SKATO', freq_col='caf', sig_col='sig_cnt'){
  if( freq_col=='caf'){data <- data[data$result_type == test, ]}
  cat('*** Correlation - CAF vs. Association Count (',test,') *** \n')
  data %>%
    na.omit() %>%
    group_by(annotation) %>%
    summarise(tidy(cor.test(get(freq_col),get(sig_col))))
}

match_gene_subset_by_caf <- function(ref_data, gene_subset, freq_col='caf', id_col='gene_id', sig_col='sig_cnt',
                                     interval=0.01, sim_times=100, seed=12345){
  ref_data = as.data.frame(ref_data)
  gene_subset = as.data.frame(gene_subset)

  ref_data = set_freq_bins(ref_data,freq_col, interval)
  sub_data = ref_data %>% filter(get(id_col) %in% gene_subset[,id_col])
  remain_data = ref_data %>% filter(!(get(id_col) %in% gene_subset[,id_col]))
  bin_sum = sub_data %>% count(bin_labels)

  set.seed(seed)
  sim_sum = data.frame()
  data_store = data.frame()
  for(k in 1:sim_times){
    n = nrow(bin_sum)
    ids = vector()
    pend = 0
    for(i in n:1){
      temp = remain_data %>% filter(bin_labels == bin_sum[i,1])
      size = bin_sum[i,2]
      if(nrow(temp)==0){
        pend = pend + size
        next
      }else{
        id = sample(row.names(temp),size+pend,replace = T)
        ids = c(ids,id)
        pend = 0
      }
    }
    sims = remain_data %>%
      filter(row.names(.) %in% ids) %>%
      mutate(sim = k)
    sums = sims %>%
      na.omit() %>%
      summarise(means=mean(get(sig_col),na.rm=TRUE),
                medians=median(get(sig_col),na.rm=TRUE),
                props=sum(get(sig_col)>0,na.rm = T)/sum(get(sig_col)>=0,na.rm=T)) %>%
      mutate(sim = k)
    sim_sum = rbind(sim_sum,sums)
    data_store = rbind(data_store,sims)
  }
  return(list(data_store, sim_sum))
}

match_gene_subset_filter_annotation <- function(annt_type, ref_data, gene_subset, freq_col='caf', id_col='gene_id', sig_col='sig_cnt',
                                     interval=0.01, sim_times=100, seed=12345){
  matched_data <- ref_data %>%
  filter(annotation == annt_type)%>%
  match_gene_subset_by_caf()
}