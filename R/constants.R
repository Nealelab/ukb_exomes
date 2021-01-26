packages <- c('GGally', 'reshape2',  'data.table', 'scales', 'Matrix', 'ggplot2', 'extrafont', 'gridExtra', 'grDevices',
              'RCurl', 'trelliscopejs', 'tidyverse', 'Hmisc', 'devtools', 'broom', 'plotly', 'slackr', 'magrittr', 'gapminder', 'readr',
              'purrr', 'skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar', 'BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges',
              'ggpubr', 'meta', 'tidygraph', 'pbapply', 'RMySQL', 'egg', 'ggwordcloud', 'patchwork', 'ggrastr', 'ggthemes', 'STRINGdb', 'ggrepel', 'dplyr')

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
names(annotation_names) <- annotation_types
annotation_color_scale <- scale_colour_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
annotation_fill_scale <- scale_fill_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
themes <-   theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 20, face = 'bold'),
                  axis.title = element_text(color = 'Black', size = 18, face = 'bold'),
                  legend.title = element_text(color = 'Black', size = 17, face = 'bold'),
                  legend.text = element_text(color = 'Black', size = 17),
                  legend.position = 'top', legend.box = 'vertical',
                  strip.text = element_text(color = 'Black', size = 17))
label_type <- labeller(trait_type2=trait_type_names, annotation = annotation_names)


get_ukb_data_url <- function() {
  return(paste0('https://storage.googleapis.com/ukbb-exome-public/summary_statistics_analysis/'))
}

get_freq_interval <- function(freq){
  interval = case_when(
    freq <= 1e-4  ~ '[0, 0.0001]',
    freq <= 1e-3 & freq > 1e-4 ~ '(0.0001, 0.001]',
    freq <= 1e-2 & freq > 1e-3 ~ '(0.001, 0.01]',
    freq <= 1e-1 & freq > 1e-2 ~ '(0.01, 0.1]',
    freq > 1e-1  ~ '(0.1, )',
  )
  return(interval)
}

set_freq_bins <- function(data, freq_name='AF', interval=0.001){
  br = seq(0,max(data[,freq_name]),by=interval)
  data$bins = cut(unlist(data[,freq_name]), breaks=br, labels=paste(head(br,-1), br[-1], sep=" - "))
  data$bin_labels = cut(unlist(data[, freq_name]), breaks = br, labels=1:(length(br)-1))
  return(data)
}

check_annotation <- function(annotation){
  if(tolower(annotation)%in%c('plof','lof')){
    annotation='pLoF'
  }else if(tolower(annotation)%in%c('synonymous','syn')){
    annotation='synonymous'
  }else if(tolower(annotation)%in%c('missense','mis','missense|lc')){
    annotation='missense|LC'
  }else{
    stop('Invalid Annotation Type')
  }
  return(annotation)
}

match_annotation_by_freq <- function(data, annt1, annt2, freq_name='AF', interval=0.001, seed=12345){
  d0 = set_freq_bins(data = data, freq_name = freq_name, interval = interval)
  annt1 = check_annotation(annt1)
  annt2 = check_annotation(annt2)
  d1 = d0[d0$annotation == annt1, ]
  d2 = d0[d0$annotation == annt2, ]
  bin1 = d1%>%dplyr::count(bin_labels)
  bin2 = d2%>%dplyr::count(bin_labels)
  n <- nrow(bin1)
  ids <- vector()
  pend <- 0
  set.seed(seed)
  for(i in n:1){
    temp <- d2[d2$bin_labels == unlist(bin1[i,1]),]
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
  sims <- d2[ids,]
 return(list(d1, sims))
}

sig_cnt_summary <- function(data,sig_cnt_col='sig_cnt'){
  sums = data[complete.cases(data),]%>%
    group_by(annotation)%>%
    summarise(means=mean(get(sig_cnt_col),na.rm=TRUE),
              prop=sum(get(sig_cnt_col)>0,na.rm = T)/sum(get(sig_cnt_col)>=0,na.rm=T))
  return(sums)
}

print_annotation_wilcoxon_test <- function(data, test_col, annt_col){
  lof <- data[data[[annt_col]]=='pLoF',]
  mis <- data[data[[annt_col]]=='missense|LC',]
  syn <- data[data[[annt_col]]=='synonymous',]
  test <- wilcox.test(lof[[test_col]], mis[[test_col]], alternative='greater')
  print(paste0('pLoF vs. Missense: p-value = ',test$p.value, '; test-statistic = ', test$statistic ))
  test <- wilcox.test(lof[[test_col]], syn[[test_col]], alternative='greater')
  print(paste0('pLoF vs. Synonymous: p-value = ',test$p.value, '; test-statistic = ', test$statistic))
  test <- wilcox.test(mis[[test_col]], syn[[test_col]], alternative='greater')
  print(paste0('Missense vs. Synonymous: p-value = ',test$p.value, '; test-statistic = ', test$statistic))
}

print_freq_sig_corr <- function(data=gene_sig, test='SKATO', freq_col='caf', sig_col='sig_cnt'){
  if( freq_col=='caf'){data <- data[data$result_type == test, ]}
  pcor <- cor.test(unlist(data[data$annotation == 'pLoF', freq_col]), unlist(data[data$annotation == 'pLoF', sig_col]))
  scor <- cor.test(unlist(data[data$annotation == 'synonymous', freq_col]), unlist(data[data$annotation == 'synonymous', sig_col]))
  mcor <- cor.test(unlist(data[data$annotation == 'missense|LC', freq_col]), unlist(data[data$annotation == 'missense|LC', sig_col]))
  cat('Correlation - CAF vs. Association Count (',test,') \n')
  cat('Missense:', mcor$estimate, '(p-value:', mcor$p.value, ')', '\n')
  cat('pLoF:', pcor$estimate, '(p-value:', pcor$p.value, ')', '\n')
  cat('Synonymous:', scor$estimate, '(p-value:', scor$p.value, ')', '\n')
}