source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
output = '~/Desktop/final_figures/'
test = 'skato'

gene_sig_after = load_ukb_file(paste0('gene_sig_cnt_filtered_', test, '_300k.txt.bgz'), subfolder = 'analysis/')
gene_info = load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', subfolder = 'analysis/')
gene_DD = load_ukb_file('forKonrad_sig31kDNM_consensus_genes_2021_01_12.txt', subfolder = 'analysis/')
## https://github.com/macarthur-lab/gene_lists.git
setwd('~/gene_lists/lists/')
figure4 = function(save_plot = F, output_path){
  universe = as.data.frame(fread('universe.tsv', quote="", header = F))
  gene_universe = gene_sig_after %>% filter(gene_symbol %in% universe[, 1])

  constrained = gene_sig_after %>%
    merge(gene_info[, c('gene_id', 'gene', 'oe_lof_upper_bin')], ., by.x = c('gene_id', 'gene'), by.y =c('gene_id', 'gene_symbol')) %>%
    filter(oe_lof_upper_bin == 0) %>% distinct(gene_id)
  matched_constrained_sum = get_subset_matched_data_summary(gene_sig_after, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = 'all_sig_pheno_cnt', oversample = 1000, 
                                                            ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Constrained')
  matched_gene470_sum = get_subset_matched_data_summary(gene_sig_after, subset = gene_DD, freq_col = 'CAF', id_col = 'gene_id', sig_col = 'all_sig_pheno_cnt', oversample = 1000, 
                                                        ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = 'Developmental Delay')
  matched = rbind(matched_constrained_sum, matched_gene470_sum)

  sub_strs = names(gene_list_names)[-(1:2)]
  for(i in 1: length(sub_strs)){
    gene_list = as.data.frame(fread(paste0(sub_strs[i], '.tsv'), quote="", header = T))
    colnames(gene_list)[1] = 'gene_symbol'
    matched_sum = get_subset_matched_data_summary(gene_universe, subset = gene_list, freq_col = 'CAF', id_col = 'gene_symbol', sig_col = 'all_sig_pheno_cnt', oversample = 1000, 
                                                  ref_label = 'Background', sub_label = 'Test Set') %>% mutate(gene_set_name = sub_strs[i])
    matched = rbind(matched, matched_sum)
  }

  matched= matched %>%
    mutate(annotation = factor(annotation, levels = annotation_types), 
           gene_set_name = factor(gene_set_name, levels = names(gene_list_names)), 
           group = factor(group, levels = c('Background', 'Test Set')))
  matched$panel = c(rep(rep(c(1, 2), each = 6), 4))

  matched_test_fisher = matched %>%
    mutate(no_sig_cnt = cnt - sig_cnt) %>%
    select(sig_cnt, no_sig_cnt, gene_set_name, annotation, panel) %>%
    group_by(annotation, gene_set_name, panel)%>%
    group_modify(~ broom::tidy(fisher.test(as.matrix(.)))) %>%
    merge(.,matched %>% group_by(annotation, gene_set_name) %>% summarise(prop = max(prop)), by = c('annotation', 'gene_set_name')) %>%
    mutate(sig_label = if_else(p.value < 0.001, '*', ''),
           annotation = factor(annotation, levels = annotation_types), ) %>%
    filter(p.value< 0.001)

  figure4_p1 = matched %>% filter(panel == 1) %>% save_subset_matched_figure2(., matched_test_fisher%>% filter(panel == 1))
  figure4_p2 = matched %>% filter(panel == 2) %>% save_subset_matched_figure2(., matched_test_fisher%>% filter(panel == 2))
  figure = ggpubr::ggarrange(figure4_p1, figure4_p2, ncol = 2, common.legend = TRUE)
  if(save_plot){
    png(output_path, height = 10, width = 7.5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure4(save_plot = T, output_path = paste0(output, 'figure4_gene_list_matching_', test, '.png'))


