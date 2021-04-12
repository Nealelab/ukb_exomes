source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/'

gene_skato = load_ukb_file(paste0('sig_cnt_after_gene_300k_', test, '.txt.bgz'))
var_sig = load_ukb_file(paste0('sig_cnt_after_SErm_var_300k_', test, '.txt.bgz'),
                        force_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(),
                                          most_severe_consequence = col_character(), mean_proportion = col_double()))

# 37 Gene Lists
## https://github.com/macarthur-lab/gene_lists.git
setwd('~OneDrive//gene_lists/lists/')
universe = as.data.frame(fread('universe.tsv',quote="", header = F))
gene_universe = gene_sig_after %>% filter(gene_symbol %in% universe[,1])
output_path1 = '~/gene_lists/sampling/'
output_path2 = '~/gene_lists/matched_plots/'
sim_times = 100
save = T
oversample = 1000

files = list.files()
files = files[grepl(".tsv$",files)]
files = files[files != 'universe.tsv']
strs = unlist(strsplit(files,".tsv"))
for(i in 1:length(strs)){
  print(i)
  gene_list =  as.data.frame(fread(files[i],quote="", header = T))
  colnames(gene_list)[1] = 'gene_symbol'
  matched_sum100 = match_gene_set_multi_sample(gene_list = gene_list, id_col = 'gene_symbol', sig_col = 'all_sig_pheno_cnt', oversample = oversample, sim_times = sim_times,save_plot = save, output_path = output_path1, filename = strs[i])
  matched_sum = get_subset_matched_data_summary(gene_universe, subset = gene_list, freq_col = 'CAF', id_col = 'gene_symbol', sig_col = 'all_sig_pheno_cnt', oversample = 1000, ref_label = 'Remainder', sub_label = strs[i])
  save_subset_matched_figure(matched_sum, save_plot = T, paste0(output_path2, strs[i], '.png'))
  data = paste0('matched_sum100_',strs[i])
  write.csv(matched_sum100, paste0(output_path1, data, '.csv'), row.names = F)
  data = paste0('matched_',strs[i])
  write.csv(matched_sum, paste0(output_path2, data, '.csv'), row.names = F)
  print('DONE')
}

# Constrained Genes Matching
match_constrained_genes = function(seed = 1024, save_plot = T, output_path){
  gene_info = load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
  constrained = gene_skato %>%
    merge(gene_info[,c('gene_id', 'gene', 'oe_lof_upper_bin')],.,by.x = c('gene_id','gene'),by.y =c('gene_id','gene_symbol')) %>%
    filter(oe_lof_upper_bin == 0) %>% distinct(gene_id)
  matched_constrained_sum = get_subset_matched_data_summary(gene_skato, subset = constrained, freq_col = 'CAF', id_col = 'gene_id', sig_col = 'all_sig_pheno_cnt',
                                                            oversample = 1000, ref_label = 'Unconstrained', sub_label = 'Constrained', seed = seed)
  save_subset_matched_figure(matched_constrained_sum, save_plot = save_plot, output_path = output_path)
  return(matched_constrained_sum)
}

matched_constrained_sum = match_constrained_genes(save_plot = T, output_path = paste0(output, 'constrained_matched.png'))

# Developmental Delay Genes Matching
match_DD_genes = function(seed = 1024, save_plot = T, output_path){
  gene_DD = load_ukb_file('forKonrad_sig31kDNM_consensus_genes_2021_01_12.txt')
  matched_geneDD_sum = get_subset_matched_data_summary(gene_skato, subset = gene_DD, freq_col = 'CAF', id_col = 'gene_id', sig_col = 'all_sig_pheno_cnt',
                                                       oversample = 1000, ref_label = 'Other', sub_label = 'Developmental Delay Genes', seed = seed)
  if(save_plot){
    save_subset_matched_figure(matched_geneDD_sum, output_path = output_path)
  }
  return(matched_geneDD_sum)
}

matched_geneDD_sum = match_DD_genes(save_plot = T, output_path = paste0(output, 'geneDD_matched.png'))

# CAF/MAF Matching
match_annotation_by_freq = function(seed = 1024, data, freq_col, sig_col, ref_group, interval = 0.01, save_plot = T, output_path){
  ref_group = check_annotation(ref_group)
  if(ref_group == 'pLoF'){
    matched_sum = get_group_matched_data_summary(data, ref_group = 'pLoF', group1 = 'missense|LC', group2 = 'synonymous', group_col = 'annotation', freq_col = freq_col, sig_col = sig_col, interval = interval, seed = seed)
  }else if(ref_group == 'missense|LC'){
    matched_sum = get_group_matched_data_summary(data, ref_group = 'missense|LC', group1 = 'pLoF', group2 = 'synonymous', group_col = 'annotation', freq_col = freq_col, sig_col = sig_col, interval = interval, seed = seed)
  }else if(ref_group == 'synonymous'){
    matched_sum = get_group_matched_data_summary(data, ref_group = 'synonymous', group1 = 'pLoF', group2 = 'missense|LC', group_col = 'annotation', freq_col = freq_col, sig_col = sig_col, interval = interval, seed = seed)
  }else{
    stop("Invalid group name!")
  }
  save_group_matched_figure(matched_sum, save_plot = save_plot, output_path = output_path)
  return(matched_sum)
}

# CAF Matching
matched_gene_sum_lof = match_annotation_by_freq(data = gene_skato, freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt', ref_group = 'pLoF', save_plot = T, output_path = paste0(output, 'gene_matched_to_lof_skato.png'))
matched_gene_sum_mis = match_annotation_by_freq(data = gene_skato, freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt', ref_group = 'missense|LC', save_plot = T, output_path = paste0(output, 'gene_matched_to_mis_skato.png'))
matched_gene_sum_syn = match_annotation_by_freq(data = gene_skato, freq_col = 'CAF', sig_col = 'all_sig_pheno_cnt', ref_group = 'synonymous', save_plot = T, output_path = paste0(output, 'gene_matched_to_syn_skato.png'))

# MAF Matching
matched_var_sum_lof = match_annotation_by_freq(data = var_sig, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt', ref_group = 'pLoF', interval = 0.01, save_plot = T, output_path = paste0(output, 'var_matched_to_lof.png'))
matched_var_sum_mis = match_annotation_by_freq(data = var_sig, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt', ref_group = 'missense|LC', interval = 0.01, save_plot = T, output_path = paste0(output, 'var_matched_to_mis.png'))
matched_var_sum_syn = match_annotation_by_freq(data = var_sig, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt', ref_group = 'synonymous', interval = 0.01, save_plot = T, output_path = paste0(output, 'var_matched_to_syn.png'))

# PolyPhen2 Variants Matching
match_polyphen_variants = function(seed = 1024, ref_group, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt', interval = 0.01, by_freq_interval = FALSE, save_plot = T, output_path){
  var_mis = var_sig %>% filter(annotation=='missense|LC')
  ref_group = tolower(ref_group)
  if(ref_group == 'probably_damaging'){
    matched_sum = get_group_matched_data_summary(var_mis, ref_group = 'probably_damaging', group1 = 'possibly_damaging', group2 = 'benign', group_col = 'polyphen2',
                                                 freq_col = freq_col, sig_col = sig_col, interval = interval, by_freq_interval = by_freq_interval, seed = seed)
  }else if(ref_group == 'possibly_damaging'){
    matched_sum = get_group_matched_data_summary(var_mis, ref_group = 'possibly_damaging', group1 = 'probably_damaging', group2 = 'benign', group_col = 'polyphen2',
                                                 freq_col = freq_col, sig_col = sig_col, interval = interval, by_freq_interval = by_freq_interval, seed = seed)
  }else if(ref_group == 'benign'){
    matched_sum = get_group_matched_data_summary(var_mis, ref_group = 'benign', group1 = 'probably_damaging', group2 = 'possibly_damaging', group_col = 'polyphen2',
                                                 freq_col = freq_col, sig_col = sig_col, interval = interval, by_freq_interval = by_freq_interval, seed = seed)
  }else{
    stop("Invalid group name!")
  }
  if(by_freq_interval){
    save_group_matched_figure2(matched_sum, levels = c('probably_damaging', 'possibly_damaging', 'benign'),
                               labels = c('Probably Damaging', 'Possibly Damaging', 'Benign'), group_col = 'polyphen2', x_lab = 'PolyPhen2', save_plot = save_plot,
                               output_path = output_path) + facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type)
  }else{
    save_group_matched_figure2(matched_sum, levels = c('probably_damaging', 'possibly_damaging', 'benign'),
                               labels = c('Probably Damaging', 'Possibly Damaging', 'Benign'), group_col = 'polyphen2', x_lab = 'PolyPhen2',
                               save_plot = save_plot, output_path = output_path)
  }
  return(matched_sum)
}

matched_polyphen_sum_probably = match_polyphen_variants(ref_group = 'probably_damaging', freq_col = 'AF', interval = 0.01, by_freq_interval = T, save_plot = T,
                                                        output_path = paste0(output, 'PolyPhen2_matched_to_probably.png'))
matched_polyphen_sum_possibly = match_polyphen_variants(ref_group = 'possibly_damaging', freq_col = 'AF', interval = 0.01, by_freq_interval = T, save_plot = T,
                                                        output_path = paste0(output, 'PolyPhen2_matched_to_possibly.png'))
matched_polyphen_sum_benign = match_polyphen_variants(ref_group = 'benign', freq_col = 'AF', interval = 0.01, by_freq_interval = T, save_plot = T,
                                                      output_path = paste0(output, 'PolyPhen2_matched_to_benign.png'))

# ClinVar Variants Matching
match_clinvar_variants = function(seed = 1024, ref_group, freq_col = 'AF', sig_col = 'all_sig_pheno_cnt', interval = 0.01, by_freq_interval = FALSE, save_plot = T, output_path){
  if(ref_group == 'B/LB'){
    matched_sum = get_group_matched_data_summary(var_sig, ref_group = 'B/LB', group1 = 'P/LP', group2 = 'VUS', group_col = 'pathogenicity',
                                                 freq_col = freq_col, sig_col = sig_col, interval = interval, by_freq_interval = by_freq_interval, seed = seed)
  }else if(ref_group == 'P/LP'){
    matched_sum = get_group_matched_data_summary(var_sig, ref_group = 'P/LP', group1 = 'B/LB', group2 = 'VUS', group_col = 'pathogenicity',
                                                 freq_col = freq_col, sig_col = sig_col, interval = interval, by_freq_interval = by_freq_interval, seed = seed)
  }else if(ref_group == 'VUS'){
    matched_sum = get_group_matched_data_summary(var_sig, ref_group = 'VUS', group1 = 'P/LP', group2 = 'B/LB', group_col = 'pathogenicity',
                                                 freq_col = freq_col, sig_col = sig_col, interval = interval, by_freq_interval = by_freq_interval, seed = seed)
  }else{
    stop("Invalid group name!")
  }
  print(matched_sum)
  if(by_freq_interval){
    save_group_matched_figure2(matched_sum, levels = c('P/LP', 'VUS', 'B/LB'),
                               labels = c( 'Pathogenic/Likely Pathogenic', 'Uncertain significance', 'Benign/Likely Benign'),
                               group_col = 'pathogenicity', x_lab = 'ClinVar Pathogenicity', save_plot = save_plot, output_path = output_path) +
      facet_grid(~interval, switch = "x", scales = "free_x", space = "free_x", labeller = label_type)
  }else{
    save_group_matched_figure2(matched_sum, levels = c('P/LP', 'VUS', 'B/LB'),
                               labels = c( 'Pathogenic/Likely Pathogenic', 'Uncertain significance', 'Benign/Likely Benign'),
                               group_col = 'pathogenicity', x_lab = 'ClinVar Pathogenicity', save_plot = save_plot, output_path = output_path)
  }
  return(matched_sum)
}

matched_clinvar_sum_B = match_clinvar_variants(ref_group = 'B/LB', freq_col = 'AF', interval = 0.01, by_freq_interval = T, save_plot = T, output_path = paste0(output, 'clinvar_matched_to_B_by_AF_filtered2.png'))
matched_clinvar_sum_P = match_clinvar_variants(ref_group = 'P/LP', freq_col = 'AF', interval = 0.01, by_freq_interval = T, save_plot = T, output_path = paste0(output, 'clinvar_matched_to_P_by_AF_filtered2.png'))
matched_clinvar_sum_V = match_clinvar_variants(ref_group = 'VUS', freq_col = 'AF', interval = 0.01, by_freq_interval = T, save_plot = T, output_path = paste0(output, 'clinvar_matched_to_V_by_AF_filtered2.png'))

# pext Matching
match_pext_variants = function(seed = 1024, ref_group, save_plot = T, output_path){
  pext = var_sig %>%
  filter(tx_annotation_csq == most_severe_consequence & most_severe_consequence %in% c('splice_acceptor_variant', 'splice_donor_variant') & lof == 'HC') %>%
  mutate(mean_prop_bin=factor(get_mean_prop_interval(mean_proportion_expressed), levels = c('[0, 0.2]','(0.2, 0.8]','(0.8, 1]')))
  acceptor = pext %>% filter(most_severe_consequence == 'splice_acceptor_variant')
  donor = pext %>% filter(most_severe_consequence == 'splice_donor_variant')

  if(ref_group == '[0, 0.2]'){
    matched_acceptor_sum = get_group_matched_data_summary(acceptor, ref_group = '[0, 0.2]', group1 = '(0.2, 0.8]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_donor_sum = get_group_matched_data_summary(donor, ref_group = '[0, 0.2]', group1 = '(0.2, 0.8]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_sum = rbind(matched_acceptor_sum, matched_donor_sum) %>% mutate(splice_variant = rep(c('Acceptor', 'Donor'), each = 3))
  }else if(ref_group == '(0.2, 0.8]'){
    matched_acceptor_sum = get_group_matched_data_summary(acceptor, ref_group = '(0.2, 0.8]', group1 = '[0, 0.2]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_donor_sum = get_group_matched_data_summary(donor, ref_group = '(0.2, 0.8]', group1 = '[0, 0.2]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_sum = rbind(matched_acceptor_sum, matched_donor_sum) %>% mutate(splice_variant = rep(c('Acceptor', 'Donor'), each = 3))
  }else if(ref_group == '(0.8, 1]'){
    matched_acceptor_sum = get_group_matched_data_summary(acceptor, ref_group = '(0.8, 1]', group1 = '[0, 0.2]', group2 = '(0.2, 0.8]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_donor_sum = get_group_matched_data_summary(donor, ref_group = '(0.8, 1]', group1 = '[0, 0.2]', group2 = '(0.2, 0.8]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_sum = rbind(matched_acceptor_sum, matched_donor_sum) %>% mutate(splice_variant = rep(c('Acceptor', 'Donor'), each = 3))
      }else{
    stop("Invalid group name!")
  }
  if(save_plot){
    plt = matched_sum %>%
      ggplot + aes(x = splice_variant, y = prop, ymin = prop-sd, ymax = prop+sd, color = mean_prop_bin, fill = mean_prop_bin) +
      geom_pointrange(stat = "identity", position = position_dodge(width = 0.4), size = 0.4) +
      labs(y = 'Proportion', x = 'Splice Variant')  +
      scale_y_continuous(label = label_percent(accuracy = 2)) +
      scale_x_discrete(labels = c('Acceptor', 'Donor')) +
      scale_color_manual(name = 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
      scale_fill_manual(name = 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
      theme_classic() + themes

    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(plt)
    dev.off()
  }
  return(matched_sum)
}

matched_pext_sum1 = match_pext_variants( ref_group = '[0, 0.2]', save_plot = T, output_path = paste0(output, 'pext_matched_to_1st.png'))
matched_pext_sum2 = match_pext_variants( ref_group = '(0.2, 0.8]', save_plot = T, output_path = paste0(output, 'pext_matched_to_2nd.png'))
matched_pext_sum3 = match_pext_variants( ref_group = '(0.8, 1]', save_plot = T, output_path = paste0(output, 'pext_matched_to_3rd.png'))