source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
gene_sig  = load_ukb_file('sig_cnt_filtered_gene_300k.txt.bgz')
var_sig = load_ukb_file('sig_cnt_filtered_var_300k.txt.bgz',
                        force_cols = cols(pathogenicity = col_character(), polyphen2 = col_character(),
                                          most_severe_consequence = col_character(),mean_proportion = col_double()))
gene_skato = gene_sig %>% filter(result_type == 'skato')
colnames(var_sig)[8] = 'sig_cnt'

# Constrained Genes Matching
match_constrained_genes = function(seed = 1024, save_plot = F, output_path){
  gene_info = load_ukb_file('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
  constrained = gene_skato %>%
    merge(gene_info[,c('gene_id', 'gene', 'oe_lof_upper_bin')],.,by.x = c('gene_id','gene'),by.y =c('gene_id','gene_symbol')) %>%
    filter(oe_lof_upper_bin==0) %>% distinct(gene_id)
  matched_constrained_sum = get_subset_matched_data_summary(gene_skato, subset = constrained, freq_col = 'CAF', id_col = 'gene_id',
                                                            ref_label = 'Unconstrained', sub_label = 'Constrained', seed= seed)
  if(save_plot){
    save_subset_matched_figure(matched_constrained_sum, output_path = output_path)
  }
  return(matched_constrained_sum)
}

matched_constrained_sum = match_constrained_genes(save_plot = F, output_path = '~/Desktop/Figures/freq_matching/constrained/constrained_matched.png')

# 470 Developmental Delay Genes Matching
match_470genes = function(seed = 1024, save_plot = F, output_path){
  gene470 = load_ukb_file('forKonrad_sig31kDNM_consensus_genes_2021_01_12.txt')
  matched_gene470_sum = get_subset_matched_data_summary(gene_skato, subset = gene470, freq_col = 'CAF', id_col = 'gene_id',
                                                        ref_label = 'Other', sub_label = 'Developmental Delay Genes', seed = seed)
  if(save_plot){
    save_subset_matched_figure(matched_gene470_sum, output_path = output_path)
  }
  return(matched_gene470_sum)
}

matched_gene470_sum = match_470genes(save_plot = F, output_path = '~/Desktop/Figures/freq_matching/gene470/gene470_matched.png')

# CAF/MAF Matching
match_annotation_by_freq = function(seed = 1024, data, freq_col, ref_group, save_plot = F, output_path){
  ref_group = check_annotation(ref_group)
  if(ref_group == 'pLoF'){
    matched_sum = get_group_matched_data_summary(data, ref_group = 'pLoF', group1 = 'missense|LC', group2 = 'synonymous', group_col = 'annotation', freq_col = freq_col, seed = seed)
  }else if(ref_group == 'missense|LC'){
    matched_sum = get_group_matched_data_summary(data, ref_group = 'missense|LC', group1 = 'pLoF', group2 = 'synonymous', group_col = 'annotation', freq_col = freq_col, seed = seed)
  }else if(ref_group == 'synonymous'){
    matched_sum = get_group_matched_data_summary(data, ref_group = 'synonymous', group1 = 'pLoF', group2 = 'missense|LC', group_col = 'annotation', freq_col = freq_col, seed = seed)
  }else{
    stop("Invalid group name!")
  }
  if(save_plot){
    save_group_matched_figure(matched_sum,output_path = output_path)
  }
  return(matched_sum)
}

# CAF Matching
matched_gene_sum_lof = match_annotation_by_freq(data = gene_skato, freq_col = 'CAF', ref_group = 'pLoF', save_plot = T, output_path ='~/Desktop/Figures/freq_matching/CAF_matching/gene_matched_to_lof_skato.png')
matched_gene_sum_mis = match_annotation_by_freq(data = gene_skato, freq_col = 'CAF', ref_group = 'missense|LC', save_plot = T, output_path ='~/Desktop/Figures/freq_matching/CAF_matching/gene_matched_to_mis_skato.png')
matched_gene_sum_syn = match_annotation_by_freq(data = gene_skato, freq_col = 'CAF', ref_group = 'synonymous', save_plot = T, output_path ='~/Desktop/Figures/freq_matching/CAF_matching/gene_matched_to_syn_skato.png')

# MAF Matching
matched_var_sum_lof = match_annotation_by_freq(data = var_sig, freq_col = 'AF', ref_group = 'pLoF', save_plot = F, output_path = '~/Desktop/Figures/freq_matching/MAF_matching/var_matched_to_lof.png')
matched_var_sum_mis = match_annotation_by_freq(data = var_sig, freq_col = 'AF', ref_group = 'missense|LC', save_plot = F, output_path = '~/Desktop/Figures/freq_matching/MAF_matching/var_matched_to_mis.png')
matched_var_sum_syn = match_annotation_by_freq(data = var_sig, freq_col = 'AF', ref_group = 'synonymous', save_plot = F, output_path = '~/Desktop/Figures/freq_matching/MAF_matching/var_matched_to_syn.png')

# PolyPhen2 Variants Matching
match_polyphen_variants = function(seed = 1024, ref_group, save_plot = F, output_path){
  var_mis = var_sig %>% filter(annotation=='missense|LC')
  ref_group = tolower(ref_group)
  if(ref_group == 'probably_damaging'){
    matched_sum = get_group_matched_data_summary(var_mis, ref_group = 'probably_damaging', group1 = 'possibly_damaging', group2 = 'benign', group_col = 'polyphen2', freq_col = 'AF', seed = seed)
  }else if(ref_group == 'possibly_damaging'){
    matched_sum = get_group_matched_data_summary(var_mis, ref_group = 'possibly_damaging', group1 = 'probably_damaging', group2 = 'benign', group_col = 'polyphen2', freq_col = 'AF', seed = seed)
  }else if(ref_group == 'benign'){
    matched_sum = get_group_matched_data_summary(var_mis, ref_group = 'benign', group1 = 'probably_damaging', group2 = 'possibly_damaging', group_col = 'polyphen2', freq_col = 'AF', seed = seed)
  }else{
    stop("Invalid group name!")
  }
  if(save_plot){
    save_group_matched_figure2(matched_sum, levels = c('probably_damaging', 'possibly_damaging', 'benign'),
                               labels = c('Probably Damaging', 'Possibly Damaging', 'Benign'), group_col = 'polyphen2', x_lab = 'PolyPhen2',
                               output_path = output_path)
  }
  return(matched_sum)
}

matched_polyphen_sum_probably = match_polyphen_variants(ref_group = 'probably_damaging', save_plot = F,
                                                        output_path = '~/Desktop/Figures/freq_matching/polyphen2/PolyPhen2_matched_to_probably.png')
matched_polyphen_sum_possibly = match_polyphen_variants(ref_group = 'possibly_damaging', save_plot = F,
                                                        output_path= '~/Desktop/Figures/freq_matching/polyphen2/PolyPhen2_matched_to_possibly.png')
matched_polyphen_sum_benign = match_polyphen_variants(ref_group = 'benign', save_plot = F,
                                                      output_path = '~/Desktop/Figures/freq_matching/polyphen2/PolyPhen2_matched_to_benign.png')

# ClinVar Variants Matching
match_clinvar_variants = function(seed = 1024, ref_group, save_plot = F, output_path){
  if(ref_group == 'B/LB'){
    matched_sum = get_group_matched_data_summary(var_sig, ref_group = 'B/LB', group1 = 'P/LP', group2 = 'VUS', group_col = 'pathogenicity', freq_col = 'AF', seed = seed)
  }else if(ref_group == 'P/LP'){
    matched_sum = get_group_matched_data_summary(var_sig, ref_group = 'P/LP', group1 = 'B/LB', group2 = 'VUS', group_col = 'pathogenicity', freq_col = 'AF', seed = seed)
  }else if(ref_group == 'VUS'){
    matched_sum = get_group_matched_data_summary(var_sig, ref_group = 'VUS', group1 = 'P/LP', group2 = 'B/LB', group_col = 'pathogenicity', freq_col = 'AF', seed = seed)
  }else{
    stop("Invalid group name!")
  }
  if(save_plot){
    save_freq_matched_figure2(matched_clinvar_sum_B, levels = c('P/LP', 'VUS', 'B/LB'),
                              labels = c( 'Pathogenic/Likely Pathogenic', 'Uncertain significance', 'Benign/Likely Benign'),
                              group_col = 'pathogenicity', x_lab = 'ClinVar Pathogenicity', output_path = output_path)
  }
  return(matched_sum)
}
matched_clinvar_sum_B = match_clinvar_variants(ref_group = 'B/LB', save_plot = F, output_path = '~/Desktop/Figures/freq_matching/clinvar/clinvar_matched_to_B.png')
matched_clinvar_sum_P = match_clinvar_variants(ref_group = 'P/LP', save_plot = F, output_path = '~/Desktop/Figures/freq_matching/clinvar/clinvar_matched_to_P.png')
matched_clinvar_sum_V = match_clinvar_variants(ref_group = 'VUS', save_plot = F, output_path = '~/Desktop/Figures/freq_matching/clinvar/clinvar_matched_to_V.png')

########################### pext ##################################
match_pext_variants = function(seed = 1024, ref_group, save_plot = F, output_path){
  pext = var_sig %>%
  filter(tx_annotation_csq == most_severe_consequence & most_severe_consequence %in% c('splice_acceptor_variant', 'splice_donor_variant') & lof == 'HC') %>%
  mutate(mean_prop_bin=factor(get_mean_prop_interval(mean_proportion_expressed), levels = c('[0, 0.2]','(0.2, 0.8]','(0.8, 1]')))
  acceptor = pext %>% filter(most_severe_consequence == 'splice_acceptor_variant')
  donor = pext %>% filter(most_severe_consequence == 'splice_donor_variant')

  if(ref_group == '[0, 0.2]'){
    matched_acceptor_sum = get_group_matched_data_summary(acceptor, ref_group = '[0, 0.2]', group1 = '(0.2, 0.8]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_donor_sum = get_group_matched_data_summary(donor, ref_group = '[0, 0.2]', group1 = '(0.2, 0.8]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_sum = rbind(matched_acceptor_sum, matched_donor_sum) %>% mutate(splice_variant = rep(c('Acceptor', 'Donor'), each=3))
  }else if(ref_group == '(0.2, 0.8]'){
    matched_acceptor_sum = get_group_matched_data_summary(acceptor, ref_group = '(0.2, 0.8]', group1 = '[0, 0.2]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_donor_sum = get_group_matched_data_summary(donor, ref_group = '(0.2, 0.8]', group1 = '[0, 0.2]', group2 = '(0.8, 1]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_sum = rbind(matched_acceptor_sum, matched_donor_sum) %>% mutate(splice_variant = rep(c('Acceptor', 'Donor'), each=3))
  }else if(ref_group == '(0.8, 1]'){
    matched_acceptor_sum = get_group_matched_data_summary(acceptor, ref_group = '(0.8, 1]', group1 = '[0, 0.2]', group2 = '(0.2, 0.8]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_donor_sum = get_group_matched_data_summary(donor, ref_group = '(0.8, 1]', group1 = '[0, 0.2]', group2 = '(0.2, 0.8]', group_col = 'mean_prop_bin', freq_col = 'AF', seed = seed)
    matched_sum = rbind(matched_acceptor_sum, matched_donor_sum) %>% mutate(splice_variant = rep(c('Acceptor', 'Donor'), each=3))
      }else{
    stop("Invalid group name!")
  }
  if(save_plot){
    plt = matched_sum %>%
      ggplot + aes(x=splice_variant, y=prop, ymin=prop-sd, ymax=prop+sd, color=mean_prop_bin, fill=mean_prop_bin) +
      geom_pointrange(stat="identity", position=position_dodge(width=0.4), size=0.4) +
      labs(y = 'Proportion', x = 'Splice Variant')  +
      scale_y_continuous(label=label_percent(accuracy=1)) +
      scale_x_discrete(labels= c('Acceptor', 'Donor')) +
      scale_color_manual(name= 'Mean proportion expressed', values = c("#9ECAE1", "#4292C6", "#084594")) +
      theme_classic() + themes

    png(output_path, height=4, width=6, units = 'in', res=300)
    print(plt)
    dev.off()
  }
  return(matched_sum)
}

matched_pext_sum1 = match_pext_variants( ref_group = '[0, 0.2]', save_plot = F, output_path = '/Users/wlu/Desktop/Figures/freq_matching/pext/pext_matched_to_1st.png')
matched_pext_sum2 = match_pext_variants( ref_group = '(0.2, 0.8]', save_plot = F, output_path = '/Users/wlu/Desktop/Figures/freq_matching/pext/pext_matched_to_2nd.png')
matched_pext_sum3 = match_pext_variants( ref_group = '(0.8, 1]', save_plot = F, output_path = '/Users/wlu/Desktop/Figures/freq_matching/pext/pext_matched_to_3rd.png')