source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

icd_var_after = load_ukb_file('icd_min_p_var_after_SErm.txt.bgz')
icd_gene_skato_after  = load_ukb_file('icd_min_p_gene_after_skato.txt.bgz')
output_path = '~/Desktop/'

# figure 2 (A)
figure2a = function(save_plot = F, output_path){
  icd_long = icd_gene_skato_after %>%
    select(-min_p) %>%
    mutate(chr = str_split(interval, ':') %>% map_chr(., 1) %>% str_replace_all(.,'\\[chr', ''),
           position = str_split(interval, ':') %>% map_chr(., 2) %>% str_split(., '-') %>% map_chr(., 1)) %>%
    mutate(chr = if_else(chr == 'X', '23', chr)) %>%
    mutate(chr = as.integer(if_else(chr == 'Y', '24', chr))) %>%
    pivot_longer(cols = contains('_min_p'), names_to = 'icd10', values_to = 'min_p') %>%
    mutate(icd10 = str_split(icd10, '_min_p') %>% map_chr(., 1)) %>%
    mutate(icd10 = if_else(icd10 == 'B', 'A', icd10))

  figure = save_icd_manhattan_figure(icd_long, p_filter = 1e-2, width = 10, spacing = 3, sig_level = 2.5e-8, save_plot = save_plot, output_path = output_path)
  return(figure)
}

# figure 2 (B)
figure2b = function(save_plot = F, output_path){
  icd_long = icd_var_after %>%
    select(-min_p) %>%
    pivot_longer(cols = contains('_min_p'), names_to = 'icd10', values_to = 'min_p') %>%
    mutate(icd10 = str_split(icd10, '_min_p') %>% map_chr(., 1),
           chr = str_split(locus, ':') %>% map_chr(., 1) %>% str_replace_all(.,'chr', ''),
           position = str_split(locus, ':') %>% map_chr(., 2)) %>%
    mutate(chr = if_else(chr == 'X', '23', chr)) %>%
    mutate(chr = as.integer(if_else(chr == 'Y', '24', chr))) %>%
    mutate(icd10 = if_else(icd10 == 'B', 'A', icd10))
  figure = save_icd_manhattan_figure(icd_long, p_filter = 1e-2, width = 10, spacing = 3, sig_level = 8e-9, save_plot = save_plot, output_path = output_path)
  return(figure)
}

# figure 2 (C)
figure2c = function(save_plot = F, output_path){
  pheno_sig_after = load_ukb_file('pheno_sig_cnt_after_gene_300k.txt.bgz')
  pheno_var_sig_after = load_ukb_file('pheno_sig_cnt_after_SErm_var_300k.txt.bgz')

  var_pheno_sum = pheno_var_sig_after %>%
    group_by(trait_type2) %>%
    sig_cnt_summary('all_sig_var_cnt') %>%
    mutate(test = "Single-Variant Test")
  gene_pheno_sum = pheno_sig_after %>%
    group_by(trait_type2) %>%
    sig_cnt_summary('all_sig_gene_cnt_skato') %>%
    mutate(test = "SKAT-O")
  pheno_sig_cnt_sum = rbind(var_pheno_sum, gene_pheno_sum) %>%
    mutate(trait_type2 = factor(trait_type2, levels = trait_types))

  figure = pheno_sig_cnt_sum %>%
    ggplot + aes(x = trait_type2, y = mean, color = trait_type2, fill = trait_type2) +
    geom_bar(stat = "identity",position = position_dodge(), width = 0.5, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .2, position = position_dodge(.9)) +
    labs(y = 'Significant Association', x = 'Phenotype')  +
    scale_x_discrete(labels = trait_type_names) +
    trait_color_scale + trait_fill_scale +
    theme_classic() + themes +
    facet_wrap(~test, scales = "free", labeller = label_type) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.placement = "top",
          strip.text = element_text(face = 'bold'),
          plot.margin = margin(0.5,0.1,0.1,0.1, "cm"),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 8),)
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

# figure 2 (D)
figure2d = function(filter = T, normalize = T, save_plot = F, output_path){
  if(filter){
    var_gene = load_ukb_file('var_gene_comparison_by_pheno_after_300k.txt.bgz')
  }else{
    var_gene = load_ukb_file('var_gene_comparison_by_pheno_before_300k.txt.bgz')
  }
  var_gene_summary = rbind(get_var_gene_overlap_count(data = var_gene, pheno_group = 'all', normalize = normalize, print = F),
                           get_var_gene_overlap_count(data = var_gene, pheno_group = 'icd_first_occurrence', normalize = normalize, print = F),
                           get_var_gene_overlap_count(data = var_gene, pheno_group = 'categorical', normalize = normalize, print = F),
                           get_var_gene_overlap_count(data = var_gene, pheno_group = 'continuous', normalize = normalize, print = F)) %>%
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
          legend.text = element_text(size = 8),)
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

# figure 2
figure2a_icd_manhattan_skato = figure2a(save_plot = T, output_path = paste0(output_path, 'figure2a_icd_manhattan_skato_filtered.png'))
figure2b_icd_manhattan_var = figure2b(save_plot = T, output_path = paste0(output_path, 'figure2b_icd_manhattan_variant_filtered.png'))
figure2c_pheno_sig_cnt = figure2c(save_plot = T, output_path = paste0(output_path, 'figure2c_pheno_sig_cnt.png'))
figure2d_var_gene = figure2d(filter = T, normalize = T, save_plot = T, output_path = paste0(output_path, 'figure2d_var_gene_comparison.png'))

figure2 = function(save_plot = F, output_path){
  if(save_plot){
    png(paste0(output_path, 'figure2.png'), height = 10, width = 12, units = 'in', res=300)
    print(ggarrange(figure2a_icd_manhattan_skato,
                    figure2b_icd_manhattan_var,
                    figure2c_pheno_sig_cnt,
                    figure2d_var_gene, ncol = 2, nrow = 2, labels = c('(A)', '(B)', '(C)', '(D)'),
                  label.args = list(gp = gpar(font = 2, cex = 1), vjust = 1)))
    dev.off()
  }
  return(figure)
}

figure2(save_plot = T, output_path = output_path)





