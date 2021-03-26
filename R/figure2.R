source('~/ukb_exomes/R/constants.R')
detach(package:plyr)
test = 'skato'
output = '~/skato/'

icd_var_after = load_ukb_file('icd_min_p_var_after_SErm.txt.bgz')
icd_gene_after  = load_ukb_file(paste0('icd_min_p_gene_after_', test, '.txt.bgz'))

# figure 2 (A)
figure2a = function(save_plot = F, output_path){
  icd_long = icd_gene_after %>%
    select(-min_p) %>%
    mutate(chr = str_split(interval, ':') %>% map_chr(., 1) %>% str_replace_all(., '\\[chr', ''), 
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
           chr = str_split(locus, ':') %>% map_chr(., 1) %>% str_replace_all(., 'chr', ''), 
           position = str_split(locus, ':') %>% map_chr(., 2)) %>%
    mutate(chr = if_else(chr == 'X', '23', chr)) %>%
    mutate(chr = as.integer(if_else(chr == 'Y', '24', chr))) %>%
    mutate(icd10 = if_else(icd10 == 'B', 'A', icd10))
  figure = save_icd_manhattan_figure(icd_long, p_filter = 1e-2, width = 10, spacing = 3, sig_level = 8e-9, save_plot = save_plot, output_path = output_path)
  return(figure)
}

# figure 2 (C)
figure2c = function(save_plot = F, output_path){
  var_gene = load_ukb_file(paste0('var_gene_comparison_by_pheno_after_300k_', test, '.txt.bgz'))
  var_gene_summary = rbind(get_var_gene_overlap_count(data = var_gene, pheno_group = 'all', normalize = T, print = F), 
                           get_var_gene_overlap_count(data = var_gene, pheno_group = 'icd_first_occurrence', normalize = T, print = F),
                           get_var_gene_overlap_count(data = var_gene, pheno_group = 'categorical', normalize = T, print = F),
                           get_var_gene_overlap_count(data = var_gene, pheno_group = 'continuous', normalize = T, print = F)) %>%
    filter(label != 'None' & category != 'all') %>%
    mutate(label = factor(label, levels = c('Significant by variant', 'Significant by both', 'Significant by gene')),
           category = factor(category, levels = trait_types)) %>%
    group_by(category) %>%
    arrange(desc(label)) %>%
    mutate(value2 = value/sum(value)) %>%
    mutate(pos = cumsum(value2) - 0.5*value2)
  figure2_p1 = var_gene_summary %>%
    ggplot + aes(x = category, y = value) +
    geom_bar(stat="identity", aes(fill = category), width = 0.5) +
    labs(y = 'Significant Association', x = '')  +
    geom_bar(stat = "identity", fill = "black", aes(alpha = label), width = 0.5) +
    scale_alpha_manual(values = c(0.6, 0.8, 0.2), name = '', labels = c( 'Significant by variant', 'Significant by both', 'Significant by gene'))  +
    trait_color_scale + trait_fill_scale + scale_x_discrete(labels = trait_type_names) +
    theme_classic() + themes +
    theme(panel.spacing = unit(1, "lines"), 
          strip.background = element_blank(), 
          strip.placement = "top", 
          strip.text = element_text(face = 'bold'), 
          plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"), 
          axis.text.x = element_text(size = 8, face = 'bold'), 
          legend.title = element_text(size = 8, face = 'bold'), 
          legend.text = element_text(size = 8), )
  figure2_p2 = var_gene_summary %>%
    ggplot + aes(x = category, y = value2, label = round(value, 2)) +
    geom_bar(stat="identity", aes(fill = category), width = 0.5) +
    labs(y = 'Significant Association', x = '')  +
    geom_bar(stat = "identity", fill = "black", aes(alpha = label), width = 0.5) +
    scale_alpha_manual(values = c(0.6, 0.8, 0.2), name = '', labels = c( 'Significant by variant', 'Significant by both', 'Significant by gene'))  +
    trait_color_scale + trait_fill_scale + scale_x_discrete(labels = trait_type_names) +
    theme_classic() + themes +
    geom_text(aes(x = category, y = round(pos, 2)), size = 2, color = 'grey') +
    theme(panel.spacing = unit(1, "lines"), 
          strip.background = element_blank(), 
          strip.placement = "top", 
          strip.text = element_text(face = 'bold'), 
          plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"), 
          axis.text.x = element_text(size = 8, face = 'bold'), 
          legend.title = element_text(size = 8, face = 'bold'), 
          legend.text = element_text(size = 8), )
  figure = ggpubr::ggarrange(figure2_p1, figure2_p2, ncol = 2, labels = c('Associations per phenotype', 'Normalized'),
                             common.legend = TRUE, vjust = 1, hjust = 0, font.label = list(size = 10, color = "black", face = "bold", family = NULL))
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

# figure 2
figure2a_icd_manhattan_skato = figure2a(save_plot = T, output_path = paste0(output, 'figure2a_icd_manhattan_filtered_skato.png'))
figure2b_icd_manhattan_var = figure2b(save_plot = T, output_path = paste0(output, 'figure2b_icd_manhattan_variant_filtered.png'))
figure2c_var_gene_comparison = figure2c(save_plot = T, output_path = paste0(output, 'figure2c_var_gene_comparison.png'))

figure2 = function(save_plot = F, output_path){
  figure = ggarrange(figure2a_icd_manhattan_skato, figure2b_icd_manhattan_var, figure2c_var_gene_comparison,
                     ncol = 1, nrow = 3, labels = c('(A) SKAT-O', '(B) Single-Variant', '(C)'), label.args = list(gp = gpar(font = 2, cex = 1), vjust = 1))
  if(save_plot){
    png(output_path, height = 12, width = 6, units = 'in', res=300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure2(save_plot = T, output_path = paste0(output, 'figure2_', test, '.png'))




