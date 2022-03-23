source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

icd_var_after = load_ukb_file(paste0('icd_min_p_var_filtered_',test,'_', tranche,'.txt.bgz'), subfolder = 'analysis/')
icd_gene_after = load_ukb_file(paste0('icd_min_p_gene_filtered_', test, '_', tranche, '.txt.bgz'), subfolder = 'analysis/')

# figure 2 (A)
figure2a = function(save_plot = F, output_path){
  icd_long = icd_var_after %>%
    # select(-min_p) %>%
    rowwise() %>%
    pivot_longer(cols = contains('_min_p'), names_to = 'icd10', values_to = 'min_p') %>%
    mutate(icd10 = strsplit(icd10, '_min_p') %>% map_chr(., 1),
           chr = strsplit(locus, ':') %>% map_chr(., 1) %>% str_replace_all(., 'chr', ''),
           position = strsplit(locus, ':') %>% map_chr(., 2)) %>%
    mutate(chr = if_else(chr == 'X', '23', chr)) %>%
    mutate(chr = as.integer(if_else(chr == 'Y', '24', chr)))
  figure = save_icd_manhattan_figure(icd_long, p_filter = 1e-2, width = 10, spacing = 3, sig_level = 8e-9, save_plot = save_plot, output_path = output_path)
  return(figure)
}

# figure 2 (B)
figure2b = function(save_plot = F, output_path){
  icd_long = icd_gene_after %>%
    # select(-min_p) %>%
    mutate(chr = str_split(interval, ':') %>% map_chr(., 1) %>% str_replace_all(., '\\[chr', ''),
           position = str_split(interval, ':') %>% map_chr(., 2) %>% str_split(., '-') %>% map_chr(., 1)) %>%
    mutate(chr = if_else(chr == 'X', '23', chr)) %>%
    mutate(chr = as.integer(if_else(chr == 'Y', '24', chr))) %>%
    # rowwise() %>%
    # mutate(A_min_p = min(A_min_p, B_min_p)) %>%
    # select(-B_min_p) %>%
    pivot_longer(cols = contains('_min_p'), names_to = 'icd10', values_to = 'min_p') %>%
    mutate(icd10 = str_split(icd10, '_min_p') %>% map_chr(., 1))
  figure = save_icd_manhattan_figure(icd_long, p_filter = 1e-2, width = 10, spacing = 3, sig_level = 2.5e-7, save_plot = save_plot, output_path = output_path)
  return(figure)
}

# figure 2 (C/D)
figure2cd = function(type = 'pheno', save_plot = F, output_path){
  var_gene = load_ukb_file(paste0('var_gene_comparison_by_', type,'_filtered_', test, '_3annt_', tranche,'.txt.bgz'), subfolder = 'analysis/')
  if(type == 'pheno'){
    set = c('categorical', 'continuous', 'icd10')
    group_type = trait_types
  }else{
    set = c('missense|LC', 'pLoF', 'synonymous')
    group_type = annotation_types
  }
  var_gene_summary = rbind(get_var_gene_overlap_count(data = var_gene, type = type, group = 'all', normalize = T, print = F),
                           get_var_gene_overlap_count(data = var_gene, type = type, group = set[1], normalize = T, print = F),
                           get_var_gene_overlap_count(data = var_gene, type = type, group = set[2], normalize = T, print = F),
                           get_var_gene_overlap_count(data = var_gene, type = type, group = set[3], normalize = T, print = F)) %>%
    filter(label != 'None' & category != 'all') %>%
    mutate(label = factor(label, levels = c('Significant by variant', 'Significant by both', 'Significant by gene')),
           category = factor(category, levels = group_type)) %>%
    group_by(category) %>%
    arrange(desc(label)) %>%
    mutate(value2 = value/sum(value)) %>%
    mutate(pos = cumsum(value2) - 0.5*value2)
  if(type == 'pheno'){
    color_scale = trait_color_scale
    axis_scale = scale_x_discrete(labels = trait_type_names)
  }else{
    color_scale = annotation_color_scale
    axis_scale = scale_x_discrete(labels = annotation_names)
  }
  figure2_p1 = var_gene_summary %>%
    ggplot + aes(x = category, y = value) +
    geom_bar(stat="identity", aes(color = category), fill='grey80', width = 0.5) +
    labs(y = 'Significant Association', x = '')  +
    geom_bar(stat = "identity", aes(alpha = label), width = 0.5) +
    scale_alpha_manual(values = c(0.6, 1, 0.1), name = '', labels = c( 'Significant by variant', 'Significant by both', 'Significant by gene'),
                       guide = guide_legend(override.aes = list(fill = 'grey30') ))  + color_scale + axis_scale +
    theme_classic() + themes +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.placement = "top",
          strip.text = element_text(face = 'bold'),
          plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"),
          axis.text.x = element_text(size = 13, face = 'bold'),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13, face = 'bold'),
          legend.title = element_text(size = 11, face = 'bold'),
          legend.text = element_text(size = 11), )

  figure2_p2 = var_gene_summary %>%
    ggplot + aes(x = category, y = value2, label = round(value, 2)) +
    geom_bar(stat="identity", aes(color = category), fill='grey80', width = 0.5) +
    labs(y = 'Significant Association', x = '')  +
    geom_bar(stat = "identity", aes(alpha = label), width = 0.5) +
    scale_alpha_manual(values = c(0.6, 1, 0.1), name = NULL, labels = c( 'Significant by variant', 'Significant by both', 'Significant by gene'),
                       guide = guide_legend(override.aes = list(fill = 'grey30') ))  + color_scale + axis_scale +
    theme_classic() + themes +
    geom_text(aes(x = category, y = round(pos, 2)), size = 3.8, color = 'black') +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.placement = "top",
          strip.text = element_text(face = 'bold'),
          plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"),
          axis.text.x = element_text(size = 13, face = 'bold'),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13, face = 'bold'),
          legend.title = element_text(size = 11, face = 'bold'),
          legend.text = element_text(size = 11), )
  figure = ggpubr::ggarrange(figure2_p1, figure2_p2, ncol = 2, labels = c(paste0('Associations per ', if_else(type=='pheno', 'phenotype', 'group')), 'Normalized'),
                             common.legend = TRUE, vjust = 1, hjust = 0, font.label = list(size = 12, color = "black", face = "bold", family = NULL))
  if(save_plot){
    png(output_path, height = 4, width = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}



# figure 2
figure2a_icd_manhattan_var = figure2a(save_plot = T, output_path = paste0(output, 'figure2a_icd_manhattan_filtered_var.png'))
figure2b_icd_manhattan_gene = figure2b(save_plot = T, output_path = paste0(output, paste0('figure2b_icd_manhattan_filtered_', test, '.png')))
figure2c_var_gene_comparison_pheno = figure2cd(type = 'pheno', save_plot = T, output_path = paste0(output, 'figure2c_var_gene_comparison_pheno.png'))
figure2d_var_gene_comparison_annotation = figure2cd(type = 'gene', save_plot = T, output_path = paste0(output, 'figure2d_var_gene_comparison_annotation.png'))

figure2 = function(save_plot = F, output_path){
  figure = ggarrange(figure2a_icd_manhattan_var, figure2b_icd_manhattan_gene, figure2c_var_gene_comparison_pheno, figure2d_var_gene_comparison_annotation,
                     ncol = 1, nrow = 4, labels = c('(A) Single-Variant', '(B) SKAT-O', '(C) Phenotypes', '(D) Groups'), label.args = list(gp = gpar(font = 2, cex = 1.5), vjust = 1))
  if(save_plot){
    png(output_path, height = 20, width = 10, units = 'in', res=300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure2(save_plot = T, output_path = paste0(output, 'figure2_', test, '.png'))





