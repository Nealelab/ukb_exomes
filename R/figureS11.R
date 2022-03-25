source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figureS11 = function(result_type, save_plot = F, output_path){
  if(result_type == 'gene'){
    gene_lambda = load_ukb_file(paste0('lambda_by_pheno_expected_ac_filtered_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
    lambda = gene_lambda %>%
      select(-(6:8)) %>%
      pivot_longer(cols = contains('n_genes_'), names_to = 'labels', names_prefix = 'n_genes_', names_repair = 'unique', values_to = 'n_genes') %>%
      merge(., gene_lambda %>%
              select(-(9:11)) %>%
              pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_prefix = 'lambda_gc_', names_repair = 'unique', values_to = 'lambda_gc'),
            by = c(colnames(gene_lambda)[-(6:11)], 'labels')) %>%
      mutate(labels = factor(labels, levels = result_types, labels = result_names)) %>%
      mutate(trait_type2 = factor(trait_type2, levels = trait_types),
             Range = factor(str_split(expected_AC_range, ':') %>% map_chr(., 2),
                            levels = ac_types, labels = ac_names))

    figure = lambda %>%
      filter(n_genes > 1000) %>%
      ggplot + aes(x = n_cases_defined, y = lambda_gc, color = trait_type2, label = phenocode) +
      geom_point(alpha = 0.5, size = 0.8) + ylim(0, 2) + theme_classic() +
      labs(y = 'Lambda GC', x = 'Expected Number of Cases/Individuals') +
      geom_hline(yintercept = 1, lty = 2) +
      scale_x_log10(label = comma, limits = c(1, NA)) +
      trait_color_scale + trait_fill_scale +
      facet_grid(labels~Range, labeller = label_type) + themes +
      theme(axis.text.x = element_text(size = 9, angle = 45, vjust = 0.8))
    if(save_plot){
      png(output_path, height = 5, width = 8, units = 'in', res = 300)
      print(figure)
      dev.off()
    }
    }else{
      var_lambda = load_ukb_file(paste0('lambda_by_pheno_expected_ac_var_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
      lambda = var_lambda %>%
        mutate(trait_type2 = factor(trait_type2, levels = trait_types),
               Range = factor(str_split(expected_AC_range, ':') %>% map_chr(., 2),
                              levels = ac_types, labels = ac_names))
      figure = lambda %>%
        filter(n_vars > 1000) %>%
        ggplot + aes(x = n_cases_defined, y = lambda_gc, color = trait_type2, label = phenocode) +
        geom_point(alpha = 0.5, size = 0.8) + ylim(0, 2) + theme_classic() +
        labs(y = 'Lambda GC', x = 'Expected Number of Cases/Individuals') +
        geom_hline(yintercept = 1, lty = 2) +
        scale_x_log10(label = comma, limits = c(1, NA)) +
        trait_color_scale + trait_fill_scale +
        facet_grid(~Range, labeller = label_type) + themes +
        theme(axis.text.x = element_text(size = 9, angle = 45, vjust = 0.8))
      if(save_plot){
        png(output_path, height = 3, width = 15, units = 'in', res = 300)
        print(figure)
        dev.off()
      }

  }


  return(figure)
}

figureS11('gene', save_plot = T, output_path = paste0(output_path, 'figureS11_lambda_expected_ac_gene.png'))



# figureS11 = function(save_plot = F, output_path){
#   lambda_gene_caf = load_ukb_file(paste0('lambda_by_pheno_freq_filtered_gene_', tranche,'.txt.bgz'), subfolder = 'qc/lambda_gc/')
#   lambda_gene_caf = lambda_gene_caf %>%
#     pivot_longer_lambda_data() %>%
#     mutate(trait_type2 = factor(trait_type2, levels = trait_types),
#            CAF_range = factor(CAF_range, levels = caf_types))
#
#   figure = lambda_gene_caf %>%
#     ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type2, label = phenocode) +
#     geom_point(alpha = 0.5, size = 0.8) + ylim(0, 2) + theme_classic() +
#     labs(y = 'Lambda GC', x = 'Number of Cases') +
#     geom_hline(yintercept = 1, lty = 2) +
#     scale_x_log10(label = comma, limits = c(2, NA)) +
#     trait_color_scale + trait_fill_scale +
#     facet_grid(result_type~CAF_range, labeller = label_type) + themes +
#     theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.8))
#
#   if(save_plot){
#     png(output_path, height = 5, width = 7.5, units = 'in', res = 300)
#     print(figure)
#     dev.off()
#   }
#   return(figure)
# }
#
# figureS11(save_plot = T, output_path = paste0(output_path, 'figureS11_lambda_gene_by_caf.png'))


