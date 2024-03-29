source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figureS9 = function(save_plot = F, save_pdf = F, output_path){
  # Filtered n_var<2 & coverage<20
  rp_lambda_gene_caf = load_ukb_file('rp_lambda_gene_caf_300k.txt.bgz', subfolder='analysis/')
  rp_lambda_gene_caf = rp_lambda_gene_caf %>%
    pivot_longer(cols = contains('lambda_gc_'), names_to = 'labels', names_repair = 'unique', values_to = 'lambda_gc') %>%
    mutate(result_type = str_split(labels, "lambda_gc_") %>% map_chr(., 2)) %>%
    mutate(trait_type2 = factor(trait_type2, levels=trait_types),
          result_type = factor(result_type, levels = result_types),
          caf = factor(CAF_range, levels = caf_types, labels = caf_names),
          heritability = as.numeric(replace(modifier, is.na(modifier), '1')))

  figure =  rp_lambda_gene_caf %>%
    filter(result_type == 'skato') %>%
    ggplot + aes(x = n_cases, y = lambda_gc, color = trait_type) +
    geom_point(alpha = 0.5) + theme_classic() +
    labs(y = 'Lambda GC \n (SKAT-O)', x = 'Number of cases') +
    geom_hline(yintercept = 1, lty = 2, size = 0.25) +
    scale_x_log10(label = comma) +
    trait_color_scale + trait_fill_scale +
    facet_grid(heritability ~ caf, scales = 'free', labeller = label_type) + themes +
    theme(strip.text = element_text(size = 8),
          plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"),
          axis.text.x = element_text(size = 8, angle = 30, vjust = 0.8),
          axis.text.y = element_text(size = 8),
          legend.key.size = unit(0.5, 'cm'),
          axis.title.y = element_text(size = 9, face = 'bold'),
          legend.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 7), )

  if(save_plot){
    png(output_path, height = 5, width = 9, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  if(save_pdf){
    ggsave(filename=output_path, plot=figure, device=cairo_pdf,
       width = 174/in2mm, height = 116/in2mm, units = "in", dpi = 300)
    # pdf(output_path, height = 116/in2mm, width = 174/in2mm)
    # print(figure)
    # dev.off()
  }
  return(figure)
}

figureS9(save_plot = F, save_pdf = T, output_path = paste0(output_path, 'final_pdf_figures/figureS9_random_pheno_lambda_gene_by_caf_h2.pdf'))
# figureS9(save_plot = T, output_path = paste0(output_path, 'figureS9_random_pheno_lambda_gene_by_caf_h2_skato.png'))

