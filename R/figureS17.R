source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figureS17 = function(type = 'full', save_plot = F, save_pdf = F, output_path){
  height_UKB  = load_ukb_file('height_beta_500k.txt.bgz', subfolder = 'analysis/', force_cols = cols(locus = col_character())) %>%
    separate(markerID, c("chrom", "position", "REF", "ALT"), sep = "[^[:alnum:]]") %>%
    mutate(locus = new_locus)

  if(type == 'full'){
    download.file("https://portals.broadinstitute.org/collaboration/giant/images/5/59/Height_All_add_SV.txt.gz", 'giant_height.txt.gz')
    height_GIANT  = read_delim(gzfile('giant_height.txt.gz'), delim = '\t') %>%
      mutate(locus = paste0(CHR, ':', `POS`))
  }else{
    download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature21039/MediaObjects/41586_2017_BFnature21039_MOESM49_ESM.xlsx", 'giant_height.xlsx')
    height_GIANT = read_excel("giant_height.xlsx", sheet = 6, skip = 2) %>%
      mutate(CHR = as.character(as.integer(Chr)), POS = as.character(`Pos (hg19)`), locus = paste0(CHR, ':', POS), REF = Ref, ALT = Alt, beta = Beta)
  }
  height = merge(height_UKB, height_GIANT, by = c('locus', 'REF', 'ALT')) %>%
    mutate(annotation = factor(if_else(annotation == 'missense', 'missense|LC', annotation), levels = annotation_types),
           beta_ukb = if_else(beta < 0, -BETA, BETA),
           beta_giant = abs(beta))

  figure = height %>%
    filter(AF > 2e-5) %>%
    ggplot + aes(x = beta_giant, y = beta_ukb, color = -log10(AF)) +
    geom_point(alpha = 0.8, position = 'jitter') +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    labs(y = "UK Biobank 450k", x = "GIANT", color = expression(bold(paste('-log'[10], '(AF)')))) +
    theme_classic() + themes +
    guides(color = guide_legend(override.aes = list(size = 5, pch = 19, alpha = 1) ) )

  if(save_plot){
    png(output_path, height = 5, width = 5, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  if(save_pdf){
    pdf(output_path, height = 85/in2mm, width = 85/in2mm)
    print(figure)
    dev.off()
  }
  return(figure)
}

figureS17(type = 'sub', save_plot = T, save_pdf = T, output_path = paste0(output_path, 'final_pdf_figures/figureS17_giant_height_comparison_sub.pdf'))
# figureS17(type = 'full', save_plot = T, save_pdf = T, output_path = paste0(output_path, 'final_pdf_figures/figureS17_giant_height_comparison_full.pdf'))
# figureS17(type = 'sub', save_plot = T, output_path = paste0(output_path, 'figureS17_giant_height_comparison_sub.png'))