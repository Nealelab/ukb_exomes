source('~/ukb_exomes/R/constants.R')
detach(package:plyr)

figure5 = function(save_plot = F, output_path){
  dti_label <- read_csv('~/OneDrive/2020Work/UKBB/Analysis300K/scrib_dti/DTI.csv')[1:110, c(1,4)]
  dti_sig8 <- read_csv('~/OneDrive/2020Work/UKBB/Analysis300K/scrib_dti/ukbiobank_dti_var_chr8_mean.csv')
  dti_sig_FA <- dti_sig8 %>%
    filter(grepl('mean_pheno', label)) %>%
    mutate(label = as.numeric(str_replace(label, pattern = 'mean_pheno', replacement = '') ))%>%
    filter(as.numeric(label)<=61) %>%
    merge(., dti_label, by = 'label', all.x=T)
  figure = dti_sig_FA %>%
    filter(phenocode %in% c('Average_FA', 'BCC_FA','SCC_FA')) %>%
    ggplot + aes(x = POS, y = -log10(P)) +
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(5e-8), lty = 2, color = 'red')+
    geom_rect(aes(xmin=143790920, xmax=143815773, ymin=0,ymax=Inf), alpha=0.2, fill="red") +
    scale_x_continuous(label = label_bytes()) +
    scale_y_continuous(trans = gwas_loglog_trans(), breaks = loglog_breaks, name = expression(bold(paste('-log'[10], '(', italic(p), ')'))))+
    labs(x = 'Chromosome 8', y = expression(bold(paste('-log'[10], '(', italic(p), ')')))) +
    theme_classic() + themes +
    facet_wrap(~phenocode, ncol=1)
  if(save_plot){
    png(output_path, height = 6, width = 4, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figure5(save_plot = T, output_path = paste0(output, 'figure5_ukbiobank_dti_chr8_manhattan_plot_FA.png'))