---
title: "Figure 1"
author: "Ilaria Cherchi"
date: "2025-02-05"
output: 
  html_document:
    keep_md: yes
---





Code to reproduce Figure 1, Supplementary Figure 1 and Supplementary Figure 2 from (add doi..).

# Figure 1

## Figure 1A
Sketch produced with Inkscape saved in Sketches/...

## Figure 1B

### Prepare data
Load samtools statistics from supplementary data:


``` r
full_combined <- fread(paste0(supp_table_path, "Supplementary_table_1.csv"))
ont_samtools <- fread(paste0(tables_path, "1kgp_long_summary_samtools.tsv"))
ont_cramino <- fread(paste0(tables_path, "1kgp_long_summary_cramino.tsv"))

ont_samtools <- left_join(ont_samtools %>% dplyr::select(-c(`average quality`)), 
                          ont_cramino %>% distinct(sample, `Average quality`, assembly) %>% 
                            dplyr::rename(`average quality` = `Average quality`), 
                          by = c("sample", "assembly"))

full_combined <- bind_rows(full_combined, ont_samtools)

full_combined$design <- factor(full_combined$design, 
                               levels = c("1kGP WGS", "CCLE WGS", 
                                          "Breast cancer\nWES", 
                                          "Prostate cancer\n3Mb panel", 
                                          "1kGP WGS\nlong-read"))
```

Filter and format statistics of interest (average quality and error rate):


``` r
full_combined <- full_combined %>% 
  mutate_at(vars(-one_of(c("assembly", "design", "sample"))) , as.numeric) %>% 
  pivot_longer(-c("assembly", "design", "sample"), names_to = "metric", values_to = "value")

rcd = c("average quality" = "Average quality", "error rate" = "Error rate (%)")

combined <- full_combined %>% mutate(metric = gsub(":", "", metric)) %>%
  filter(metric %in% c("average quality", "error rate")) %>%
  mutate(value = case_when(metric == "error rate" ~ value*100, 
                           .default = value)) %>%
  mutate(metric = recode(metric, !!!rcd))
```

Obtain statistics increment or decrement from hg38 to T2T-CHM13:


``` r
combined <- combined %>% filter(!is.na(value)) %>% group_by(sample, metric) %>% filter(n() == 2) %>% ungroup()
combined %>% pivot_wider(names_from = "assembly", values_from = "value") %>%
  mutate(change = - ((hg38 - `T2T-CHM13`) / hg38) * 100) %>%
  dplyr::summarize(median = median(change[!is.na(change)]), 
                   mean = mean(change[!is.na(change)]), 
                   min = min(change[!is.na(change)]), 
                   max = max(change[!is.na(change)]), .by = c("metric"))
```

```
## # A tibble: 2 × 5
##   metric          median   mean    min    max
##   <chr>            <dbl>  <dbl>  <dbl>  <dbl>
## 1 Error rate (%)  -10.5  -10.6  -23.3  -0.876
## 2 Average quality   1.70   1.97  -1.50 14.9
```

Format data and cohort labels:


``` r
counts <- combined %>% filter(metric == "Error rate (%)", assembly == "hg38") %>% 
  summarize(n=n(), value = min(value), metric = unique(metric), .by = "design")

levels(combined$design) <- levels(counts$design) <- c("1kGP WGS", "CCLE WGS", 
                                  "Breast cancer\nWES 33Mb", "Prostate cancer\n3Mb panel",
                                  "1kGP WGS\nlong-read")
levels(combined$design) <- paste0(levels(combined$design), "\nN=", counts %>% arrange(design) %>% pull(n))
```

Plot pairwise comparison and wilcoxon paired test of selected statistics values for hg38 versus T2T-CHM13:

![](Figure_1_S1_S2/figure1B.png)<!-- -->

## Figure 1C
Sketch ...

## Figure 1D

Load AF data and perform classification of SNP status:


``` r
# Load data
load(paste0(rdata_path, "targeted_af_het.RData"))

load(paste0(rdata_path, "wes_af_het_all.RData"))

load(paste0(rdata_path, "vcf_af_het.RData"))
vcf_comp <- vcf_comp %>% filter(cov.hg38 > 60 | cov.T2T > 60)
annotations <- fread(paste0(utilities_path, "igsr_samples.tsv"))

ann <- annotations %>% distinct(`Sample name`, `Superpopulation code`, "Sex")
colnames(ann) <- c("sample","population","sex")

vcf_comp <- left_join(vcf_comp, ann, by = "sample") %>% filter(population != "EUR,AFR")

levels(vcf_comp$cl) <- rev(c("VAF consistent", "VAF het for hg38 only", 
                                 "VAF het for T2T-CHM13 only", "VAF inconsistent", "Allele swapping"))
```

Transform to long format and compute peaks:


``` r
targeted_af_long <- compute_af_long(targeted_comp)
wes_af_long <- compute_af_long(wes_comp)
vcf_af_long <- compute_af_long(vcf_comp)

targeted_peaks <- peaks_af(targeted_af_long)
wes_peaks <- peaks_af(wes_af_long)
vcf_peaks <- peaks_af(vcf_af_long)

wes_density <- create_density(wes_af_long, wes_peaks)
wes_density <- create_density(wes_af_long, wes_peaks)
vcf_density <- create_density(vcf_af_long, vcf_peaks)

af_long <- bind_rows(list("Prostate cancer 3Mb panel (200-1000x)" = targeted_af_long, 
                          "Breast cancer WES 33Mb (96-500x)" = wes_af_long,
                          "1kGP VCF (cov > 60)" = vcf_af_long), .id = "design")
peaks <- bind_rows(list("1kGP VCF (cov > 60)" = vcf_peaks, 
                        "Breast cancer WES 33Mb (96-500x)" = wes_peaks,
                        "Prostate cancer 3Mb panel (200-1000x)" = targeted_peaks), .id = "design")
```

Compute difference with respect to expected 0.5 threshold and prepare plots:

``` r
dt <- peaks %>% group_by(assembly) %>% summarize(median = median(peak), diff = 0.5 - median,
                                                 y = median, yend = 0.5)

a <- (ggplot(peaks, aes(x = assembly))
  + geom_hline(yintercept = 0.5, linetype = "dotted")
  + geom_textsegment(data = dt, aes(x = as.numeric(factor(assembly)) - 0.2, 
                                    xend = as.numeric(factor(assembly)) - 0.2, color = assembly, 
                                    y = y, yend = yend, 
                                    label = round(diff, 3)),
                 show.legend = F)
  + geom_line(aes(x = as.numeric(as.factor(assembly)), y = as.numeric(peak), group = design),
               colour = "gray50", alpha = 0.4, position = position_dodge(0.2))
  + geom_point(aes(x =  as.numeric(as.factor(assembly)), y =as.numeric(peak), 
                   group = design, shape = design, color = assembly),
                alpha = 1, position = position_dodge(0.2), size = 3) #color = "darkgray")
  + coord_flip()
  + theme_publication()
  + scale_color_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                               labels = c("hg38", "T2T-CHM13"))
  + labs(y = "VAF", x = "", fill = "Assembly", shape = "Design", color = "Assembly version")
  + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top",
          strip.text.y.left = element_text(angle = 0),
           legend.box = "vertical",
    legend.title.align = 0.5)
)

shape_map <- c("Circle" = "1kGP VCF (cov > 60)", 
               "Triangle" = "Breast cancer WES 33Mb (96-500x)", 
               "Square" = "Prostate cancer 3Mb panel (200-1000x)")
symbol_map <- c("Circle" = "\u25CF", "Triangle" = "\u25CF", "Square" = "\u25A0")

levels(af_long$design) <- c("1kGP VCF (cov > 60)", "Breast cancer WES 33Mb (96-500x)", "Prostate cancer 3Mb panel (200-1000x)")
facet_labels <- paste0(symbol_map, " ", levels(af_long$design))

b <- (ggplot(af_long, aes(x = af, color = assembly))
          + geom_vline(xintercept = 0.5, color = "gray40", lwd = 0.3, alpha = 0.5)
          + geom_vline(data = peaks, aes(xintercept = peak, color = assembly), 
                       linetype = "dotted", lwd = 0.5, show.legend = F)
          + labs(x = "VAF", y = "", colour = "Assembly version")
          + theme(legend.position = "bottom", panel.grid = element_blank())
          + geom_text(data = peaks %>% filter(assembly == "hg38"), 
                      aes(label = round(peak, 3) , x = peak, y = -0.15, colour = assembly), 
                      hjust = 1.2, show.legend  =  FALSE, size = 4)
          + geom_text(data = peaks %>% filter(assembly == "T2T"), 
                      aes(label = round(peak, 3), x = peak, y = -0.15, colour = assembly), 
                      hjust = -0.25, show.legend  =  FALSE, size = 4)
          + stat_density(
            geom = "line",position = "identity", bw = 0.01, show.legend = F)
          + theme_publication()
          + facet_wrap(design ~ ., scales = "free_y", nrow = 1, strip.position = "top",
                       labeller = labeller(design = facet_labels))
          + theme(legend.position = "bottom")
          + scale_x_continuous(limits = c(0.2, 0.8), breaks = c(0.2,0.5,0.8))
          + scale_color_manual(values = c("hg38" = "#7f061b","T2T" = "#ED7D31"),
                               labels = c("hg38", "T2T-CHM13"))
          + theme(legend.position = c(.85, .85), legend.text = element_text(size = 8),
                  legend.title = element_blank(), legend.box.background = element_rect(colour = "gray"),
                  legend.margin = margin(1.75, 1.75, 1.75, 1.75),
                   #strip.background = element_blank(),strip.text = element_blank(),
                  #axis.text.y = element_blank(), axis.ticks.y = element_blank())
          )
          
)
```


``` r
patchwork::wrap_plots(a , b , ncol= 1, heights = c(3, 6)) & plot_annotation(theme = theme(legend.margin  = margin(t = -4)))
```

![](Figure_1_S1_S2/figure1D.png)<!-- -->

## Figure 1E

Assign VAF status to SNPs:

``` r
targeted_comp <- assign_snp_status(targeted_comp)
wes_comp <- assign_snp_status(wes_comp)

af_stats <- bind_rows(list("1kGP VCF (cov > 60)" = vcf_comp, 
                        "Breast cancer WES\n33Mb (96-500x)" = wes_comp,
                        "Prostate cancer\n3Mb panel (200-1000x)" = targeted_comp), .id = "design")
```

Compare VAF:

``` r
targeted_df_comp_diff <- extract_snp_diff(targeted_comp)
targeted_af_status <- assign_af_status(targeted_df_comp_diff)

wes_df_comp_diff <- extract_snp_diff(wes_comp)
wes_af_status <- assign_af_status(wes_df_comp_diff)

vcf_df_comp_diff <- extract_snp_diff(vcf_comp)
vcf_af_status <- assign_af_status(vcf_df_comp_diff)

status <- bind_rows(list("Prostate cancer\n3Mb panel (200-1000x)" = targeted_af_status, 
                          "Breast cancer WES\n33Mb (96-500x)" = wes_af_status,
                          "1kGP VCF (cov > 60)" = vcf_af_status), .id = "design")

af_stats$cl <- ifelse(af_stats$cl %in% 
                       c("VAF het for hg38 only", 
                         "VAF het for T2T-CHM13 only", 
                         "VAF inconsistent"),
                     "VAF inconsistent (VAF delta > 0.005)",
                     ifelse(af_stats$cl == "Allele swapping", 
                            "Allele swapping", "VAF consistent")
)
```

Compute difference in VAF and pvalue:

``` r
pval_df <- bind_rows(list("Prostate cancer\n3Mb panel (200-1000x)" = 
                           data.frame(pval = wilcox.test(targeted_df_comp_diff$af.T2T , 
                                       targeted_df_comp_diff$af.hg38, 
                                       alternative = "greater", paired = T )$p.value), 
                          "Breast cancer WES\n33Mb (96-500x)" =  
                           data.frame(pval = wilcox.test(wes_df_comp_diff$af.T2T , 
                                       wes_df_comp_diff$af.hg38, 
                                       alternative = "greater", paired = T )$p.value),
                          "1kGP VCF (cov > 60)" =  
                           data.frame(pval = wilcox.test(vcf_df_comp_diff$af.T2T , 
                                       vcf_df_comp_diff$af.hg38, alternative = "greater", paired = T )$p.value)), 
                    .id = "design")

pval_df$p <- ifelse(pval_df$pval == 0, "p<2.22e-16", as.character(format(pval_df$pval)))

pval_df <- left_join(pval_df, status)
pval_df <- pval_df %>% group_by(design, p) %>% summarize(y.position = max(n))
pval_df$group1 <- "> hg38"
pval_df$group2 <- "> T2T-CHM13"
```

Produce figure:

``` r
a <- (ggplot(af_stats,
        aes( x = 1))
   + facet_wrap(design ~ ., 
                scales = "free_y", nrow= 1, strip.position = "top",
                       labeller = labeller(design = facet_labels))
 
  + geom_bar(width = 0.5, aes(fill = cl))
  + theme_publication()
  + labs(x = "", y = "N SNPs", fill = "VAF status")
  + scale_fill_manual(values = c("VAF consistent" = "gray", 
                                 "VAF inconsistent (VAF delta > 0.005)" = "#8bad0d", 
                                 "Allele swapping" = "#86bfd5"),
                      breaks = c("Allele swapping", "VAF consistent", "VAF inconsistent (VAF delta > 0.005)"),
                      labels = c("VAF consistent" = "VAF coherent",
                                 "VAF inconsistent (VAF delta > 0.005)" = "VAF incoherent (VAF delta > 0.005)", 
                                 "Allele swapping" = "Allele swapping")
  )
  + geom_label_repel(
    aes(label= unlist(tapply(after_stat(count), list(..x.., ..PANEL..), 
                                        function(a) paste(round(100*a/sum(a), 2), '%'))),
                    
                    y = after_stat(count), fill = cl),
    stat = "count",  direction = "y", 
    #nudge_x = 0.2, 
    #min.segment.length = 0.2, 
    position=position_stack(vjust = 0.5), 
    color = "white", vjust = 0.5,  size = 5,
    show.legend = F) 
  + theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.subtitle = element_text(size = 18), legend.position = "bottom",
          legend.title.position = "top")
  + scale_x_discrete(expand = c(0,0.1))
 
)

b <- (ggplot(status %>% mutate(af_status = ifelse(af_status == "VAF increase",
                                               "> T2T-CHM13",
                                               "> hg38")))
  + geom_col( aes(x = af_status, y = n, fill = af_status), show.legend = F)
  + theme_publication()
  + facet_wrap(design ~ ., scales = "free_y", nrow = 1, strip.position = "top")
  + geom_text(data = pval_df, aes(x = 1.5, y = y.position, label = p), vjust = -0.5,
              size = 4)
  + scale_fill_manual(values = c("> T2T-CHM13" = "#ED7D31" ,
                                 "> hg38" = "#7f061b"))
  + scale_y_continuous(expand=expansion(mult=c(0,0.15)))
  + theme(legend.position = "bottom",
          legend.title.position = "left",
          strip.background = element_blank(),
          strip.text = element_blank())
  + labs(x = "VAF", y  = "N SNPs")
  
  )

patchwork::wrap_plots(a, b , nrow = 2, heights = c(4,2), guides = "collect") +
  plot_annotation(theme = theme(legend.position = "top", legend.title.position = "left"))
```

![](Figure_1_S1_S2/figure1E.png)<!-- -->

# Supplementary Figure 1

## Figure S1A

Compare more metrics from samtools and format data:

``` r
ve <- full_combined %>% pivot_wider(names_from = "assembly", values_from = "value") %>%
  mutate(change = - ((hg38 - `T2T-CHM13`) / hg38) * 100) %>% filter(abs(change) < 70 ) %>% 
  dplyr::summarize(median = median(change), mean = mean(change), min = min(change), max = max(change), .by = c("metric", "design")) 

ve <- ve %>%  mutate(metric = Hmisc::capitalize(metric)) %>% 
  filter(!is.na(min) & !is.nan(min), abs(min) > 0, metric != "Reads unmapped", !grepl("pair|length", metric))

order <- ve %>% ungroup() %>% dplyr::summarize(min = min(min), .by = "metric") %>% arrange(min) %>% pull(metric)

ve$metric <- factor(ve$metric, levels = order)
```

Plot percentage change from hg38 to T2T-CHM13 over multiple samtools metrics:


``` r
(ggplot(ve)
  + geom_linerange(
    aes(xmin = min, xmax = max, y = metric, color = design), 
    position = position_dodge(0.4))
  + geom_point(
    aes(x = median, y = metric, color = design), 
    position = position_dodge(0.4))
  + theme_publication()
  + scale_y_discrete(expand = c(.07, .07)) 
  + geom_label(data = . %>% ungroup() %>% dplyr::summarize(min = min(min), .by = "metric"),
               aes(y = metric, x = min, label = metric), size = 5, hjust = 1, nudge_x = -1,
              label.size = 0, fill = "#FFFCFC") 
  + labs(x = "Change from hg38 to T2T-CHM13 (%)", y = "", color = "Design" )
  + scale_x_continuous(expand = c(0.3,0.1,0.05,0.1), breaks = c( -25, 0 , 25, 50))
  + theme(
    legend.position = "bottom",
  plot.margin = margin(15, 0, 15, 15),
  panel.border = element_blank(),
  #plot.background = element_rect(color = "#FFFCFC", fill = "#FFFCFC"),
  panel.grid.major.x = element_line(color = "grey94"),
  panel.grid.major.y = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(color = "grey40"),
  axis.text.y = element_blank(),
  #axis.title = element_blank(),
  axis.ticks = element_blank()
    )
)
```

![](Figure_1_S1_S2/figureS1A.png)<!-- -->

## Figure S1B

Load Picard statistics data and combine all designs into a single data frame:


``` r
picard_df <- fread(paste0(supp_table_path, "Supplementary_table_1_Picard.csv"))
picard_df$design <- factor(picard_df$design, levels = c("1kGP WGS", "CCLE WGS", "Breast cancer\nWES", "Prostate cancer\n3Mb panel"))
levels(picard_df$design) <- c("1kGP WGS", "CCLE WGS", "Breast cancer\nWES 33Mb", "Prostate cancer\n3Mb panel")
```

Plot pairwise comparison and t-test of selected statistics values for hg38 versus T2T:

![](Figure_1_S1_S2/figureS1B.png)<!-- -->


# Supplementary figure 2

Load sample peaks from germline data (prostate cancer gDNA, breast cancer gDNA and 1kGP WGS data) from supplementary Table 2:

``` r
summary_peaks <- fread(paste0(supp_table_path, "Supplementary_table_2.csv"))

summary_peaks <- summary_peaks %>% pivot_longer(
    cols = -c(design, sample),
    names_to = c("variable", "assembly"),
    names_sep = "\\.",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  )


summary_peaks <- summary_peaks %>% filter(!sample %like% c("cfDNA|CFDNA"), design != "CCLE WGS (30x)") %>% filter(assembly == "T2T-CHM13") 
summary_peaks$design[summary_peaks$design == "Breast cancer WES (96-500x)"] <- "Germline WES (96x)"
summary_peaks$design[summary_peaks$design == "Prostate cancer 3Mb panel (200-1000x)"] <- "Germline 3Mb panel (200-1000x)"

palette <- c("#b66dff", "#009999", "#004949")
ord <- c("1kGP WGS (30x)",
         "Germline WES (96x)",
         "Germline 3Mb panel (200-1000x)")
setNames(palette, ord)
```

```
##                 1kGP WGS (30x)             Germline WES (96x) 
##                      "#b66dff"                      "#009999" 
## Germline 3Mb panel (200-1000x) 
##                      "#004949"
```


``` r
(ggplot(summary_peaks %>% mutate(design = factor(design, levels = ord)), 
        aes(x = sd, y = log(mean_cov), color = design))
  + geom_point(size = 3, alpha = 0.5)
   + theme_publication()
  + labs(x = "VAF standard deviation", y = "Log(Mean SNP coverage)", color = "Design")
  + theme(legend.position = "bottom", panel.grid = element_blank())
  + scale_color_manual(values = palette)
          
  + stat_cor(show.legend = F, size = 6, label.x.npc = "right",
             label.y.npc = "top", hjust = 1)
  #+ guides(color =  guide_legend(nrow = 2))
)
```

![](Figure_1_S1_S2/figureS2A.png)<!-- -->

## Figure S2B

Bin coverage from beast WES data:

``` r
comp_bin <- wes_af_long %>% filter(assembly == "hg38") %>%
  mutate(cov_bin = cut(cov.hg38, breaks = c(5, 20, 30, 50, 100, 200))) %>% 
  filter(!is.na(cov_bin)) %>% mutate(type = "Observed (hg38 WES 96-500x)", assembly = "Observed (hg38 WES 96-500x)")

n_snps_df <- comp_bin %>% group_by(cov_bin) %>% 
  dplyr::count() %>% 
  mutate(cov = sapply(cov_bin, function(x) as.numeric(gsub("\\(","", strsplit(as.character(x), ",")[[1]][1])))
         )
```


Compute a theoretical distribution of VAF at increasing coverage:

``` r
# Set alleles
universe_size <- 100000
truth_vect <- c(rep("A",(universe_size / 2)), rep("B", (universe_size / 2)))

# Set coverage bins
coverages <- c(5, 20, 30, 50, 100, 200)

set.seed(123)
all_AF_vect<-c()
all_cov_vect<-c()
for(j in seq(from = 1, to = length(coverages)-1)){
    # Extract number of SNPs per bin
    n_SNPs <- n_snps_df %>% filter(cov == coverages[j]) %>% pull(n)
    AF_vect <- rep(NA,n_SNPs)
    cov_vect <- rep(NA,n_SNPs)
    for(i in seq(from = 1, to = n_SNPs)){
        # Sample
        tmp_cov <- sample(seq(from = coverages[j], to = coverages[j+1]), 1)
        tmp_extr <- sample(truth_vect, tmp_cov, replace = FALSE)
        AF_vect[i] <- sum(tmp_extr == "A") / length(tmp_extr) 
        cov_vect[i] <- tmp_cov
    }
    all_AF_vect <- c(all_AF_vect,AF_vect)
    all_cov_vect <- c(all_cov_vect,cov_vect)
}

AF_df <- data.frame(all_AF_vect,all_cov_vect)
colnames(AF_df) <- c("af","coverage")

AF_df <- AF_df %>%
  mutate(cov_bin = cut(coverage, breaks = c(5, 20, 30, 50, 100, 200))) %>% 
  filter(!is.na(cov_bin)) %>% mutate(type = "Theoretical")
```


Compare theoretical and observed VAF distribution at increasing coverage:


``` r
a <- (ggplot(comp_bin, aes(x = round(af, 4), y = cov_bin))
  + geom_vline(xintercept = 0.5, color = "gray40", lwd = 0.3, alpha = 0.5)
  + geom_density_ridges2(bandwidth = 0.005, alpha = 0.7, scale = 1) 
  + facet_wrap(.~type, scales = "free")
  + theme_publication()
  + scale_x_continuous(limits = c(0.2, 0.8), breaks = c(0.2,0.5,0.8))
  + labs(x = "VAF", y = "Coverage bin", fill = "Classification")
  + theme(legend.position = "bottom", panel.grid = element_blank())
)

b <- (ggplot(AF_df, aes(x = round(af, 4), y = cov_bin))
  + geom_vline(xintercept = 0.5, color = "gray40", lwd = 0.3, alpha = 0.5)
  + geom_density_ridges2(bandwidth = 0.005, alpha = 0.7, scale = 0.9) 
  + facet_wrap(.~type, scales = "free")
  + theme_publication()
  + scale_x_continuous(limits = c(0.2, 0.8), breaks = c(0.2,0.5,0.8))
  + labs(x = "VAF", y = "Coverage bin", fill = "Classification")
  + theme(legend.position = "bottom", panel.grid = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
)

patchwork::wrap_plots(a,  b, nrow = 1, axis_titles = "collect")
```

![](Figure_1_S1_S2/figureS2B.png)<!-- -->

## Figure S2C

Plot SNP direct comparison and coverage for data presented in Figure 1D and 1E:


``` r
targeted_histogram <- create_histogram(targeted_df_comp_diff, targeted_af_status)
targeted_coverage <- create_cov_plot(targeted_df_comp_diff)

wes_histogram <- create_histogram(wes_df_comp_diff, wes_af_status)
wes_coverage <- create_cov_plot(wes_df_comp_diff)


vcf_histogram <- (ggplot()
      + geom_col(data = vcf_df_comp_diff %>% group_by(value) %>% dplyr::count(),
                 aes(x = value, y = n), fill = "#8bad0d", width = 0.009,
                 show.legend = F)
      + theme_publication()
      + labs(x = "", y = "N SNPs", fill = "AF status", color = "AF status")
      + scale_x_continuous(breaks = round(seq(-0.6, 0.6, 0.2),2), 
                           limits = c(-0.606, 0.606))
      + geom_segment(data = vcf_af_status,
                     aes(x = sx, xend = dx, y = -500 , yend = - 500, color = af_status), lwd = 3)
      #+ geom_segment(data = vcf_af_status, aes(x = sx, xend = dx, color = af_status, y = 0 , yend = 0), 
      #               lwd = 3, position = position_nudge(y = -5))
      + scale_color_manual(values = c("VAF increase"="#f44560",
                                      "VAF decrease" = "#44d1df"))
      + scale_fill_manual(values = c("VAF increase"="#f44560",
                                     "VAF decrease" = "#44d1df"))
      + geom_label(data = vcf_af_status,
                   aes(x = (dx + sx) / 2, y = -500 , 
                       fill = af_status, label = scales::percent(freq)), 
                   color = "white",show.legend = F, vjust = 0.5)
      + theme(axis.title.x = element_blank(),
              legend.position = "inside",
              legend.position.inside = c(.8, .3), 
              legend.text = element_text(size = 8),
              legend.title = element_blank(), 
              legend.box.background = element_rect(colour = "gray"),
              legend.margin = margin(1.75, 1.75, 1.75, 1.75))
)

vcf_coverage <- (ggplot()
      + geom_boxplot(data = vcf_df_comp_diff,
                     aes(x = value, y = cov.T2T, group = value))
      + theme_publication()
      + labs(x = "VAF T2T-CHM13 - VAF hg38", y = "Coverage", fill = "AF status", color = "AF status")
      + scale_y_reverse(labels = c("     0", "   100", "   200"),
                        breaks = c(0,100,200)
                        ) 
      + scale_x_continuous(position = 'top', 
                           breaks = round(seq(-0.6, 0.6, 0.2),2),
                           limits = c(-0.606, 0.606),
                           sec.axis = dup_axis())
      + theme( axis.text.x = element_blank(),
                axis.title.x.top = element_blank(),
               axis.ticks.x.bottom = element_blank())
)

histogram_row <- wrap_elements(patchwork::wrap_plots((vcf_histogram + labs(title = '1KGP VCF (cov > 60)',
                                      subtitle = paste0('N samples = ',
                                                        length(unique(vcf_comp$sample)))) +
                                  theme(plot.title = element_text(hjust = 0.52, size = 14, face = 'bold'),
                                              plot.subtitle = element_text(size = 12, hjust = 0.52))),
                               (wes_histogram +  guides(color="none") + 
                                    labs(title = 'Breast cancer WES (96-500x)',
                                      subtitle = paste0('N samples = ', length(unique(wes_comp$sample))))+
                                  theme(plot.title = element_text(hjust = 0.52, size = 14, face = 'bold'),
                                              plot.subtitle = element_text(size = 12, hjust = 0.52))),
                              (targeted_histogram +  guides(color="none") + 
                                    labs(title = 'Prostate cancer 3Mb panel (200-1000x)',
                                     subtitle = paste0('N samples = ', length(unique(targeted_comp$sample))))+
                                  theme(plot.title = element_text(hjust = 0.52, size = 14, face = 'bold'),
                                              plot.subtitle = element_text(size = 12, hjust = 0.52))),
                               widths = c(9,9,9), nrow = 1, axis_titles = "collect") +
                                plot_annotation(theme = theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5)))
                               )
coverage_row <- wrap_elements(patchwork::wrap_plots(vcf_coverage,
                                                    wes_coverage,
                                                    targeted_coverage,
                                                    widths = c(9,9,9), nrow = 1, axis_titles = "collect") +
                                plot_annotation(theme = theme(plot.margin = margin(-8, 5.5, 5.5, 5.5))))

histogram_row + plot_spacer() + coverage_row + plot_layout(ncol = 1, heights = c(3,-0.38,2))
```

![](Figure_1_S1_S2/figureS2C.png)<!-- -->

## Figure S2D


``` r
# Obtained from running Rscript plot_supp_2.R
load(paste0(rdata_path, "panels_kgp_ccle.RData"))
```


``` r
ggpubr::as_ggplot(panel)
```

![](Figure_1_S1_S2/figureS2D.png)<!-- -->

## Figure S2E

Compute peaks by population:


``` r
peaks <- vcf_af_long %>% filter(af >= 0.2 & af <= 0.8, population != "EUR,AFR") %>%
   group_by(assembly, population) %>% summarize(peak = extract_peaks(stats::density(af, bw  = 0.01, na.rm  =  T)))
```


``` r
a <- (ggplot(af_long %>% filter(af >= 0.2 & af <= 0.8) %>%
               filter(population != "EUR,AFR", (cov.T2T > 60 | cov.hg38 > 60 )), aes(x = af))
          + facet_grid(~ population )
          + geom_vline(xintercept = 0.5, color = "gray40", lwd = 0.3, alpha = 0.5)
          + geom_vline(data = peaks, aes(xintercept = peak, color = assembly), linetype = "dotted", lwd = 0.5, show.legend = F)
          + labs(x = "VAF", y = "Density", colour = "Assembly version")
          + theme(legend.position = "bottom", panel.grid = element_blank())
          + geom_text(data = peaks %>% group_by(population) %>% dplyr::slice_min(peak) %>% arrange(desc(assembly)) %>%
                        dplyr::slice(1), aes(label = round(peak, 3) ,
                                                                       x = peak, y = -0.15, colour = assembly),
                      hjust = 1.2, show.legend  =  FALSE, size = 5)
          + geom_text(data = peaks %>% group_by(population) %>% dplyr::slice_max(peak) %>% 
                        arrange(assembly) %>% dplyr::slice(1), aes(label = round(peak, 3),
                                                                      x = peak, y = -0.15, colour = assembly),
                      hjust = -0.25, show.legend  =  FALSE, size = 5)
          + geom_text(data = . %>% group_by(population) %>% distinct(sample) %>% dplyr::count(),
                      aes( x = 0.3, y = 5.5, label = paste0("N samples: ", n) ), size = 12/.pt)
          + stat_density(
            geom = "line",position = "identity", bw = 0.01, aes(color = assembly))
          + theme_publication()
          + theme(legend.position = "top", panel.grid = element_blank())
          + scale_x_continuous(limits = c(0.2, 0.8), breaks = c(0.2,0.5,0.8))
          + scale_color_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                               labels = c("hg38", "T2T-CHM13"))
)

b <- (ggplot(vcf_comp %>% filter(population != "EUR,AFR"),
        aes( x = "", fill = cl))
  + facet_grid( ~ population, scales = "free")
  + geom_bar(width = 0.8)
  + theme_publication()
  + labs(x = "", y = "N SNPs", fill = "VAF status")
  + scale_fill_manual(values = c("VAF consistent" = "gray", 
                                 "VAF inconsistent" = "#8bad0d", 
                                 "Allele swapping" = "#86bfd5"),
                      breaks = c("VAF consistent", "VAF inconsistent", "Allele swapping"),
                      labels = c("VAF consistent" = "VAF coherent",
                                 "VAF inconsistent" = "VAF incoherent (VAF delta > 0.005)", 
                                 "Allele swapping" = "Allele swapping"))
  + geom_label_repel(
    aes(label= unlist(tapply(after_stat(count), list(..x.., ..PANEL..),
                                        function(a) paste(round(100*a/sum(a), 2), '%'))),

                    y = after_stat(count)),
    stat = "count",  direction = "x",
    #nudge_x = 0.2,
    #min.segment.length = 0.2,
    position=position_stack(vjust = 0.5),
    color = "white", vjust = 0.5,  size = 5,
    show.legend = F)
  + theme(#axis.ticks.x = element_blank(),
          plot.subtitle = element_text(size = 18),
          legend.direction = "horizontal",
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.ticks.y = element_blank())
  + coord_flip()

)
#patchwork::wrap_plots(a, b, nrow = 2, guides = "collect", axis = "collect", heights = c(6,2)) + plot_annotation(theme = theme(legend.position = "top", legend.title.position = "left"))
a + b + plot_layout(nrow = 2, guides = "collect", axis = "collect", heights = c(6,2)) + plot_annotation(theme = theme(legend.position = "top", legend.title.position = "left"))
```

![](Figure_1_S1_S2/figureS2E.png)<!-- -->

