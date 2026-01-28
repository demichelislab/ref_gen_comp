---
title: "Figure 3"
author: "Ilaria Cherchi"
date: "2025-11-12"
output: 
  html_document:
    keep_md: yes
---





# Figure 3

## Figure 3A

Sketch

## Figure 3B

Load summary of general calls:


``` r
suppressWarnings(load(paste0(rdata_path,"clinvar.RData")))
def <- de
def$cv <- ifelse(def$cv == "Other" & def$class == "Likely_benign", "Likely benign", as.character(def$cv))
def$cv <- ifelse(is.na(def$cv), "Other", as.character(def$cv))
levels(def$cv) <- c("Benign", "Likely benign", "Other", "Conflicting interpretations", "Uncertain significance", "Likely pathogenic", "Pathogenic")

load(paste0(rdata_path, "variant_calls_ClinVar.RData"))
calls <- left_join(calls %>% dplyr::select(-c(nc)) %>% mutate(ID = as.character(ID.y)), def, by = "ID")
calls <- calls %>% mutate(identifier = paste0(sample, ID))
calls <- calls %>% group_by(sample, ID) %>% dplyr::filter(n()==2) %>% ungroup() %>% 
  mutate(seq = ifelse(nc != "Equal", "Differences\nin proximity", ifelse(nc == "Equal", "Identical\nproximity", NA)))
```

Compare calls from different callers:


``` r
vaf_thr <- 0.05
cov_thr <- 50
call_strelka <- calls %>% filter(AF.indel >= vaf_thr & cov.indel > cov_thr &
                                   nchar(REF.indel) == 1 & nchar(ALT.indel) == 1 &
                                   !is.na(AF.indel), nc != "Equal") %>% 
  distinct(identifier, assembly, AF.indel, cv, QUAL.indel, cov.indel, nc)

call_strelka$method <- "Strelka2"
call_bcf <- calls %>% filter(AF.bcf >= vaf_thr & cov.bcf > cov_thr &
                               nchar(REF.bcf) == 1 & nchar(ALT.bcf) == 1 &
                                   !is.na(AF.bcf), nc != "Equal") %>% 
  distinct(identifier, assembly, AF.bcf, cv, QUAL.bcf, cov.bcf, nc)
call_bcf$method <- "BCFtools"

colnames(call_strelka) <- colnames(call_bcf) <- c("identifier", "assembly", "VAF", "cv", "QUAL", "cov", "nc", "method")
call_comparison <- rbind(call_strelka,  
                          call_bcf)

call_comparison %>% group_by(method, cv) %>% summarize(increment = 
        (length(assembly[assembly == "T2T"]) - length(assembly[assembly == "hg38"])) / length(assembly[assembly == "hg38"]) * 100)
```

```
## # A tibble: 14 × 3
## # Groups:   method [2]
##    method   cv                          increment
##    <chr>    <fct>                           <dbl>
##  1 BCFtools Benign                         -2.70 
##  2 BCFtools Likely benign                  -0.257
##  3 BCFtools Other                           2.95 
##  4 BCFtools Conflicting interpretations    15.3  
##  5 BCFtools Uncertain significance          1.73 
##  6 BCFtools Likely pathogenic              13.0  
##  7 BCFtools Pathogenic                     24.7  
##  8 Strelka2 Benign                        -14.7  
##  9 Strelka2 Likely benign                  -5.69 
## 10 Strelka2 Other                         -20.4  
## 11 Strelka2 Conflicting interpretations    13.7  
## 12 Strelka2 Uncertain significance         -0.440
## 13 Strelka2 Likely pathogenic               8.33 
## 14 Strelka2 Pathogenic                     22.6
```


``` r
df <- call_comparison %>% group_by(cv, method) %>% 
  summarise(label = scales::percent(
  -((length(identifier[assembly == "hg38"]) - length(identifier[assembly == "T2T"])) / 
      length(identifier[assembly == "hg38"])), 2),
  y = max(length(identifier[assembly == "T2T"]), length(identifier[assembly == "hg38"])),
   x = unique(method))
df$x <- rep(c(1:length(unique(df$cv))), 2)
df$x <- rep(c(1, 2), nrow(df)/2)

a <- (ggplot(call_comparison, 
       aes(y = method))
 + theme_publication()
 + geom_bar(position = "dodge", aes(fill = assembly))
 #+ facet_wrap(cv ~ ., scales = "free", strip.position = "top", nrow = 1)
 + facet_wrap(cv ~ ., scales = "free", strip.position = "right", ncol = 1)
  + scale_x_continuous(expand = expansion(mult = c(0, .1)))
 + scale_fill_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                       labels = c("hg38", "T2T-CHM13"))
  + theme(#axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
                      legend.position = "bottom",
          strip.text.y.right = element_text(angle = 0),
          legend.margin = margin(l = -.05, unit = "npc"))
+ labs(y = "", x = "N SNVs with VAF>=5% and coverage > 50", fill = "Assembly")
+ geom_text(data = df, aes(y=x, x=y, label=ifelse(grepl("-", label),label, ifelse(label != "0%", paste0("+", label), label))), 
                           vjust = 0.5 , hjust = -0.25, size = 4)
+ geom_segment(data = df, aes(y = x-0.45, yend = x+0.45, x = y+ y*0.005, xend =  y+ y*0.005), lwd = 0.5) 
)

b <- (ggplot(call_comparison, 
       aes(x = assembly, y = as.numeric(QUAL), color = assembly))
 + theme_publication()
 + geom_boxplot(show.legend = F, outlier.shape = NA)
 + facet_wrap(method ~ ., scales = "free", ncol = 2)
 + stat_compare_means( comparison = list(c("hg38", "T2T")), alternative = "less",
                        tip.length = 0, angle = 0, show.legend =  F)
 + theme(legend.position = "bottom", 
          panel.grid = element_blank(),
          #strip.text = element_text(size = 18),
          plot.title = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.subtitle = element_text(size = 16),
         plot.margin = margin(1,1,10,1)
         #legend.margin = margin(l = -.25, unit = "npc")
         # legend.margin = margin(-5,0,-5,0)
         ) 
 + geom_text( data = . %>% ungroup() %>% filter(!is.na(QUAL)) %>% 
                distinct(identifier, method, assembly) %>% dplyr::summarize(n = n(), .by = c("method", "assembly")),
                aes (y = 0, label = paste0 ("N = ", n)), angle = 0, size = 4, vjust = 1.25, show.legend = F) 
 + labs( y = "SNV Phred-scaled quality score \n(QUAL)", color = "Assembly")
 + scale_color_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                       labels = c("hg38", "T2T-CHM13"))
  )

(a + plot_spacer() + b + plot_layout(nrow = 1, guides = "collect", widths = c(5, 0.1,4)) + plot_annotation(theme =  theme(legend.position = "bottom")))
```

![](Figure_3_S6/figure3B.png)<!-- -->

## Figure 3C

Compare VAF from Strelka2 calls:


``` r
dek <- calls %>% filter(cov.indel > cov_thr & AF.indel >= vaf_thr &
                                   !is.na(AF.indel)) %>% 
  distinct(identifier, assembly, AF.indel, cv, QUAL.indel, cov.indel, nc, seq, hm_match, aln_hg38vsT2T, chr.x, pos.hg38, pos.T2T, ID.y, stretch)

dek <- dek %>% group_by(identifier) %>% dplyr::filter(n()==2) %>% ungroup()
p1 <- (ggplot(dek,
        aes(y = assembly, x = AF.indel))
  + geom_vline(xintercept = 0.5, color ="gray40", lwd= 1)
  + theme_publication()
  #+ geom_violin( aes(fill = assembly))
  + geom_boxplot(aes(color = assembly), outlier.shape = NA)
  + stat_compare_means(paired = TRUE, size = 3.5, vjust = 0, method = "wilcox.test", 
                       label.y = 0.83, 
                        comparison = list(c("hg38", "T2T")), alternative = "less",
                        tip.length = 0, angle = 0)
  + labs( x = "VAF", color = "Assembly version", shape = "Class")
  + scale_color_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                                  labels = c("hg38", "T2T-CHM13"))
+ geom_text( data = . %>% distinct(identifier, seq) %>% dplyr::summarize(n = n(), .by = "seq"),
                aes(x = 0.25, y = 1.5, label = paste0 ("N = ", n)), angle = 0, size = 4) 
  + facet_wrap( seq ~ ., ncol = 1, strip.position = "left")
    + theme(legend.position = "bottom", 
          panel.grid = element_blank(),
          #strip.text = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          plot.subtitle = element_text(size = 16)
         # legend.margin = margin(-5,0,-5,0)
         )
  + coord_cartesian(xlim = c(0.2, 0.9))
  + scale_x_continuous(breaks = c(0.2, 0.4, 0.5, 0.6, 0.8))
)
```


``` r
df_barplot <- dek %>% mutate(hp_status = case_when((hm_match == "F" | is.na(hm_match))  & 
                                                            aln_hg38vsT2T != "." & !is.na(stretch) ~ "Difference outside homopolymer", 
                                              hm_match == "T"  & !is.na(stretch) ~ "Difference in homopolymer",
                                              #aln_hg38vsT2T != "."  & !is.na(stretch) ~ "Variation in homopolymer",
                                              is.na(stretch) & aln_hg38vsT2T != "." ~ "Difference outside homopolymer",
                                              is.na(stretch)  & aln_hg38vsT2T == "." ~ "Identical proximity",
                                              !is.na(stretch) & aln_hg38vsT2T == "." ~ "Identical proximity", 
                                              .default = NA)) %>%
  mutate(stretch = ifelse(is.na(stretch), "No homopolymer stretch", stretch),
         hp_status = factor(hp_status, levels = c( "Difference in homopolymer",  "Difference outside homopolymer", "Identical proximity")))

df_barplot <- df_barplot %>% group_by(chr.x, pos.hg38, pos.T2T, hp_status) %>% dplyr::summarise(IDs = paste0(ID.y, collapse = "|"), 
                                                                                               
                 VAF_status = case_when( all((AF.indel[assembly == "T2T"] - AF.indel[assembly == "hg38"]) != 0)  ~ "Discordant VAF", 
                                         all(AF.indel[assembly == "T2T"] - AF.indel[assembly == "hg38"] == 0)  ~ "Equal VAF",
                                         .default = "Other")) %>% filter(VAF_status != "Other") %>%
  mutate(VAF_status = factor(VAF_status, levels = c( "Discordant VAF", "Equal VAF")))
  
group_pairs <- combn(unique(df_barplot$VAF_status), 2, simplify = FALSE)

p_value_table <- bind_rows(lapply(group_pairs, perform_chisq_test, 
                                  df = df_barplot, varA = "VAF_status", varB = "hp_status"))

p2 <- (ggplot(df_barplot, aes(x = VAF_status))
  + theme_publication()
  + geom_bar(aes(fill = hp_status), position = "fill")
  + geom_label(
     aes(fill = hp_status, label = scales::percent( after_stat(count) /
        tapply(after_stat(count), after_stat(x), sum)[as.character(after_stat(x))], accuracy = 0.1)), stat ="count",
position=position_fill(vjust = 0.5), color = "white", show.legend = F)
   #+ facet_grid(~stretch)
  + scale_fill_manual(values = c("Difference in homopolymer" = "#003665", "Difference outside homopolymer" = "#29bec3", "Identical proximity" = "gray60"), 
                      breaks = c("Identical proximity","Difference outside homopolymer", "Difference in homopolymer" ))
  + labs(y = "Proportion", x = "", fill = "Proximal sequence status (±150bp) for hg38vsT2T-CHM13")
  + stat_pvalue_manual(p_value_table ,
                           label = "p.value",  y.position = c(1), 
                       size = 3.7, tip.length = 0.25, bracket.shorten = 0.08, 
                       vjust = -0.5, bracket.nudge.y = 0.02, coord.flip = T)
  + theme(legend.position = "bottom", plot.caption = element_text(size = 9),
          strip.text.y = element_text(angle = 0),
          plot.subtitle = element_text(size = 16))
  + coord_flip()
)

#p2
```


``` r
p1 + p2 + plot_layout(nrow = 1,  widths = c(4,4), guides = "collect") &
   plot_annotation(theme = theme(legend.position = "bottom",  legend.margin = margin(0,0,0,0)))
```

![](Figure_3_S6/figure3C.png)<!-- -->

## Figure 3D

### BRCA1


``` r
id <- "125860"
width <- 75

pos_t2t <- as.numeric(def$pos.T2T[def$ID == id])
pos_hg38 <- as.numeric(def$pos.hg38[def$ID == id])
chrom <- def$chr[def$ID == id]
var_base <- def$alt.T2T[def$ID == id]
gene <- def$gene_info[def$ID == id]
```



``` r
lev <- compute_lev_distance(chrom, pos_hg38, pos_t2t, width)
lev
```

```
## [1] 0.8192771
```


``` r
exons <- brca1_exons
pp <- plot_blast(gene, exons)
plot_pp <- (pp
  + annotate(geom = "rect", xmin = pos_hg38 - unique(exons$gene_start),
             xmax =  pos_hg38 + 1 - unique(exons$gene_start),
             ymin = 0.3, ymax = 0.7, color = "#a12d69", alpha = 0.1)
  + annotate(geom = "rect", xmin = pos_hg38 - unique(exons$gene_start),
             xmax =  pos_hg38 + 1 - unique(exons$gene_start),
             ymin = 1.3, ymax = 1.6, color = "#a12d69", alpha = 0.1)
  + annotate(geom = "text", label = paste0("ClinVar ID: ",id),
             x = pos_hg38 - unique(exons$gene_start), y = 0.5,
             hjust = -0.05, color = "#a12d69", size = 5)
  + annotate(geom = "text", label = paste0("Levenshtein similarity (±75 bp): ", round(lev, 2)),
             x = pos_hg38 - unique(exons$gene_start), y = 0.35,
             hjust = -0.03, color = "#a12d69", size = 4.3)
  + theme(plot.margin = margin(5,5,10,0))
  + geom_hline(yintercept = 2.55, color = "gray30")
  + guides(fill="none")
  + coord_cartesian(clip = "off")
  + annotate(geom = "text", label = paste0(round((unique(exons$gene_end) - unique(exons$gene_start)) / 1e6, 2), " Mbp"),
             x = (unique(exons$gene_end) - unique(exons$gene_start)) / 2, y = 2.45, size = 3.5, fontface = "bold")
)
plot_pp
```

![](Figure_3_S6/figure3D_msabrca1.png)<!-- -->

### TP53


``` r
id <- "485033"
width <- 75
pos_t2t <- def$pos.T2T[def$ID == id]
pos_hg38 <- def$pos.hg38[def$ID == id]
chrom <- def$chr[def$ID == id]
var_base <- def$alt.T2T[def$ID == id]
gene <- def$gene_info[def$ID == id]
```


``` r
lev <- compute_lev_distance(chrom, pos_hg38, pos_t2t, width)
lev
```

```
## [1] 0.9235669
```


``` r
exons <- tp53_exons
pp <- plot_blast(gene, exons)
plot_pp <- (pp 
  + annotate(geom = "rect", xmin = pos_hg38 - unique(exons$gene_start), 
             xmax =  pos_hg38 + 1 - unique(exons$gene_start), 
             ymin = 0.3, ymax = 0.7, color = "#a12d69", alpha = 0.1) 
  + annotate(geom = "rect", xmin = pos_hg38 - unique(exons$gene_start), 
             xmax =  pos_hg38 + 1 - unique(exons$gene_start), 
             ymin = 1.3, ymax = 1.6, color = "#a12d69", alpha = 0.1) 
  + annotate(geom = "text", label = paste0("ClinVar ID: ",id), 
             x = pos_hg38 - unique(exons$gene_start), y = 0.55, 
             hjust = 1.05, color = "#a12d69", size = 5)
  + annotate(geom = "text", label = paste0("Levenshtein similarity (±75 bp): ", round(lev, 2)), 
             x = pos_hg38 - unique(exons$gene_start), y = 0.35, 
             hjust = 1.03, color = "#a12d69", size = 4.3)
  + theme(plot.margin = margin(5,5,10,0))
  + geom_hline(yintercept = 2.55, color = "gray30")
  + guides(fill="none")
  + coord_cartesian(clip = "off")
  + annotate(geom = "text", label = paste0(round((unique(exons$gene_end) - unique(exons$gene_start)) / 1e6, 2), " Mbp"), 
             x = (unique(exons$gene_end) - unique(exons$gene_start)) / 2, y = 2.45, size = 3.5, fontface = "bold")
)
plot_pp
```

![](Figure_3_S6/figure3D_msatp53.png)<!-- -->


# Supplementary 6

## Figure S6A


``` r
vaf_thr <- 0.05
cov_thr <- 50
qual_thr <- 200

calls_venn <- calls %>% filter(grepl("athogenic", cv))

hg38 <- list("Strelka2\n" = unique(calls_venn %>% filter(AF.indel >= vaf_thr & 
                                                                !is.na(AF.indel) &
                                                               cov.indel >= cov_thr & 
                                                                nchar(REF.indel) == 1 &   nchar(ALT.indel) == 1 &
                                                           QUAL.indel > qual_thr,
                                                                assembly == "hg38") %>% pull(identifier)),
          " Strelka2 noindel\n" = unique(calls_venn %>% filter(AF.noindel >=vaf_thr & 
                                                                     !is.na(AF.noindel) &
                                                                     cov.noindel >= cov_thr & 
                                                                 QUAL.noindel > qual_thr,
                                                                      nchar(REF.noindel) == 1 &   nchar(ALT.noindel) == 1 &
                                                                     assembly == "hg38") %>% pull(identifier)), 
          "BCFtools\n" = unique(calls_venn %>% filter(AF.bcf >= vaf_thr & 
                                                             !is.na(AF.bcf) &
                                                             cov.bcf >= cov_thr & 
                                                        QUAL.bcf > qual_thr,
                                                              nchar(REF.bcf) == 1 &   nchar(ALT.bcf) == 1 &
                                                           assembly == "hg38") %>% pull(identifier))
           )
palette <- c("#3e9cbf","#226258","#f2c43d","#f17c37","#f26d50", "#226258", "#356920", "#966729")
#show total N in group
names(hg38) <- paste0(names(hg38),"N = ",lapply(hg38,length))
#show numbers and percentages

a <- (ggvenn(hg38, show_elements = F, show_percentage = TRUE, label_sep = "\n", 
       fill_color = palette[1:length(hg38)],       
       stroke_linetype = 1, stroke_size = 0.5,
       set_name_size = 5,
       set_name_color = palette[1:length(hg38)],
       auto_scale = FALSE, text_size = 3.8, padding = 0.1) 
      + ggtitle("hg38") 
      + theme(plot.title = element_text(size = 16, face = "bold", hjust= 0.5, vjust = -1, margin = margin(b = 1, unit = "cm"))))

t2t <- list("Strelka2\n" = unique(calls_venn %>% filter(AF.indel >= vaf_thr & 
                                                                !is.na(AF.indel) &
                                                               cov.indel >= cov_thr & 
                                                                nchar(REF.indel) == 1 &   nchar(ALT.indel) == 1 &
                                                                assembly == "T2T") %>% pull(identifier)),
          "Strelka2 noindel\n" = unique(calls_venn %>% filter(AF.noindel >= vaf_thr & 
                                                                     !is.na(AF.noindel) &
                                                                     cov.noindel >= cov_thr & 
                                                                      nchar(REF.noindel) == 1 &   nchar(ALT.noindel) == 1 &
                                                                     assembly == "T2T") %>% pull(identifier)), 
          "BCFtools\n" = unique(calls_venn %>% filter(AF.bcf >= vaf_thr & 
                                                             !is.na(AF.bcf) &
                                                             cov.bcf >= cov_thr & 
                                                              nchar(REF.bcf) == 1 &  nchar(ALT.bcf) == 1 &
                                                           assembly == "T2T") %>% pull(identifier))
           )
palette <- c("#3e9cbf","#226258","#f2c43d","#f17c37","#f26d50", "#226258", "#356920", "#966729")
#show total N in group
names(t2t) <- paste0(names(t2t),"N = ",lapply(t2t,length))
#show numbers and percentages

b <- (ggvenn(t2t, show_elements = F, show_percentage = TRUE, label_sep = "\n", 
       fill_color = palette[1:length(t2t)],       
       stroke_linetype = 1, stroke_size = 0.5,
       set_name_size = 5,
       set_name_color = palette[1:length(t2t)],
       auto_scale = FALSE, text_size = 3.8, padding = 0.1) 
      + ggtitle("T2T-CHM13") 
      + theme(plot.title = element_text(size = 16, face = "bold", hjust= 0.5, vjust = -1, margin = margin(b = 1, unit = "cm"))))

#consensus <- Reduce(intersect, x)

x <- list("hg38 calls\n" = Reduce(union, hg38),
          "T2T calls\n" = Reduce(union, t2t))

layout <- 'AB'
a + b + plot_layout(design = layout)
```

![](Figure_3_S6/figureS6A.png)<!-- -->

## Supplementary Figure 6B


``` r
a <- (ggplot(calls %>% filter(AF.indel >= vaf_thr, cov.indel > cov_thr) %>% 
          group_by(sample, ID, seq) %>% 
            dplyr::filter(length(AF.indel[!is.na(AF.indel)]) ==2 & length(unique(assembly)) == 2) %>% ungroup() ,
        aes(y = assembly, x = AF.indel))
  + geom_vline(xintercept = 0.5, color ="gray40", lwd= 1)
  + theme_publication()
  #+ geom_violin( aes(fill = assembly))
  + geom_boxplot(aes(color = assembly), outlier.shape = NA)
  + stat_compare_means(paired = TRUE, size = 3.5, vjust = 0, method = "wilcox.test", 
                       label.y = 0.83, 
                        comparison = list(c("hg38", "T2T")), alternative = "less",
                        tip.length = 0, angle = 0)
  + labs( x = "VAF", color = "Assembly version", shape = "Class", title = "Strelka2 indel")
  + scale_color_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                       labels = c("hg38", "T2T-CHM13"))
  + facet_wrap( seq ~ ., ncol = 1, strip.position = "left")
  + geom_text( data = . %>% distinct(sample, ID.y, seq) %>% dplyr::summarize(n = n(), .by = "seq"),
                aes(x = 0.25, y = 1.5, label = paste0 ("N = ", n)), angle = 0, size = 4)
  + theme(legend.position = "bottom", 
          panel.grid = element_blank(),
          #strip.text = element_text(size = 18),
          plot.title = element_text(size = 13),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          plot.subtitle = element_text(size = 16)
         # legend.margin = margin(-5,0,-5,0)
         )
  + coord_cartesian(xlim = c(0.2, 0.9))
  + scale_x_continuous(breaks = c(0.2, 0.4, 0.5, 0.6, 0.8))
)  

b <- (ggplot(calls %>% filter(AF.noindel >= vaf_thr, cov.noindel > cov_thr) %>% 
          group_by(sample, ID, seq) %>% dplyr::filter(length(AF.noindel[!is.na(AF.noindel)])==2 & length(unique(assembly)) == 2) %>% ungroup(),
        aes(y = assembly, x = AF.noindel))
  + geom_vline(xintercept = 0.5, color ="gray40", lwd= 1)
  + theme_publication()
  #+ geom_violin( aes(fill = assembly))
  + geom_boxplot(aes(color = assembly), outlier.shape = NA)
  + stat_compare_means(paired = TRUE, size = 3.5, vjust = 0, method = "wilcox.test", 
                       label.y = 0.83, 
                        comparison = list(c("hg38", "T2T")), alternative = "less",
                        tip.length = 0, angle = 0)
  + labs( x = "VAF", color = "Assembly version", shape = "Class", title = "Strelka2 no indel")
  + scale_color_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                       labels = c("hg38", "T2T-CHM13"))
  + facet_wrap( seq ~ ., ncol = 1, strip.position = "left")
  + geom_text( data = . %>% distinct(sample, ID.y, seq) %>% dplyr::summarize(n = n(), .by = "seq"),
                aes(x = 0.25, y = 1.5, label = paste0 ("N = ", n)), angle = 0, size = 4)
  + theme(legend.position = "bottom", 
          panel.grid = element_blank(),
          #strip.text = element_text(size = 18),
          plot.title = element_text(size = 13),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          plot.subtitle = element_text(size = 16)
         # legend.margin = margin(-5,0,-5,0)
         )
  + coord_cartesian(xlim = c(0.2, 0.9))
  + scale_x_continuous(breaks = c(0.2, 0.4, 0.5, 0.6, 0.8))
)  
 
c <- (ggplot(calls %>% filter(AF.bcf >= vaf_thr, cov.bcf > cov_thr) %>% 
          group_by(sample, ID, seq) %>% dplyr::filter(length(AF.bcf[!is.na(AF.bcf)]) == 2 & length(unique(assembly)) == 2) %>% ungroup() ,
        aes(y = assembly, x = AF.bcf))
  + theme_publication()
  + geom_vline(xintercept = 0.5, color ="gray40", lwd= 1)
  #+ geom_violin( aes(fill = assembly))
  + geom_boxplot(aes(color = assembly), outlier.shape = NA)
  + stat_compare_means(paired = TRUE, size = 3.5, vjust = 0, method = "wilcox.test", 
                       label.y = 0.83, 
                        comparison = list(c("hg38", "T2T")), alternative = "less",
                        tip.length = 0, angle = 0)
  + labs( x = "VAF", color = "Assembly version", shape = "Class", title = "BCFtools")
  + scale_color_manual(values = c( "T2T" = "#ED7D31",  "hg38" = "#7f061b"),
                       labels = c("hg38", "T2T-CHM13"))
  + facet_wrap( seq ~ ., ncol = 1, strip.position = "left")
  + geom_text( data = . %>% distinct(sample, ID.y, seq) %>% dplyr::summarize(n = n(), .by = "seq"),
                aes(x = 0.25, y = 1.5, label = paste0 ("N = ", n)), angle = 0, size = 4)
  + theme(legend.position = "bottom", 
          panel.grid = element_blank(),
          plot.title = element_text(size = 13),
          #strip.text = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
           strip.text.y.left = element_text(angle = 0),
          #strip.background = element_blank(),
          #strip.text.y = element_blank(),
          plot.subtitle = element_text(size = 16)
         # legend.margin = margin(-5,0,-5,0)
         )
  + coord_cartesian(xlim = c(0.2, 0.9))
  + scale_x_continuous(breaks = c(0.2, 0.4, 0.5, 0.6, 0.8))
)  

a + b + c + plot_layout(ncol = 1, guides = "collect", axis_titles = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))
```

![](Figure_3_S6/figureS6B.png)<!-- -->

## Supplementary Figure 6C


``` r
brca1_df <- calls %>% filter(gene == "BRCA1", !is.na(cov.bcf))

genen <- "BRCA1"

## Annotate with exons
exons_hg38 <- genes_hg38 %>% filter(hg38.kgXref.geneSymbol == genen) %>% filter(`#hg38.knownGene.name` == "ENST00000352993.7") %>%
  separate_longer_delim(cols = c(hg38.knownGene.exonStarts, hg38.knownGene.exonEnds), delim = ",") 

exons_hg38 <- exons_hg38 %>% mutate(id_exon = 1:nrow(exons_hg38)) %>%
  dplyr::rename(chr = `hg38.knownGene.chrom`, gene = hg38.kgXref.geneSymbol, gene_start = hg38.knownGene.txStart,
                gene_end = hg38.knownGene.txEnd, start = hg38.knownGene.exonStarts, end = hg38.knownGene.exonEnds) %>%
  mutate(start.hg38 = as.numeric(start), end.hg38 = as.numeric(end)) %>%
  dplyr::select(chr,gene,gene_start,gene_end,start.hg38,end.hg38,id_exon) %>% filter(!is.na(start.hg38))

inpt <- brca1_df %>% dplyr::filter(cv == "Pathogenic") %>% group_by(pos.hg38) %>% 
  summarize(med_cov = median(cov.bcf[assembly == "T2T"] - cov.bcf[assembly == "hg38"]),
            cov_perc = median((cov.bcf[assembly == "T2T"] - cov.bcf[assembly == "hg38"]) / cov.bcf[assembly == "T2T"]) * 100)

inpt$id <- NA  # Initialize

for (i in seq_len(nrow(inpt))) {
  pos_i <- inpt$pos.hg38[i]
  region <- exons_hg38$id_exon[exons_hg38$start <= pos_i & exons_hg38$end >= pos_i]
  if (length(region) == 1) {
    inpt$id[i] <- region
  } else if (length(region) > 1) {
    # If overlapping regions, assign the first one or handle accordingly
    inpt$id[i] <- region[1]
  }
}

median(inpt %>% dplyr::filter(id == 17) %>% pull(cov_perc))
```

```
## [1] 0.422833
```

``` r
a <- (ggplot(inpt, 
        aes(x = pos.hg38 / 1e6, y = cov_perc))
   + annotate("segment",x = exons_hg38$start.hg38 / 1e6, xend = exons_hg38$end.hg38 / 1e6, y = -0.45, yend = -0.45, lwd = 2, colour ="gray80") 
   + annotate("text", x = ((exons_hg38$start.hg38 + exons_hg38$end.hg38) / 2) / 1e6, y = -0.7, 
            label = exons_hg38$id_exon, color = "gray60", hjust = 0.5)
  + geom_hline(yintercept = 0)
  + geom_point()
  + facet_zoom( xy = id == "17", horizontal = F, show.area = F)
  + theme_publication()
  + labs(x = "Position hg38 (Mb)", y = "Percentage of median absolute change in\ncoverage (T2T-CHM13 - hg38)\nfor BRCA1 pathogenic variants", subtitle = "BRCA1")
)
b <- (ggplot(brca1_df %>% dplyr::filter(cv == "Pathogenic", pos.hg38 %in% inpt$pos.hg38[inpt$id == 17]), aes(x = assembly, y = cov.bcf, fill = assembly))
       + geom_boxplot()
       + theme_publication()
       + stat_compare_means(paired = T, label.x = 1.5, comparisons = list(c("hg38", "T2T")),
                            alternative = "less")
       + scale_fill_manual(values = c("T2T" = "#ED7D31",  "hg38" = "#7f061b"), labels = c("hg38", "T2T-CHM13"))
       + labs(x = "Assembly version", y = "Coverage of exon 17\npathogenic variants", fill = "")
       + theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
               legend.position = "bottom", legend.margin = margin(l = -.05, unit = "npc"))
     )
lt <- "A#
AB"
a + b + plot_layout(design = lt, widths = c(10,1), heights = c(5,7))
```

![](Figure_3_S6/figureS6C.png)<!-- -->

## Supplementary Figure 6D


``` r
actionable_genes <- fread(paste0(utilities_path, "Actionability_AllData_v14_GRCh37.tsv")) %>% distinct(GENE) %>% pull(GENE)

sum_calls <- calls %>% group_by(sample, ID) %>% filter(n() == 2) %>% group_by(ID, nc, sample, gene) %>%
  summarize(delta_cov = cov.bcf[assembly == "T2T"] - cov.bcf[assembly == "hg38"],
            delta_VAF = AF.bcf[assembly == "T2T"] - AF.bcf[assembly == "hg38"])

 
sc <- left_join(sum_calls, def %>% distinct(cv, ID)) %>% filter(cv %in% c("Pathogenic", "Likely pathogenic")) %>% 
  ungroup() %>% summarize(tot = n(),
                                      perc_cov_T2T = length(delta_cov[delta_cov >0]) / tot,
                                      perc_VAF_T2T = length(delta_VAF[delta_VAF >0]) / tot,
                                      N = length(unique(ID)), .by = "gene")

(ggplot(sc %>% filter(N > 100) %>% mutate(act = ifelse(gene %in% actionable_genes, "Yes", "No")), 
        aes(x = perc_cov_T2T, y = perc_VAF_T2T))
  + geom_point(aes( size = N, color = act))
  + geom_label_repel(aes(label = gene))
  + theme_publication()
  + scale_color_manual(values = c("Yes" = "#a2315d", "No" = "gray50"))
  + labs(x = "Proportion of pathogenic / likely pathogenic IDs with\ncoverage increase in T2T-CHM13", 
         y = "Propotion of pathogenic / likely pathogenic IDs\nwith VAF increase in T2T-CHM13", 
         size = "N patho IDs",
         color = "Actionable gene")
  + stat_cor(size = 5)
)
```

![](Figure_3_S6/figureS6D.png)<!-- -->

## Supplementary Figure 6E


``` r
load("/shares/CIBIO-Storage/CO/SPICE/personal/ilaria/long_1kgp/dfplot_cons_1kgp_hg38_t2t_patho.RData")

count_df <- df_plot %>% filter(ID == 100326, assembly == "hg38") %>% 
  group_by(`Superpopulation code`, ID, AS, assembly) %>%
  dplyr::count() %>%
  mutate(y = 0)  # set y-position for label

(ggplot(df_plot, 
        #%>% filter(ID == 100326), 
        aes(y = value, x = `Superpopulation code`, fill = assembly))
  + theme_bw()
  + geom_violin()
  + geom_boxplot(width = 0.05, position = position_dodge(0.9), color = "gray90", show.legend = F,
                 outlier.shape = NA)
  + facet_grid( . ~ `Superpopulation code`, scales = "free_x")
  + scale_fill_manual(values = c("hg19" = "#a2315d", "t2t" = "#ED7D31",  "hg38" = "#7f061b"),
                      labels = c("hg38", "T2T-CHM13"))
  + theme_publication()
  + labs(y = "Levenshtein distance\nread consensus vs reference genome", 
         x = "Assembly",
         fill = "Assembly version")
  + theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
   + stat_compare_means(paired = TRUE,  vjust = -0.5, method = "wilcox.test", size = 4,
                        #comparison = list(c("hg38", "t2t")),
                        method.args = list(alternative = "less"), 
                        tip.length = 0)
  + geom_text(data = count_df, aes(x = `Superpopulation code`, y = 0.1, 
                                   label = paste0("N individuals = ", n)),
              inherit.aes = FALSE, vjust = 1, size = 4)
  + scale_y_continuous(
    trans = focus_above_threshold_trans(threshold = 0.98),
    breaks = c(0.2, 0.5, 0.98, 1),
    #limits = c(0.7, 1)
  ) # + scale_fill_manual(values = palette_population)
)
```

![](Figure_3_S6/figureS6E.png)<!-- -->
