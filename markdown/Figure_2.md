---
title: "Figure 2"
author: "Ilaria Cherchi"
date: "2025-02-05"
output: 
  html_document:
    keep_md: yes
---






# Figure 2

## Figure 2A
Left side: sketch

For right side load data:


``` r
all_tracks <- fread(paste0(supp_table_path, "Supplementary_table_4.csv"))
all_tracks <- all_tracks %>% 
  mutate(Symbol = factor(Symbol, levels = rev(c("Match", "Mismatch", "Deletion", "Insertion")))) %>% 
  filter(Symbol != "Match") 
```



``` r
fig_f_df <- all_tracks %>% filter(comparison == "hg38vsT2T") %>%
          mutate(total_bases = (EVENT/prop_bases_event),
                 perc_event_num = as.numeric(sub("%", "", perc_peaks_event)))
fig_f_df$comparison <- "hg38vsT2T-CHM13"
```


``` r
(ggplot(fig_f_df, 
        aes( y = Symbol, x = perc_event_num, fill = Symbol))
  + geom_col()
  + scale_fill_manual(values = palette_pairwise,
                      breaks = c("Match", "Mismatch", "Deletion", "Insertion"))
  + theme_publication()
  + facet_grid(track ~ comparison,  switch = "y", space = "free_x")
  + labs (y = "", x = "Sequences with\nat least one modification (%)", fill = "Alignment outcome")
  + geom_label(data = . %>% filter(Symbol == "Mismatch" ), 
               aes(label = perc_peaks_event ), hjust = 1.1, size = 10/.pt, show.legend = F)
  + geom_label(data = . %>% filter(Symbol != "Mismatch" ), 
               aes(label = perc_peaks_event ), hjust = -0.1, size = 10/.pt, show.legend = F)
  + theme(legend.position = "left", legend.justification = "bottom", 
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
)
```

![](Figure_2_S5/figure2A.png)<!-- -->


## Figure 2B
Load data:


``` r
load(paste0(rdata_path, "clinvar.RData"))
def <- def %>% mutate(hp_status = case_when(
  (hm_match == "F" | is.na(hm_match)) & aln_hg38vsT2T != "." & !is.na(stretch) ~ 
    "Difference outside homopolymer", 
  hm_match == "T"  & !is.na(stretch) ~ "Difference in homopolymer",
  is.na(stretch) & aln_hg38vsT2T != "." ~ "Difference outside homopolymer", 
  is.na(stretch)  & aln_hg38vsT2T == "." ~ "Identical proximity",
  !is.na(stretch) & aln_hg38vsT2T == "." ~ "Identical proximity", 
  .default = NA)) %>%
  mutate(stretch = ifelse(is.na(stretch), "No homopolymer stretch", stretch),
         hp_status = factor(hp_status, 
  levels = c( "Difference in homopolymer",  
              "Difference outside homopolymer", 
              "Identical proximity")))

def_barplot <- def %>% distinct(ID, nc, cv, stretch, hm_match, 
                                aln_hg38vsT2T, pos.hg38, pos.T2T, chr, gene_info, hp_status,
                                disease, st_dbsnp, type,
                                var_type, CLNREVSTAT) %>% mutate(ID.y = as.numeric(ID)) 

def_barplot$dep <- sub(".*\\((\\d+)\\).*", "\\1", def_barplot$aln_hg38vsT2T)
```


``` r
def_barplot$cv <- factor(def_barplot$cv, levels =  
              c("Benign", "Likely benign", "Other", "Conflicting interpretations", 
                "Uncertain significance", "Likely pathogenic", "Pathogenic"))

dat <- def_barplot %>% 
  mutate(nc = ifelse(hp_status == "Difference in homopolymer", 
                                          "Difference in homopolymer", as.character(nc))) %>%
  group_by(cv) %>% mutate(tot = n()) %>% ungroup() %>% 
  summarize(perc = n() / unique(tot), n= n(), tot = unique(tot), .by = c("cv", "nc"))

dat$nc <- factor(dat$nc, levels = rev(names(palette_seq_mod)))

(ggplot(dat %>% filter(nc != "Equal"), aes( y = cv, x = perc))
              + geom_col(aes(fill = nc))
              + labs(y = "ClinVar classification", 
                     x = "Proportion of ClinVar IDs with differences\nin their proximal sequence", 
                     fill = "Proximal sequence status (±150bp) for hg38vsT2T-CHM13")
              + theme_publication()
              + theme(panel.grid = element_blank(), legend.position = "bottom")
              + scale_fill_manual(values = palette_seq_mod, breaks = names(palette_seq_mod))
              + theme(legend.title = element_text(size = 14),
                      legend.margin = margin(0,10,0,-120))
              + guides(fill = guide_legend(nrow = 1, title.position = "top"))
              + scale_x_continuous(labels = scales::percent, limits = c(0,0.4))
              + geom_text(data = . %>% group_by(cv) %>% 
                            summarize(perc = sum(perc), tot = unique(tot), n = sum(n)),
                          aes(x = perc, 
                              label = paste0(scales::percent(round(perc, 2)),
                              " (", n, "/", tot, ")")
                                             ),
                          vjust = 0.5, hjust = -0.15)
              + scale_x_continuous(expand = expansion(mult = c(0, .45)))
              + geom_label(data = . %>% filter(round(perc, 2)*100 > 0), 
                           aes(label = paste0(#"N=",n, "\n", 
                             scales::percent(round(perc, 2))) , x = perc, fill = nc), 
                                 position=position_stack(vjust = 0.5), 
                                 color = "white", size = 3, 
                           show.legend = F,  segment.color = "white")
)
```

![](Figure_2_S5/figure2B.png)<!-- -->

## Figure 2C

Extract percentage of difference in pathogenic:

``` r
def %>% mutate(ab = ifelse(nc == "Equal", "Equal", "Sequence modification")) %>%
  group_by(cv,ab) %>%  
  summarise(n = n(), .groups = 'drop') %>% 
  group_by(cv) %>%  
  mutate(percentage =  scales::percent((n / sum(n)))) %>% 
  ungroup() %>% filter(cv=="Pathogenic")
```

```
## # A tibble: 2 × 4
##   cv         ab                        n percentage
##   <chr>      <chr>                 <int> <chr>     
## 1 Pathogenic Equal                 91984 85%       
## 2 Pathogenic Sequence modification 15789 15%
```

``` r
patho_diff <- def %>% filter(cv %in% c("Pathogenic")) 
```

Load COSMIC actionable genes:


``` r
actionable_genes <- fread(paste0(utilities_path, "Actionability_AllData_v14_GRCh37.tsv")) %>% 
  distinct(GENE) %>% pull(GENE)

length(actionable_genes[actionable_genes %in% patho_diff$gene_info])
```

```
## [1] 195
```

Compute pathogenic percentages per gene:


``` r
patho_focus <- patho_diff %>% filter(#var_type == "SNV", 
  gene_info %in% actionable_genes)
patho_actio <- patho_focus %>% group_by(gene_info) %>% 
                        summarize(n_IDs = n(), 
                                  n_mod = length(nc[nc != "Equal" & nc != "Allele swapping" & nc != "Reference mismatch"]),
                        prop_diff = n_mod / n_IDs)

patho_actio %>% filter(n_IDs > 10) %>% arrange(desc(prop_diff))
```

```
## # A tibble: 110 × 4
##    gene_info n_IDs n_mod prop_diff
##    <chr>     <int> <int>     <dbl>
##  1 HRAS         29    17     0.586
##  2 CSF1R        50    28     0.56 
##  3 RECQL4      144    78     0.542
##  4 FLT4         15     8     0.533
##  5 NRAS         23    11     0.478
##  6 FANCA       371   152     0.410
##  7 BRCA1      2960  1183     0.400
##  8 STK11       168    62     0.369
##  9 TP53        451   166     0.368
## 10 PAX3         51    18     0.353
## # ℹ 100 more rows
```

Build and annotate barplot:

``` r
df_barplot_patho <- patho_actio %>% 
          filter(n_IDs > 10, prop_diff > 0.05, n_mod > 20) %>%
          mutate(div = case_when(prop_diff > 0.2 ~ "a",
                                 #prop_mism > 0.30 ~ "b",
                                 .default = "c"))
gene_annot <- annotate_gene_list(df_barplot_patho %>% pull(gene_info))
df_barplot_patho <- left_join(df_barplot_patho, gene_annot %>% distinct(SYMBOL, PathwayName) %>%
                                dplyr::rename(gene_info = SYMBOL, pathway = PathwayName))

df_barplot_patho <- df_barplot_patho %>% ungroup() %>% summarize(n=n(), 
path = case_when(any(grepl("cancer", pathway)) ~ "Pathway in cancer", 
                 .default = head(gsub(" -.*", "", pathway), 1)),
                         .by = c("gene_info", "prop_diff" , "n_mod", "div"))
```


``` r
(ggplot(df_barplot_patho, 
 aes(x = reorder(gene_info, - prop_diff), 
     y = round(prop_diff, 2)*100))
 + geom_col(aes(fill = path))
 + theme_publication()
 + geom_label(aes(label = paste0("N=", n_mod) , y = 0), size = 3,
              angle = 90, hjust = -0.1, fill = "transparent", show.legend = F, color = "white")     
 + facet_wrap(div ~ ., scales = "free", ncol = 1 )
 + labs(x = "Actionable genes included in ClinVar", 
        y = "Pathogenic SNVs with differences\nin their proximal sequence (%)",
        fill = "KEGG pathway hits")
  + scale_fill_manual(values = 
                        c("#f7b6d3", "#d6a9c8", "#9467bd",
                          "#c5b0d5", "#aec7e8", "#4a6d88", "#1f77b4"), 
                    na.value = "gray80")
 + scale_color_manual(values =
                        c("#f7b6d3", "#d6a9c8", "#9467bd", 
                          "#c5b0d5", "#aec7e8", "#4a6d88", "#1f77b4"),
                      na.value = "gray80")
  + geom_text(aes(label = scales::percent(prop_diff, 2)), size = 4, vjust = -0.5)
  + scale_y_continuous(expand = expansion(mult = c(0, .1)))
 + theme(
         strip.background = element_blank(),
         strip.text = element_blank(),
         legend.position = "bottom",
         legend.title.position = "top",
         legend.margin = margin(0,0,0,-20),
         axis.text.x = element_text(size = 11, face = "italic"),
        axis.ticks.x = element_blank()
         )
  + guides(fill = guide_legend(nrow = 3))
)
```

![](Figure_2_S5/figure2C.png)<!-- -->

## Figure 2D


``` r
def_ds <- def_barplot %>% filter(!is.na(disease), 
                         !is.na(nc), 
                         disease != "not specified", cv %in% c("Pathogenic", "Likely pathogenic")) 
tp <- def_ds %>% 
   distinct(hp_status, disease, ID, cv) %>% group_by(disease) %>% 
  mutate(tot = n()) %>% ungroup() %>% group_by(hp_status, disease) %>% 
  summarize(N = n(), perc = n()/unique(tot), tot = unique(tot))
tp$disease <- ifelse(grepl("Neurodevelopmental", tp$disease),
                     "Neurodevelopmental disorder with hypotonia,\nfacial dysmorphism, and brain abnormalities",
                     tp$disease)
tp$disease <- ifelse(grepl("17-alpha", tp$disease),
                     "17-alpha-hydroxylase/17,20-lyase\ndeficiency, combined complete",
                     tp$disease)

tp <- tp %>% filter(N>1)
tp <- tp %>% filter(!(hp_status == "Identical proximity" & perc == 1) & tot > 5 & hp_status != "Identical proximity")
order <- tp %>% ungroup() %>% group_by(disease) %>% 
  summarize(
   perc_var = sum(perc[hp_status != "Identical proximity"]),
   perc_equal = ifelse( "Identical proximity" %in% hp_status,
                        perc[hp_status == "Identical proximity"], 0)) %>%
   arrange(perc_var) %>% pull(disease)

tp <- tp %>% mutate(disease = factor(disease, levels = order))
```


``` r
(ggplot(tp, 
        aes( x = perc, y = disease))
  + geom_col( aes(fill = hp_status))
  + theme_publication()
  + scale_fill_manual(values = palette_homopolymer, 
                      breaks = c("Identical proximity","Difference outside homopolymer", "Difference in homopolymer" ))
  + scale_x_continuous(limits = c(-0.1, 0.63),expand = c(0, 0))
  + labs(x = "Proportion of pathogenic and likely pathogenic variants\nwith differences in their proximal sequence", 
         y = "", fill = "Proximal sequence status (±150bp) for hg38vsT2T-CHM13")
   + geom_label(aes(label = paste0("N=",N) , x = perc, fill = hp_status), 
                                 position=position_stack(vjust = 0.5), 
                                 color = "white", size = 3, 
                                  show.legend = F,  segment.color = "white")
 + geom_text(data = . %>% group_by(disease) %>% summarize(perc = sum(perc)),
                          aes(x = perc, 
                              label = scales::percent(round(perc, 2), 1)),
                          vjust = 0.5, hjust = -0.5, size = 4)
 + scale_x_continuous(expand =  expansion(mult = c(0, .08)))
   + theme(
          legend.position = "bottom",
          legend.margin = margin(l = -.45, unit = "npc"))
  )
```

![](Figure_2_S5/figure2D.png)<!-- -->

## Figure 2E
Load mutation calls from cBioportal for MSK-IMPACT 50k and GDC TCGA cohorts overlapping ClinVar annotated sites:

``` r
load("../cbioportal/tcga_clinvar.RData")
clinical_tcga <- fread("../cbioportal/tcga_gdc_2025_clinical_sample.txt", fill = T)

load("../cbioportal/msk_clinvar.RData")
clinical_msk <- fread("../cbioportal/msk_impact_50k_2026_clinical_sample.txt", 
                      skip = 4 )
```


``` r
tot <- length(unique(clinical_tcga$PATIENT_ID))
a <- tcga %>% 
  group_by(ID, nc, Hugo_Symbol) %>% dplyr::count()
a <- a %>% mutate(freq = (n/tot) * 100)
#length(unique(a$ID))
tot <- length(unique(clinical_msk$PATIENT_ID))
#tot <- length(unique(msk$Tumor_Sample_Barcode))
b <- msk %>% 
  group_by(ID, nc, Hugo_Symbol) %>% dplyr::count()
b <- b %>% mutate(freq = (n/tot) * 100)

comp <- bind_rows(list("TCGA GDC 2025" = a,
                       "MSK-IMPACT 50k 2026" = b),
                  .id = "proj")

(ggplot(comp, 
        aes(x = proj, fill = nc))
  + geom_bar(data = . %>% filter(nc != "Equal"))
  + geom_label_repel(data = . %>% group_by(proj) %>% mutate(tot = n()) %>%
                 ungroup() %>% group_by(proj, nc) %>% dplyr::summarize(
                   n=n(),
                   perc = n()/unique(tot)) %>% filter(nc != "Equal"),
               position = position_stack(0.5),
               aes(y=n, label = scales::percent(perc, 0.01)),
               show.legend = F,
               direction = "y", color = "white")
  + theme_publication()
  + labs(x = "Cohort", 
         y = "N ClinVar IDs\noverlapping with mutation calls",
         fill = "Proximal sequence status\n(±150bp) for hg38vsT2T-CHM13")
  + theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
          legend.position = "bottom",
          legend.title.position = "top",
          legend.margin = margin(l = -.35, unit = "npc")
          )
   + scale_fill_manual(values = palette_seq_mod)
  + guides(fill = guide_legend(nrow = 3))
)
```

![](Figure_2_S5/figure2E.png)<!-- -->


``` r
comp %>% group_by(proj) %>% mutate(tot = n()) %>%
                 ungroup() %>% group_by(proj,nc) %>% dplyr::summarize(
                   n=n(), perc = n()/unique(tot)) %>% filter(nc != "Equal") %>%
  ungroup() %>% group_by(proj) %>% dplyr::summarize(perc = sum(perc))
```

```
## # A tibble: 2 × 2
##   proj                 perc
##   <chr>               <dbl>
## 1 MSK-IMPACT 50k 2026 0.142
## 2 TCGA GDC 2025       0.142
```


``` r
tb <- comp %>% filter(nc != "Equal") %>% dplyr::select(-c(n)) %>% mutate(freq = round(freq,3)) %>% 
  arrange(desc(freq)) %>% pivot_wider(names_from = proj, values_from = freq) 
write.table(tb, paste0(supp_table_path, "Supplementary_table_7.csv"),
            sep = ",", row.names = F)
```



``` r
tot <- length(unique(clinical_tcga$PATIENT_ID))
#tot <- length(unique(tcga$case_id))
a <- tcga %>% 
  group_by(ID, nc, cv, Hugo_Symbol) %>% dplyr::count()
a <- a %>% mutate(freq = (n/tot) * 100)

tot <- length(unique(clinical_msk$PATIENT_ID))
#tot <- length(unique(msk$Tumor_Sample_Barcode))
b <- msk %>% 
  group_by(ID, nc, cv, Hugo_Symbol) %>% dplyr::count()
b <- b %>% mutate(freq = (n/tot) * 100)

comp <- bind_rows(list("TCGA GDC 2025" = a,
                       "MSK-IMPACT 50k 2026" = b),
                  .id = "proj")

comp$cv <- factor(comp$cv,levels =  c("Benign", "Likely benign", "Other",
                                      "Conflicting interpretations", "Uncertain significance", 
                                      "Likely pathogenic", "Pathogenic"))

(ggplot(comp, 
        aes(x = cv, fill = nc))
  + geom_bar(data = . %>% filter(nc != "Equal"))
  + geom_label_repel(data = . %>% group_by(proj, cv) %>% mutate(tot = n()) %>%
                 ungroup() %>% group_by(proj, cv, nc) %>% dplyr::summarize(
                   n=n(),
                   perc = n()/unique(tot)) %>% filter(nc != "Equal"),
               position = position_stack(0.5),
               aes(y=n, label = scales::percent(perc, 0.01)),
               show.legend = F,
               direction = "y", color = "white")
  + theme_publication()
  + labs(x = "ClinVar classification", 
         y = "N ClinVar IDs\noverlapping with mutation calls",
         fill = "Proximal sequence status (±150bp) for hg38vsT2T-CHM13")
  + facet_grid(proj ~ ., scales = "free_y")
  + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "bottom",
          legend.title.position = "top"
          )
   + scale_fill_manual(values = palette_seq_mod)
)
```

![](Figure_2_S5/figure_rev3A.png)<!-- -->

Check overlap among mutation calls between the two cohorts:

``` r
int <- full_join(a %>% ungroup() %>% distinct(ID) %>% mutate(proj = "TCGA"), 
          b %>% ungroup() %>% distinct(ID) %>% mutate(proj = "MSK"), by = "ID") %>%
  mutate(cl = ifelse(!is.na(proj.x) & is.na(proj.y), "TCGA",
                      ifelse(is.na(proj.x) & !is.na(proj.y), "MSK", "Both")) )
table(int$cl)
```

```
## 
##  Both   MSK  TCGA 
##  1358 21010  7462
```


``` r
e <- comp %>% filter(nc != "Equal") %>% filter(n>1) %>%
  filter(proj == "TCGA GDC 2025")

e <- e %>% 
          group_by(Hugo_Symbol) %>% 
          filter(n()>1) %>% 
  filter(!grepl("enign", cv))
(ggplot(e, 
        aes(x = freq, y = cv))
  + geom_violin(color = "gray70")
   + geom_beeswarm(aes( color = Hugo_Symbol))
   + theme_publication()
   + facet_wrap(proj~ .)
   + labs(x = "Frequency of mutations overlapping\nClinVar IDs (%)", 
         y = "ClinVar classification",
         color = "Gene")
   + theme(legend.position = "bottom",
           legend.text= element_text(face = "italic"),
            legend.justification = c(0,1),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -100, unit = "pt"))
  + geom_label_repel(data = . %>% filter(freq > 0.05),
                     aes(label = ID), show.legend = F)
   + scale_color_manual(values = paletteer::paletteer_d("palettetown::blastoise"))
   + guides(color = guide_legend(nrow = 2, title.position = "left"))
)
```

![](Figure_2_S5/figure_rev3B.png)<!-- -->


# Supplementary 5

## Supplementary 5A


``` r
dat <- def_barplot %>% group_by(cv) %>% mutate(hp_class = 
          case_when(hp_status == "Difference in homopolymer" ~ "Difference in homopolymer",
                    stretch == "Homopolymer stretch" ~ "Homopolymer",
                    .default = "No homopolymer")) %>% mutate(tot = n()) %>% ungroup() %>% 
  summarize(perc = n() / unique(tot), n= n(), tot = unique(tot), .by = c("cv", "hp_class"))

(ggplot(dat, aes( y = cv, x = perc))
              + geom_col(aes(fill = hp_class))
              + labs(y = "ClinVar classification", 
                     x = "Proportion of ClinVar IDs with differences\nin their proximal sequence for hg38vsT2T-CHM13", 
                     fill = "")
              + theme_publication()
              + theme(panel.grid = element_blank(), legend.position = "bottom")
              + scale_fill_manual(values = c("Difference in homopolymer" = "#003665", 
                                             "Homopolymer" = "#8bad0d", "No homopolymer" = "gray60"), 
                      breaks = c("No homopolymer","Homopolymer", "Difference in homopolymer" ))
 
              + theme(#axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
                       legend.title = element_text(size = 14),
                      legend.margin = margin(0,10,0,-120))
              + guides(fill = guide_legend(nrow = 1, title.position = "left"))
              + scale_x_continuous(expand =  expansion(mult = c(0, .15)))
          
  + geom_label(aes(label = 
                     paste0("N=",n, "\n", scales::percent(round(perc, 2), 1)) , 
                   x = perc, fill = hp_class), 
                                 position=position_stack(vjust = 0.5), 
                                 color = "white", size = 3, 
                           #direction = "x",
                            #      max.overlap = 0,
                                  show.legend = F,  segment.color = "white")
              
)
```

![](Figure_2_S5/figureS5A.png)<!-- -->


## Supplementary 5B

Load ClinVar annotation for hg19vsT2T and hg19vsgh38



``` r
def_hg19vsT2T$cv <- ifelse(def_hg19vsT2T$cv == "Other" & def_hg19vsT2T$class == "Likely_benign", "Likely benign", as.character(def_hg19vsT2T$cv))
def_hg19vsT2T$cv <- ifelse(is.na(def_hg19vsT2T$cv), "Other", as.character(def_hg19vsT2T$cv))
def_hg19vsT2T$cv <- factor(def_hg19vsT2T$cv, levels =  c("Benign", "Likely benign", "Other", "Conflicting interpretations", "Uncertain significance", "Likely pathogenic", "Pathogenic"))

def_hg19vshg38$cv <- ifelse(def_hg19vshg38$cv == "Other" & def_hg19vshg38$class == "Likely_benign", "Likely benign", as.character(def_hg19vshg38$cv))
def_hg19vshg38$cv <- ifelse(is.na(def_hg19vshg38$cv), "Other", as.character(def_hg19vshg38$cv))
def_hg19vshg38$cv <- factor(def_hg19vshg38$cv, levels =  c("Benign", "Likely benign", "Other", "Conflicting interpretations", "Uncertain significance", "Likely pathogenic", "Pathogenic"))

a <- (ggplot(def %>% filter(!is.na(nc)), 
                     aes( x = cv, fill = nc))
              + geom_bar()
              + labs(x = "ClinVar classification", y = "N ClinVar IDs", fill = "Sequence alignment")
              + theme_publication()
              + theme(panel.grid = element_blank(), legend.position = "bottom",
                         axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
                      plot.title = element_text(size = 15)
                      )
               + scale_fill_manual(values = palette_seq_mod)
)
b <- (ggplot(def_hg19vsT2T %>% filter(!is.na(nc)), 
                     aes( x = cv, fill = nc))
              + geom_bar()
              + labs(x = "ClinVar classification", y = "N ClinVar IDs", fill = "Sequence alignment")
              + theme_publication()
              + theme(panel.grid = element_blank(), legend.position = "bottom",
                          axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
                       plot.title = element_text(size = 15)
                      )
              + scale_fill_manual(values = palette_seq_mod)
)
c <- (ggplot(def_hg19vshg38, 
                     aes( x = cv, fill = nc))
              + geom_bar(show.legend = F)
              + labs(x = "ClinVar classification", y = "N ClinVar IDs", fill = "Sequence alignment")
              + theme_publication()
              + theme(panel.grid = element_blank(), legend.position = "bottom",
                        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
                       plot.title = element_text(size = 15)
                      )
              + scale_fill_manual(values = palette_seq_mod)
             )

patchwork::wrap_plots(a + labs(title = "hg38vsT2T-CHM13", x = "") , b + labs(title = "hg19vsT2T-CHM13", y = ""), 
                      c + labs(title = "hg19vshg38", y = "", x = ""),
                      nrow= 1, widths = c(4,4,4), guides = "collect", axes = "collect") & plot_annotation(theme = theme(legend.position = 'bottom'))
```

![](Figure_2_S5/figureS5B.png)<!-- -->

## Supplementary 5C


``` r
selected <- def %>% filter(gene_info %in% df_barplot_patho$gene_info)
```


``` r
(ggplot(selected, 
        aes( x = reorder(gene_info, pos.hg38), fill = nc))
              + geom_bar(position = "fill")
              + labs(x = "Gene", y = "Proportion of ClinVar IDs", fill = "Sequence alignment hg38vsT2T-CHM13")
              + theme_publication()
              + theme(panel.grid = element_blank(), legend.position = "bottom",
                       axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, face = "italic")
             )
              + facet_grid(.~ chr, scales = "free", space = "free")
              + scale_fill_manual(values = palette_seq_mod)
)
```

![](Figure_2_S5/figureS5C.png)<!-- -->

## Supplementary 5D


``` r
id <-  405391
pos_t2t <- def$pos.T2T[def$ID == id]
pos_hg38 <- def$pos.hg38[def$ID == id]
pos_hg19 <- def$pos.hg19[def$ID == id]
chrom <- def$chr[def$ID == id]
gene <- def$gene_info[def$ID == id]

width <- 150

def[def$ID == id]
```

```
##        ID    chr  pos.T2T ref.T2T alt.T2T pos.hg19 ref.hg19 alt.hg19 pos.hg38
##    <char> <fctr>    <int>  <char>  <char>    <int>   <char>   <char>    <int>
## 1: 405391   chr3 37026979       A       C 37067126        A        C 37025635
##    ref.hg38 alt.hg38      class                    type aln_hg19vshg38
##      <char>   <char>     <char>                  <char>         <char>
## 1:        A        C Pathogenic splice_acceptor_variant              .
##    match_hg19vshg38 mismatch_hg19vshg38 ins_hg19vshg38 del_hg19vshg38
##               <int>               <int>          <int>          <num>
## 1:              301                   0              0              0
##    score_hg19vshg38 idhg19                                       aln_hg19vsT2T
##               <int> <char>                                              <char>
## 1:              301 405391 ?37067094(1)|?37067096(1)|?37067098(1)|?37067119(1)
##    match_hg19vsT2T mismatch_hg19vsT2T ins_hg19vsT2T del_hg19vsT2T
##              <int>              <int>         <num>         <int>
## 1:             297                  4             0             0
##    score_hg19vsT2T idhg38                                       aln_hg38vsT2T
##              <int> <char>                                              <char>
## 1:             289 405391 ?37025603(1)|?37025605(1)|?37025607(1)|?37025628(1)
##    match_hg38vsT2T mismatch_hg38vsT2T ins_hg38vsT2T del_hg38vsT2T
##              <int>              <int>         <int>         <int>
## 1:             297                  4             0             0
##    score_hg38vsT2T              cl         cv                            tp
##              <int>          <fctr>     <char>                        <char>
## 1:             289 Equal reference Pathogenic Splice donor/acceptor variant
##         rsid             st_dbsnp      CUI
##       <char>               <char>   <char>
## 1: 267607815 Shared with dbSNP155 C0009405
##                                         disease gene_info
##                                          <char>    <char>
## 1: Hereditary nonpolyposis colorectal neoplasms      MLH1
##                         diff                nc var_type
##                       <char>            <fctr>   <char>
## 1: At least one modification [1, 5] mismatches      SNV
##                                                                                                                                                                                                                                                                                                                           reps_t2t
##                                                                                                                                                                                                                                                                                                                             <char>
## 1: T|T|T|T|T|T|T|T|T|T|T|T|T|T|T|T:37026947|37026949|37026950|37026951|37026953|37026955|37026956|37026957|37026958|37026959|37026960|37026961|37026962|37026963|37026964|37026967-37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972|37026972
##                                                                                                                                                                                                                                                              reps_hg38
##                                                                                                                                                                                                                                                                 <char>
## 1: T:37025609-37025629|T:37025611-37025629|T:37025612-37025629|T:37025613-37025629|T:37025614-37025629|T:37025615-37025629|T:37025616-37025629|T:37025617-37025629|T:37025618-37025629|T:37025619-37025629|T:37025620-37025629|T:37025623-37025629|T:37025625-37025629
##                stretch hm_match CLNREVSTAT                 hp_status
##                 <char>   <char>     <char>                    <fctr>
## 1: Homopolymer stretch        T       <NA> Difference in homopolymer
```


``` r
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1, heights=c(0.2, 0.5,0.5))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text(
  bquote("ClinVar ID: " ~ .(id) ~ " - Gene: " ~ italic(.(gene))),
  gp = gpar(fontsize = 16, fontface = "bold")
)
grid.legend(labels = c("A","C","T","G"), nrow = 1,
           pch = 15,  # square symbols
           gp = gpar(col = c("#9cd0e3", "#ea9249", "#2466a8", "#cd0016"), cex = 0.8),
           hgap = unit(0.5, "lines"),
           vp = viewport(x = 0.85, y = 0.5))
popViewport(1)  
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
plot_sequence_top(BSgenome.Hsapiens.UCSC.hg38, chrom, pos_hg38, width)
grid.lines(x = unit(c(0.014, 0.014), "npc"),
           y = unit(c(0.03, 0.97), "npc"),
           gp = gpar(col = "black", lwd = 2))
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
plot_sequence_bottom(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0, chrom, pos_t2t, width)
grid.lines(x = unit(c(0.014, 0.014), "npc"),
           y = unit(c(0.03, 0.97), "npc"),
           gp = gpar(col = "black", lwd = 2))

popViewport(1)
```

![](Figure_2_S5/figureS5D.png)<!-- -->

