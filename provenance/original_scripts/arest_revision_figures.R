#!/usr/bin/env Rscript
# AREST revised manuscript figures, 2026-10-31.
# Standalone companion to arestAnalysis_revised.R; reads its exported tables/RDS.
# Run: Rscript arest_revision_figures.R [revision_directory] [output_directory]
# No raw-read filtering or assembly is performed. Keep technical filters pooled.
# Dependencies (install once in R):
# install.packages(c('dplyr','tidyr','readr','ggplot2','patchwork','vegan','ggh4x','circlize','BiocManager'))
# BiocManager::install('ComplexHeatmap', ask=FALSE, update=FALSE)
# Default source is pinned to the validated run, NOT automatically the newest run.
args <- commandArgs(trailingOnly=TRUE)
RUN <- if(length(args)>=1) args[1] else '/home/sersancar/Documents/AREST/results/revision_20261030_214057'
OUT <- if(length(args)>=2) args[2] else file.path(RUN,paste0('manuscript_figures_',format(Sys.time(),'%Y%m%d_%H%M%S')))
MAKE_SITE_ONLY_DESCRIPTIVE <- FALSE # optional extra S8; not a replacement for site+date tests
TOP_TEMPORAL_GENES <- 20L
BUBBLE_SELECTION <- 'original_list' # 'all' recommended; 'original_list' reproduces the declared subset
ELLIPSES <- TRUE # descriptive 95% normal data ellipses, not confidence regions of means
GENE_LABEL_PAGES <- TRUE # S2 full gene labels in a multipage PDF
pkgs <- c('dplyr','tidyr','readr','ggplot2','patchwork','vegan','ggh4x','ComplexHeatmap','circlize')
missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install missing packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2);library(patchwork);library(vegan);library(ggh4x);library(ComplexHeatmap);library(circlize)})
if(!dir.exists(RUN)) stop('Revision folder not found: ',RUN)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUT,'plot_data'),showWarnings=FALSE)
# NPG palette used in the original script: A=3, B=1, C=4, freshwater=2.
COL <- c(A='#00A087',B='#E64B35',C='#3C5488',D='#4DBBD5')
THEME <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face='bold'),strip.background=element_rect(fill='grey95'),
  legend.position='bottom',plot.caption=element_text(hjust=0,size=9))
theme_set(THEME)
read_table <- function(name,needed) {
 f <- file.path(RUN,'tables',paste0(name,'.tsv'))
 if(!file.exists(f)) stop('Missing input: ',f)
 z <- read_tsv(f,show_col_types=FALSE,progress=FALSE)
 if(!all(needed %in% names(z))) stop(name,': missing columns ',paste(setdiff(needed,names(z)),collapse=', '))
 if('Date' %in% names(z)) z$Date <- as.Date(z$Date)
 if('Site' %in% names(z)) z$Site <- as.character(z$Site)
 z
}
export <- function(x,name) write_tsv(as.data.frame(x),file.path(OUT,'plot_data',paste0(name,'.tsv')))
save_plot <- function(p,name,w=10,h=7) {
 ggsave(file.path(OUT,paste0(name,'.pdf')),p,width=w,height=h,device='pdf',useDingbats=FALSE,limitsize=FALSE)
 ggsave(file.path(OUT,paste0(name,'.png')),p,width=w,height=h,dpi=600,bg='white',limitsize=FALSE)
 message('Written: ',name)
}
md <- read_table('S1_biological_environment',c('BioSample','Site','Date','TotalReads','nFilters','Salinity','pH','WaterTemperature')) %>% arrange(Site,Date)
if(anyDuplicated(md$BioSample)||nrow(md)!=20||any(md$nFilters!=3)) stop('Expected 20 unique biological samples with three filters each.')
if(!setequal(unique(md$Site),names(COL))) stop('Expected site labels A/B/C/D.')
if(any(!is.finite(md$TotalReads)|md$TotalReads<=0)||anyNA(md$Date)) stop('Invalid biological sample metadata.')
DATES <- sort(unique(md$Date)); DATE_LABEL <- format(DATES,'%d %b %Y')
SHAPES <- setNames(c(16,17,15,3,7,8)[seq_along(DATES)],DATE_LABEL)
md$DateLabel <- factor(format(md$Date,'%d %b %Y'),levels=DATE_LABEL)
md$Site <- factor(md$Site,levels=names(COL))
ref <- read_table('S5_MEGARes_reference_biological',c('BioSample','gene_accession','GeneId','Type','Class','Counts','RPM','length','RPKM'))
ref$Type <- ifelse(ref$Type=='Multi-compount','Multi-compound',ref$Type)
if(anyDuplicated(ref[c('BioSample','gene_accession')])) stop('Duplicated reference/sample rows.')
if(!setequal(ref$BioSample,md$BioSample)) stop('Reference/metadata sample mismatch.')
if(any(!is.finite(ref$RPKM)|ref$RPKM<0)||any(!is.finite(ref$RPM)|ref$RPM<0)) stop('Invalid abundance values.')
if(any(!is.finite(ref$length)|ref$length<=0)) stop('Reference lengths are required for RPKM sensitivity.')
cb <- read_table('S3_CARPDM_biological_thresholds',c('BioSample','gene','GeneId','Class','Screen','nPassFilters','RPM','RPKM')) %>% filter(Screen=='strict')
if(!setequal(cb$BioSample,md$BioSample)) stop('CARPDM/metadata sample mismatch.')
if(anyDuplicated(cb[c('BioSample','gene')])||any(!cb$nPassFilters %in% 0:3)) stop('Invalid CARPDM support rows.')
# Require explicit zero rows from the completed revised run; never invent missing samples.
if(nrow(cb)!=nrow(md)*n_distinct(cb$gene)) stop('CARPDM table must be complete across samples and references.')
make_matrix <- function(metric) {
 z <- ref %>% group_by(BioSample,GeneId) %>% summarise(Value=sum(.data[[metric]]),.groups='drop') %>%
  pivot_wider(names_from=GeneId,values_from=Value,values_fill=0)
 x <- as.matrix(z[,-1]); rownames(x) <- z$BioSample
 x <- x[md$BioSample,,drop=FALSE]; x <- x[,colSums(x)>0,drop=FALSE]
 if(any(rowSums(x)<=0)) stop('Zero-sum sample in ',metric)
 x
}
MATS <- lapply(c('RPM','RPKM'),make_matrix); names(MATS)<-c('RPM','RPKM')
# One sample order used by the class/gene heatmaps, clustered WITHIN each site.
X <- MATS$RPM
order_ids <- unlist(lapply(names(COL),function(s) {
 ids<-md$BioSample[md$Site==s]
 if(length(ids)>1) ids[hclust(vegdist(decostand(X[ids,,drop=FALSE],'hellinger'),'bray'),method='complete')$order] else ids
}),use.names=FALSE)
short_ids <- setNames(paste(as.character(md$Site),format(md$Date,'%d-%b')),md$BioSample)
# Colour values remain log10(abundance+1), as in the original script.
tile_plot <- function(z,metric,title,row_levels=NULL,label_rows=TRUE) {
 if(is.null(row_levels)) row_levels <- sort(unique(z$Feature))
 z <- z %>% mutate(BioSample=factor(BioSample,levels=order_ids),Feature=factor(Feature,levels=rev(row_levels)),
   Site=factor(Site,levels=names(COL)))
 ggplot(z,aes(BioSample,Feature,fill=log10(Value+1)))+geom_tile()+
  facet_grid(.~Site,scales='free_x',space='free_x')+
  scale_fill_gradientn(colours=c('#FFFFFF','#FFF7BC','#FEC44F','#D95F0E','#7F0000'),name=paste0('log10(',metric,' + 1)'))+
  scale_x_discrete(labels=short_ids)+labs(x=NULL,y=NULL,title=title)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8),panel.grid=element_blank(),
    axis.text.y=if(label_rows) element_text(size=7) else element_blank(),axis.ticks.y=element_blank())
}
# Figure 1A: previous theme and typography with revised biological-sample means.
# Each site mean gives equal weight to dates after pooling the three filters.

cg <- cb %>%
  group_by(BioSample, GeneId, Class) %>%
  summarise(RPM = sum(RPM), RPKM = sum(RPKM), .groups = 'drop') %>%
  left_join(md %>% select(BioSample, Site, Date), by = 'BioSample')

cmean <- cg %>%
  group_by(Site, GeneId, Class) %>%
  summarise(
    MeanRPM = mean(RPM),
    n_dates = n(),
    SD = sd(RPM),
    .groups = 'drop'
  ) %>%
  mutate(
    SE = SD / sqrt(n_dates),
    CI95 = qt(.975, n_dates - 1) * SE
  )

export(cmean, 'Figure1A_CARPDM_means')

card <- cmean %>% filter(Site %in% c('B', 'D'))
lev <- sort(unique(card$Class))
ng <- card %>% distinct(Class, GeneId) %>% count(Class)
angles <- ifelse(ng$n[match(lev, ng$Class)] <= 3, 90, 0)
card$Class <- factor(card$Class, levels = lev)

dodge <- position_dodge(width = .75)

class_strip_text <- ggh4x::elem_list_text(
  angle = angles,
  face = rep("bold", length(lev)),
  size = rep(14, length(lev))
)

p1a <- ggplot(card, aes(GeneId, MeanRPM, fill = Site)) +
  geom_col(position = dodge, width = .7) +
  geom_errorbar(
    aes(
      ymin = pmax(0, MeanRPM - CI95),
      ymax = MeanRPM + CI95
    ),
    position = dodge,
    width = .25
  ) +
  ggh4x::facet_grid2(
    . ~ Class,
    scales = 'free_x',
    space = 'free_x',
    labeller = labeller(
      Class = function(x) {
        vapply(x, function(z) paste(strwrap(z, width = 18), collapse = '\n'),
               character(1))
      }
    ),
    strip = ggh4x::strip_themed(
      text_x = class_strip_text
    )
  ) +
  scale_fill_manual(values = COL, breaks = c('B', 'D')) +
  labs(
    x = 'ARG',
    y = 'RPM',
    fill = 'Site'
    ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = 'outside',
    strip.clip = 'off',
    axis.text.x = element_text(
      size = 13, face = 'bold', angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(size = 20, face = 'bold'),
    axis.title = element_text(size = 25, face = 'bold'),
    legend.key.size = grid::unit(1.2, 'cm'),
    legend.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 15),
    panel.spacing.x = grid::unit(1, 'cm'),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

save_plot(p1a, 'Figure1A_CARPDM', 23, 10)

p1a_sqrt <- p1a +
  scale_y_sqrt() +
  labs(
    y = "RPM (square-root scale)"
)

save_plot(p1a_sqrt, "Figure1A_CARPDM_sqrt", 23, 10)

# Figure 1B: descriptive read-level taxonomic classifications of CARPDM-matching reads.
# Keep all observed genera by default; highlight a prespecified clinical/environmental set.
# Red labels indicate genera containing recognized priority or aquatic/enteric pathogens;
# they do not indicate species identification or confirmed ARG hosts.
original_list <- c('Escherichia','Salmonella','Vibrio','Clostridium','Bacteroides','Shigella',
                   'Campylobacter','Yersinia','Legionella','Pseudomonas','Mycobacterium','Microcystis',
                   'Dolichospermum','Nodularia','Klebsiella','Enterobacter','Citrobacter','Planktothrix',
                   'Cylindrospermopsis','Aphanizomenon','Acinetobacter','Burkholderia','Stenotrophomonas',
                   'Sphingomonas','Ralstonia','Aeromonas','Arcobacter','Brucella','Coxiella','Edwardsiella',
                   'Enterococcus','Erysipelothrix','Francisella','Helicobacter','Lawsonia','Leptospira',
                   'Listeria','Pasteurella','Proteus','Providencia','Streptococcus','Bordetella','Haemophilus',
                   'Moraxella','Nocardia','Staphylococcus','Bacillus','Corynebacterium','Trueperella',
                   'Anaplasma','Bartonella','Chlamydia','Chlamydophila','Ehrlichia','Mycoplasma','Rickettsia','Tropheryma')
highlight_genera <- c('Acinetobacter','Aeromonas','Campylobacter','Citrobacter',
                      'Enterococcus','Escherichia','Klebsiella','Proteus','Pseudomonas','Salmonella',
                      'Staphylococcus','Streptococcus')
tx <- read_table('S4_CARPDM_taxonomy_all_rows',c('BioSample','ARG_ID','genus','read_count'))
if(any(!is.finite(tx$read_count)|tx$read_count<0)) stop('Invalid taxonomy counts.')
lookup <- cb %>% distinct(gene,GeneId,Class)
if(anyDuplicated(lookup$gene)) stop('Ambiguous CARPDM reference annotation.')
tax <- tx %>% select(BioSample,ARG_ID,genus,read_count) %>%
  left_join(lookup,by=c('ARG_ID'='gene')) %>%
  filter(!is.na(genus),!tolower(genus) %in% c('unclassified','unassigned',''),!is.na(GeneId)) %>%
  group_by(BioSample,GeneId,Class,genus) %>% summarise(Reads=sum(read_count),.groups='drop') %>%
  left_join(md %>% select(BioSample,Site,TotalReads),by='BioSample') %>% mutate(RPM=1e6*Reads/TotalReads)
# Sum observed biological RPM then divide by ALL sampled dates, including zeros.
# Unlike the old code, multiple species within a genus are summed before averaging.
den <- md %>% count(Site,name='n_dates')
bubble_all <- tax %>% group_by(Site,GeneId,Class,genus) %>% summarise(SumRPM=sum(RPM),.groups='drop') %>%
  left_join(den,by='Site') %>% mutate(RPM=SumRPM/n_dates) %>% filter(Site %in% c('B','D'),RPM>0)
export(bubble_all,'Figure1B_all_genera_means')
bubble <- bubble_all
if(BUBBLE_SELECTION=='original_list') {
  bubble <- filter(bubble,genus %in% original_list)
  genera <- bubble %>% group_by(genus) %>% summarise(x=sum(RPM),.groups='drop') %>% arrange(desc(x),genus) %>% pull(genus)
  bubble <- filter(bubble,genus %in% genera)
  genes <- bubble %>% group_by(GeneId) %>% summarise(x=sum(RPM),.groups='drop') %>% arrange(desc(x),GeneId) %>% pull(GeneId)
  bubble <- filter(bubble,GeneId %in% genes)
} else if(BUBBLE_SELECTION!='all') stop('BUBBLE_SELECTION must be all or original_list.')
bl_order <- bubble %>% group_by(Class,GeneId) %>% summarise(x=sum(RPM),.groups='drop') %>% arrange(Class,desc(x),GeneId) %>% pull(GeneId) %>% unique()
gen_order <- bubble %>% group_by(genus) %>% summarise(x=sum(RPM),.groups='drop') %>% arrange(desc(x),genus) %>% pull(genus)
bubble <- bubble %>% mutate(GeneId=factor(GeneId,levels=bl_order),genus=factor(genus,levels=gen_order))
export(bubble,'Figure1B_plotted_genera')
selection_text <- if(BUBBLE_SELECTION=='all') 'All assigned genera shown; no top-N filtering.' else 'Selected genus list, with top taxa/genes; full annotations exported separately.'
genus_labeller <- function(x) {
  ifelse(x %in% highlight_genera,
         paste0("<span style='color:#B2182B'><b><i>",x,"</i></b></span>"),
         paste0('<i>',x,'</i>'))
}

CARPDM_CLASS_COLS <- c(
  'Aminoglycosides'   = '#0072B2',
  'Betalactams'       = '#D55E00',
  'Cat. AM peptides'  = '#009E73',
  'D & B res.'        = '#CC79A7',
  'Lipopeptides'      = '#B28E00',
  'MLS'               = '#56B4E9',
  'Multi-D res.'      = '#882255',
  'Nucleosides'       = '#4D4D4D',
  'Phenicol'          = '#332288',
  'Sulfonamides'      = '#8C564B',
  'Tetracyclines'     = '#117733',
  'Trimethoprim'      = '#AA4499'
)


p1b <- ggplot(bubble,aes(GeneId,genus))+
  geom_point(aes(size=RPM,colour=Class,shape=Site),alpha=.85,
             position=position_dodge2(width=.6,preserve='single'))+
  scale_colour_manual(values=CARPDM_CLASS_COLS,breaks=names(CARPDM_CLASS_COLS),drop=FALSE)+ 
  scale_shape_manual(values=c(B=16,D=17))+
  scale_size_continuous(range=c(2.5,15),breaks=scales::pretty_breaks(n=4))+
  scale_y_discrete(labels=genus_labeller)+
  labs(x='ARG',y='Kraken2 genus assignment',colour='AMR class',
       size='Mean RPM',shape='Site')+
  theme_bw(base_size=12)+
  theme(panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,size=17,face='bold'),
        axis.title.x=element_text(size=25,face='bold'),
        axis.text.y=ggtext::element_markdown(size=18),
        axis.title.y=element_text(size=25,face='bold'),
        legend.key.size=grid::unit(1.2,'cm'),
        legend.text=element_text(size=15),
        legend.title=element_text(size=16,face='bold'),
        plot.margin=margin(t=15,r=20,b=20,l=15))+
  guides(colour=guide_legend(override.aes=list(size=5)),
         shape=guide_legend(override.aes=list(size=5)))
save_plot(p1b,'Figure1B_taxonomy_bubbles',23,max(10,.34*length(gen_order)))

save_plot(p1a/p1b+plot_layout(heights=c(1,1.3)),'Figure1_combined',23,max(20,10+.28*length(gen_order)))

# Figure 2A: class heatmap with Type stripe and site annotation.
cl <- ref %>% group_by(BioSample,Type,Class) %>% summarise(Value=sum(RPM),.groups='drop') %>%
  mutate(Feature=paste(Type,Class,sep=' | '))
classanno <- cl %>% distinct(Feature,Type,Class) %>% arrange(Type,Class)
w <- cl %>% select(Feature,BioSample,Value) %>% pivot_wider(names_from=BioSample,values_from=Value,values_fill=0)
hm <- as.matrix(w[,-1]); rownames(hm)<-w$Feature
hm <- log10(hm[classanno$Feature,order_ids,drop=FALSE]+1)
class_short_map <- c(
  'Acetate_resistance'='Acetate res.',
  'Acid_resistance'='Acid res.',
  'Aldehyde_resistance'='Aldehyde res.',
  'Biguanide_resistance'='Biguanide res.',
  'Multi-biocide_resistance'='Multi-B res.',
  'Peroxide_resistance'='Peroxide res.',
  'Phenolic_compound_resistance'='Phenolic res.',
  'Quaternary_Ammonium_Compounds_(QACs)_resistance'='QACs res.',
  'Cationic_antimicrobial_peptides'='Cat. AM peptides',
  'Aminocoumarins'='Aminocoumarins',
  'Aminoglycosides'='Aminoglycosides',
  'Elfamycins'='Elfamycins',
  'Fluoroquinolones'='Fluoroquinolones',
  'Fosfomycin'='Fosfomycin',
  'Fusidic_acid'='Fusidic acid',
  'Glycopeptides'='Glycopeptides',
  'Lipopeptides'='Lipopeptides',
  'MLS'='MLS',
  'Metronidazole'='Metronidazole',
  'Multi-drug_resistance'='Multi-D res.',
  'Mupirocin'='Mupirocin',
  'Mycobacterium_tuberculosis-specific_Drug'='M. tuberculosis spec. D',
  'Nucleosides'='Nucleosides',
  'Oxazolidinone'='Oxazolidinone',
  'Phenicol'='Phenicol',
  'Pleuromutilin'='Pleuromutilin',
  'Rifampin'='Rifampin',
  'Spiropyrimidinetriones'='Spiropyrimidinetriones',
  'Sulfonamides'='Sulfonamides',
  'Tetracenomycin'='Tetracenomycin',
  'Tetracyclines'='Tetracyclines',
  'Thiopeptides'='Thiopeptides',
  'Trimethoprim'='Trimethoprim',
  'betalactams'='Betalactams',
  'Drug_and_biocide_resistance'='D and B res.',
  'Biocide_and_metal_resistance'='B & M res.',
  'Drug_and_biocide_and_metal_resistance'='D, B & M res.',
  'Multi-metal_resistance'='Multi-M res.',
  'Arsenic_resistance'='As res.',
  'Chromium_resistance'='Cr res.',
  'Copper_resistance'='Cu res.',
  'Iron_resistance'='Fe res.',
  'Lead_resistance'='Pb res.',
  'Mercury_resistance'='Hg res.',
  'Nickel_resistance'='Ni res.',
  'Sodium_resistance'='Na res.',
  'Tellurium_resistance'='Te res.',
  'Zinc_resistance'='Zn res.', 
  "Bacitracin" = "Bacitracin",
  "Cadmium_resistance" = "Cd res."
)
heatmap_class_labels <- unname(class_short_map[classanno$Class])
missing_labels <- is.na(heatmap_class_labels)
heatmap_class_labels[missing_labels] <- gsub('_',' ',classanno$Class[missing_labels],fixed=TRUE)
names(heatmap_class_labels) <- classanno$Feature
short_ids_fig2 <- setNames(
  paste(as.character(md$Site),format(md$Date,'%d-%b-%Y')),
  md$BioSample)
anno_legend_gp <- list(
  title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 14),
  grid_width = grid::unit(0.8, "cm"),
  grid_height = grid::unit(0.8, "cm")
)
hm_legend <- list(
  title = expression(log[10](RPM + 1)),
  title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 14),
  legend_height = grid::unit(6, "cm")
)
site_annotation_fig2 <- function(ids) {
  ComplexHeatmap::HeatmapAnnotation(
    Site = factor(
      md$Site[match(ids, md$BioSample)],
      levels = names(COL)
    ),
    col = list(Site = COL),
    annotation_name_gp = grid::gpar(fontsize = 16, fontface = "bold"),
    annotation_legend_param = list(Site = anno_legend_gp)
  )
}
p2a <- ComplexHeatmap::rowAnnotation(Type=classanno$Type,col=list(Type=TYPECOL),show_annotation_name=FALSE,
                                     annotation_legend_param=list(Type=anno_legend_gp))+
  ComplexHeatmap::Heatmap(hm,name='log10(RPM+1)',col=heat_col(hm),cluster_rows=FALSE,
                          cluster_columns=FALSE,column_split=factor(md$Site[match(order_ids,md$BioSample)],levels=names(COL)),
                          top_annotation=site_annotation_fig2(order_ids),row_split=classanno$Type,row_title=NULL,
                          row_labels=heatmap_class_labels[classanno$Feature],column_labels=short_ids_fig2[order_ids],column_names_rot=90,
                          row_names_gp=grid::gpar(fontsize=12),column_names_gp=grid::gpar(fontsize=11),
                          row_title_gp=grid::gpar(fontsize=13,fontface='bold'),column_title_gp=grid::gpar(fontsize=16,fontface='bold'),
                          heatmap_legend_param=hm_legend)
export(cl,'Figure2A_class_abundance'); save_heat(p2a,'Figure2A_class_heatmap',15,11)

# Figure 2B and S7: classical PCoA, Cailliez additive correction as in revised script.
# Axis percentages use the positive eigenvalues of the corrected distance matrix.
pcoa <- function(metric,ids,title,key) {
  H <- decostand(MATS[[metric]][ids,,drop=FALSE],'hellinger')
  d <- vegdist(H,'bray'); ord <- cmdscale(d,k=2,eig=TRUE,add=TRUE)
  vals <- pmax(ord$eig,0); pct <- 100*vals[1:2]/sum(vals)
  pts <- tibble(BioSample=rownames(ord$points),Axis1=ord$points[,1],Axis2=ord$points[,2]) %>% left_join(md,by='BioSample')
  export(pts,paste0(key,'_coordinates'))
  write_tsv(tibble(Axis=seq_along(ord$eig),Eigenvalue=ord$eig),file.path(OUT,'plot_data',paste0(key,'_eigenvalues.tsv')))
  saveRDS(list(ordination=ord,distance=d,ids=ids),file.path(OUT,'plot_data',paste0(key,'.rds')))
  ellipse_data <- pts %>% group_by(Site) %>% filter(n()>=4) %>% ungroup()
  # Skip degenerate covariance groups as well as sites with too few dates.
  if(nrow(ellipse_data)) ellipse_data <- ellipse_data %>% group_by(Site) %>% filter(det(cov(cbind(Axis1,Axis2)))>1e-12) %>% ungroup()
  ggplot(pts,aes(Axis1,Axis2,colour=Site,shape=DateLabel))+
    {if(ELLIPSES && nrow(ellipse_data)) stat_ellipse(data=ellipse_data,aes(group=Site),type='norm',level=.95,linewidth=.7)}+
    geom_point(size=4,alpha=.9)+
    scale_colour_manual(values=COL,drop=FALSE)+scale_shape_manual(values=SHAPES,drop=FALSE,name='Date')+
    coord_equal()+
    labs(x=sprintf('PCoA1 (%.1f%%)',pct[1]),y=sprintf('PCoA2 (%.1f%%)',pct[2])) +
    theme(
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(size = 14), 
      legend.title = element_text(face = "bold", size = 15), 
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.2, "cm"),
      legend.position = "right"
    )
}

p2b <- pcoa('RPM',md$BioSample,'B  Resistome composition: all 20 water samples','Figure2B')
save_plot(p2b,'Figure2B_PCoA',10,7)

p2a_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(p2a,heatmap_legend_side='right',annotation_legend_side='right',merge_legends=TRUE),width=15,height=11)
save_plot(patchwork::wrap_elements(full=p2a_grob)+p2b+plot_layout(widths=c(1.5,1)),'Figure2_combined',24,11)

# Figure S2: broad-resistome gene heatmap and site-overlap Venn diagram.
# Figure S2: original layout with the revised 20 biological samples.
# Requires your complete class_short_map to be defined.

local({
  
  RUN_S2 <- "/home/sersancar/Documents/AREST/results/revision_20261030_214057"
  OUT_S2 <- file.path(RUN_S2, "manuscript_figures_20261101_183051")
  dir.create(OUT_S2, recursive = TRUE, showWarnings = FALSE)
  
  ComplexHeatmap::ht_opt(RESET = TRUE)
  
  # ------------------------------------------------------------
  # 1) Read revised biological-sample tables
  # ------------------------------------------------------------
  dat <- readr::read_tsv(
    file.path(RUN_S2, "tables", "S5_MEGARes_reference_biological.tsv"),
    show_col_types = FALSE
  )
  
  meta <- readr::read_tsv(
    file.path(RUN_S2, "tables", "S1_biological_samples.tsv"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      BioSample = as.character(BioSample),
      Site = as.character(Site),
      Date = as.Date(Date)
    ) %>%
    dplyr::arrange(Site, Date)
  
  stopifnot(
    nrow(meta) == 20L,
    !anyDuplicated(meta$BioSample),
    !anyNA(meta$Site),
    !anyNA(meta$Date)
  )
  
  site_cols <- c(
    A = "#00A087",
    B = "#E64B35",
    C = "#3C5488",
    D = "#4DBBD5"
  )
  
  # ------------------------------------------------------------
  # 2) Ordinary numeric genes × samples RPM matrix
  # ------------------------------------------------------------
  abund <- dat %>%
    dplyr::transmute(
      GeneId = as.character(GeneId),
      BioSample = as.character(BioSample),
      RPM = as.numeric(RPM)
    ) %>%
    dplyr::group_by(GeneId, BioSample) %>%
    dplyr::summarise(RPM = sum(RPM), .groups = "drop")
  
  stopifnot(
    all(is.finite(abund$RPM)),
    all(abund$RPM >= 0),
    setequal(abund$BioSample, meta$BioSample)
  )
  
  rpm_table <- stats::xtabs(
    RPM ~ GeneId + BioSample,
    data = abund
  )
  
  # Remove the contingency-table class explicitly.
  rpm <- matrix(
    as.numeric(rpm_table),
    nrow = nrow(rpm_table),
    ncol = ncol(rpm_table),
    dimnames = dimnames(rpm_table)
  )
  
  rpm <- rpm[
    rowSums(rpm) > 0,
    meta$BioSample,
    drop = FALSE
  ]
  
  stopifnot(
    is.matrix(rpm),
    is.numeric(rpm),
    !inherits(rpm, "table"),
    ncol(rpm) == 20L,
    identical(colnames(rpm), meta$BioSample),
    all(is.finite(rpm))
  )
  
  # ------------------------------------------------------------
  # 3) Original Type/Class/GeneId row ordering
  # ------------------------------------------------------------
  annotation <- dat %>%
    dplyr::transmute(
      GeneId = as.character(GeneId),
      Type = gsub(
        "Multi-compount", "Multi-compound",
        as.character(Type), fixed = TRUE
      ),
      Class = as.character(Class)
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Type, Class, GeneId) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE) %>%
    dplyr::filter(GeneId %in% rownames(rpm)) %>%
    dplyr::arrange(Type, Class, GeneId)
  
  rpm <- rpm[annotation$GeneId, , drop = FALSE]
  amr_data <- log10(rpm + 1)
  
  stopifnot(
    identical(rownames(amr_data), annotation$GeneId),
    identical(colnames(amr_data), meta$BioSample),
    ncol(amr_data) == 20L,
    all(is.finite(amr_data))
  )
  
  type_levels <- unique(annotation$Type)
  type_split <- factor(annotation$Type, levels = type_levels)
  
  type_palette <- c(
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
  )
  
  type_cols <- setNames(
    type_palette[seq_along(type_levels)],
    type_levels
  )
  
  # ------------------------------------------------------------
  # 4) Original black/white class-membership panel
  # ------------------------------------------------------------
  if (!exists("class_short_map")) {
    stop("Define your complete class_short_map before running this block.")
  }
  
  class_map <- class_short_map
  class_map["Quaternary_Ammonium_Compounds_(QACs)_resistance"] <- "QACs res."
  class_map["Bacitracin"] <- "Bacitracin"
  class_map["Cadmium_resistance"] <- "Cd res."
  
  sub_levels <- unique(annotation$Class)
  sub_labels <- unname(class_map[sub_levels])
  
  missing_labels <- is.na(sub_labels)
  sub_labels[missing_labels] <- gsub(
    "_", " ", sub_levels[missing_labels], fixed = TRUE
  )
  
  sub_mat <- matrix(
    0L,
    nrow = nrow(amr_data),
    ncol = length(sub_levels),
    dimnames = list(rownames(amr_data), sub_levels)
  )
  
  sub_mat[cbind(
    seq_len(nrow(amr_data)),
    match(annotation$Class, sub_levels)
  )] <- 1L
  
  # ------------------------------------------------------------
  # 5) Bray–Curtis / complete-linkage clustering
  # ------------------------------------------------------------
  # For custom column clustering, ComplexHeatmap supplies:
  # rows = samples; columns = genes.
  column_cluster <- function(m) {
    
    ids <- rownames(m)
    
    if (is.null(ids) || anyNA(match(ids, colnames(rpm)))) {
      stop(
        "Sample row names supplied to column_cluster do not match ",
        "the biological-sample RPM matrix."
      )
    }
    
    profiles <- t(rpm[, ids, drop = FALSE])
    
    if (any(rowSums(profiles) <= 0)) {
      stop(
        "Zero-total RPM profile for: ",
        paste(
          rownames(profiles)[rowSums(profiles) <= 0],
          collapse = ", "
        )
      )
    }
    
    distances <- vegan::vegdist(
      profiles,
      method = "bray"
    )
    
    if (any(!is.finite(distances))) {
      stop("Undefined Bray–Curtis distances.")
    }
    
    stats::hclust(
      distances,
      method = "complete"
    )
  }
  
  # ------------------------------------------------------------
  # 6) Legends and annotations
  # ------------------------------------------------------------
  legend_style <- list(
    title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 17),#
    grid_width = grid::unit(1, "cm"),
    grid_height = grid::unit(1, "cm")
  )
  
  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    Site = factor(meta$Site, levels = names(site_cols)),
    col = list(Site = site_cols),
    annotation_name_gp = grid::gpar(
      fontsize = 25,#
      fontface = "bold"
    ),
    annotation_legend_param = list(Site = legend_style)
  )
  
  left_ha <- ComplexHeatmap::rowAnnotation(
    Type = type_split,
    col = list(Type = type_cols),
    show_annotation_name = FALSE,
    annotation_width = grid::unit(1.5, "mm"),
    annotation_legend_param = list(Type = legend_style)
  )
  
  # ------------------------------------------------------------
  # 7) Original abundance palette
  # ------------------------------------------------------------
  lims <- as.numeric(stats::quantile(
    amr_data,
    probs = c(0.02, 0.98),
    names = FALSE
  ))
  
  lims <- c(floor(lims[1]), ceiling(lims[2]))
  
  if (lims[2] <= lims[1]) {
    lims[2] <- lims[1] + 1
  }
  
  col_fun <- circlize::colorRamp2(
    c(lims[1], mean(lims), lims[2]),
    c("#2c7bb6", "#ffffbf", "#d7191c")
  )
  
  # ------------------------------------------------------------
  # 8) Class-membership heatmap
  # ------------------------------------------------------------
  ht_sub <- ComplexHeatmap::Heatmap(
    sub_mat,
    name = "Class membership",
    col = circlize::colorRamp2(
      c(0, 1),
      c("white", "black")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_split = type_split,
    row_title = NULL,
    row_gap = grid::unit(1, "mm"),
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_labels = sub_labels,
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 18),#
    rect_gp = grid::gpar(col = NA),
    width = grid::unit(30, "cm"),
    show_heatmap_legend = FALSE
  )
  
  # ------------------------------------------------------------
  # 9) Main abundance heatmap: 20 biological samples
  # ------------------------------------------------------------
  ht_main <- ComplexHeatmap::Heatmap(
    amr_data,
    name = "log10(RPM+1)",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = column_cluster,
    cluster_column_slices = FALSE,
    column_dend_reorder = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_labels = paste(
      meta$Site,
      format(meta$Date, "%d-%b-%Y")
    ),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 18),#
    top_annotation = top_ha,
    row_split = type_split,
    row_title = NULL,
    row_gap = grid::unit(1, "mm"),
    # column_split = factor(
    #   meta$Site,
    #   levels = names(site_cols)
    # ),
    column_title = NULL,
    heatmap_legend_param = list(
      title = expression(bold(log[10](RPM + 1))),
      title_gp = grid::gpar(
        fontsize = 18,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(fontsize = 16),
      legend_height = grid::unit(6, "cm")
    ), 
    width = grid::unit(14, "cm")
  )
  
  ht <- ht_sub + left_ha + ht_main
  
  message(
    "Abundance matrix: ", nrow(amr_data),
    " genes × ", ncol(amr_data), " biological samples."
  )
  message(
    "Class-membership panel: ",
    ncol(sub_mat), " resistance classes."
  )
  
  # ------------------------------------------------------------
  # 10) Export plotted values
  # ------------------------------------------------------------
  plot_data <- data.frame(
    GeneId = rownames(amr_data),
    as.data.frame.matrix(amr_data),
    row.names = NULL,
    check.names = FALSE
  )
  
  readr::write_tsv(
    plot_data,
    file.path(OUT_S2, "FigureS2_plotted_matrix.tsv")
  )
  
  # ------------------------------------------------------------
  # 11) Save PDF and PNG
  # ------------------------------------------------------------
  save_s2 <- function(filename, format) {
    
    if (format == "pdf") {
      grDevices::pdf(
        filename,
        width = 25,
        height = 17,
        useDingbats = FALSE
      )
    } else {
      grDevices::png(
        filename,
        width = 25,
        height = 17,
        units = "in",
        res = 600
      )
    }
    
    on.exit(grDevices::dev.off(), add = TRUE)
    
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      main_heatmap = "log10(RPM+1)", 
      padding = grid::unit(c(25, 5, 5, 5), "mm"),
    )
  }
  
  save_s2(
    file.path(OUT_S2, "FigureS2_original_layout_20_samples.pdf"),
    "pdf"
  )
  
  save_s2(
    file.path(OUT_S2, "FigureS2_original_layout_20_samples.png"),
    "png"
  )
  
  message("Files saved in: ", OUT_S2)
})

# Separate Venn diagram: genes detected at each site on dates shared by all sites.
required_sites_s2 <- c("A", "B", "C", "D")
dates_by_site_s2 <- split(sample_meta_s2$Date, sample_meta_s2$Site)

if (!all(required_sites_s2 %in% names(dates_by_site_s2))) {
  stop("Metadata must include sites A, B, C and D for the Venn diagram.")
}

shared_dates_s2 <- Reduce(
  intersect,
  dates_by_site_s2[required_sites_s2]
)
if (!length(shared_dates_s2)) {
  stop("No sampling dates are shared across sites A, B, C and D.")
}

gene_rpm_s2 <- ref %>%
  dplyr::mutate(
    BioSample = as.character(BioSample),
    GeneId = as.character(GeneId)
  ) %>%
  dplyr::group_by(BioSample, GeneId) %>%
  dplyr::summarise(RPM = sum(RPM, na.rm = TRUE), .groups = "drop")

venn_data_s2 <- gene_rpm_s2 %>%
  dplyr::left_join(
    md_s2 %>% dplyr::select(BioSample, Site, Date),
    by = "BioSample"
  ) %>%
  dplyr::filter(
    Date %in% shared_dates_s2,
    Site %in% required_sites_s2,
    RPM > 0
  )

site_sets_s2 <- lapply(required_sites_s2, function(site_i) {
  venn_data_s2 %>%
    dplyr::filter(Site == site_i) %>%
    dplyr::pull(GeneId) %>%
    unique()
})
names(site_sets_s2) <- required_sites_s2

venn_grob_s2 <- VennDiagram::venn.diagram(
  x = site_sets_s2,
  filename = NULL,
  fill = unname(COL[required_sites_s2]),
  alpha = 0.45,
  lwd = 1.5,
  cex = 1.2,
  cat.cex = 1.3,
  fontface = "bold",
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  main = "Broad resistome genes\n(shared sampling dates)"
)

grDevices::pdf(
  file.path(OUT, "FigureS2_broad_resistome_venn.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

grDevices::png(
  file.path(OUT, "FigureS2_broad_resistome_venn.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 600
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()



# Separate Venn diagram using dates sampled at all four sites.
required_sites_s2 <- c("A", "B", "C", "D")
dates_by_site_s2 <- split(md_s2$Date, md_s2$Site)

if (!all(required_sites_s2 %in% names(dates_by_site_s2))) {
  stop("Metadata must include sites A, B, C and D for the Venn diagram.")
}

shared_dates_s2 <- Reduce(
  intersect,
  dates_by_site_s2[required_sites_s2]
)

if (!length(shared_dates_s2)) {
  stop("No sampling dates are shared across sites A, B, C and D.")
}

ref_for_venn_s2 <- gene_rpm_s2 %>%
  dplyr::left_join(
    md_s2 %>% dplyr::select(BioSample, Site, Date),
    by = "BioSample"
  ) %>%
  dplyr::filter(
    Date %in% shared_dates_s2,
    Site %in% required_sites_s2,
    RPM > 0
  )

site_sets_s2 <- lapply(required_sites_s2, function(site_i) {
  ref_for_venn_s2 %>%
    dplyr::filter(Site == site_i) %>%
    dplyr::pull(GeneId) %>%
    unique()
})
names(site_sets_s2) <- required_sites_s2

venn_grob_s2 <- VennDiagram::venn.diagram(
  x = site_sets_s2,
  filename = NULL,
  fill = unname(COL[required_sites_s2]),
  alpha = 0.45,
  lwd = 1.5,
  cex = 1.2,
  cat.cex = 1.3,
  fontface = "bold",
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  main = "Broad resistome genes\n(shared sampling dates)"
)

grDevices::pdf(
  file.path(OUT, "FigureS2_broad_resistome_venn.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

grDevices::png(
  file.path(OUT, "FigureS2_broad_resistome_venn.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 400
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()


# Separate Venn diagram based on dates sampled at all four sites.
md_s2 <- md %>%
  dplyr::mutate(
    Date = as.Date(Date),
    Site = as.character(Site),
    BioSample = as.character(BioSample)
  )

required_sites_s2 <- c("A", "B", "C", "D")
dates_by_site_s2 <- split(md_s2$Date, md_s2$Site)

if (!all(required_sites_s2 %in% names(dates_by_site_s2))) {
  stop("Metadata must include sites A, B, C and D for the Venn diagram.")
}

shared_dates_s2 <- Reduce(intersect, dates_by_site_s2[required_sites_s2])
if (!length(shared_dates_s2)) {
  stop("No sampling dates are shared across sites A, B, C and D.")
}

ref_for_venn_s2 <- ref %>%
  dplyr::select(BioSample, GeneId, RPM) %>%
  dplyr::group_by(BioSample, GeneId) %>%
  dplyr::summarise(RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(
    md_s2 %>% dplyr::select(BioSample, Site, Date),
    by = "BioSample"
  ) %>%
  dplyr::filter(
    Date %in% shared_dates_s2,
    Site %in% required_sites_s2,
    RPM > 0
  )

site_sets_s2 <- lapply(required_sites_s2, function(site_i) {
  ref_for_venn_s2 %>%
    dplyr::filter(Site == site_i) %>%
    dplyr::pull(GeneId) %>%
    unique()
})
names(site_sets_s2) <- required_sites_s2

venn_grob_s2 <- VennDiagram::venn.diagram(
  x = site_sets_s2,
  filename = NULL,
  fill = unname(COL[required_sites_s2]),
  alpha = 0.45,
  lwd = 1.5,
  cex = 1.2,
  cat.cex = 1.3,
  fontface = "bold",
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  main = "Broad resistome genes\n(shared sampling dates)"
)

grDevices::pdf(
  file.path(OUT, "FigureS2_broad_resistome_venn.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

grDevices::png(
  file.path(OUT, "FigureS2_broad_resistome_venn.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 400
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

# Separate Venn diagram using sampling dates represented at all four sites.
md_s2 <- md %>%
  dplyr::mutate(
    Date = as.Date(Date),
    Site = as.character(Site),
    BioSample = as.character(BioSample)
  )

required_sites_s2 <- c("A", "B", "C", "D")
dates_by_site_s2 <- split(md_s2$Date, md_s2$Site)

if (!all(required_sites_s2 %in% names(dates_by_site_s2))) {
  stop("Metadata must include sites A, B, C and D for the Venn diagram.")
}

shared_dates_s2 <- Reduce(intersect, dates_by_site_s2[required_sites_s2])
if (!length(shared_dates_s2)) {
  stop("No sampling dates are shared across sites A, B, C and D.")
}

ref_for_venn_s2 <- ref %>%
  dplyr::select(BioSample, GeneId, RPM) %>%
  dplyr::group_by(BioSample, GeneId) %>%
  dplyr::summarise(RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(
    md_s2 %>% dplyr::select(BioSample, Site, Date),
    by = "BioSample"
  ) %>%
  dplyr::filter(
    Date %in% shared_dates_s2,
    Site %in% required_sites_s2,
    RPM > 0
  )

site_sets_s2 <- lapply(required_sites_s2, function(site_i) {
  ref_for_venn_s2 %>%
    dplyr::filter(Site == site_i) %>%
    dplyr::pull(GeneId) %>%
    unique()
})
names(site_sets_s2) <- required_sites_s2

venn_grob_s2 <- VennDiagram::venn.diagram(
  x = site_sets_s2,
  filename = NULL,
  fill = unname(COL[required_sites_s2]),
  alpha = 0.45,
  lwd = 1.5,
  cex = 1.2,
  cat.cex = 1.3,
  fontface = "bold",
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  main = "Broad resistome genes\n(shared sampling dates)"
)

grDevices::pdf(
  file.path(OUT, "FigureS2_broad_resistome_venn.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

grDevices::png(
  file.path(OUT, "FigureS2_broad_resistome_venn.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 400
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

# Venn diagram: gene presence at sites on dates sampled at all four sites.
md_s2 <- md %>%
  dplyr::mutate(
    Date = as.Date(Date),
    Site = as.character(Site),
    BioSample = as.character(BioSample)
  )

required_sites_s2 <- c("A", "B", "C", "D")
dates_by_site_s2 <- split(md_s2$Date, md_s2$Site)

if (!all(required_sites_s2 %in% names(dates_by_site_s2))) {
  stop("Metadata must include sites A, B, C and D for the Venn diagram.")
}

shared_dates_s2 <- Reduce(
  intersect,
  dates_by_site_s2[required_sites_s2]
)

if (!length(shared_dates_s2)) {
  stop("No sampling dates are shared across sites A, B, C and D.")
}

ref_for_venn_s2 <- ref %>%
  dplyr::select(BioSample, GeneId, RPM) %>%
  dplyr::group_by(BioSample, GeneId) %>%
  dplyr::summarise(RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(
    md_s2 %>% dplyr::select(BioSample, Site, Date),
    by = "BioSample"
  ) %>%
  dplyr::filter(
    Date %in% shared_dates_s2,
    Site %in% required_sites_s2,
    RPM > 0
  )

site_sets_s2 <- lapply(required_sites_s2, function(site_i) {
  ref_for_venn_s2 %>%
    dplyr::filter(Site == site_i) %>%
    dplyr::pull(GeneId) %>%
    unique()
})
names(site_sets_s2) <- required_sites_s2

venn_grob_s2 <- VennDiagram::venn.diagram(
  x = site_sets_s2,
  filename = NULL,
  fill = unname(COL[required_sites_s2]),
  alpha = 0.75,
  lwd = 1.5,
  cex = 2.5,
  cat.cex = 2.5,
  fontface = "bold",
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08
)

# Save the Venn diagram as its own figure.
grDevices::pdf(
  file.path(OUT, "FigureS2_broad_resistome_venn.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

grDevices::png(
  file.path(OUT, "FigureS2_broad_resistome_venn.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 400
)
grid::grid.newpage()
grid::grid.draw(venn_grob_s2)
grDevices::dev.off()

# S10: gene-by-date support, B and D only; a filter contributes at most once
# per gene even if several catalogue reference sequences pass in that filter.
fh_path <- paste0(
  "/home/sersancar/Documents/AREST/results/",
  "revision_20261030_214057/tables/",
  "S3_CARPDM_filter_all_available_hits.tsv"
)

required <- c(
  "BioSample", "SampleID", "GeneId", "Class", "Pass_strict"
)

fh <- readr::read_tsv(
  fh_path,
  col_types = readr::cols(.default = readr::col_character()),
  na = c("", "NA"),
  trim_ws = TRUE,
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(readr::problems(fh)) > 0L) {
  stop("Parsing problems in CARPDM input; inspect readr::problems(fh).")
}

missing_columns <- setdiff(required, names(fh))

if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

fh <- fh %>%
  dplyr::select(dplyr::all_of(required))

# Convert strict-pass flags explicitly; reject unexpected values.
pass_text <- toupper(trimws(fh$Pass_strict))

invalid <- !is.na(pass_text) &
  !pass_text %in% c("TRUE", "FALSE", "1", "0", "PASS", "FAIL", "T", "F")

if (any(invalid)) {
  stop(
    "Unexpected Pass_strict values: ",
    paste(unique(pass_text[invalid]), collapse = ", ")
  )
}

fh$Pass_strict <- dplyr::case_when(
  pass_text %in% c("TRUE", "1", "PASS", "T") ~ TRUE,
  pass_text %in% c("FALSE", "0", "FAIL", "F") ~ FALSE,
  TRUE ~ NA
)

message(
  "Loaded ", nrow(fh), " rows; ",
  sum(fh$Pass_strict, na.rm = TRUE), " strict-passing reference hits."
)
if(!is.logical(fh$Pass_strict)||anyNA(fh$Pass_strict)) stop('Pass_strict must contain TRUE/FALSE without missing values.')
support_hits <- fh %>% filter(Pass_strict) %>% distinct(BioSample,SampleID,GeneId) %>%
 count(BioSample,GeneId,name='NFilters')
if(any(support_hits$NFilters>3)) stop('More than three filters support a gene: inspect sample IDs.')
bd <- md %>% filter(Site %in% c('B','D')) %>% select(BioSample,Site,Date,DateLabel)
keep_genes <- support_hits %>% filter(BioSample %in% bd$BioSample) %>% distinct(GeneId)

carpdm_class_labels <- function(x) {
  dplyr::recode(
    as.character(x),
    "D & B res." = "Drug and biocide resistance",
    "Multi-D res." = "Multidrug resistance",
    "Cat. AM peptides" = "Cationic antimicrobial peptides",
    .default = as.character(x)
  )
}

gene_classes <- fh %>%
  dplyr::filter(GeneId %in% keep_genes$GeneId) %>%
  dplyr::mutate(Class = carpdm_class_labels(Class)) %>%
  dplyr::distinct(GeneId, Class) %>%
  dplyr::group_by(GeneId) %>%
  dplyr::summarise(
    Class = paste(
      sort(unique(Class[!is.na(Class) & nzchar(Class)])),
      collapse = " / "
    ),
    .groups = "drop"
  )

support <- tidyr::crossing(BioSample=bd$BioSample,GeneId=keep_genes$GeneId) %>%
 left_join(support_hits,by=c('BioSample','GeneId')) %>% mutate(NFilters=replace_na(NFilters,0L)) %>%
 left_join(bd,by='BioSample') %>% left_join(gene_classes,by='GeneId')
# Same genes in both panels, ordered by class then gene name for easy comparison.
support_order <- gene_classes %>% arrange(Class,GeneId) %>% pull(GeneId)
support <- support %>% mutate(GeneId=factor(GeneId,levels=rev(support_order)),
 Support=factor(NFilters,levels=0:3)) %>% filter(Site == B | Site == D)
export(support,'FigureS10_gene_filter_support')

SITE_COLS <- c(
  B = "#E64B35",
  D = "#4DBBD5"
)

psup <- ggplot(support,aes(DateLabel,GeneId,fill=Support)) +
 geom_tile(colour='grey80',linewidth=.25) +
  ggh4x::facet_grid2(
    . ~ Site,
    drop = TRUE,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = unname(SITE_COLS[c("B", "D")]),
        colour = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = c("white", "black"),
        face = "bold",
        size = 18
      )
    )
  ) +
 scale_fill_manual(values=c('0'='white','1'='#C6DBEF','2'='#6BAED6','3'='#08519C'),
 breaks=c('0','1','2','3'),labels=c('0 replicates','1 replicate','2 replicates','3 replicates'),
 drop=FALSE,name='Replicates supporting gene') +
 labs(x=NULL,y='CARPD genes') +
 theme(axis.text.x=element_text(angle=45, hjust=1, size=16),
       axis.text.y=element_text(size=16),
       axis.title = element_text(face='bold',size=18),
       panel.grid=element_blank(),
       strip.text=element_text(face='bold',size=18),
       legend.position='right', 
       legend.text = element_text(size=16), 
       legend.title = element_text(size=18, face = 'bold'), 
       legend.key.size = grid::unit(1.2, 'cm'))

save_plot(psup,'FigureS10_gene_filter_support',12,max(8,.23*length(support_order)+2))

# ComplexHeatmap layouts restored from original code, with aligned biological data.
save_heat <- function(ht,name,w,h) {
 pdf(file.path(OUT,paste0(name,'.pdf')),width=w,height=h,useDingbats=FALSE)
 ComplexHeatmap::draw(ht,heatmap_legend_side='right',annotation_legend_side='right',merge_legends=TRUE); dev.off()
 png(file.path(OUT,paste0(name,'.png')),width=w,height=h,units='in',res=400)
 ComplexHeatmap::draw(ht,heatmap_legend_side='right',annotation_legend_side='right',merge_legends=TRUE); dev.off()
 message('Written: ',name)
}
heat_col <- function(m) {
 lim <- as.numeric(quantile(m,c(.02,.98),na.rm=TRUE)); lim<-c(floor(lim[1]),ceiling(lim[2]))
 if(lim[2]<=lim[1]) lim[2]<-lim[1]+1
 circlize::colorRamp2(c(lim[1],mean(lim),lim[2]),c('#2c7bb6','#ffffbf','#d7191c'))
}
site_annotation <- function(ids) ComplexHeatmap::HeatmapAnnotation(
 Site=factor(md$Site[match(ids,md$BioSample)],levels=names(COL)),col=list(Site=COL),
 annotation_name_gp=grid::gpar(fontsize=12,fontface='bold'))
types <- sort(unique(ref$Type)); TYPECOL <- setNames(grDevices::hcl.colors(length(types),'Set 2'),types)

# Figure 2A: class heatmap with Type stripe and site annotation.
cl <- ref %>% group_by(BioSample,Type,Class) %>% summarise(Value=sum(RPM),.groups='drop') %>%
 mutate(Feature=paste(Type,Class,sep=' | '))
classanno <- cl %>% distinct(Feature,Type,Class) %>% arrange(Type,Class)
w <- cl %>% select(Feature,BioSample,Value) %>% pivot_wider(names_from=BioSample,values_from=Value,values_fill=0)
hm <- as.matrix(w[,-1]); rownames(hm)<-w$Feature
hm <- log10(hm[classanno$Feature,order_ids,drop=FALSE]+1)
p2a <- ComplexHeatmap::rowAnnotation(Type=classanno$Type,col=list(Type=TYPECOL),show_annotation_name=FALSE)+
 ComplexHeatmap::Heatmap(hm,name='log10(RPM+1)',col=heat_col(hm),cluster_rows=FALSE,
 cluster_columns=FALSE,column_split=factor(md$Site[match(order_ids,md$BioSample)],levels=names(COL)),
 top_annotation=site_annotation(order_ids),row_split=classanno$Type,row_title=NULL,
 row_labels=gsub('_',' ',classanno$Class),column_labels=short_ids[order_ids],column_names_rot=90,
 row_names_gp=grid::gpar(fontsize=10),column_title_gp=grid::gpar(fontsize=16,fontface='bold'))
export(cl,'Figure2A_class_abundance'); save_heat(p2a,'Figure2A_class_heatmap',15,11)




# S5: categorical dates, no lines across the unsampled April-December interval.
ct <- cg %>% group_by(BioSample,Site,Date) %>% summarise(TotalRPM=sum(RPM),.groups='drop') %>%
 mutate(DateLabel=factor(format(Date,'%d %b %Y'),levels=DATE_LABEL))
export(ct,'FigureS5_totals')
p5a <- ggplot(
  filter(ct, Site %in% c("B", "D")),
  aes(DateLabel, log10(TotalRPM + 1), colour = Site)
) +
  geom_line(
    aes(group = interaction(Site, format(Date, "%m") == "12")),
    linewidth = .8
  ) +
  geom_point(size = 3) +
  scale_colour_manual(values = COL) +
  labs(
    x = NULL,
    y = expression(bold(log[10](Sigma~RPM + 1)))
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
    axis.text.y = element_text(size = 15), 
    axis.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size=16), 
    legend.title = element_text(size=18, face = 'bold'), 
  )

top <- cg %>% group_by(GeneId) %>% summarise(Total=sum(RPM),.groups='drop') %>% arrange(desc(Total),GeneId) %>% slice_head(n=TOP_TEMPORAL_GENES) %>% pull(GeneId)
temporal <- cg %>% filter(GeneId %in% top,Site %in% c('B','D')) %>% group_by(BioSample,Site,Date,GeneId,Class) %>%
 summarise(RPM=sum(RPM),.groups='drop') %>% mutate(DateLabel=factor(format(Date,'%d %b %Y'),levels=DATE_LABEL),GeneId=factor(GeneId,levels=rev(top)))
export(temporal,'FigureS5_top_genes')
pseudo <- min(temporal$RPM[temporal$RPM>0])/2
 temporal$Fill <- ifelse(temporal$RPM>0,log10(temporal$RPM+pseudo),NA_real_)
 temporal$GeneLabel <- paste(temporal$GeneId,temporal$Class,sep=' | ')
 temporal$GeneLabel <- factor(temporal$GeneLabel,levels=rev(unique(temporal$GeneLabel[order(match(temporal$GeneId,top))])))
 export(tibble(Pseudocount=pseudo),'FigureS5_pseudocount')
 temporal$Site <- factor(temporal$Site, levels = c("B", "D"))
 
 fill_values <- temporal$Fill[is.finite(temporal$Fill)]
 
 if (!length(fill_values)) {
   stop("No finite abundance values available for Figure 5b.")
 }
 
 # Robust colour limits; values outside these limits are saturated.
 lims <- as.numeric(
   quantile(fill_values, probs = c(.02, .98), names = FALSE)
 )
 
 if (lims[1] == lims[2]) {
   lims <- lims + c(-0.5, 0.5)
 }
 
 mid <- mean(lims)
 
 p5b <- ggplot(
   temporal,
   aes(DateLabel, GeneLabel, fill = Fill)
 ) +
   geom_tile(colour = "white", linewidth = .25) +
   
   ggh4x::facet_grid2(
     Site ~ .,
     strip = ggh4x::strip_themed(
       background_y = ggh4x::elem_list_rect(
         fill = unname(COL[c("B", "D")]),
         colour = NA
       ),
       text_y = ggh4x::elem_list_text(
         colour = c("black", "black"),
         face = c("bold", "bold"),
         size = c(16, 16)
       )
     )
   ) +
   scale_fill_gradientn(
     colours = c("#2c7bb6", "#ffffbf", "#d7191c"),
     values = scales::rescale(
       c(lims[1], mid, lims[2]),
       from = lims
     ),
     limits = lims,
     oob = scales::squish,
     na.value = "white",
     name = expression(log[10](RPM + epsilon))
   )  +
   labs(x = NULL, y = NULL) +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
     axis.text.y = element_text(size = 15), 
     legend.text = element_text(size=16), 
     legend.title = element_text(size=18, face = 'bold'), 
     legend.key.size = grid::unit(1.2, 'cm')
   )
 

save_plot(p5a/p5b+plot_layout(heights=c(1,2)),'FigureS5_temporal_CARPDM',12,15)

# S6: sequencing depth and diversity. QC is descriptive and cannot prove equal sensitivity.
alpha <- bind_rows(lapply(names(MATS),function(metric) {
 x<-MATS[[metric]]
 md %>% mutate(Metric=metric,Total=rowSums(x),Richness=specnumber(x),q1=exp(diversity(x,'shannon')),q2=diversity(x,'invsimpson'))
}))
export(alpha,'FigureS6_alpha_all_metrics')
a <- filter(alpha,Metric=='RPM')
qcplot <- function(v,ylab,title) {
 z<-cor.test(a$TotalReads,a[[v]],method='spearman',exact=FALSE)
 ggplot(a,aes(TotalReads,if(v=='Total') .data[[v]]+1 else .data[[v]],colour=Site))+
 geom_point(size=3)+scale_colour_manual(values=COL)+scale_x_log10()+
 {if(v=='Total') scale_y_log10()}+
 labs(title=title,x='Pooled sequencing reads (log10 axis)',y=if(v=='Total') 'Total RPM + 1 (log10 axis)' else ylab,
 subtitle=sprintf('Descriptive Spearman rho = %.2f; P = %.3g; n = %d',unname(z$estimate),z$p.value,nrow(a)))
}
p6a<-qcplot('Total','Total RPM','Sequencing depth and abundance') + 
  theme(legend.position = 'none')
p6b<-qcplot('Richness','Detected gene richness','Sequencing depth and richness') + 
  theme(legend.position = 'none')
div<-a %>% select(BioSample,Site,Date,Richness,q1,q2) %>% pivot_longer(c(Richness,q1,q2),names_to='Metric',values_to='Value') %>%
 mutate(Metric=factor(Metric,levels=c('Richness','q1','q2'),labels=c('Richness (q = 0)','Exponential Shannon (q = 1)','Inverse Simpson (q = 2)')))
p6c<-ggplot(div,aes(Site,Value,colour=Site,fill=Site))+
 geom_boxplot(data=filter(div,Site!='C'),outlier.shape=NA,alpha=.35,linewidth=.5)+
 scale_fill_manual(values=COL)+geom_point(position=position_jitter(width=.12,height=0,seed=1),size=2.5)+
 facet_wrap(~Metric,scales='free_y')+scale_colour_manual(values=COL)+
  labs(title='Biological-sample diversity',y='Effective number of genes',x='Site')
save_plot((p6a+p6b)/p6c,'FigureS6_QC_diversity',13,9)


# S7: balanced design and length-normalisation sensitivity (four panels).
shared <- md %>% count(Date) %>% filter(n==4) %>% pull(Date)
subsets <- list(ABD_six_dates=md$BioSample[md$Site %in% c('A','B','D')],ABCD_shared_dates=md$BioSample[md$Date %in% shared])
plots7<-list(); stat7<-list()
for(metric in names(MATS)) for(tag in names(subsets)) {
 key<-paste0('FigureS7_',metric,'_',tag)
 pp<-pcoa(metric,subsets[[tag]],paste(metric,if(tag=='ABD_six_dates') 'A/B/D: six dates' else 'A/B/C/D: two December dates'),key)
 st<-read_table(paste0('PERMANOVA_',metric,'_',tag),c('R2','F','Pr(>F)'))
 pp<-pp+labs(subtitle=sprintf('Site PERMANOVA: R2 = %.3f; F = %.2f; P = %.4g',st$R2[1],st$F[1],st[['Pr(>F)']][1]))
 plots7[[key]]<-pp; stat7[[key]]<-mutate(st,Metric=metric,Subset=tag)
}
export(bind_rows(stat7),'FigureS7_site_statistics')
save_plot(wrap_plots(plots7,ncol=2,guides='collect'),'FigureS7_design_RPM_RPKM',15,12)

# S8: use exact saved RPM models/tests from the validated run. No fresh P-values
# or selective significant-vector plotting. Rank-one fits plotted as RDA1 vs PC1;
# the vertical coordinate is UNCONSTRAINED, not a second environmental effect.
variables<-c('WaterTemperature','pH','Salinity'); s8<-list(); s8stats<-list()
pretty<-c(WaterTemperature='Water temperature',pH='pH',Salinity='Salinity')
for(v in variables) {
 f<-file.path(RUN,'audit',paste0('pRDA_',v,'.rds'))
 if(!file.exists(f)) stop('Missing saved pRDA model: ',f)
 obj<-readRDS(f); fit<-obj$model; test<-obj$test
 if(is.null(fit$CCA)||fit$CCA$rank!=1) stop('Expected a rank-one saved partial RDA for ',v)
 ids<-as.character(obj$metadata$BioSample)
 # Detect accidentally loading a model from another version of the data.
 if(!setequal(ids,md$BioSample)) stop('Saved pRDA metadata does not match biological samples.')
 ss<-as.data.frame(scores(fit,display='sites',choices=1:2,scaling=1))
 if(nrow(ss)!=length(ids)) stop('pRDA score/metadata mismatch.')
 z<-tibble(BioSample=ids,RDA1=ss[[1]],PC1=ss[[2]]) %>% left_join(md,by='BioSample')
 export(z,paste0('FigureS8_',v,'_scores'))
 pval<-test[['Pr(>F)']][1]; fval<-test$F[1]; adj<-RsquareAdj(fit)$adj.r.squared
 st<-tibble(Variable=v,F=fval,P=pval,Adjusted_R2=adj,Model=paste(v,'+ Condition(Site + Date)'),Metric='RPM')
 s8stats[[v]]<-st
 s8[[v]]<-ggplot(z,aes(RDA1,PC1,colour=Site,shape=DateLabel))+geom_point(size=3)+
 scale_colour_manual(values=COL)+scale_shape_manual(values=SHAPES,name='Date')+
 labs(title=pretty[[v]],subtitle=sprintf('F = %.3f; raw P = %.4g; adjusted R2 = %.4f',fval,pval,adj),
 x='RDA1 (constrained by tested variable)',y='PC1 (unconstrained residual axis)') +
   theme(axis.title = element_text(face = "bold"))
}
s8tab<-bind_rows(s8stats) %>% mutate(Q_BH_three_tests=p.adjust(P,'BH'))
export(s8tab,'FigureS8_environment_tests')
save_plot(wrap_plots(s8,ncol=3,guides='collect'),'FigureS8_pRDA_site_date',18,6)

# Optional exploratory site-only comparison. No significance tests: shared date
# effects are not removed, and date-only weather cannot be tested within dates.
if(MAKE_SITE_ONLY_DESCRIPTIVE) {
 H<-decostand(MATS$RPM,'hellinger'); dd<-as.data.frame(md); dd$Site<-factor(dd$Site)
 pp<-list()
 for(v in variables) {
 dd$Predictor<-as.numeric(scale(dd[[v]])); fit<-rda(H~Predictor+Condition(Site),data=dd)
 ss<-scores(fit,display='sites',choices=1:2,scaling=1)
 z<-dd %>% mutate(RDA1=ss[,1],PC1=ss[,2])
 pp[[v]]<-ggplot(z,aes(RDA1,PC1,colour=Site,shape=DateLabel))+geom_point(size=3)+
 scale_colour_manual(values=COL)+scale_shape_manual(values=SHAPES,name='Date')+
 labs(title=pretty[[v]],subtitle='Site conditioned; date NOT removed; descriptive only',x='RDA1 (constrained)',y='PC1 (unconstrained)')
 }
 save_plot(wrap_plots(pp,ncol=3,guides='collect'),'FigureS8_extra_site_only_DESCRIPTIVE',18,6)
}
# Optional below-depth sensitivity only if the upstream analysis actually supplied it.
all_cb<-read_table('S3_CARPDM_biological_thresholds',c('Screen','BioSample','RPM'))
if('breadth' %in% all_cb$Screen) {
 z<-all_cb %>% group_by(Screen,BioSample) %>% summarise(RPM=sum(RPM),.groups='drop') %>% left_join(md,by='BioSample')
 p<-ggplot(z,aes(DateLabel,log10(RPM+1),colour=Screen,shape=Screen))+geom_point(position=position_dodge(.4),size=2.5)+facet_wrap(~Site)+
 labs(x=NULL,y='log10(total RPM + 1)',title='Additional CARPDM depth-threshold sensitivity')+theme(axis.text.x=element_text(angle=45,hjust=1))
 save_plot(p,'CARPDM_depth_sensitivity_UPSTREAM_AVAILABLE',12,8)
}
writeLines(c('Figure mapping and interpretation',
 '1A: original single-row class facets; B/D mean and descriptive t interval across six dates, including zeros.',
 '1B: restored taxonomy bubble plot; all genera by default; original-list subset optional and explicitly labelled. No host inference.',
 'S10: B/D gene-by-date heatmap; three blue levels for 1/2/3 distinct filters, white for no passing call; not depth-cutoff sensitivity.',
 '2A: summed reference RPM by Type and Class. Includes drug, metal and biocide resistance categories.',
 '2B: pooled samples; RPM PCoA with descriptive data ellipses where n>=4, no ellipse for C.',
 'S2: full gene heatmap plus labelled multipage version; no unequal-effort Venn inference.',
 'S5: original temporal line/heatmap layout; no within-date biological CI (n=1); lines broken between campaigns.',
 'S6: original log-depth scatter and boxplot layout; C shown as points (n=2); correlations descriptive.',
 'S7: balanced site contrasts and RPM/RPKM sensitivity; site dispersion differs for A/B/D.',
 'S8: saved RPM partial RDAs, one predictor each, conditioned on site and date; raw P values displayed.',
 'S8 table adds BH across all three water-chemistry tests; this differs from old two-test/single-test grouping.',
 'The optional site-only panel is descriptive and disabled by default.',
 'S8 does not implement RPKM environmental sensitivity.',
 'S9 is intentionally omitted: old class-environment regressions are not reproduced.',
 'No assembly, host assignment, below-threshold read recovery or absolute-load correction is performed.',
 paste('Inputs:',normalizePath(RUN)),paste('Created:',Sys.time())),file.path(OUT,'FIGURE_NOTES.txt'))
write_tsv(tibble(File=list.files(file.path(RUN,'tables'),full.names=TRUE),MD5=unname(tools::md5sum(list.files(file.path(RUN,'tables'),full.names=TRUE)))),file.path(OUT,'input_checksums.tsv'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
message('COMPLETE. Figures and plotting data: ',OUT)
