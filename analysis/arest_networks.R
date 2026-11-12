#!/usr/bin/env Rscript
# AREST S3: descriptive resistance-class / taxonomic-affiliation networks.
# Run while Cytoscape is open on this same computer:
# Rscript arest_FigureS3_networks.R [taxonomy_directory] [revision_directory] [output_directory]
# install.packages(c('dplyr','tidyr','readr','ggplot2'))
# BiocManager::install('RCy3')  # optional Cytoscape integration

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop('Supply taxonomy and revision directories.')
TAX_DIR <- args[1]
RUN <- args[2]
OUT <- if(length(args)>=3) args[3] else file.path(RUN,'FigureS3_networks_core_filtered_mean')
TAXON_RANK <- 'class'              # input column; its names may represent mixed taxonomic ranks
NODE_FILTER <- 'above_median'            # 'top_n', 'above_median', or 'none'
TOP_N <- 20L                     # separately per site and node type, after persistence
DOMAIN_MAP <- file.path(RUN,'taxon_domain_map.tsv') # optional Taxon,Domain columns
DOMAIN_COLS <- c(Bacteria='#0072B2',Archaea='#E69F00',Eukarya='#009E73',Viruses='#CC79A7')
CYTOSCAPE <- identical(toupper(Sys.getenv('AREST_CYTOSCAPE', unset='FALSE')), 'TRUE')                 # static R plots are always produced
MIN_MEAN_RPM <- 0                 # minimum edge mean RPM
MIN_BIO_SAMPLES <- 1L             # edge detection rule; node persistence is applied separately
BASE_URL <- "http://127.0.0.1:1234/v1"
SITE_ORDER <- c('A','B','C','D')
needed <- c('dplyr','tidyr','readr','ggplot2')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2)})
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

CLASS_SHORT <- c(
 Acetate_resistance='Acetate res.',Acid_resistance='Acid res.',
 Aldehyde_resistance='Aldehyde res.',Biguanide_resistance='Biguanide res.',
 'Multi-biocide_resistance'='Multi-B res.',Peroxide_resistance='Peroxide res.',
 Phenolic_compound_resistance='Phenolic res.',
 'Quaternary_Ammonium_Compounds_(QACs)_resistance'='QACs res.',
 Quaternary_ammonium_compound_resistance='QACs res.',
 Cationic_antimicrobial_peptides='Cat. AM peptides',
 Aminocoumarins='Aminocoumarins',Aminoglycosides='Aminoglycosides',
 Bacitracin='Bacitracin',Elfamycins='Elfamycins',Fluoroquinolones='Fluoroquinolones',
 Fosfomycin='Fosfomycin',Fusidic_acid='Fusidic acid',Glycopeptides='Glycopeptides',
 Lipopeptides='Lipopeptides',MLS='MLS',Metronidazole='Metronidazole',
 'Multi-drug_resistance'='Multi-D res.',Mupirocin='Mupirocin',
 'Mycobacterium_tuberculosis-specific_Drug'='M. tuberculosis spec. D',
 Nucleosides='Nucleosides',Oxazolidinone='Oxazolidinone',Phenicol='Phenicol',
 Pleuromutilin='Pleuromutilin',Rifampin='Rifampin',
 Spiropyrimidinetriones='Spiropyrimidinetriones',Sulfonamides='Sulfonamides',
 Tetracenomycin='Tetracenomycin',Tetracyclines='Tetracyclines',Thiopeptides='Thiopeptides',
 Trimethoprim='Trimethoprim',betalactams='Betalactams',
 Drug_and_biocide_resistance='D and B res.',Biocide_and_metal_resistance='B & M res.',
 Drug_and_biocide_and_metal_resistance='D, B & M res.',
 'Multi-metal_resistance'='Multi-M res.',Arsenic_resistance='As res.',
 Cadmium_resistance='Cd res.',Chromium_resistance='Cr res.',Copper_resistance='Cu res.',
 Iron_resistance='Fe res.',Lead_resistance='Pb res.',Mercury_resistance='Hg res.',
 Nickel_resistance='Ni res.',Sodium_resistance='Na res.',Tellurium_resistance='Te res.',
 Zinc_resistance='Zn res.')
# Stable shared colours across all sites; core colours match the revised ARG figures.
CORE_COLS <- c(Aminoglycosides='#0072B2',Betalactams='#D55E00',
 'Cat. AM peptides'='#009E73','D and B res.'='#CC79A7',Lipopeptides='#B28E00',
 MLS='#56B4E9','Multi-D res.'='#882255',Nucleosides='#4D4D4D',Phenicol='#332288',
 Sulfonamides='#8C564B',Tetracyclines='#117733',Trimethoprim='#AA4499')
canon <- function(x) {
 x <- toupper(trimws(x)); x <- gsub('FRESHWATER','FW',x,fixed=TRUE)
 x <- gsub('SEAWATER','SW',x,fixed=TRUE)
 x <- sub('^([123])FW','F\\1FW',x)
 x
}
read_checked <- function(path,cols) {
 if(!file.exists(path)) stop('Missing file: ',path)
 z <- read_tsv(path,show_col_types=FALSE,progress=FALSE)
 if(!all(cols %in% names(z))) stop('Missing columns in ',path,': ',paste(setdiff(cols,names(z)),collapse=', '))
 z
}
manifest <- read_checked(file.path(RUN,'tables','S1_library_manifest.tsv'),
 c('SampleID','BioSample','Site','TotalReads')) %>%
 mutate(Key=canon(SampleID),Site=as.character(Site),TotalReads=as.numeric(TotalReads))
if(anyDuplicated(manifest$Key)) stop('Non-unique canonical sample IDs in manifest.')
if(any(!is.finite(manifest$TotalReads)|manifest$TotalReads<=0)) stop('Invalid total-read denominators.')
md <- manifest %>% group_by(BioSample,Site) %>%
 summarise(TotalReads=sum(TotalReads),nFilters=n(),.groups='drop')
if(nrow(manifest)!=60 || nrow(md)!=20 || any(md$nFilters!=3)) stop('Expected 60 filters / 20 biological samples / three filters each.')
if(!setequal(md$Site,SITE_ORDER) || anyDuplicated(md$BioSample)) stop('Invalid biological sample/site mapping.')
files <- list.files(TAX_DIR,pattern='\\.ARG_taxonomy\\.tsv$',recursive=TRUE,full.names=TRUE)
if(!length(files)) stop('No broad-resistome .ARG_taxonomy.tsv files in ',TAX_DIR)
file_map <- tibble(Path=files,Key=canon(sub('\\.ARG_taxonomy\\.tsv$','',basename(files))))
write_tsv(anti_join(file_map,manifest,by='Key'),file.path(OUT,'unused_taxonomy_files.tsv'))
file_map <- inner_join(manifest,file_map,by='Key')
if(anyDuplicated(file_map$Key)) stop('More than one taxonomy file matches a filter. Remove duplicate copies or narrow TAX_DIR.')
missing_files <- anti_join(manifest,file_map,by='Key')
if(nrow(missing_files)) {write_tsv(missing_files,file.path(OUT,'missing_filters.tsv'));stop('Missing filter taxonomy files; see missing_filters.tsv. Missing is not zero.')}
write_tsv(file_map,file.path(OUT,'input_manifest.tsv'))

parts <- vector('list',nrow(file_map))
for(i in seq_len(nrow(file_map))) {
 z <- read_checked(file_map$Path[i],c('ARG_ID',TAXON_RANK,'read_count'))
 counts <- suppressWarnings(as.numeric(z$read_count))
 if(any(!is.finite(counts)|counts<0)) stop('Invalid read_count in ',file_map$Path[i])
 # Parse MEGARes correctly: accession | Type | Class | Mechanism | GeneId | optional note.
 fields <- strsplit(as.character(z$ARG_ID),'|',fixed=TRUE)
 if(any(lengths(fields)<5)) stop('Non-MEGARes ARG_ID in ',file_map$Path[i],'; do not use CARPDM taxonomy files.')
 parts[[i]] <- tibble(SampleID=file_map$SampleID[i],
  Type=vapply(fields,`[`,character(1),2),Class=vapply(fields,`[`,character(1),3),
  Taxon=trimws(as.character(z[[TAXON_RANK]])),Reads=counts) %>%
  filter(!is.na(Taxon),!tolower(Taxon) %in% c('','unclassified','unassigned','na'),Reads>0) %>%
  group_by(SampleID,Type,Class,Taxon) %>% summarise(Reads=sum(Reads),.groups='drop')
}
filter_counts <- bind_rows(parts) %>%
 mutate(Type=ifelse(Type=='Multi-compount','Multi-compound',Type))
if(!nrow(filter_counts)) stop('No classified positive associations at the requested taxonomic rank.')
# These counts must represent non-overlapping assignments in the upstream taxonomy tables.
# Multiple rows at the chosen rank are summed (e.g. species sharing one class).
filter_counts <- filter_counts %>% left_join(manifest %>% select(SampleID,BioSample,Site),by='SampleID')
write_tsv(filter_counts,file.path(OUT,'filter_association_counts.tsv'))
pairs <- filter_counts %>% distinct(Type,Class,Taxon) %>% arrange(Type,Class,Taxon) %>% mutate(PairID=row_number())
observed <- filter_counts %>% left_join(pairs,by=c('Type','Class','Taxon')) %>%
 group_by(BioSample,PairID) %>% summarise(Reads=sum(Reads),.groups='drop')
# Include all dates, including genuine zero observations, before site summarisation.
bio <- tidyr::crossing(BioSample=md$BioSample,PairID=pairs$PairID) %>%
 left_join(observed,by=c('BioSample','PairID')) %>% mutate(Reads=replace_na(Reads,0)) %>%
 left_join(md,by='BioSample') %>% left_join(pairs,by='PairID') %>% mutate(RPM=1e6*Reads/TotalReads)
site_edges <- bio %>% group_by(Site,Type,Class,Taxon) %>%
 summarise(MeanRPM=mean(RPM),MedianRPM=median(RPM),nBio=n(),nPositive=sum(Reads>0),
  Prevalence=nPositive/nBio,TotalMappedReadEnds=sum(Reads),.groups='drop')
write_tsv(bio,file.path(OUT,'biological_association_RPM_all_zeros.tsv'))
write_tsv(site_edges,file.path(OUT,'site_associations_all.tsv'))
# NODE persistence is evaluated independently of individual edge recurrence.
# A taxon is present if ANY of its ARG-class affiliations is positive in that sample.
# An AMR class is present if ANY of its classified taxonomic affiliations is positive.
node_bio <- bind_rows(
 bio %>% group_by(Site,BioSample,Taxon) %>% summarise(NodeRPM=sum(RPM),.groups='drop') %>%
  transmute(Site,BioSample,id=paste0('taxon::',Taxon),RawLabel=Taxon,Kind='Taxon',NodeRPM),
 bio %>% group_by(Site,BioSample,Class) %>% summarise(NodeRPM=sum(RPM),.groups='drop') %>%
  transmute(Site,BioSample,id=paste0('class::',Class),RawLabel=Class,Kind='AMR class',NodeRPM))
node_stats <- node_bio %>% group_by(Site,id,RawLabel,Kind) %>%
 summarise(SumRPM=sum(NodeRPM),MeanRPM=mean(NodeRPM),nBio=n(),
  nPositive=sum(NodeRPM>0),Prevalence=nPositive/nBio,.groups='drop') %>%
 mutate(CoreNode=nPositive==nBio & nBio>0,
  Label=ifelse(Kind=='Taxon',RawLabel,unname(CLASS_SHORT[RawLabel])))

node_stats$Label[is.na(node_stats$Label)] <- gsub('_',' ',node_stats$RawLabel[is.na(node_stats$Label)],fixed=TRUE)
write_tsv(node_bio,file.path(OUT,'node_RPM_each_biological_sample.tsv'))
write_tsv(node_stats,file.path(OUT,'node_persistence_all.tsv'))
# Domain assignments: exact NCBI scientific-name queries, cached; manual overrides take precedence.
# Viruses is a display grouping alongside the three cellular domains.
domain_cache_file <- file.path(OUT,'taxon_domain_cache.tsv')
domain_cache <- tibble(Taxon=character(),Domain=character(),TaxID=character(),Source=character())
if(file.exists(domain_cache_file)) domain_cache <- read_tsv(domain_cache_file,
 col_types=cols(.default=col_character()),show_col_types=FALSE)
if(file.exists(DOMAIN_MAP)) {
 overrides <- read_tsv(DOMAIN_MAP,col_types=cols(.default=col_character()),show_col_types=FALSE)
 if(!all(c('Taxon','Domain') %in% names(overrides))) stop('Domain map needs Taxon and Domain columns.')
 overrides <- overrides %>% filter(!is.na(Domain),Domain!='') %>%
  mutate(Domain=ifelse(Domain=='Eukaryota','Eukarya',Domain))
 if(anyDuplicated(overrides$Taxon) || any(!overrides$Domain %in% names(DOMAIN_COLS))) stop('Invalid or duplicate manual domain assignments.')
 domain_cache <- bind_rows(overrides %>% transmute(Taxon,Domain,TaxID=NA_character_,Source='manual override'),
  anti_join(domain_cache,overrides,by='Taxon'))
}
taxa_needed <- sort(unique(node_stats$RawLabel[node_stats$Kind=='Taxon']))
missing_domains <- setdiff(taxa_needed,domain_cache$Taxon[!is.na(domain_cache$Domain)])
if(length(missing_domains) && !identical(toupper(Sys.getenv('AREST_ALLOW_LIVE_TAXONOMY', unset='FALSE')), 'TRUE')) {
 stop('Frozen domain map does not cover: ', paste(missing_domains, collapse=', '),
      '. Supply taxon_domain_map.tsv; live taxonomy is disabled for reproducibility.')
}
if(length(missing_domains)) {
 for(pkg in c('rentrez','XML')) if(!requireNamespace(pkg,quietly=TRUE)) {
  stop('Install domain-lookup packages: install.packages(c("rentrez", "XML")); or provide ',DOMAIN_MAP)
 }
 resolve_domain <- function(taxon) {
  tryCatch({
   result <- rentrez::entrez_search(db='taxonomy',
                                    term = paste0(
                                      '"', gsub('"', '', taxon, fixed = TRUE),
                                      '"[Scientific Name]'
                                    ),
    retmax=2)
   if(result$count!=1) stop('No unique exact scientific-name match')
   record <- rentrez::entrez_fetch(db='taxonomy',id=result$ids,rettype='xml',parsed=TRUE)
   ids <- XML::xpathSApply(record,'//LineageEx/Taxon/TaxId',XML::xmlValue)
   ids <- unique(c(ids,result$ids))
   roots <- c('2'='Bacteria','2157'='Archaea','2759'='Eukarya','10239'='Viruses')
   groups <- unique(unname(roots[intersect(ids,names(roots))]))
   if(length(groups)!=1) stop('Unresolved broad lineage')
   tibble(Taxon=taxon,Domain=groups,TaxID=as.character(result$ids),
          Source="NCBI Taxonomy exact scientific name; lineage-based domain")
  },error=function(err) {
   message('Domain unresolved for ',taxon,': ',conditionMessage(err))
   tibble(Taxon=taxon,Domain=NA_character_,TaxID=NA_character_,Source=conditionMessage(err))
  })
 }
 for(taxon in missing_domains) {
  message('Resolving domain: ',taxon)
  resolved <- resolve_domain(taxon)
  domain_cache <- bind_rows(filter(domain_cache,Taxon!=taxon),resolved)
  write_tsv(domain_cache,domain_cache_file)
  Sys.sleep(.4)
 }
}
write_tsv(domain_cache,domain_cache_file)
unresolved <- tibble(Taxon=taxa_needed) %>% left_join(domain_cache,by='Taxon') %>%
 filter(is.na(Domain) | !Domain %in% names(DOMAIN_COLS))
if(nrow(unresolved)) {
 write_tsv(unresolved,file.path(OUT,'unresolved_taxon_domains.tsv'))
 stop('Unresolved taxonomy: fill Taxon,Domain in ',DOMAIN_MAP,
      ' using unresolved_taxon_domains.tsv, then rerun. No guessed domain colours were used.')
}
# Keep the complete statistics, then restrict taxonomic candidates before selection.
node_stats <- node_stats %>%
 left_join(domain_cache %>% select(RawLabel=Taxon,Domain),by='RawLabel') %>%
 mutate(Domain=ifelse(Kind=='Taxon',Domain,NA_character_),
        EligibleDomain=Kind!='Taxon' | Domain %in% c('Bacteria','Archaea'))
node_stats_all <- node_stats
write_tsv(node_stats_all,file.path(OUT,'node_statistics_all_domains.tsv'))
node_stats <- node_stats %>% filter(EligibleDomain)
write_tsv(node_stats,file.path(OUT,'node_statistics_eligible_domains.tsv'))
write_tsv(tibble(Domain=names(DOMAIN_COLS),Colour=unname(DOMAIN_COLS)),
 file.path(OUT,'domain_palette.tsv'))

persistent_nodes <- node_stats %>% filter(CoreNode,MeanRPM>0)
if(!NODE_FILTER %in% c('top_n','above_median','none')) stop('Invalid NODE_FILTER.')
if(length(TOP_N)!=1 || is.na(TOP_N) || TOP_N<1) stop('TOP_N must be positive.')
core_nodes <- persistent_nodes %>% group_by(Site,Kind) %>%
 arrange(desc(MeanRPM),RawLabel,.by_group=TRUE)
if(NODE_FILTER=='top_n') core_nodes <- core_nodes %>% slice_head(n=TOP_N)
if(NODE_FILTER=='above_median') core_nodes <- core_nodes %>% filter(MeanRPM>median(MeanRPM))
core_nodes <- ungroup(core_nodes)
selection_audit <- node_stats_all %>% left_join(
 core_nodes %>% transmute(Site,id,Selected=TRUE),by=c('Site','id')) %>%
 mutate(Selected=replace_na(Selected,FALSE),Filter=NODE_FILTER,
        ExclusionReason=case_when(!EligibleDomain~'Excluded domain',!CoreNode~'Not persistent',
         Selected~'Selected',TRUE~'Abundance display filter'))
write_tsv(selection_audit,file.path(OUT,'node_selection_audit.tsv'))
if(!nrow(core_nodes)) stop('No nodes are present in all biological samples within a site.')
# Keep the previous edge criterion and original global edge-width reference.
edge_candidates <- site_edges %>% filter(MeanRPM>MIN_MEAN_RPM,nPositive>=MIN_BIO_SAMPLES) %>%
 mutate(ClassLabel=unname(CLASS_SHORT[Class]))
edge_candidates$ClassLabel[is.na(edge_candidates$ClassLabel)] <- gsub('_',' ',edge_candidates$Class[is.na(edge_candidates$ClassLabel)],fixed=TRUE)
if(!nrow(edge_candidates)) stop('No edges pass the edge display settings.')
max_rpm <- max(edge_candidates$MeanRPM)
max_node_rpm <- max(core_nodes$MeanRPM)
# Node size uses MEAN associated RPM over all dates and all partners before pruning.
core_nodes <- core_nodes %>% mutate(NodeSize=70*sqrt(MeanRPM/max_node_rpm))
plotted <- edge_candidates %>%
 semi_join(core_nodes %>% filter(Kind=='Taxon') %>% transmute(Site,Taxon=RawLabel),by=c('Site','Taxon')) %>%
 semi_join(core_nodes %>% filter(Kind=='AMR class') %>% transmute(Site,Class=RawLabel),by=c('Site','Class'))

write_tsv(core_nodes,file.path(OUT,'core_nodes_all_sites.tsv'))
# Include colours for isolated retained class nodes too.
labels <- sort(unique(c(edge_candidates$ClassLabel,core_nodes$Label[core_nodes$Kind=='AMR class'])))
palette <- setNames(grDevices::hcl.colors(length(labels),'Dark 3'),labels)
palette[intersect(names(CORE_COLS),labels)] <- CORE_COLS[intersect(names(CORE_COLS),labels)]
write_tsv(tibble(ClassLabel=names(palette),Colour=unname(palette)),file.path(OUT,'class_palette.tsv'))
write_tsv(plotted,file.path(OUT,'site_associations_plotted.tsv'))
# Global coordinate and edge-width scales: positions are comparable across sites.
taxa <- sort(unique(core_nodes$RawLabel[core_nodes$Kind=='Taxon'])); classes <- sort(unique(core_nodes$RawLabel[core_nodes$Kind=='AMR class']))
pos <- bind_rows(
 tibble(id=paste0('taxon::',taxa),x=0,y=seq(0,1,length.out=length(taxa))),
 tibble(id=paste0('class::',classes),x=1,y=seq(0,1,length.out=length(classes))))
# max_rpm remains the pre-pruning reference, preserving previous edge widths.
cy_available <- FALSE
if(CYTOSCAPE && requireNamespace('RCy3',quietly=TRUE)) {
 cy_available <- tryCatch({RCy3::cytoscapePing(base.url=BASE_URL);TRUE},error=function(e){message('Cytoscape unavailable: ',conditionMessage(e));FALSE})
}
if(CYTOSCAPE && !cy_available) message('Static R plots and Cytoscape import tables will be saved; open Cytoscape and install RCy3 to also export editable networks.')
cy_errors <- character()
for(site in SITE_ORDER) {
 e <- plotted %>% filter(Site==site) %>% mutate(
  source=paste0('taxon::',Taxon),target=paste0('class::',Class),interaction='read-level affiliation',
  Colour=unname(palette[ClassLabel]),Width=0.4+5*sqrt(MeanRPM/max_rpm))
 nodes <- core_nodes %>% filter(Site==site) %>%
  mutate(Colour=ifelse(Kind=='Taxon',unname(DOMAIN_COLS[Domain]),unname(palette[Label]))) %>%
  left_join(pos,by='id') %>%
  group_by(Kind) %>% arrange(desc(MeanRPM),RawLabel,.by_group=TRUE) %>%
  mutate(y=if(n()==1) .5 else seq(0,1,length.out=n())) %>% ungroup()
 if(!nrow(nodes)) {message('No persistent nodes for site ',site);next}
 # Keep persistent nodes even if pruning their neighbours leaves them isolated.
 incident <- bind_rows(
  e %>% transmute(id=source,RPM=MeanRPM),e %>% transmute(id=target,RPM=MeanRPM)) %>%
  group_by(id) %>% summarise(RetainedIncidentMeanRPM=sum(RPM),.groups='drop')
 nodes <- nodes %>% left_join(incident,by='id') %>%
  mutate(RetainedIncidentMeanRPM=replace_na(RetainedIncidentMeanRPM,0))
 edges <- e %>% select(source,target,interaction,MeanRPM,MedianRPM,nBio,nPositive,Prevalence,Colour,Width)
 write_tsv(nodes,file.path(OUT,paste0('S3_',site,'_nodes.tsv')))
 write_tsv(edges,file.path(OUT,paste0('S3_',site,'_edges.tsv')))
 coords <- e %>% left_join(nodes %>% select(id,x,y) %>% rename(source=id,x1=x,y1=y),by='source') %>%
  left_join(nodes %>% select(id,x,y) %>% rename(target=id,x2=x,y2=y),by='target')
 p <- ggplot() +
  geom_segment(data=coords,aes(x=x1,y=y1,xend=x2,yend=y2,colour=ClassLabel,linewidth=MeanRPM),alpha=.45)+
  scale_colour_manual(values=palette,drop=FALSE)+
  scale_linewidth_continuous(range=c(.2,2.5),limits=c(0,max_rpm),trans='sqrt',name='Mean RPM')+
  geom_point(data=nodes,aes(x=x,y=y,shape=Kind,size=MeanRPM,fill=Colour),colour='#333333',stroke=.3)+
  scale_fill_identity(name='Taxonomic domain/group',breaks=unname(DOMAIN_COLS[c('Bacteria','Archaea')]),labels=c('Bacteria','Archaea'),guide='legend')+
  scale_size_area(max_size=12,limits=c(0,max_node_rpm),name='Node mean associated RPM')+
  scale_shape_manual(values=c(Taxon=21,'AMR class'=22))+
  geom_text(data=filter(nodes,Kind=='Taxon'),aes(x=x-.025,y=y,label=Label),hjust=1,size=3)+
  geom_text(data=filter(nodes,Kind=='AMR class'),aes(x=x+.025,y=y,label=Label),hjust=0,size=3)+
  coord_cartesian(xlim=c(-.55,1.65),ylim=c(-.03,1.03),clip='off')+theme_void()+
  labs(title=paste('Site',site),subtitle=paste(unique(nodes$nBio),'biological samples | persistent nodes | bacterial/archaeal groups; input column:',TAXON_RANK),
   caption='Read-level taxonomic affiliations; edges do not establish ARG hosts or genomic linkage.')+
  theme(legend.position='bottom',plot.title=element_text(face='bold',size=18),plot.margin=margin(15,20,15,20))+
  guides(colour='none',fill=guide_legend(override.aes=list(shape=21,size=4)))
 h <- max(9,.30*max(table(nodes$Kind)))
 ggsave(file.path(OUT,paste0('FigureS3_',site,'.pdf')),p,width=15,height=h,limitsize=FALSE)
 ggsave(file.path(OUT,paste0('FigureS3_',site,'.png')),p,width=15,height=h,dpi=300,limitsize=FALSE,bg='white')
 if(cy_available) tryCatch({
  net <- RCy3::createNetworkFromDataFrames(as.data.frame(nodes),as.data.frame(edges),
   title=paste0('AREST_S3_core_filtered_',site,'_',TAXON_RANK),collection='AREST S3 persistent abundance-filtered nodes',base.url=BASE_URL)
  maps <- list(
   RCy3::mapVisualProperty('node label','Label','p',network=net,base.url=BASE_URL),
   RCy3::mapVisualProperty('node size','NodeSize','p',network=net,base.url=BASE_URL),
   RCy3::mapVisualProperty('node fill color','Colour','p',network=net,base.url=BASE_URL),
   RCy3::mapVisualProperty('node shape','Kind','d',c('Taxon','AMR class'),c('ELLIPSE','RECTANGLE'),network=net,base.url=BASE_URL),
   RCy3::mapVisualProperty('edge stroke unselected paint','Colour','p',network=net,base.url=BASE_URL),
   RCy3::mapVisualProperty('edge width','Width','p',network=net,base.url=BASE_URL))
  style <- paste0('AREST_S3_',site,'_',format(Sys.time(),'%H%M%S'))
  RCy3::createVisualStyle(style,defaults=list(NETWORK_BACKGROUND_PAINT='#FFFFFF',NODE_SIZE=25,NODE_SIZE_LOCKED=TRUE,
   NODE_LABEL_FONT_SIZE=13,EDGE_TRANSPARENCY=120,EDGE_TARGET_ARROW_SHAPE='NONE',EDGE_SOURCE_ARROW_SHAPE='NONE'),
   mappings=maps,base.url=BASE_URL)
  RCy3::setVisualStyle(style,network=net,base.url=BASE_URL)
  RCy3::setNodePositionBypass(nodes$id,new.x.locations=nodes$x*900,
   new.y.locations=nodes$y*max(900,30*max(table(nodes$Kind))),network=net,base.url=BASE_URL)
  RCy3::fitContent(network=net,base.url=BASE_URL)

  # Export a vector PDF.
  export_base <- file.path(
    OUT,
    paste0("FigureS3_", site, "_Cytoscape")
  )

  RCy3::fitContent(
    network = net,
    base.url = BASE_URL
  )

  RCy3::exportImage(
    filename = export_base,
    type = "PDF",
    network = net,
    base.url = BASE_URL,
    overwriteFile = TRUE
  )

  # Convert PDF to PNG without using Cytoscape's bitmap exporter.
  pdf_file <- paste0(export_base, ".pdf")

  if (!file.exists(pdf_file)) {
    stop("Cytoscape PDF was not created: ", pdf_file)
  }

  converter <- Sys.which("pdftoppm")
  if (!nzchar(converter)) {
    stop("Install the converter with: sudo apt install poppler-utils")
  }

  conversion_status <- system2(
    converter,
    args = c(
      "-png",
      "-r", "600",
      "-singlefile",
      shQuote(pdf_file),
      shQuote(export_base)
    )
  )

  if (conversion_status != 0L) {
    stop("PDF-to-PNG conversion failed for site ", site)
  }

 },error=function(err){cy_errors <<- c(cy_errors,paste(site,conditionMessage(err)));message('Cytoscape export failed for ',site,': ',conditionMessage(err))})
}
writeLines(cy_errors,file.path(OUT,'Cytoscape_errors.txt'))
writeLines(c('Networks describe same-read ARG/taxonomy affiliations, not proven hosts, co-occurrence correlations or physical linkage.',
 paste('Taxonomic rank:',TAXON_RANK),
 'All 60 filter files are required. Counts are summed across filters, divided by summed total-read denominators, then averaged equally across biological site-date samples including zeros.',
 'Taxonomy count rows must represent non-overlapping read assignments upstream; this script cannot establish that provenance.',
 'Taxonomic labels come from the configured input column and may have mixed ranks; describe them as taxonomic groups. No pathogen-list selection is applied.',
 paste('Display thresholds: mean RPM >',MIN_MEAN_RPM,'and at least',MIN_BIO_SAMPLES,'positive biological samples.'),
 'Only nodes present in every pooled biological sample within their own site are retained: six dates at A/B/D and two at C.',
 'Node persistence is marginal: individual edges need not recur on every date. Only edges connecting two retained nodes are displayed.',
 'Node SumRPM sums that node\'s associated RPM across all partners and all dates BEFORE pruning neighbours. It is not total organism abundance.',
 'Node size uses mean associated RPM across biological samples, including zeros, before pruning partners. AMR-class totals retain all original classified taxonomic partners, including excluded domains; exclusion is a display rule, not recalculation of AMR abundance. A common scale is used across sites and node types.',
 paste('Node abundance filter:',NODE_FILTER,'; TOP_N =',TOP_N,'; applied separately by site and node type after persistence. Ties in top-N are broken alphabetically.'),
 'Above-median filtering uses strict > median within each site and node kind, after domain restriction and persistence. A singleton pool or entirely tied pool retains no nodes; there is no automatic fallback.',
 'Only Bacteria and Archaea are eligible taxonomic nodes. Eukarya and Viruses are excluded before computing taxonomic-node medians; complete input statistics are retained.',
 'Domain assignments, provenance and manual overrides are recorded in taxon_domain_cache.tsv. AMR class colours and edge-width reference are retained.',
 'Persistent isolated nodes are retained. RetainedIncidentMeanRPM reports the sum of displayed incident edge means.',
 'Nodes are arranged within each site by decreasing mean associated RPM. Edge widths retain the original shared unpruned-network reference.',
 'Class colours are shared across sites; the exported class_palette.tsv records the exact mapping.',
 'A/B/D have six biological samples, C two: compare recurrence counts with their denominators; these are descriptive networks.'),file.path(OUT,'ANALYSIS_NOTES.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
message('COMPLETE: ',normalizePath(OUT))

# Interactive layout adjustments run only when explicitly enabled.
if(CYTOSCAPE && cy_available) {

library(RCy3)

BASE_URL <- "http://127.0.0.1:1234/v1"

X_FACTOR <- 0   # Gap between the two node columns
Y_FACTOR <- 1.8   # Vertical spacing between nodes

networks <- unlist(
  getNetworkList(base.url = BASE_URL),
  use.names = FALSE
)

networks <- networks[
  grepl("^AREST_S3_core_filtered_", networks)
]

if (!length(networks)) {
  stop("No matching S3 networks are open in Cytoscape.")
}

# Retrieve and validate all coordinates before changing anything.
positions_all <- lapply(networks, function(net) {

  positions <- getNodePosition(
    network = net,
    base.url = BASE_URL
  )

  if (ncol(positions) != 2L || !nrow(positions)) {
    stop("Unexpected coordinate table for network: ", net)
  }

  x <- suppressWarnings(
    as.numeric(as.character(positions[[1]]))
  )
  y <- suppressWarnings(
    as.numeric(as.character(positions[[2]]))
  )

  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Invalid coordinates in network: ", net)
  }

  data.frame(
    node = rownames(positions),
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
})

names(positions_all) <- networks

for (net in networks) {

  positions <- positions_all[[net]]

  centre_x <- mean(positions$x)
  centre_y <- mean(positions$y)

  setNodePositionBypass(
    node.names = positions$node,
    new.x.locations =
      centre_x + (positions$x - centre_x) * X_FACTOR,
    new.y.locations =
      centre_y + (positions$y - centre_y) * Y_FACTOR,
    network = net,
    base.url = BASE_URL
  )

  fitContent(network = net, base.url = BASE_URL)

  message("Increased spacing: ", net)
}


}
