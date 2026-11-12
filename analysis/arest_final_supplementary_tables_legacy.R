#!/usr/bin/env Rscript
# Assemble revised tables. Does not rerun inferential analyses or mapping.
# Rscript arest_final_supplementary_tables.R RUN INPUT BASE_XLSX [OUTPUT_XLSX]
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(openxlsx)})
a <- commandArgs(TRUE)
if(length(a)<4) stop('Supply RUN INPUT BASE_XLSX OUTPUT_XLSX. This legacy packager requires the original workbook schema.')
RUN <- a[1]
INPUT <- a[2]
BASE <- a[3]
DEST <- a[4]
OUT <- file.path(RUN,'final_supplementary_tables')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
if(normalizePath(BASE,mustWork=FALSE)==normalizePath(DEST,mustWork=FALSE)) stop('Output must differ from source workbook.')
if(file.exists(DEST)) stop('Output already exists; choose a new output filename to preserve it.')
need <- function(x,cols,label) {if(!all(cols %in% names(x))) stop(label,': missing ',paste(setdiff(cols,names(x)),collapse=', '));x}
readtab <- function(p) {
 if(!file.exists(p)) stop('Missing required input: ',p)
 z<-read_tsv(p,show_col_types=FALSE,progress=FALSE,guess_max=100000)
 if(nrow(problems(z))) stop('Parsing errors: ',p)
 z
}
rt <- function(n) readtab(file.path(RUN,'tables',paste0(n,'.tsv')))
canon <- function(x) {
 x<-toupper(trimws(x));x<-gsub('FRESHWATER','FW',x,fixed=TRUE);x<-gsub('SEAWATER','SW',x,fixed=TRUE)
 sub('^([123])FW','F\\1FW',x)
}
uniq <- function(x,k,label) {if(anyDuplicated(x[k])) stop('Duplicate key in ',label);invisible(x)}
flag <- function(x) {
 v<-toupper(trimws(as.character(x)))
 if(anyNA(v)||any(!v %in% c('TRUE','FALSE','T','F','1','0','PASS','FAIL'))) stop('Invalid/missing threshold flag.')
 v %in% c('TRUE','T','1','PASS')
}
manifest<-rt('S1_library_manifest') %>% mutate(Key=canon(SampleID),Date=as.Date(Date))
need(manifest,c('Key','BioSample','Site','Date','TotalReads'),'manifest');uniq(manifest,'Key','manifest')
if(nrow(manifest)!=60||n_distinct(manifest$BioSample)!=20||any(!is.finite(manifest$TotalReads)|manifest$TotalReads<=0)) stop('Invalid 60-filter manifest.')
bio<-manifest %>% group_by(BioSample,Site,Date) %>% summarise(TotalReads=sum(TotalReads),nFilters=n(),.groups='drop')
if(any(bio$nFilters!=3)) stop('Expected three filters per biological sample.')
# S1: retain original environmental measurements, explicitly add biological IDs.
old1<-read.xlsx(BASE,sheet='S1_Metadata',detectDates=TRUE,check.names=FALSE)
old1<-need(old1,'Sample_name','original S1') %>% mutate(Key=canon(Sample_name))
uniq(old1,'Key','original S1')
if(!setequal(old1$Key,manifest$Key)) stop('Original S1 and revised manifest differ.')
S1<-manifest %>% left_join(old1 %>% select(-any_of(c('SampleID','BioSample','Site','Date','TotalReads'))),by='Key') %>% select(-Key)
S1bio<-rt('S1_biological_environment')
# S2: read actual rows only; the old sheet has formatting extending to Excel's limit.
S2<-read.xlsx(BASE,sheet='S2_previous_studies',skipEmptyRows=TRUE,check.names=FALSE)
S2<-S2[rowSums(!is.na(S2))>0,,drop=FALSE]
# S3: primary calls at filter/reference level retain breadth/depth and exact accession.
car<-rt('S3_CARPDM_filter_all_available_hits')
need(car,c('SampleID','BioSample','Site','gene','GeneId','Class','breadth','mean_depth','RPM','Pass_strict'),'CARPDM')
car$Pass_strict<-flag(car$Pass_strict)
uniq(car,c('SampleID','gene'),'CARPDM filter calls')
S3filter<-car %>% filter(Pass_strict) %>% arrange(Site,BioSample,SampleID,GeneId,gene)
cb<-rt('S3_CARPDM_biological_thresholds');need(cb,c('Screen','RPM','BioSample','gene'),'CARPDM pooled calls')
S3bio<-cb %>% filter(Screen=='strict',RPM>0)
# S4: all ranks, including unclassified rows; primary calls already gated upstream.
tx<-rt('S4_CARPDM_taxonomy_all_rows');need(tx,c('SampleID','ARG_ID','read_count'),'CARPDM taxonomy')
if(any(!is.finite(tx$read_count)|tx$read_count<0)) stop('Invalid taxonomy read counts.')
lookup<-car %>% distinct(SampleID,gene,GeneId,Class,Pass_strict)
# Remove metadata/annotations supplied by authoritative manifest/lookup before joining.
tx<-tx %>% select(-any_of(c('BioSample','Site','Date','TotalReads','GeneId','Class','Pass_strict','RPM')))
S4<-tx %>% inner_join(lookup %>% filter(Pass_strict),by=c('SampleID','ARG_ID'='gene')) %>%
 left_join(manifest %>% select(SampleID,BioSample,Site,Date,TotalReads),by='SampleID') %>% mutate(RPM=1e6*read_count/TotalReads)
if(anyNA(S4$BioSample)) stop('Unmatched CARPDM taxonomy libraries.')
# Gene reference lengths are retained for reference-level sensitivity validation.
S5<-rt('S5_MEGARes_reference_biological')
need(S5,c('gene_accession','GeneId','Type','Class','Mechanism','RPM','RPKM','length'),'MEGARes revised reference table')
if(any(!is.finite(S5$length)|S5$length<=0)||any(!is.finite(S5$RPKM))) stop('RPKM or reference lengths unavailable.')
# S6: preserve original read-level taxonomy at FILTER level with explicit BioSample.
# Do not infer exact reference accessions or read counts from rounded RPM.
S6<-read.xlsx(BASE,sheet='S6_ARG_taxonomy',check.names=FALSE)
need(S6,c('SampleID','RPM'),'original broad taxonomy');S6<-S6 %>% mutate(Key=canon(SampleID)) %>%
 select(-any_of(c('BioSample','Site','Date','TotalReads'))) %>%
 left_join(manifest %>% select(Key,BioSample,Site,Date,TotalReads),by='Key') %>% select(-Key)
if(anyNA(S6$BioSample)) stop('Unmatched broad taxonomy libraries.')
# S7: independent BAM-derived sensitivity table, complete reference-by-filter grid required.
sens<-readtab(file.path(INPUT,'CARPDM_depth_sensitivity_per_filter.tsv'))
need(sens,c('sample_id','reference_id','breadth_fraction','mean_depth',
 'original_MAPQ30_breadth80_depth5','sensitivity_MAPQ30_breadth80_no_depth_cutoff'),'CARPDM sensitivity')
sens<-sens %>% mutate(Key=canon(sample_id)) %>% left_join(manifest %>% select(Key,BioSample,Site,Date,SampleID),by='Key')
if(anyNA(sens$BioSample)||!setequal(sens$Key,manifest$Key)) stop('Sensitivity must cover all 60 filters.')
uniq(sens,c('Key','reference_id'),'sensitivity')
if(nrow(sens)!=60*n_distinct(sens$reference_id)) stop('Sensitivity needs all targets in all filters, including zeros.')
sens<-sens %>% mutate(Primary=flag(original_MAPQ30_breadth80_depth5),Relaxed=flag(sensitivity_MAPQ30_breadth80_no_depth_cutoff))
if(any(sens$Primary != (sens$breadth_fraction>=.8 & sens$mean_depth>=5))||any(sens$Relaxed != (sens$breadth_fraction>=.8))) stop('Threshold flags disagree with numeric metrics.')
sl<-sens %>% pivot_longer(c(Primary,Relaxed),names_to='Screen',values_to='Pass')
S7<-sl %>% group_by(Site,Screen) %>% summarise(Filters=n_distinct(Key),Biological_samples=n_distinct(BioSample),
 Filter_target_calls=sum(Pass),Biological_sample_target_calls=n_distinct(paste(BioSample[Pass],reference_id[Pass],sep='::')),
 Unique_targets=n_distinct(reference_id[Pass]),.groups='drop')
S7all<-sl %>% group_by(Screen) %>% summarise(Site='All',Filters=n_distinct(Key),Biological_samples=n_distinct(BioSample),
 Filter_target_calls=sum(Pass),Biological_sample_target_calls=n_distinct(paste(BioSample[Pass],reference_id[Pass],sep='::')),
 Unique_targets=n_distinct(reference_id[Pass]),.groups='drop')
S7<-bind_rows(S7,S7all)
# S8: per-filter quality metrics and pooled biological-sample RPM; no ARG/marker ratios.
cr<-readtab(file.path(INPUT,'crassphage_all_filter_metrics.tsv'))
need(cr,c('sample','mapq_threshold','mapped_read_ends','breadth_fraction','mean_depth'),'crAssphage metrics')
cr<-cr %>% mutate(Key=canon(sample)) %>% left_join(manifest %>% select(Key,BioSample,Site,Date,TotalReads),by='Key')
uniq(cr,'Key','crAssphage metrics')
if(nrow(cr)!=60||anyNA(cr$BioSample)||!setequal(cr$Key,manifest$Key)||any(cr$mapq_threshold!=30)) stop('Invalid crAssphage library set/MAPQ.')
if(any(!is.finite(cr$mapped_read_ends)|cr$mapped_read_ends<0)||any(!is.finite(cr$breadth_fraction)|cr$breadth_fraction<0|cr$breadth_fraction>1)||any(!is.finite(cr$mean_depth)|cr$mean_depth<0)) stop('Invalid marker metrics.')
cr<-cr %>% mutate(Pass=breadth_fraction>=.8 & mean_depth>=5)
S8bio<-cr %>% group_by(BioSample,Site,Date) %>% summarise(Filters=n(),Mapped_read_ends=sum(mapped_read_ends),
 TotalReads=sum(TotalReads),Passing_filters=sum(Pass),RPM=1e6*Mapped_read_ends/TotalReads,.groups='drop')
S8<-S8bio %>% group_by(Site) %>% summarise(Biological_samples=n(),Total_filters=sum(Filters),
 Passing_filters=sum(Passing_filters),Samples_any_passing=sum(Passing_filters>=1),Samples_all_three_passing=sum(Passing_filters==3),
 Mean_RPM=mean(RPM),Median_RPM=median(RPM),Minimum_pooled_read_ends=min(Mapped_read_ends),Maximum_pooled_read_ends=max(Mapped_read_ends),.groups='drop') %>%
 left_join(cr %>% group_by(Site) %>% summarise(Filters_with_alignments=sum(mapped_read_ends>0),
 Maximum_filter_breadth=max(breadth_fraction),Maximum_filter_depth=max(mean_depth),.groups='drop'),by='Site')
# S9: retain heterogeneous tests in separate component sheets with their original columns.
collect <- function(pattern) {
 f<-list.files(file.path(RUN,'tables'),pattern=pattern,full.names=TRUE)
 if(!length(f)) stop('No test tables matching ',pattern)
 # Numeric identifiers/metadata are not silently converted to text by mixing schemas.
 bind_rows(lapply(f,function(p) readtab(p) %>% mutate(Source_table=basename(p))))
}
S9perm<-collect('^PERMANOVA_.*\\.tsv$')
S9fried<-collect('^friedman_.*\\.tsv$')
S9pair<-collect('^paired_contrasts_.*\\.tsv$')
S9qc<-rt('sequencing_depth_QC')
# Dispersion diagnostics with n=2 stay descriptive, as specified by revised pipeline.
fdisp<-list.files(file.path(RUN,'tables'),pattern='^dispersion_R.*\\.tsv$',full.names=TRUE)
S9disp<-if(length(fdisp)) bind_rows(lapply(fdisp,function(p) readtab(p) %>% mutate(Source_table=basename(p)))) else NULL
S10<-rt('partial_RDA')
need(S10,c('Variable','p'),'partial RDA');S10<-S10 %>% mutate(Q_BH_all_three=p.adjust(p,'BH'))
S10audit<-rt('environment_identifiability')
# S11 is supplied by the author: do not invent historical software versions.
provpath<-file.path(INPUT,'software_reference_provenance.tsv')
if(!file.exists(provpath)) {
 write_tsv(tibble(Component=c('fastp','Kraken2','Kraken2 database','Bracken','AMR++','MEGARes','Bowtie2','CARPDM catalogue','samtools','R analysis'),
 Version=NA_character_,Reference_or_accession=NA_character_,Command_or_function=NA_character_,Source_log=NA_character_),
 file.path(OUT,'software_reference_provenance_TEMPLATE.tsv'))
 stop('Fill software_reference_provenance_TEMPLATE.tsv with actual run provenance, save as ',provpath,', and rerun. No final workbook written.')
}
S11<-readtab(provpath);need(S11,c('Component','Version','Reference_or_accession','Command_or_function','Source_log'),'provenance')
if(anyNA(S11$Version)||any(!nzchar(S11$Version))) stop('Missing actual software/reference versions.')
# Useful tidal transparency table: prediction is not an observed synchronised tidal stage.
tidecols<-grep('tid',names(S1bio),ignore.case=TRUE,value=TRUE)
S1tide<-if(length(tidecols)) S1bio %>% select(any_of(c('BioSample','Site','Date')),all_of(tidecols)) else NULL
sheets<-list(S1_Filter_metadata=S1,S1_Biological_samples=S1bio,S2_Previous_studies=S2,
 S3_CARPDM_filter_calls=S3filter,S3_CARPDM_biological=S3bio,S4_CARPDM_taxonomy=S4,
 S5_MEGARes_biological=S5,S6_Broad_taxonomy_filters=S6,S7_Threshold_summary=S7,S7_Threshold_filter_detail=sens,
 S8_Marker_site_summary=S8,S8_Marker_biological=S8bio,S8_Marker_filter_detail=cr,
 S9_PERMANOVA=S9perm,S9_Friedman=S9fried,S9_Paired_contrasts=S9pair,S9_Sequencing_depth=S9qc,
 S10_Environmental_tests=S10,S10_Covariate_audit=S10audit,S11_Software_references=S11)
if(!is.null(S9disp)) sheets$S9_Dispersion<-S9disp
if(!is.null(S1tide)) sheets$S1_Tidal_context<-S1tide
# Separate raw taxonomy affiliations from host identification; preserve all annotations.
notes<-tibble(Item=c('Sampling unit','RPM pooling','S3 primary calls','S4 and S6 taxonomy','S6 provenance','S7 targets',
 'S8 interpretation','S9 pairwise multiplicity','S10 multiplicity','Tidal interpretation','Historical S2','Gene length'),
 Definition=c('60 filter libraries represent 20 original site-date samples; filters are not biological replicates.',
 'Biological RPM = 1e6 * sum eligible counts / sum total reads across all three filters.',
 'Thresholds applied per filter. Breadth is a fraction. No quantitative limit of detection calibrated.',
 'Read-level affiliations do not identify ARG hosts, viability or physical linkage. Mixed taxonomic ranks may occur.',
 'S6 retains original filter-level RPM and adds biological IDs; it is not silently recalculated or assigned exact accessions.',
 'Targets are exact reference sequences; gene labels and references are different counting units. All-site union is deduplicated.',
 'Single-reference screen NC_024711.1. Low-coverage mapping does not establish absence or support reliable dilution inference.',
 'BH applied separately within each outcome and sampling subset, as in the saved tests. Friedman P values are asymptotic.',
 'Q_BH_all_three corrects the three displayed environmental tests together; original analysis-specific correction is also retained.',
 'Shared date-only tide values cannot demonstrate matched sampling-time tidal stages; sampling times/datum are needed.',
 'Measurements and isolates from previous studies are historical context, not same-sample validation.',
 'RPKM corrects reference length, not bacterial biomass, cells or absolute copies per litre.'))
sheets$Definitions<-notes
# Build a NEW workbook; retain the source unchanged. Also export every sheet as TSV.
wb<-createWorkbook()
header <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#234E70",
  textDecoration = "bold",
  wrapText = TRUE,
  valign = "top"
)
for(n in names(sheets)) {
 z<-as.data.frame(sheets[[n]])
 if(nrow(z)>1048575||ncol(z)>16384) stop('Excel limits exceeded: ',n,'; use TSV instead.')
 write_tsv(as_tibble(z),file.path(OUT,paste0(n,'.tsv')))
 addWorksheet(wb,n)
 writeData(wb,n,z,withFilter=TRUE,keepNA=FALSE,headerStyle=header)
 freezePane(wb,n,firstRow=TRUE)
 setColWidths(wb,n,cols=seq_len(ncol(z)),widths=22)
 setRowHeights(wb,n,rows=1,heights=38)
}
saveWorkbook(wb,DEST,overwrite=FALSE)
message('Written: ',DEST,'\nTSV tables: ',OUT)
