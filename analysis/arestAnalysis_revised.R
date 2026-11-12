#!/usr/bin/env Rscript
# AREST reviewer revision pipeline. Adapted from Sergio's arestAnalysis.R.
# Run: AREST_INPUT=/path/to/inputs AREST_RUN=/path/to/output Rscript analysis/arestAnalysis_revised.R
# No automatic package installation; install missing CRAN packages once.
# Biological unit = one original water sample per Site x Date, split over 3 filters.
# All inferential analyses use biological units, never filters as independent n.
# Scope: existing read/mapping tables. Assembly, MGE detection, faecal-marker
# mapping, MAPQ filtering and reference-length extraction must occur upstream.
# Configurable scientific choices are documented here; do not hide changes.

INPUT <- Sys.getenv('AREST_INPUT', unset='input')
OUTPUT <- Sys.getenv('AREST_OUTPUT', unset='analysis_results')
SEED <- 20261030L
NPERM <- 9999L
# TRUE only after verifying CARD_genes contains UNFILTERED gene summaries
# from MAPQ >=30 alignments, including genes failing breadth/depth thresholds.
CARPDM_UNFILTERED <- identical(toupper(Sys.getenv('AREST_CARPDM_UNFILTERED', unset='FALSE')), 'TRUE')
# Counts from AMR++ are assumed already subject to MAPQ>=30 and breadth>=0.8.
# CARPDM: primary is >=0.8 breadth AND >=5 mean depth on each filter.
# Counts passing those filters are summed, then divided by ALL 3 library totals.
# A sample is positive if >=1 filter passes; >=2 and all-3 are sensitivities.
# This preserves original library criteria; it is NOT a pooled-BAM depth call.
CARPDM_BREADTH <- 0.80
CARPDM_DEPTH <- 5

packages <- c('tidyverse','vegan','permute','ggsci','patchwork','openxlsx')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install required packages: install.packages(c(',
  paste(sprintf('"%s"',missing),collapse=', '), '))')
suppressPackageStartupMessages({library(tidyverse); library(vegan); library(permute)
  library(ggsci); library(patchwork)})
set.seed(SEED)
RUN <- Sys.getenv('AREST_RUN', unset=file.path(OUTPUT,paste0('revision_',format(Sys.time(),'%Y%m%d_%H%M%S'))))
for(d in c('tables','figures','audit','templates')) dir.create(file.path(RUN,d),recursive=TRUE)
status <- tibble(Module=character(),Status=character(),Detail=character())
record <- function(module,state,detail) {
  status <<- bind_rows(status,tibble(Module=module,Status=state,Detail=detail))
  write_tsv(status,file.path(RUN,'audit','module_status.tsv'))
  message(module,': ',state,' — ',detail)
}
write_tab <- function(x,name) write_tsv(as_tibble(x),file.path(RUN,'tables',paste0(name,'.tsv')))
save_plot <- function(p,name,w=10,h=7) {
  ggsave(file.path(RUN,'figures',paste0(name,'.pdf')),p,width=w,height=h)
  ggsave(file.path(RUN,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white')
}
cols <- setNames(ggsci::pal_npg()(4)[c(3,1,4,2)],c('A','B','C','D'))
theme_set(theme_bw(base_size=12))
read_t <- function(p) read_tsv(p,show_col_types=FALSE,progress=FALSE)
need <- function(x,n,context) {
  z <- setdiff(n,names(x)); if(length(z)) stop(context,': missing ',paste(z,collapse=', '))
}
unique_key <- function(x,key,context) {
  if(anyDuplicated(x[key])) stop(context,': duplicated ',paste(key,collapse=' + '))
}
locate <- function(name,required=TRUE) {
  z <- list.files(INPUT,recursive=TRUE,full.names=TRUE)
  z <- z[basename(z)==name]
  if(length(z)>1) stop('Ambiguous input: ',name,': ',paste(z,collapse='; '))
  if(!length(z)) {if(required) stop('Missing input: ',file.path(INPUT,name)); return(NULL)}
  z
}
# INPUT CONTRACT (based on original template):
# metadata.tsv: Sample_name; environment columns as in original script/workbook.
# total_reads.tsv: sample,total_reads (individual reads, consistent with Counts).
# dedup_AMR_analytic_matrix2.csv: gene_accession plus library count columns.
# CARD_genes/<SampleID>.<suffix>.tsv: gene,length,cov_bases,breadth,mean_depth,num_reads.
# Optional megares_lengths.tsv: gene_accession,length (EXACT reference ID).
# Optional sample_manifest.tsv: SampleID,Site,Date (YYYY-MM-DD),BioSample.
# Optional CARPDM taxonomy: *.CARD_ARG_taxonomy.tsv from original template.
# Dates/site can otherwise be parsed from original AREST IDs (F becomes D).
# README output records all omissions. Never use missing modules as negative results.

megares_class <- c(
  "AAC(3)-Ia"       = "Aminoglycosides",
  "AAC6_Ie_APH2_Ia" = "Aminoglycosides",
  "ANT(2'')-Ia"     = "Aminoglycosides",
  "ANT(3'')-IIa"    = "Aminoglycosides",
  "ANT(6)-Ia"       = "Aminoglycosides",
  "ANT(9)-Ia"       = "Aminoglycosides",
  "ANT3Ii_AAC6_IID" = "Aminoglycosides",
  "APH(3'')-Ib"     = "Aminoglycosides",
  "APH(3')-IIIa"    = "Aminoglycosides",
  "APH(3')-Ia"      = "Aminoglycosides",
  "APH(6)-Id"       = "Aminoglycosides",
  "aad(6)"          = "Aminoglycosides",
  "aadA"            = "Aminoglycosides",
  "aadA15"          = "Aminoglycosides",
  "aadA2"           = "Aminoglycosides",
  "aadA22"          = "Aminoglycosides",
  "aadA24"          = "Aminoglycosides",
  "aadA25"          = "Aminoglycosides",
  "aadA3"           = "Aminoglycosides",
  "aadA5"           = "Aminoglycosides",

  "EreA2"           = "MLS",
  "ErmB"            = "MLS",
  "ErmF"            = "MLS",
  "ErmT"            = "MLS",
  "ErmX"            = "MLS",
  "lnuG"            = "MLS",
  "mef(B)"          = "MLS",
  "mel"             = "MLS",
  "mphA"            = "MLS",
  "mphE"            = "MLS",
  "msrE"            = "MLS",

  "OXA-1"           = "Betalactams",
  "OXA-10"          = "Betalactams",
  "OXA-2"           = "Betalactams",
  "OXA-4"           = "Betalactams",
  "OXA-58"          = "Betalactams",
  "TEM-160"         = "Betalactams",
  "TEM-2"           = "Betalactams",
  "VEB-1"           = "Betalactams",
  "VEB-3"           = "Betalactams",

  "SAT-4"           = "Nucleosides",

  "cmlA5"           = "Phenicol",
  "cmlA9"           = "Phenicol",
  "catB3"           = "Phenicol",
  "floR"            = "Phenicol",

  "dfrA1"           = "Trimethoprim",
  "dfrA14"          = "Trimethoprim",
  "dfrF"            = "Trimethoprim",

  "sul1"            = "Sulfonamides",
  "sul2"            = "Sulfonamides",

  "tet(A)"          = "Tetracyclines",
  "tet(B)"          = "Tetracyclines",
  "tet(C)"          = "Tetracyclines",
  "tet(L)"          = "Tetracyclines",
  "tet(M)"          = "Tetracyclines",
  "tet(O)"          = "Tetracyclines",
  "tet(T)"          = "Tetracyclines",
  "tet(X3)"         = "Tetracyclines",
  "tet(X4)"         = "Tetracyclines",
  "tet(X6)"         = "Tetracyclines",

  "qacE"            = "D & B res.",
  "qacEdelta1"      = "D & B res.",
  "qacL"            = "D & B res.",

  "ICR-Mo"          = "Lipopeptides",
  "MCR-3.1"         = "Lipopeptides",
  "ugd"             = "Cat. AM peptides",

  "optrA"           = "Multi-D res."
)

# Standardise only the two documented naming conventions. Site F remains
# in library IDs; reporting Site D is a separate field. No input files are edited.
id_audit <- tibble(Source=character(),OriginalID=character(),SampleID=character())
normalise_ids <- function(x, source, strict=TRUE) {
  old <- as.character(x)
  y <- trimws(old)
  y <- gsub("Seawater", "SW", y, ignore.case=TRUE)
  y <- gsub("Freshwater", "FW", y, ignore.case=TRUE)
  y <- toupper(y)
  valid <- !is.na(y) & grepl("^(?:[ABC][123]SW|F[123]FW)[0-9]{6}$",y,perl=TRUE)
  if(strict && any(!valid)) {
    print(data.frame(Source=source,OriginalID=old[!valid],Normalised=y[!valid]))
    stop('Unexpected sample-ID format; inspect identifiers above.',call.=FALSE)
  }
  if(any(valid)) {
    id_audit <<- bind_rows(id_audit,tibble(Source=source,OriginalID=old[valid],SampleID=y[valid])) %>% distinct()
    write_tsv(id_audit,file.path(RUN,'audit','sample_ID_crosswalk.tsv'))
  }
  if(strict) y else ifelse(valid,y,old)
}
assert_read_join <- function(manifest,tr) {
  absent <- anti_join(manifest,tr,by='SampleID')
  if(nrow(absent)) {
    print(as.data.frame(absent),row.names=FALSE)
    write_tab(absent,'ERROR_missing_read_counts')
    stop('Identifiers above have no sequencing totals after ID normalisation.',call.=FALSE)
  }
  invalid <- tr %>% filter(SampleID %in% manifest$SampleID,
                          !is.finite(TotalReads) | TotalReads<=0)
  if(nrow(invalid)) {
    print(as.data.frame(invalid),row.names=FALSE)
    stop('Non-positive or missing sequencing totals.',call.=FALSE)
  }
}

main <- function() {
  meta <- read_t(locate('metadata.tsv'))
  need(meta,'Sample_name','metadata.tsv')
  names(meta)[names(meta)=='Sample_name'] <- 'SampleID'
  meta$SampleID <- normalise_ids(meta$SampleID,'metadata.tsv')
  unique_key(meta,'SampleID','metadata')
  tr <- read_t(locate('total_reads.tsv')); need(tr,c('sample','total_reads'),'total_reads')
  tr <- tr %>% transmute(SampleID=normalise_ids(sample,'total_reads.tsv'),TotalReads=as.numeric(total_reads))
  unique_key(tr,'SampleID','total_reads')
  manifest_file <- locate('sample_manifest.tsv',FALSE)
  if(is.null(manifest_file)) {
    manifest <- meta %>% transmute(SampleID,
      Site=recode(substr(SampleID,1,1),F='D'),
      Date=as.Date(str_extract(SampleID,'[0-9]{6}$'),format='%d%m%y')) %>%
      mutate(BioSample=paste(Site,Date,sep='_'))
  } else {
    manifest <- read_t(manifest_file); need(manifest,c('SampleID','Site','Date','BioSample'),'manifest')
    manifest <- manifest %>% mutate(SampleID=normalise_ids(SampleID,'sample_manifest.tsv'),
      Site=recode(as.character(Site),F='D'),Date=as.Date(Date))
  }
  unique_key(manifest,'SampleID','manifest')
  if(!setequal(meta$SampleID,manifest$SampleID)) stop('Manifest and metadata library sets differ.')
  if(anyNA(manifest) || any(!manifest$Site %in% c('A','B','C','D'))) stop('Invalid site/date/sample manifest.')
  assert_read_join(manifest,tr)
  lib <- manifest %>% left_join(tr,by='SampleID')
  if(any(!is.finite(lib$TotalReads) | lib$TotalReads<=0)) stop('Missing/invalid total reads.')
  counts <- lib %>% count(BioSample,name='nFilters')
  if(any(counts$nFilters!=3)) stop('Expected 3 filters per biological sample; inspect manifest.')
  if(nrow(lib)!=60 || n_distinct(lib$BioSample)!=20) stop('Expected 60 filters / 20 biological samples.')
  if(any((lib %>% distinct(BioSample,Site,Date) %>% count(BioSample))$n!=1)) stop('BioSample has conflicting site/date.')
  bio <- lib %>% group_by(BioSample,Site,Date) %>%
    summarise(TotalReads=sum(TotalReads),nFilters=n(),.groups='drop') %>% arrange(Date,Site)
  if(anyDuplicated(bio[c('Site','Date')])) stop('This design expects one biological sample per site/date.')
  write_tab(lib,'S1_library_manifest'); write_tab(bio,'S1_biological_samples')
  write_tab(bio %>% count(Site),'sample_counts')
  record('replication','DONE','60 filter libraries collapsed to 20 original site-date samples.')

  # Environment: names accepted from both original script and manuscript workbook.
  aliases <- list(Salinity=c('Salinity','Salinity (g/L)'),pH='pH',
    WaterTemperature=c('Sample_temperature','Sample temperature (°C)'),
    AirTemperature=c('Temperature','Temperature  (°C)'),Wind=c('Wind','Wind (Km/h)'),
    UV=c('Max_UV_index','Max UV index'),Humidity=c('Rel_humidity','Rel. humidity (%)'),
    Rainfall=c('Rainfall','Rainfall (mm)'),Tide=c('Tidal_level','Tidal_level (m)'))
  env <- lib %>% select(SampleID,BioSample,Site,Date)
  for(v in names(aliases)) {
    k <- intersect(aliases[[v]],names(meta))
    env[[v]] <- if(length(k)) suppressWarnings(as.numeric(meta[[k[1]]][match(env$SampleID,meta$SampleID)])) else NA_real_
  }
  env_audit <- env %>% pivot_longer(all_of(names(aliases)),names_to='Variable') %>%
    group_by(BioSample,Site,Date,Variable) %>% summarise(n_values=n_distinct(value,na.rm=TRUE),
      missing=sum(is.na(value)),.groups='drop')
  write_tab(env_audit,'environment_filter_agreement')
  if(any(env_audit$n_values>1)) stop('Environmental values differ between filters of one original sample. Resolve metadata; do not average silently.')
  env_bio <- env %>% group_by(BioSample) %>% summarise(across(all_of(names(aliases)),
    ~if(all(is.na(.x))) NA_real_ else first(na.omit(.x))),.groups='drop')
  bio <- bio %>% left_join(env_bio,by='BioSample')
  write_tab(bio,'S1_biological_environment')

  # MEGARes reference-level counts; preserve alleles until length normalization.
  raw <- read_csv(locate('dedup_AMR_analytic_matrix2.csv'),show_col_types=FALSE)
  names(raw) <- normalise_ids(names(raw),'AMR matrix columns',strict=FALSE)
  if(anyDuplicated(names(raw))) stop('Duplicate AMR columns after sample-ID normalisation.')
  need(raw,c('gene_accession',lib$SampleID),'AMR matrix')
  unique_key(raw,'gene_accession','AMR matrix')
  anno <- raw %>% select(gene_accession) %>% mutate(parts=str_split(gene_accession,fixed('|')))
  if(any(lengths(anno$parts)<5)) stop('Expected MEGAResID|Type|Class|Mechanism|GeneId[|Note].')
  anno <- anno %>% transmute(gene_accession,Type=map_chr(parts,~.x[2]),
    Class=map_chr(parts,~.x[3]),Mechanism=map_chr(parts,~.x[4]),GeneId=map_chr(parts,~.x[5]))
  long <- raw %>% select(gene_accession,all_of(lib$SampleID)) %>%
    pivot_longer(-gene_accession,names_to='SampleID',values_to='Counts')
  if(any(!is.finite(long$Counts) | long$Counts<0)) stop('Invalid AMR counts.')
  long <- long %>% left_join(lib,by='SampleID') %>% left_join(anno,by='gene_accession') %>%
    mutate(RPM=1e6*Counts/TotalReads)
  pooled <- long %>% group_by(BioSample,gene_accession,Type,Class,Mechanism,GeneId) %>%
    summarise(nPositiveFilters=sum(Counts>0),Counts=sum(Counts),.groups='drop') %>%
    left_join(bio %>% select(BioSample,Site,Date,TotalReads),by='BioSample') %>%
    mutate(RPM=1e6*Counts/TotalReads)
  lf <- locate('megares_lengths.tsv',FALSE)
  have_lengths <- !is.null(lf)
  if(have_lengths) {
    le <- read_t(lf); need(le,c('gene_accession','length'),'lengths'); unique_key(le,'gene_accession','lengths')
    pooled <- pooled %>% left_join(le %>% select(gene_accession,length),by='gene_accession')
    if(any(!is.finite(pooled$length)|pooled$length<=0)) stop('Missing/invalid MEGARes reference lengths.')
    pooled <- pooled %>% mutate(RPKM=RPM/(length/1000))
    record('length_normalization','DONE','Reference-level RPKM before gene/class aggregation; not biomass correction.')
  } else {
    write_tsv(anno %>% select(gene_accession) %>% mutate(length=NA_real_),file.path(RUN,'templates','megares_lengths.tsv'))
    record('length_normalization','MISSING INPUT','Provide megares_lengths.tsv; RPM results remain provisional for between-gene/class rankings.')
  }
  write_tab(pooled,'S5_MEGARes_reference_biological')
  genes <- pooled %>% group_by(BioSample,Site,Date,GeneId) %>%
    summarise(RPM=sum(RPM),.groups='drop')
  write_tab(genes,'S5_MEGARes_gene_biological')
  metrics <- c('RPM',if(have_lengths) 'RPKM')
  make_matrix <- function(metric) {
    x <- pooled %>% group_by(BioSample,GeneId) %>% summarise(Value=sum(.data[[metric]]),.groups='drop') %>%
      pivot_wider(names_from=GeneId,values_from=Value,values_fill=0) %>% column_to_rownames('BioSample') %>% as.matrix()
    x <- x[bio$BioSample,,drop=FALSE]; x[,colSums(x)>0,drop=FALSE]
  }
  # Permutations exchange site labels WITHIN a date, at biological sample level.
  # Date-only weather predictors cannot be tested under this scheme.
  # Tests describe these locations, not replicated habitat types.
  pmatrix <- function(md) shuffleSet(nrow(md),nset=NPERM,
    control=how(blocks=factor(md$Date),within=Within(type='free')))
  alpha_all <- list()
  for(metric in metrics) {
    X <- make_matrix(metric)
    if(any(rowSums(X)==0)) stop('Zero-total resistome sample: inspect before Hellinger analyses.')
    H <- decostand(X,'hellinger')
    alpha <- bio %>% select(BioSample,Site,Date,TotalReads) %>% mutate(
      Total=rowSums(X),Richness=specnumber(X),q1=exp(diversity(X,'shannon')),
      q2=diversity(X,'invsimpson'),Metric=metric)
    alpha_all[[metric]] <- alpha
    write_tab(alpha,paste0('alpha_',metric))
    p <- alpha %>% pivot_longer(c(Total,Richness,q1,q2),names_to='Outcome') %>%
      ggplot(aes(Site,value,color=Site))+geom_boxplot(outlier.shape=NA)+
      geom_point(position=position_jitter(width=.08,height=0),size=2)+
      facet_wrap(~Outcome,scales='free_y')+scale_color_manual(values=cols)+
      labs(y=metric,title='One point per original water sample')
    save_plot(p,paste0('alpha_',metric))
    # Balanced site tests: three sites on six dates; four sites on two common dates.
    common_dates <- bio %>% count(Date) %>% filter(n==4) %>% pull(Date)
    subsets <- list(ABD_six_dates=which(bio$Site %in% c('A','B','D')),
                    ABCD_shared_dates=which(bio$Date %in% common_dates))
    for(tag in names(subsets)) {
      ix <- subsets[[tag]]; md <- as.data.frame(bio[ix,]); md$Site <- factor(md$Site); md$Date <- factor(md$Date)
      Y <- H[ix,,drop=FALSE]; P <- pmatrix(md)
      d <- vegdist(Y,'bray')
      fit <- adonis2(d~Date+Site,data=md,permutations=P,by='margin')
      # Export SITE inference only: date is blocked and cannot be tested by P.
      write_tab(as.data.frame(fit)['Site',,drop=FALSE] %>% rownames_to_column('Term'),paste0('PERMANOVA_',metric,'_',tag))
      saveRDS(list(model=fit,permutations=P,metadata=md),file.path(RUN,'audit',paste0('PERMANOVA_',metric,'_',tag,'.rds')))
      # Descriptive dispersion; blocked permutation diagnostic (few n in winter).
      bd <- betadisper(d,md$Site,bias.adjust=TRUE)
      write_tab(tibble(BioSample=md$BioSample,Site=md$Site,Distance=bd$distances),
                paste0('dispersion_distances_',metric,'_',tag))
      if(any(table(md$Site)<3)) {
        record(paste0('dispersion_',metric,'_',tag),'DESCRIPTIVE ONLY',
               'Two biological samples per site: within-site distances are degenerate; no dispersion F-test reported.')
      } else {
        bt <- permutest(bd,permutations=P)
        write_tab(as.data.frame(bt$tab) %>% rownames_to_column('Term'),paste0('dispersion_',metric,'_',tag))
      }
      ord <- cmdscale(d,k=2,eig=TRUE,add=TRUE)
      pts <- cbind(md,as.data.frame(ord$points)); names(pts)[(ncol(pts)-1):ncol(pts)] <- c('Axis1','Axis2')
      save_plot(ggplot(pts,aes(Axis1,Axis2,color=Site,shape=Date))+geom_point(size=3)+
        scale_color_manual(values=cols)+labs(title=paste(metric,tag),subtitle='Bray-Curtis on Hellinger profiles; additive PCoA correction'),
        paste0('PCoA_',metric,'_',tag))
      # Friedman and EXACT paired sign-flip tests: no assumption of 3 independent filters.
      af <- alpha[ix,]; res <- list()
      for(out in c('Total','Richness','q1','q2')) {
        aa <- af %>% transmute(Site=factor(Site),Date=factor(Date),value=.data[[out]])
        ft <- tryCatch(friedman.test(value~Site|Date,data=aa),error=function(e) NULL)
        write_tab(tibble(Outcome=out,ChiSquared=if(is.null(ft)) NA_real_ else unname(ft$statistic),
          p_asymptotic=if(is.null(ft)) NA_real_ else ft$p.value,
          Note='Small number of dates: asymptotic Friedman p is descriptive; use exact paired contrasts.'),
          paste0('friedman_',metric,'_',tag,'_',out))
        for(pair in combn(sort(unique(af$Site)),2,simplify=FALSE)) {
          wide <- af %>% filter(Site %in% pair) %>% select(Date,Site,all_of(out)) %>%
            pivot_wider(names_from=Site,values_from=all_of(out))
          dif <- wide[[pair[1]]]-wide[[pair[2]]]
          signs <- as.matrix(expand.grid(rep(list(c(-1,1)),length(dif))))
          null <- as.vector(signs %*% dif / length(dif))
          pval <- mean(abs(null)>=abs(mean(dif))-1e-12)
          res[[length(res)+1]] <- tibble(Outcome=out,Site1=pair[1],Site2=pair[2],n_dates=length(dif),
            MeanDifference=mean(dif),MedianDifference=median(dif),p_exact=pval)
        }
      }
      write_tab(bind_rows(res) %>% group_by(Outcome) %>% mutate(q_BH=p.adjust(p_exact,'BH')) %>% ungroup(),
                paste0('paired_contrasts_',metric,'_',tag))
    }
    # Descriptive gene heatmap and classes, with all zero observations preserved.
    top <- names(sort(colSums(X),decreasing=TRUE))[seq_len(min(60,ncol(X)))]
    heat <- as.data.frame(X[,top,drop=FALSE]) %>% rownames_to_column('BioSample') %>%
      pivot_longer(-BioSample,names_to='Gene',values_to='Value') %>% left_join(bio %>% select(BioSample,Site),by='BioSample')
    save_plot(ggplot(heat,aes(BioSample,Gene,fill=log10(Value+1)))+geom_tile()+
      facet_grid(~Site,scales='free_x',space='free_x')+scale_fill_viridis_c()+
      theme(axis.text.x=element_text(angle=90,hjust=1))+labs(title=paste('Top 60 genes:',metric)),paste0('gene_heatmap_',metric),14,11)
    cl <- pooled %>% group_by(BioSample,Site,Date,Type,Class) %>% summarise(Value=sum(.data[[metric]]),.groups='drop')
    write_tab(cl,paste0('class_abundance_',metric))
  }
  alpha <- alpha_all[['RPM']]
  # Sequencing-depth QC is descriptive at n=20; it does not establish equal sensitivity.
  qc <- bind_rows(lapply(c('Total','Richness'),function(v) {
    z <- cor.test(alpha$TotalReads,alpha[[v]],method='spearman',exact=FALSE)
    tibble(Outcome=v,rho=unname(z$estimate),p=z$p.value,n=nrow(alpha))
  }))
  write_tab(qc,'sequencing_depth_QC')
  # Recurrence: 1 biological sample/date. No comparisons of 2-date and 6-date cores.
  for(tag in c('ABD_six_dates','ABCD_shared_dates')) {
    eligible <- if(tag=='ABD_six_dates') bio %>% filter(Site!='C') else bio %>% filter(Date %in% common_dates)
    recur <- genes %>% filter(BioSample %in% eligible$BioSample) %>% group_by(Site,GeneId) %>%
      summarise(n_detected_dates=sum(RPM>0),n_dates=n(),fraction=n_detected_dates/n_dates,
        recurrent_all_dates=all(RPM>0),.groups='drop')
    write_tab(recur,paste0('recurrence_',tag))
  }

  # Environmental analyses: prespecified primary pH and water temperature,
  # site + date conditioned. Each tested separately; no significance-driven selection.
  # Meteorological predictors date-only: descriptive audit, no pseudo-replicated tests.
  H <- decostand(make_matrix('RPM'),'hellinger')
  audit <- bind_rows(lapply(names(aliases),function(v) {
    x <- bio[[v]]; k <- split(x,bio$Date)
    tibble(Variable=v,Missing=sum(is.na(x)),WithinDateVariable=any(vapply(k,function(y)n_distinct(y,na.rm=TRUE)>1,logical(1))))
  }))
  write_tab(audit,'environment_identifiability')
  env_results <- list()
  for(v in c('WaterTemperature','pH','Salinity')) {
    md <- as.data.frame(bio); md$Site <- factor(md$Site); md$Date <- factor(md$Date)
    x <- md[[v]]
    if(anyNA(x)||sd(x)==0) {record(paste0('pRDA_',v),'SKIPPED','Missing or constant predictor.');next}
    C <- model.matrix(~Site+Date,md)
    residual <- qr.resid(qr(C),x)
    if(sum(residual^2)<=1e-10*sum((x-mean(x))^2)) {
      record(paste0('pRDA_',v),'NON-ESTIMABLE','Predictor lies in site/date covariate space.');next
    }
    md$Predictor <- as.numeric(scale(x))
    fit <- rda(H~Predictor+Condition(Site+Date),data=md)
    P <- pmatrix(md); test <- anova(fit,permutations=P)
    env_results[[v]] <- tibble(Variable=v,Analysis=if(v=='Salinity') 'exploratory' else 'primary',
      residual_fraction=sum(residual^2)/sum((x-mean(x))^2),
      adjusted_R2=RsquareAdj(fit)$adj.r.squared,F=test$F[1],p=test$`Pr(>F)`[1])
    saveRDS(list(model=fit,test=test,metadata=md,permutations=P),file.path(RUN,'audit',paste0('pRDA_',v,'.rds')))
  }
  # The final supplement defines one BH family across all three predictors.
  if(length(env_results)) write_tab(bind_rows(env_results) %>%
    mutate(q_BH=p.adjust(p,'BH')),'partial_RDA')
  record('environment','DONE','No date-only meteorological hypothesis tests; six dates insufficient for the previous multivariable weather model. No legacy S9 singular regressions retained.')

  # CARPDM reference summaries. Original nomenclature map copied below from template.
  files <- list.files(file.path(INPUT,'CARD_genes'),pattern='\\.tsv$',full.names=TRUE)
  files <- files[!grepl('taxonomy',basename(files),ignore.case=TRUE)]
  if(!length(files)) {
    record('CARPDM','MISSING INPUT','Provide CARD_genes per-library mapping summaries. No absence inference.')
  } else {
    file_ids <- normalise_ids(str_split_fixed(basename(files),fixed('.'),2)[,1],'CARPDM filenames')
    if(anyDuplicated(file_ids)) stop('Multiple CARPDM summary files for the same library; remove duplicate input versions.')
    if(!setequal(file_ids,lib$SampleID)) {
      print(list(Missing=setdiff(lib$SampleID,file_ids),Unexpected=setdiff(file_ids,lib$SampleID)))
      stop('Require a CARPDM summary file for each library; header-only files are allowed for zero strict hits.')
    }
    car <- map_dfr(seq_along(files),function(i) {
      f <- files[i]
      # Fixed import types prevent per-file guessing (especially empty files).
      x <- readr::read_tsv(f,col_types=readr::cols(.default=readr::col_character()),
                           show_col_types=FALSE,progress=FALSE)
      required <- c('gene','length','cov_bases','breadth','mean_depth','num_reads')
      need(x,required,f)
      x <- x %>% select(all_of(required))
      # Remove only exact repeated header records from concatenated summaries.
      header <- rep(TRUE,nrow(x))
      for(nm in required) header <- header & !is.na(x[[nm]]) & trimws(x[[nm]])==nm
      if(any(header)) {
        record(paste0('CARPDM_header_',basename(f)),'CLEANED',
               paste('Removed',sum(header),'exact repeated header rows.'))
        x <- x[!header,,drop=FALSE]
      }
      for(nm in setdiff(required,'gene')) {
        original <- x[[nm]]
        parsed <- suppressWarnings(readr::parse_double(original,na=c('','NA')))
        bad <- !is.finite(parsed)
        if(any(bad)) {
          detail <- tibble(File=f,DataRow=which(bad),Column=nm,
                           Gene=x$gene[bad],Value=original[bad])
          write_tab(detail,'ERROR_CARPDM_numeric_values')
          print(as.data.frame(detail),row.names=FALSE)
          stop('Invalid CARPDM numeric values; see ERROR_CARPDM_numeric_values.tsv. No rows silently discarded.',call.=FALSE)
        }
        x[[nm]] <- parsed
      }
      x$SampleID <- rep(file_ids[i],nrow(x)); x
    })
    unique_key(car,c('SampleID','gene'),'CARPDM')
    if(!nrow(car)) stop('All CARPDM summaries are empty; no gene catalogue available to analyse.')
    if(any(!is.finite(car$breadth)|car$breadth<0|car$breadth>1)) stop('CARPDM breadth must be fractions 0-1.')
    for(v in c('length','cov_bases','mean_depth','num_reads')) if(any(!is.finite(car[[v]])|car[[v]]<0)) stop('Invalid CARPDM ',v)
    if(any(car$length==0)) stop('Zero reference length.')
    ck <- car %>% group_by(gene) %>% summarise(n=n_distinct(length),.groups='drop')
    if(any(ck$n!=1)) stop('CARPDM reference lengths inconsistent.')
    if(any(abs(car$breadth-car$cov_bases/car$length)>0.011)) stop('CARPDM breadth inconsistent with covered bases / reference length.')
    car <- car %>% mutate(GeneId=map_chr(str_split(gene,fixed('|')),~if(length(.x)>=6) .x[6] else .x[1]),
      Class=unname(megares_class[GeneId]),Class=replace_na(Class,'Unmapped; inspect')) %>%
      left_join(lib,by='SampleID') %>% mutate(RPM=1e6*num_reads/TotalReads,RPKM=RPM/(length/1000),
        Pass_strict=breadth>=CARPDM_BREADTH & mean_depth>=CARPDM_DEPTH,
        Pass_breadth=breadth>=CARPDM_BREADTH)
    write_tab(car,'S3_CARPDM_filter_all_available_hits')
    screens <- if(isTRUE(CARPDM_UNFILTERED)) c('strict','breadth') else 'strict'
    if(!isTRUE(CARPDM_UNFILTERED)) record('CARPDM_depth_sensitivity','BLOCKED',
      'Set CARPDM_UNFILTERED=TRUE only after checking input provenance; discarded mappings cannot be recovered by R.')
    if(isTRUE(CARPDM_UNFILTERED)) {
      # Require explicit all-reference rows: no silent zero imputation for missing calls.
      expected <- n_distinct(car$gene)*nrow(lib)
      if(nrow(car)!=expected) stop('Unfiltered CARPDM summaries must contain every catalogue reference in every library (including zero rows).')
    }
    all_cb <- list()
    for(screen in screens) {
      x <- car %>% mutate(Pass=if(screen=='strict') Pass_strict else Pass_breadth) %>%
        group_by(BioSample,gene,GeneId,Class,length) %>%
        summarise(PassingReads=sum(if_else(Pass,num_reads,0)),nPassFilters=sum(Pass),
          nRows=n(),.groups='drop')
      # Complete missing BIOLOGICAL combinations only for strict filtered tables;
      # zeros mean no passing call, not proven biological absence.
      refs <- car %>% distinct(gene,GeneId,Class,length)
      x <- crossing(BioSample=bio$BioSample,gene=refs$gene) %>%
        left_join(x,by=c('BioSample','gene')) %>% select(-GeneId,-Class,-length) %>%
        left_join(refs,by='gene') %>% mutate(across(c(PassingReads,nPassFilters,nRows),~replace_na(.x,0))) %>%
        left_join(bio %>% select(BioSample,Site,Date,TotalReads),by='BioSample') %>%
        mutate(RPM=1e6*PassingReads/TotalReads,RPKM=RPM/(length/1000),
          Positive_any=nPassFilters>=1,Positive_two=nPassFilters>=2,Positive_all=nPassFilters==3,
          Screen=screen)
      all_cb[[screen]] <- x
    }
    cb <- bind_rows(all_cb); write_tab(cb,'S3_CARPDM_biological_thresholds')
    write_tab(cb %>% group_by(Screen,Site) %>% summarise(n_biosamples=n_distinct(BioSample),
      n_positive_samples=n_distinct(BioSample[Positive_any]),n_refs=n_distinct(gene[Positive_any]),
      .groups='drop'),'CARPDM_screen_summary')
    totals <- cb %>% group_by(Screen,BioSample,Site,Date) %>% summarise(RPM=sum(RPM),RPKM=sum(RPKM),.groups='drop')
    write_tab(totals,'CARPDM_sample_totals')
    save_plot(ggplot(totals,aes(Date,RPKM,color=Site,group=Site))+geom_line()+geom_point(size=2)+
      facet_wrap(~Screen)+scale_color_manual(values=cols)+scale_y_continuous(trans='log1p')+
      labs(y='Length-normalised CARPDM abundance (RPKM; log1p axis)',subtitle='Technical filters combined; no independent filter error bars'),
      'CARPDM_threshold_timeseries',12,6)
    # Complete date prevalence and abundance table, inclusive of zeros.
    write_tab(cb %>% group_by(Screen,Site,GeneId) %>% summarise(n_dates=n_distinct(Date),
      n_positive_dates=n_distinct(Date[Positive_any]),Mean_RPKM=sum(RPKM)/n_distinct(Date),.groups='drop'),
      'CARPDM_site_gene_summary')
    record('CARPDM','DONE','Filter-level thresholds; read-weighted biological abundance. >=1, >=2, all-three support flags exported.')
  }
  # Taxonomy: preserve all ranks and genera; no pathogen-selected host inference.
  taxfiles <- list.files(INPUT,pattern='\\.CARD_ARG_taxonomy.tsv$',recursive=TRUE,full.names=TRUE)
  if(length(taxfiles) && exists('car')) {
    tx <- map_dfr(taxfiles,function(f) {
      # taxid and taxonomic ranks are identifiers, not numeric measurements.
      # Fix every column to character before binding files; only counts are parsed.
      z <- readr::read_tsv(f,col_types=readr::cols(.default=readr::col_character()),
                           show_col_types=FALSE,progress=FALSE)
      need(z,c('ARG_ID','read_count','genus'),f)
      header <- rep(TRUE,nrow(z))
      for(nm in names(z)) header <- header & !is.na(z[[nm]]) & trimws(z[[nm]])==nm
      if(any(header)) {
        record(paste0('taxonomy_header_',basename(f)),'CLEANED',
               paste('Removed',sum(header),'exact repeated header rows.'))
        z <- z[!header,,drop=FALSE]
      }
      original <- z$read_count
      parsed <- suppressWarnings(readr::parse_double(original,na=c('','NA')))
      bad <- !is.finite(parsed) | parsed<0 | abs(parsed-round(parsed))>1e-8
      if(any(bad)) {
        detail <- tibble(File=f,DataRow=which(bad),ARG_ID=z$ARG_ID[bad],Value=original[bad])
        write_tab(detail,'ERROR_taxonomy_read_counts')
        print(as.data.frame(detail),row.names=FALSE)
        stop('Invalid taxonomy read_count values; see ERROR_taxonomy_read_counts.tsv.',call.=FALSE)
      }
      z$read_count <- parsed
      sid <- normalise_ids(str_split_fixed(basename(f),fixed('.'),2)[,1],
                           'CARPDM taxonomy filenames')
      if(!sid %in% lib$SampleID) stop('Unknown taxonomy library: ',sid,' in ',f)
      z$SampleID <- rep(sid,nrow(z))
      z
    })
    # Must have every library represented; zero-call libraries may lack taxonomic files.
    tx <- tx %>% inner_join(car %>% filter(Pass_strict) %>% select(SampleID,gene),by=c('SampleID','ARG_ID'='gene')) %>%
      left_join(lib %>% select(SampleID,BioSample,Site,Date),by='SampleID')
    if(anyNA(tx$BioSample)) stop('Unknown taxonomy samples.')
    tx <- tx %>% mutate(genus=if_else(is.na(genus)|genus %in% c('','Unassigned','unclassified'),'Unclassified',genus))
    write_tab(tx,'S4_CARPDM_taxonomy_all_rows')
    tg <- tx %>% group_by(BioSample,Site,Date,genus) %>% summarise(Reads=sum(read_count),.groups='drop')
    write_tab(tg,'S4_all_genera_biological')
    record('taxonomy','DESCRIPTIVE','All genera exported. Classification fractions require verified non-overlapping read assignments; no automatic host inference.')
  } else record('taxonomy','MISSING INPUT','No CARPDM taxonomy summaries available.')

  # Explicit crosswalk of output tables; no silent renumbering of prior-studies S2.
  writeLines(c('Main analyses use n=20 original water samples, not n=60 filters.',
    'A/B/D: six dates; C: two dates. Only these locations are replicated over time.',
    'Counts pooled across technical filters. CARPDM calls remain filter-level.',
    'Site tests exchange biological units within date; date effects are not tested by this permutation.',
    'Paired exact contrasts have limited power: 2 dates cannot support conventional significance.',
    'Weather/date-only effects are descriptive. pRDA tests only residual within-date information.',
    'RPKM is not copies per cell or per litre. MAPQ/breadth enforcement upstream must be verified.',
    'No conclusions about viability, pathogen carriage, HGT events or dilution follow from this script.',
    'No assemblies, faecal markers, coordinates, filter specifications or tidal timestamps were supplied.',
    'Retain the previous-studies S2 table separately; join output tables into the final supplement after review.',
    'Do not interpret missing/blocked modules as negative findings.',
    'Check module_status.tsv before editing the manuscript.'),file.path(RUN,'ANALYSIS_NOTES.txt'))
  record('assembly_MGE_faecal_markers','UPSTREAM REQUIRED','R cannot infer missing assemblies, marker mappings or physical linkage from RPM tables.')
  # Workbook contains the concise result tables; full long tables stay in TSV.
  tabfiles <- list.files(file.path(RUN,'tables'),pattern='\\.tsv$',full.names=TRUE)
  tabfiles <- tabfiles[!grepl('reference_biological|filter_all_available|taxonomy_all_rows|all_genera',basename(tabfiles))]
  wb <- openxlsx::createWorkbook()
  index <- tibble(Sheet=sprintf('T%02d',seq_along(tabfiles)),File=basename(tabfiles))
  openxlsx::addWorksheet(wb,'Index');openxlsx::writeData(wb,'Index',index)
  for(i in seq_along(tabfiles)) {
    z <- read_t(tabfiles[i]); sh <- index$Sheet[i]
    openxlsx::addWorksheet(wb,sh);openxlsx::writeData(wb,sh,z)
    openxlsx::freezePane(wb,sh,firstRow=TRUE)
  }
  openxlsx::saveWorkbook(wb,file.path(RUN,'Revision_analysis_tables.xlsx'),overwrite=TRUE)
  record('pipeline','COMPLETE','Review audit and blocked modules; output is not automatically publication-final.')
}

tryCatch(main(), error=function(e) {
  record("pipeline","FAILED",conditionMessage(e))
  stop(e)
}, finally={
  capture.output(sessionInfo(),file=file.path(RUN,"audit","sessionInfo.txt"))
  writeLines(c(paste("Seed",SEED),paste("Requested permutations",NPERM),
    paste("CARPDM_UNFILTERED",CARPDM_UNFILTERED)),file.path(RUN,"audit","configuration.txt"))
})
message("Results: ",RUN)
