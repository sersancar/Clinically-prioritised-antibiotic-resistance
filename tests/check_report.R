a<-commandArgs(TRUE);root<-a[1]
read<-function(name)read.delim(file.path(root,'tables',paste0(name,'.tsv')),check.names=FALSE)
stopifnot(nrow(read('TableS10_PERMANOVA'))==6,nrow(read('TableS10_Friedman'))==16)
pair<-read('TableS10_paired_contrasts');stopifnot(nrow(pair)==72,all(pair$p_BH>=pair$p_exact-1e-12))
env<-read('TableS11_partial_RDA');stopifnot(nrow(env)==3,isTRUE(all.equal(env$p_BH,p.adjust(env$p,'BH'))))
files<-list.files(file.path(root,'plot_data'),pattern='^permutations_.*rds$',full.names=TRUE)
stopifnot(length(files)==6)
for(f in files){z<-readRDS(f);P<-z$permutations;dates<-as.character(z$metadata$date)
 stopifnot(all(vapply(seq_len(nrow(P)),function(i)all(dates[P[i,]]==dates),logical(1))))}
stopifnot(!file.exists(file.path(root,'figures','Figure1A.pdf')))
message('Report methods and date-blocked permutation invariants passed.')
