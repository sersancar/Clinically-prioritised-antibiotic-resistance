#!/usr/bin/env Rscript
# Method-based analysis and batch figures. No historical output values are used.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2)
 library(vegan);library(permute);library(patchwork)})
a<-commandArgs(TRUE);if(length(a)<2)stop('Usage: report.R DATA_DIRECTORY OUTPUT_DIRECTORY [SEED] [NPERM]')
input<-a[1];out<-a[2];seed<-if(length(a)>2)as.integer(a[3])else 20261112L
nperm<-if(length(a)>3)as.integer(a[4])else 9999L
set.seed(seed);dir.create(out,recursive=TRUE,showWarnings=FALSE)
file.copy(file.path(input,'tables'),out,recursive=TRUE)
dir.create(file.path(out,'figures'),showWarnings=FALSE);dir.create(file.path(out,'plot_data'),showWarnings=FALSE)
tab<-function(name){
 z<-read_tsv(file.path(out,'tables',name),col_types=cols(.default=col_character()),show_col_types=FALSE,progress=FALSE)
 numeric_fields<-c('total_read_ends','filters','salinity','pH','water_temperature','air_temperature','wind','uv','humidity','rainfall','tide',
  'length','breadth','mean_depth','count','raw_count','positive_filters','RPM','RPKM','association_count','filter_RPM')
 for(k in intersect(numeric_fields,names(z)))z[[k]]<-parse_double(z[[k]])
 for(k in intersect(c('primary','relaxed'),names(z)))z[[k]]<-tolower(z[[k]])=='true'
 z
}
save_tab<-function(x,name)write_tsv(as_tibble(x),file.path(out,'tables',paste0(name,'.tsv')))
plot_data<-function(x,name)write_tsv(as_tibble(x),file.path(out,'plot_data',paste0(name,'.tsv')))
save_plot<-function(p,name,w=12,h=8){
 ggsave(file.path(out,'figures',paste0(name,'.pdf')),p,width=w,height=h,limitsize=FALSE)
 ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=180,bg='white',limitsize=FALSE)
}
empty_plot<-function(title,text='No qualifying observations')ggplot()+annotate('text',x=0,y=0,label=text)+labs(title=title)+theme_void()
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
colours<-c(A='#00A087',B='#E64B35',C='#3C5488',D='#4DBBD5')
md<-tab('TableS1_biological.tsv')%>%mutate(date=as.Date(date),site=factor(site,levels=names(colours)))%>%arrange(site,date)
cp<-tab('TableS5_CARPDM_biological.tsv')%>%mutate(date=as.Date(date))
car<-tab('TableS4_CARPDM_filter_screening.tsv')%>%mutate(date=as.Date(date))
bp<-tab('TableS8_broad_biological.tsv')
ctax<-tab('TableS7_CARPDM_taxonomy.tsv');btax<-tab('TableS9_broad_taxonomy.tsv')
if(anyDuplicated(md$biosample)||nrow(md)!=20||any(md$filters!=3))stop('Expected 20 biological samples with three filters each')
date_levels<-sort(unique(md$date));shape_values<-setNames(c(16,17,15,3,7,8),as.character(date_levels))
md$date_label<-factor(as.character(md$date),levels=names(shape_values))
date_scale<-function()scale_shape_manual(values=shape_values,name='Sampling date')
site_scale<-function()scale_colour_manual(values=colours,drop=FALSE,name='Site')
primary<-cp%>%filter(screen=='primary')
gene_car<-primary%>%group_by(biosample,site,date,gene,amr_class)%>%summarise(RPM=sum(RPM),.groups='drop')
mean_car<-gene_car%>%group_by(site,gene,amr_class)%>%summarise(mean_RPM=mean(RPM),n=n(),
 lower=pmax(0,mean(RPM)-qt(.975,n()-1)*sd(RPM)/sqrt(n())),
 upper=mean(RPM)+qt(.975,n()-1)*sd(RPM)/sqrt(n()),.groups='drop')%>%filter(mean_RPM>0)
save_tab(mean_car,'CARPDM_site_mean_CI');plot_data(mean_car,'Figure1B')
p1<-if(nrow(mean_car))ggplot(mean_car,aes(x=reorder(gene,mean_RPM),y=mean_RPM,colour=site))+
 geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge(width=.6))+
 facet_wrap(~amr_class,scales='free_y')+coord_flip()+scale_y_sqrt()+site_scale()+
 labs(x=NULL,y='Mean RPM (square-root scale)',title='CARPDM-prioritised resistance signals',
 caption='95% t-based intervals describe variation across sampled dates; affiliations do not establish hosts.') else empty_plot('CARPDM-prioritised resistance signals')
save_plot(p1,'Figure1B_CARPDM',14,10)

# Gene matrices preserve exact references until length normalisation.
gene_broad<-bp%>%group_by(biosample,gene)%>%summarise(RPM=sum(RPM),RPKM=sum(RPKM),.groups='drop')
matrix_for<-function(metric){
 z<-gene_broad%>%select(biosample,gene,all_of(metric))%>%pivot_wider(names_from=gene,values_from=all_of(metric),values_fill=0)
 x<-as.matrix(z[,-1]);rownames(x)<-z$biosample;x<-x[md$biosample,,drop=FALSE]
 x[,colSums(x)>0,drop=FALSE]
}
mats<-list(RPM=matrix_for('RPM'),RPKM=matrix_for('RPKM'))
if(!ncol(mats$RPM)||any(rowSums(mats$RPM)<=0))stop('Broad resistome has an empty biological profile; ordination is not estimable')
cls<-bp%>%group_by(biosample,amr_type,amr_class)%>%summarise(RPM=sum(RPM),.groups='drop')%>%left_join(md%>%select(biosample,site),by='biosample')
class_levels<-cls%>%distinct(amr_type,amr_class)%>%arrange(amr_type,amr_class)%>%pull(amr_class)%>%unique()
p2a<-ggplot(cls,aes(factor(biosample,levels=md$biosample),factor(amr_class,levels=rev(class_levels)),fill=log10(RPM+1)))+
 geom_tile()+facet_grid(amr_type~site,scales='free',space='free')+
 scale_fill_viridis_c(name='log10(RPM + 1)')+labs(x=NULL,y=NULL,title='A  Broad resistance classes')+
 theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6),axis.text.y=element_text(size=7))
plot_data(cls,'Figure2A')

# Date-blocked permutations; raw distances for tests, Cailliez only for display.
perms<-function(d)shuffleSet(nrow(d),nset=nperm,control=how(blocks=factor(d$date)))
ordination<-function(metric,ids,title){
 d<-md[match(ids,md$biosample),,drop=FALSE];x<-mats[[metric]][ids,,drop=FALSE]
 h<-decostand(x,'hellinger');distance<-vegdist(h,'bray')
 fit<-wcmdscale(distance,k=2,eig=TRUE,add='cailliez')
 points<-as.data.frame(fit$points);names(points)<-c('Axis1','Axis2');points<-bind_cols(d,points)
 positive<-fit$eig[fit$eig>0];pct<-100*positive[1:2]/sum(positive)
 ell<-points%>%group_by(site)%>%filter(n()>=4)%>%ungroup()
 if(nrow(ell)){
  valid<-ell%>%group_split(site)%>%Filter(function(z)det(cov(z[c('Axis1','Axis2')]))>1e-12,.)
  ell<-bind_rows(valid)
 }
 p<-ggplot(points,aes(Axis1,Axis2,colour=site,shape=date_label))+geom_point(size=2.6)+site_scale()+date_scale()+
  labs(title=title,x=sprintf('PCoA1 (%.1f%%)',pct[1]),y=sprintf('PCoA2 (%.1f%%)',pct[2]))
 if(nrow(ell))p<-p+stat_ellipse(data=ell,aes(Axis1,Axis2,colour=site,group=site),inherit.aes=FALSE,type='norm',level=.95,show.legend=FALSE)
 list(plot=p,points=points,distance=distance,metadata=d,hellinger=h)
}
all_ord<-ordination('RPM',md$biosample,'B  Broad resistome composition')
save_plot(p2a+all_ord$plot+plot_layout(widths=c(1.5,1)),'Figure2',18,10)
plot_data(all_ord$points,'Figure2B')

test_rows<-list();disp_rows<-list();fried_rows<-list();paired_rows<-list();alpha_rows<-list()
for(metric in names(mats)){
 x<-mats[[metric]];p<-x/rowSums(x)
 alpha<-md%>%select(biosample,site,date,total_read_ends)%>%mutate(metric=metric,
  summed_abundance=rowSums(x),richness=rowSums(x>0),
  exponential_shannon=exp(-rowSums(ifelse(p>0,p*log(p),0))),inverse_simpson=1/rowSums(p^2))
 alpha_rows[[metric]]<-alpha
}
alpha<-bind_rows(alpha_rows);save_tab(alpha,'Alpha_diversity')
common_dates<-Reduce(intersect,lapply(split(md$date,md$site),as.character))
subsets<-list(ABD_six_dates=md$biosample[md$site%in%c('A','B','D')],
 ABCD_shared_winter=md$biosample[as.character(md$date)%in%common_dates],All=md$biosample)
for(metric in names(mats))for(subset in names(subsets)){
 ids<-subsets[[subset]];d<-as.data.frame(md[match(ids,md$biosample),]);d$site<-droplevels(factor(d$site));d$date<-factor(d$date)
 distance<-vegdist(decostand(mats[[metric]][ids,,drop=FALSE],'hellinger'),'bray');P<-perms(d)
 fit<-adonis2(distance~date+site,data=d,permutations=P,by='margin')
 test_rows[[paste(metric,subset)]]<-tibble(metric,subset,term='site',R2=fit['site','R2'],F=fit['site','F'],p=fit['site','Pr(>F)'],permutations=nrow(P))
 bd<-betadisper(distance,d$site,bias.adjust=TRUE)
 if(all(table(d$site)==2)){
  disp_rows[[paste(metric,subset)]]<-tibble(metric,subset,F=NA_real_,p=NA_real_,interpretation='Not tested: two biological samples per site')
 }else{
  bt<-permutest(bd,permutations=P)
  disp_rows[[paste(metric,subset)]]<-tibble(metric,subset,F=bt$tab[1,'F'],p=bt$tab[1,'Pr(>F)'],
   interpretation=if(any(table(d$site)<3))'Small site group: interpret cautiously'else'Date-restricted permutations')
 }
 saveRDS(list(permutations=P,metadata=d),file.path(out,'plot_data',paste0('permutations_',metric,'_',subset,'.rds')))
 if(subset=='All')next
 z<-alpha%>%filter(.data$metric==.env$metric,biosample%in%ids)
 z$site<-droplevels(z$site)
 for(outcome in c('summed_abundance','richness','exponential_shannon','inverse_simpson')){
  ft<-friedman.test(as.formula(paste(outcome,'~site|date')),data=z)
  fried_rows[[paste(metric,subset,outcome)]]<-tibble(metric,subset,outcome,chi_squared=unname(ft$statistic),df=unname(ft$parameter),p=ft$p.value)
  for(pair in combn(sort(unique(as.character(z$site))),2,simplify=FALSE)){
   zz<-z%>%filter(site%in%pair)%>%select(date,site,all_of(outcome))%>%pivot_wider(names_from=site,values_from=all_of(outcome))%>%arrange(date)
   differences<-zz[[pair[1]]]-zz[[pair[2]]]
   signs<-as.matrix(expand.grid(rep(list(c(-1,1)),length(differences))))
   null<-abs(as.numeric(signs%*%differences)/length(differences));observed<-abs(mean(differences))
   paired_rows[[paste(metric,subset,outcome,paste(pair,collapse='_'))]]<-tibble(metric,subset,outcome,
    site1=pair[1],site2=pair[2],mean_difference=mean(differences),p_exact=mean(null>=observed-1e-12),matched_dates=length(differences))
  }
 }
}
save_tab(bind_rows(test_rows),'TableS10_PERMANOVA');save_tab(bind_rows(disp_rows),'TableS10_dispersion')
save_tab(bind_rows(fried_rows),'TableS10_Friedman')
save_tab(bind_rows(paired_rows)%>%group_by(metric,subset,outcome)%>%mutate(p_BH=p.adjust(p_exact,'BH'))%>%ungroup(),'TableS10_paired_contrasts')
correlations<-bind_rows(lapply(c('summed_abundance','richness'),function(v){
 z<-alpha%>%filter(metric=='RPM');test<-cor.test(z$total_read_ends,z[[v]],method='spearman',exact=FALSE)
 tibble(outcome=v,rho=unname(test$estimate),p=test$p.value,interpretation='Descriptive sequencing-depth correlation')}))
save_tab(correlations,'TableS10_sequencing_depth')

# Separate residual environmental tests; no date-only weather hypothesis tests.
environment<-list();env_plots<-list();H<-decostand(mats$RPM,'hellinger')
for(v in c('water_temperature','pH','salinity')){
 d<-as.data.frame(md);d$site<-factor(d$site);d$date<-factor(d$date);value<-as.numeric(d[[v]])
 C<-model.matrix(~site+date,d);residual<-qr.resid(qr(C),value)
 if(anyNA(value)||sd(value)==0||sum(residual^2)<1e-10*sum((value-mean(value))^2)){
  environment[[v]]<-tibble(variable=v,F=NA_real_,p=NA_real_,adjusted_R2=NA_real_,residual_variance_fraction=NA_real_,status='Not estimable')
  env_plots[[v]]<-empty_plot(v,'Predictor not estimable after site/date conditioning');next
 }
 d$predictor<-as.numeric(scale(value));fit<-rda(H~predictor+Condition(site+date),data=d)
 P<-perms(d);test<-anova(fit,permutations=P);s<-as.data.frame(scores(fit,display='sites',choices=1:2,scaling=1));names(s)<-c('Constrained','Residual')
 s<-bind_cols(md,s);environment[[v]]<-tibble(variable=v,F=test$F[1],p=test$`Pr(>F)`[1],adjusted_R2=RsquareAdj(fit)$adj.r.squared,
  residual_variance_fraction=sum(residual^2)/sum((value-mean(value))^2),status='Tested')
 env_plots[[v]]<-ggplot(s,aes(Constrained,Residual,colour=site,shape=date_label))+geom_point(size=2.5)+site_scale()+date_scale()+
  labs(title=v,subtitle=sprintf('pseudo-F %.3f; p %.4f',test$F[1],test$`Pr(>F)`[1]),x='Constrained axis',y='First residual axis')
 plot_data(s,paste0('FigureS8_',v));saveRDS(list(model=fit,permutations=P),file.path(out,'plot_data',paste0('pRDA_',v,'.rds')))
}
save_tab(bind_rows(environment)%>%mutate(p_BH=p.adjust(p,'BH')),'TableS11_partial_RDA')
audit<-bind_rows(lapply(c('salinity','pH','water_temperature','air_temperature','wind','uv','humidity','rainfall','tide'),function(v){
 tibble(variable=v,missing=sum(is.na(md[[v]])),within_date_variable=any(vapply(split(md[[v]],md$date),function(x)n_distinct(x)>1,logical(1))))}))
save_tab(audit,'TableS11_covariate_audit');save_tab(md%>%select(biosample,site,date,tide),'TableS11_tidal_context')
save_plot(wrap_plots(env_plots,ncol=3,guides='collect'),'FigureS8_partial_RDA',16,6)

# S1: no interpolation across the unsampled spring/winter interval.
temporal<-gene_car%>%group_by(biosample,site,date)%>%summarise(total_RPM=sum(RPM),.groups='drop')%>%
 mutate(campaign=if_else(as.integer(format(date,'%m'))<7,'Spring','Winter'))
topgenes<-gene_car%>%group_by(gene)%>%summarise(total=sum(RPM),.groups='drop')%>%arrange(desc(total))%>%filter(total>0)%>%slice_head(n=20)%>%pull(gene)
t<-gene_car%>%filter(gene%in%topgenes);positive<-t$RPM[t$RPM>0];epsilon<-if(length(positive))min(positive)/2 else 1
pa<-ggplot(temporal,aes(date,log10(total_RPM+1),colour=site,group=interaction(site,campaign)))+geom_line()+geom_point()+site_scale()+labs(y='log10(summed RPM + 1)',x=NULL)
pb<-if(nrow(t))ggplot(t,aes(factor(date),gene,fill=ifelse(RPM>0,log10(RPM+epsilon),NA_real_)))+geom_tile()+facet_grid(.~site)+
 scale_fill_viridis_c(na.value='white',name='log10(RPM + epsilon)')+theme(axis.text.x=element_text(angle=90,hjust=1))+labs(x=NULL,y=NULL)else empty_plot('Temporal CARPDM targets')
plot_data(t,'FigureS1');save_tab(tibble(pseudocount=epsilon),'Temporal_heatmap_pseudocount')
save_plot(pa/pb+plot_layout(heights=c(1,2)),'FigureS1_temporal',14,12)

# S2: each filter contributes at most once to a gene, regardless of alleles.
support<-car%>%group_by(filter_id,biosample,site,date,gene)%>%summarise(positive=any(primary),.groups='drop')%>%
 group_by(biosample,site,date,gene)%>%summarise(positive_filters=sum(positive),.groups='drop')
showgenes<-support%>%group_by(gene)%>%summarise(n=sum(positive_filters),.groups='drop')%>%filter(n>0)%>%pull(gene)
sp<-support%>%filter(gene%in%showgenes);plot_data(sp,'FigureS2')
p<-if(nrow(sp))ggplot(sp,aes(factor(date),gene,fill=factor(positive_filters,levels=0:3)))+geom_tile(colour='grey90')+facet_grid(.~site)+
 scale_fill_manual(values=c('white','#9ecae1','#4292c6','#08519c'),drop=FALSE,name='Supporting filters')+
 theme(axis.text.x=element_text(angle=90,hjust=1))+labs(x=NULL,y=NULL)else empty_plot('CARPDM filter support')
save_plot(p,'FigureS2_filter_support',14,max(6,min(18,length(showgenes)*.18)))

# S3: selected genera, with zero-inclusive site means of pooled associations.
candidate<-readLines(file.path(input,'candidate_genera.txt'));genus_plot<-ctax%>%filter(!is.na(genus),genus!='',genus%in%candidate)%>%
 group_by(biosample,site,gene,genus,amr_class)%>%summarise(count=sum(association_count),.groups='drop')%>%
 left_join(md%>%select(biosample,total_read_ends),by='biosample')%>%mutate(RPM=1e6*count/total_read_ends)%>%
 group_by(site,gene,genus,amr_class)%>%summarise(sum_RPM=sum(RPM),.groups='drop')%>%
 left_join(md%>%count(site,name='n_dates'),by='site')%>%mutate(mean_RPM=sum_RPM/n_dates)
plot_data(genus_plot,'FigureS3')
p<-if(nrow(genus_plot))ggplot(genus_plot,aes(gene,genus,size=mean_RPM,colour=amr_class,shape=site))+geom_point()+
 scale_size_area(max_size=9)+labs(x=NULL,y=NULL,size='Mean association RPM',colour='Resistance class',title='Read-level affiliations; genomic hosts are not established')+
 theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))else empty_plot('Genus-level CARPDM affiliations')
save_plot(p,'FigureS3_genus_affiliations',15,max(7,min(18,n_distinct(genus_plot$genus)*.3)))

# S4: full gene heatmap plus shared-date overlap. All underlying values exported.
order_ids<-rownames(mats$RPM)[hclust(vegdist(mats$RPM,'bray'),method='complete')$order]
genes_meta<-bp%>%distinct(gene,amr_type,amr_class)%>%arrange(amr_type,amr_class,gene)
# A gene label can occur in multiple reference classes: retain each class separately here.
heat<-bp%>%group_by(biosample,gene,amr_type,amr_class)%>%summarise(RPM=sum(RPM),.groups='drop')%>%
 mutate(feature=paste(amr_type,amr_class,gene,sep=' | '))
features<-heat%>%distinct(feature,amr_type,amr_class,gene)%>%arrange(amr_type,amr_class,gene)%>%pull(feature)
p<-ggplot(heat,aes(factor(biosample,levels=order_ids),factor(feature,levels=rev(features)),fill=log10(RPM+1)))+
 geom_tile()+facet_grid(amr_type~.,scales='free_y',space='free_y')+scale_fill_viridis_c()+labs(x=NULL,y='Resistance genes (ordered by type and class)')+
 theme(axis.text.x=element_text(angle=90,hjust=1,size=6),axis.text.y=if(length(features)>100)element_blank()else element_text(size=5),axis.ticks.y=element_blank())
plot_data(heat,'FigureS4_heatmap');save_plot(p,'FigureS4_gene_heatmap',13,14)
sets<-lapply(names(colours),function(site_i){ids<-md$biosample[md$site==site_i&as.character(md$date)%in%common_dates];colnames(mats$RPM)[colSums(mats$RPM[ids,,drop=FALSE])>0]});names(sets)<-names(colours)
overlap<-tibble(gene=sort(unique(unlist(sets))));for(s in names(sets))overlap[[s]]<-overlap$gene%in%sets[[s]]
save_tab(overlap,'TableS8_shared_date_gene_overlap')
pdf(file.path(out,'figures','FigureS4_gene_overlap.pdf'),9,9)
grid::grid.newpage();if(length(unique(unlist(sets))))grid::grid.draw(VennDiagram::venn.diagram(sets,filename=NULL,fill=unname(colours),alpha=.5))
dev.off()

# S5: depth diagnostics and all three Hill numbers.
z<-alpha%>%filter(metric=='RPM');pa<-ggplot(z,aes(total_read_ends,summed_abundance+1,colour=site))+geom_point()+scale_x_log10()+scale_y_log10()+site_scale()+labs(x='Pooled post-QC read ends',y='Summed RPM + 1')
pb<-ggplot(z,aes(total_read_ends,richness,colour=site))+geom_point()+scale_x_log10()+site_scale()+labs(x='Pooled post-QC read ends',y='Gene richness')
za<-z%>%pivot_longer(c(richness,exponential_shannon,inverse_simpson),names_to='diversity',values_to='value')
pc<-ggplot(za,aes(site,value,colour=site))+geom_boxplot(data=filter(za,site!='C'),outlier.shape=NA)+geom_point(position=position_jitter(width=.08,seed=seed))+facet_wrap(~diversity,scales='free_y')+site_scale()+labs(x=NULL,y='Hill diversity')
save_plot((pa+pb)/pc,'FigureS5_depth_diversity',13,10)

# S6: persistent above-median nodes; edge persistence is not required.
edge_bio<-btax%>%group_by(biosample,site,amr_type,amr_class,taxon,taxon_id,taxon_rank,domain)%>%summarise(count=sum(association_count),.groups='drop')%>%
 left_join(md%>%select(biosample,total_read_ends),by='biosample')%>%mutate(RPM=1e6*count/total_read_ends)
edges<-edge_bio%>%group_by(site,amr_type,amr_class,taxon,taxon_id,taxon_rank,domain)%>%summarise(sum_RPM=sum(RPM),.groups='drop')%>%
 left_join(md%>%count(site,name='n_dates'),by='site')%>%mutate(mean_RPM=sum_RPM/n_dates)
class_bio<-edge_bio%>%group_by(biosample,site,amr_type,amr_class)%>%summarise(RPM=sum(RPM),.groups='drop')%>%
 mutate(id=paste0('class:',amr_type,':',amr_class),label=amr_class,kind='Resistance class',domain=NA_character_)
tax_bio<-edge_bio%>%group_by(biosample,site,taxon,taxon_id,domain)%>%summarise(RPM=sum(RPM),.groups='drop')%>%
 mutate(id=paste0('taxon:',taxon_id),label=taxon,kind='Taxonomic affiliation',amr_type=NA_character_)
nodes<-bind_rows(class_bio,tax_bio)%>%group_by(site,id,label,kind,domain,amr_type)%>%summarise(sum_RPM=sum(RPM),positive_dates=sum(RPM>0),.groups='drop')%>%
 left_join(md%>%count(site,name='n_dates'),by='site')%>%mutate(mean_RPM=sum_RPM/n_dates,persistent=positive_dates==n_dates)
eligible<-nodes%>%filter(persistent,kind=='Resistance class'|domain%in%c('Bacteria','Archaea'))%>%
 group_by(site,kind)%>%mutate(selected=mean_RPM>median(mean_RPM))%>%ungroup()
keep<-eligible%>%filter(selected)%>%group_by(site,kind)%>%arrange(desc(mean_RPM),label,.by_group=TRUE)%>%
 mutate(x=if_else(kind=='Resistance class',1,0),y=if(n()==1).5 else seq(0,1,length.out=n()))%>%ungroup()
edges<-edges%>%mutate(source=paste0('taxon:',taxon_id),target=paste0('class:',amr_type,':',amr_class))
plotted<-edges%>%inner_join(keep%>%select(site,source=id,x,y),by=c('site','source'))%>%
 inner_join(keep%>%select(site,target=id,xend=x,yend=y),by=c('site','target'))
save_tab(nodes,'Network_node_statistics');save_tab(eligible,'Network_node_selection');save_tab(keep,'Network_nodes');save_tab(plotted,'Network_edges')
plot_data(edge_bio,'FigureS6_biological_associations')
p<-if(nrow(keep))ggplot()+geom_segment(data=plotted,aes(x,y,xend=xend,yend=yend,linewidth=mean_RPM),colour='grey65',alpha=.6)+
 geom_point(data=keep,aes(x,y,size=mean_RPM,colour=if_else(kind=='Resistance class',amr_type,domain)))+
 geom_text(data=keep,aes(x=if_else(x==0,-.04,1.04),y=y,label=gsub('_',' ',label),hjust=if_else(x==0,1,0)),size=2.2)+
 facet_wrap(~site,ncol=2)+scale_size_area(max_size=13,name='Mean association RPM')+scale_linewidth(range=c(.2,3),name='Mean edge RPM')+
 scale_x_continuous(limits=c(-1.3,2.3))+labs(colour='Domain / resistance type')+theme_void()+theme(legend.position='bottom') else empty_plot('Persistent resistance–taxonomic affiliation networks')
save_plot(p,'FigureS6_networks',17,14)

plots<-list();for(metric in c('RPM','RPKM'))for(subset in c('ABD_six_dates','ABCD_shared_winter')){
 o<-ordination(metric,subsets[[subset]],paste(metric,subset));plots[[paste(metric,subset)]]<-o$plot;plot_data(o$points,paste0('FigureS7_',metric,'_',subset))
}
save_plot(wrap_plots(plots,ncol=2,guides='collect'),'FigureS7_sensitivity',14,11)
marker<-tab('TableS12_marker_biological.tsv')
save_tab(marker%>%group_by(site)%>%summarise(biological_samples=n(),mean_RPM=mean(RPM),median_RPM=median(RPM),
 passing_filters=sum(positive_filters),samples_any_qualifying=sum(positive_filters>0),.groups='drop'),'TableS12_marker_sites')
writeLines(c('All statistics use 20 biological site–date samples, not 60 independent filters.',
 'No date-only weather hypothesis tests are performed. Environmental associations are not climate-change effects.',
 'Read-level taxonomic affiliations do not establish genomic hosts, viability or mobile-element context.',
 'Marker qualification is operational. No ARG-to-marker dilution comparison is inferred from low coverage.',
 'Figure 1A is intentionally omitted. Figure 1B and Figures 2/S1–S8 are generated.',
 paste('Seed:',seed,'; requested date-restricted permutations:',nperm)),file.path(out,'INTERPRETATION.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
packages<-c('vegan','permute','ggplot2','patchwork','readr','dplyr','tidyr','VennDiagram')
save_tab(tibble(component=c('R',packages),version=c(as.character(getRversion()),vapply(packages,function(p)as.character(packageVersion(p)),character(1)))),'TableS3_R_packages')
