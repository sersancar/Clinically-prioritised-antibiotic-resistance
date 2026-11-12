#!/usr/bin/env python3
"""Assemble filter measurements and pooled biological tables for the report."""
import argparse
import itertools
import json
import shutil
from collections import defaultdict
from pathlib import Path
from common import read_tsv,reference_annotation,write_tsv
from ena_manifest import validate_filters


def assemble(filters_path,inputs,metadata,outdir):
    filters=read_tsv(filters_path);groups=validate_filters(filters,require_runs=False)
    byid={r['filter_id']:r for r in filters};inputs=Path(inputs);metadata=Path(metadata)
    out=Path(outdir);tables=out/'tables';tables.mkdir(parents=True,exist_ok=True)
    qc={}
    for path in inputs.rglob('*.fastp.json'):
        fid=path.name.removesuffix('.fastp.json')
        if fid in qc:raise ValueError(f'Duplicate QC file {fid}')
        qc[fid]=int(json.loads(path.read_text())['summary']['after_filtering']['total_reads'])
    if set(qc)!=set(byid) or min(qc.values())<=0:raise ValueError('Require positive runtime totals for all filters')
    totals={bio:sum(qc[r['filter_id']] for r in rows) for bio,rows in groups.items()}
    meta=[{k:v for k,v in r.items() if k!='expected_post_qc_read_ends'}|{'total_read_ends':qc[r['filter_id']]} for r in filters]
    write_tsv(tables/'TableS1_filters.tsv',meta)
    bio_meta=[]
    for bio,rows in sorted(groups.items()):
        row={k:v for k,v in rows[0].items() if k not in ['filter_id','replicate','run_accession','expected_post_qc_read_ends']}
        row.update(total_read_ends=totals[bio],filters=len(rows));bio_meta.append(row)
    write_tsv(tables/'TableS1_biological.tsv',bio_meta)
    shutil.copyfile(metadata/'historical_context.tsv',tables/'TableS2_historical_context.tsv')
    classes={r['gene']:r['class'] for r in read_tsv(metadata/'carpdm_gene_classes.tsv')}
    carp=[];broad=[];marker=[];seen=set()
    for kind,pattern,target in [('carpdm','*.carpdm.metrics.tsv',carp),('broad','*.broad.counts.tsv',broad),('crass','*.crass.metrics.tsv',marker)]:
        refs=defaultdict(set)
        for path in inputs.rglob(pattern):
            for r in read_tsv(path):
                fid=r['filter_id'];ref=r['reference'];key=(kind,fid,ref)
                if fid not in byid or key in seen:raise ValueError('Unknown or duplicated filter/reference measurement')
                seen.add(key);refs[fid].add(ref);m=byid[fid]
                row=dict(filter_id=fid,biosample=m['biosample'],site=m['site'],date=m['date'],reference=ref,
                    length=int(r['length']),breadth=float(r['breadth']),mean_depth=float(r['mean_depth']),
                    total_read_ends=qc[fid])
                if kind=='broad':row.update(count=int(r['count']),raw_count=int(r['raw_count']),requires_snp=r['requires_snp'])
                else:row.update(count=int(r['mapped_read_ends']),primary=r['primary']=='True',relaxed=r['relaxed']=='True')
                if kind!='crass':
                    annotation=reference_annotation(ref,kind,classes)
                    row.update(gene=annotation['gene'],amr_type=annotation['type'],amr_class=annotation['class'],mechanism=annotation['mechanism'])
                target.append(row)
        if set(refs)!=set(byid) or len({frozenset(v) for v in refs.values()})!=1:
            raise ValueError(f'{kind}: all references, including zeros, must be measured in every filter')
    write_tsv(tables/'TableS4_CARPDM_filter_screening.tsv',carp)
    write_tsv(tables/'Broad_filter_counts.tsv',broad)
    write_tsv(tables/'TableS12_marker_filters.tsv',marker)

    def pool(rows,screens):
        sums=defaultdict(int);support=defaultdict(int);annotations={}
        for r in rows:
            annotations[r['reference']]={k:r[k] for k in ['reference','length','gene','amr_type','amr_class','mechanism']}
            for screen in screens:
                keep=r['count']>0 if screen=='broad' else r[screen]
                key=(r['biosample'],r['reference'],screen)
                sums[key]+=r['count'] if keep else 0;support[key]+=keep
        result=[]
        for bio,ref,screen in itertools.product(sorted(groups),sorted(annotations),screens):
            count=sums[(bio,ref,screen)];n=support[(bio,ref,screen)];m=groups[bio][0]
            result.append(dict(biosample=bio,site=m['site'],date=m['date'],screen=screen,
               **annotations[ref],count=count,positive_filters=n,total_read_ends=totals[bio],
               RPM=1e6*count/totals[bio],RPKM=1e9*count/totals[bio]/annotations[ref]['length']))
        return result
    cp=pool(carp,['primary','relaxed']);bp=pool(broad,['broad'])
    write_tsv(tables/'TableS5_CARPDM_biological.tsv',cp)
    write_tsv(tables/'TableS8_broad_biological.tsv',bp)
    sensitivity=[]
    for site in ['A','B','C','D','All']:
        rows=[r for r in carp if site=='All' or r['site']==site]
        for screen in ['primary','relaxed']:
            hit=[r for r in rows if r[screen]]
            sensitivity.append(dict(site=site,screen=screen,filters=len({r['filter_id'] for r in rows}),
                biological_samples=len({r['biosample'] for r in rows}),unique_targets=len({r['reference'] for r in hit}),
                filter_target_calls=len(hit),biological_sample_target_calls=len({(r['biosample'],r['reference']) for r in hit})))
    write_tsv(tables/'TableS6_threshold_sensitivity.tsv',sensitivity)
    tax=[];completeness=[]
    for path in inputs.rglob('taxonomy.tsv'):
        for r in read_tsv(path):
            fid=r['filter_id'];m=byid[fid];annotation=reference_annotation(r['reference'],r['catalogue'],classes)
            row=r|dict(biosample=m['biosample'],site=m['site'],date=m['date'],association_count=int(r['association_count']),
                       gene=annotation['gene'],amr_type=annotation['type'],amr_class=annotation['class'],
                       filter_RPM=1e6*int(r['association_count'])/qc[fid])
            tax.append(row)
    tax_fields=['filter_id','biosample','site','date','catalogue','reference','gene','amr_type','amr_class','association_count',
                'filter_RPM','taxid','name','rank','domain','lineage','phylum','class','order','family','genus','species','taxon','taxon_id','taxon_rank']
    for catalogue,name in [('carpdm','TableS7_CARPDM_taxonomy.tsv'),('broad','TableS9_broad_taxonomy.tsv')]:
        write_tsv(tables/name,[r for r in tax if r['catalogue']==catalogue],tax_fields)
    for path in inputs.rglob('completeness.tsv'):completeness.extend(read_tsv(path))
    expected={(fid,c) for fid in byid for c in ['carpdm','broad']}
    keys=[(r['filter_id'],r['catalogue']) for r in completeness]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:raise ValueError('Incomplete/duplicate taxonomy linkage results')
    comp=[];categories=['genus_resolved','above_genus','root_only','no_retained_taxid','eligible_associations']
    for catalogue in ['carpdm','broad']:
        selected=[r for r in completeness if r['catalogue']==catalogue]
        for level in ['filter','site','study']:
            partitions=defaultdict(list)
            for r in selected:partitions[r['filter_id'] if level=='filter' else byid[r['filter_id']]['site'] if level=='site' else 'All'].append(r)
            for key,items in sorted(partitions.items()):
                row=dict(catalogue=catalogue,level=level,identifier=key)
                row.update({c:sum(int(r[c]) for r in items) for c in categories})
                n=row['eligible_associations'];annotated=row['genus_resolved']+row['above_genus']
                row.update(genus_percent=100*row['genus_resolved']/n if n else '',
                    annotated_percent=100*annotated/n if n else '',
                    genus_percent_of_annotated=100*row['genus_resolved']/annotated if annotated else '')
                comp.append(row)
    write_tsv(tables/'TableS7_completeness.tsv',comp)
    marker_bio=[]
    for bio,group in groups.items():
        rows=[r for r in marker if r['biosample']==bio];count=sum(r['count'] for r in rows)
        marker_bio.append(dict(biosample=bio,site=group[0]['site'],date=group[0]['date'],count=count,
            total_read_ends=totals[bio],RPM=1e6*count/totals[bio],positive_filters=sum(r['primary'] for r in rows)))
    write_tsv(tables/'TableS12_marker_biological.tsv',marker_bio)
    bracken=[];community=[]
    for path in inputs.rglob('*.bracken.tsv'):
        fid=path.name.removesuffix('.bracken.tsv');m=byid[fid]
        for r in read_tsv(path):bracken.append(dict(filter_id=fid,biosample=m['biosample'],site=m['site'],**r))
    if bracken:
        write_tsv(tables/'Community_species_filter.tsv',bracken)
        species=defaultdict(float);names={}
        for r in bracken:
            taxid=r['taxonomy_id'];key=(r['biosample'],taxid)
            species[key]+=float(r['new_est_reads']);names[taxid]=r['name']
        species_totals=defaultdict(float)
        for (bio,taxid),n in species.items():species_totals[bio]+=n
        rows=[]
        for (bio,taxid),n in sorted(species.items()):
            rows.append(dict(biosample=bio,site=groups[bio][0]['site'],taxid=taxid,name=names[taxid],
                estimated_read_pairs=n,total_post_QC_read_pairs=totals[bio]/2,
                RPM_read_pairs=2e6*n/totals[bio],fraction_of_estimated_species=n/species_totals[bio] if species_totals[bio] else 0))
        write_tsv(tables/'Community_species_biological.tsv',rows)
    for path in inputs.rglob('community.json'):community.append(json.loads(path.read_text()))
    write_tsv(tables/'Community_classification_QC.tsv',community)
    versions=[]
    for path in inputs.rglob('*.versions.json'):
        record=json.loads(path.read_text())
        for tool,version in record.items():versions.append(dict(component=tool,version=version,source=path.name))
    write_tsv(tables/'TableS3_software.tsv',versions,['component','version','source'])
    (out/'provenance').mkdir(exist_ok=True)
    statuses=[json.loads(p.read_text()) for p in inputs.rglob('*.snp_status.json')]
    if statuses:write_tsv(out/'provenance/SNP_verification_status.tsv',statuses)
    for path in inputs.rglob('resource_provenance.json'):
        shutil.copyfile(path,out/'provenance/reference_resources.json')
        record=json.loads(path.read_text())
        write_tsv(tables/'TableS3_reference_checksums.tsv',
            [dict(resource=k,sha256=v) for k,v in record.get('sha256',{}).items()],['resource','sha256'])
    shutil.copyfile(metadata/'candidate_genera.txt',out/'candidate_genera.txt')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['filters','inputs','metadata','out']:p.add_argument('--'+key,required=True)
    a=p.parse_args();assemble(a.filters,a.inputs,a.metadata,a.out)
