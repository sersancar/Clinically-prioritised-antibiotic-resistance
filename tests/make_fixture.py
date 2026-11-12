#!/usr/bin/env python3
"""Synthetic 60-filter fixture for method and report validation; not study data."""
import json
import math
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from common import read_tsv,write_tsv


def make(root,out,empty_carpdm=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);filters=read_tsv(root/'metadata/filters.tsv')
    refs=[f'MEG_{i}|Drugs|'+('MLS' if i%2 else 'Tetracyclines')+f'|Mechanism|gene{i}' for i in range(1,13)]
    carp=['gb|example|+|1-1000|ARO:1|sul2','gb|example|+|1-1000|ARO:2|tet(O)']
    tax_fields=['filter_id','catalogue','reference','association_count','taxid','name','rank','domain','lineage','phylum','class','order','family','genus','species','taxon','taxon_id','taxon_rank']
    dates=sorted({r['date'] for r in filters})
    for f in filters:
        fid=f['filter_id'];site=f['site'];rep=int(f['replicate']);date=dates.index(f['date']);si='ABCD'.index(site)
        (out/(fid+'.fastp.json')).write_text(json.dumps({'summary':{'after_filtering':{'total_reads':1000000+si*1000+date*1100+rep*100}}}))
        broad=[];car=[];marker=[];taxonomy=[]
        for i,ref in enumerate(refs):
            count=round(20+10*si+date*3+rep+6*math.sin((i+1)*(date+1)+si))
            if i>4+si and (date+i)%3!=0:count=0
            broad.append(dict(filter_id=fid,reference=ref,length=1000,breadth=1 if count else 0,mean_depth=8 if count else 0,raw_count=count,count=count,requires_snp=False))
            if count:
                genus='Escherichia' if i%2 else 'Aeromonas';name='Gammaproteobacteria'
                taxonomy.append(dict(filter_id=fid,catalogue='broad',reference=ref,association_count=max(1,count//2),
                    taxid=11,name=genus,rank='genus',domain='Bacteria',lineage='Bacteria;'+name+';'+genus,
                    phylum='Pseudomonadota',**{'class':name},order='',family='',genus=genus,species='',taxon=name,taxon_id=10,taxon_rank='class'))
        for i,ref in enumerate(carp):
            positive=not empty_carpdm and site in ['B','D'] and not(site=='D' and i==1 and date%2==0)
            shallow=not empty_carpdm and site=='A' and date%2==0 and rep==1
            breadth=1 if positive or shallow else 0;depth=6 if positive else 2 if shallow else 0;count=40+date+rep if positive else 4 if shallow else 0
            car.append(dict(filter_id=fid,catalogue='carpdm',reference=ref,length=1000,covered_bases=int(breadth*1000),breadth=breadth,mean_depth=depth,mapped_read_ends=count,primary=positive,relaxed=bool(breadth)))
            if positive:
                taxonomy.append(dict(filter_id=fid,catalogue='carpdm',reference=ref,association_count=max(1,count//2),taxid=11,name='Escherichia',rank='genus',domain='Bacteria',lineage='Bacteria;Escherichia',
                    phylum='Pseudomonadota',**{'class':'Gammaproteobacteria'},order='',family='',genus='Escherichia',species='',taxon='Gammaproteobacteria',taxon_id=10,taxon_rank='class'))
        marker.append(dict(filter_id=fid,catalogue='crass',reference='NC_024711.1',length=1000,covered_bases=1000 if site=='B' else 100,breadth=1 if site=='B' else .1,mean_depth=8 if site=='B' else .2,mapped_read_ends=100 if site=='B' else 2,primary=site=='B',relaxed=site=='B'))
        write_tsv(out/(fid+'.broad.counts.tsv'),broad);write_tsv(out/(fid+'.carpdm.metrics.tsv'),car);write_tsv(out/(fid+'.crass.metrics.tsv'),marker)
        links=out/(fid+'.links');links.mkdir();write_tsv(links/'taxonomy.tsv',taxonomy,tax_fields)
        comp=[]
        for cat in ['broad','carpdm']:
            n=sum(r['association_count'] for r in taxonomy if r['catalogue']==cat)
            comp.append(dict(filter_id=fid,catalogue=cat,genus_resolved=n,above_genus=0,root_only=0,no_retained_taxid=0,eligible_associations=n))
        write_tsv(links/'completeness.tsv',comp)
        (links/'community.json').write_text(json.dumps(dict(filter_id=fid,total_pairs=1000,classified_pairs=800,retained_pairs=700,custom_score_threshold=.2)))
    (out/'fixture.versions.json').write_text(json.dumps({'test_fixture':'Synthetic data only'}))


if __name__=='__main__':make(Path(__file__).resolve().parents[1],sys.argv[1],empty_carpdm='--empty-carpdm' in sys.argv)
