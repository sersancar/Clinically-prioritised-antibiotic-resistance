#!/usr/bin/env python3
"""Link unique query/reference associations and rebuild community reports.

Only ARG-associated query identifiers are kept in memory. SQLite deduplicates
mates without a BAM-sized Python set; the full Kraken output is scanned once.
"""
import argparse
import json
import re
import sqlite3
from collections import Counter
from pathlib import Path
import pysam
from common import Taxonomy,query_id,read_tsv,write_tsv


def taxid_from_label(label):
    if label.isdigit():return int(label)
    m=re.search(r'\(taxid\s+(\d+)\)',label)
    if not m:raise ValueError(f'Unrecognised Kraken taxid field: {label}')
    return int(m.group(1))


def support_score(taxid,field):
    assigned=total=0
    for token in field.split():
        if ':' not in token:continue
        key,value=token.split(':',1)
        if key in ['A','0']:continue
        try:n=int(value)
        except ValueError:continue
        if n<0:raise ValueError('Negative k-mer count')
        total+=n
        if key==str(taxid):assigned+=n
    return assigned/total if total else None


def kraken_report(direct,taxonomy,path):
    clade=Counter();children={};active=set()
    for taxid,n in direct.items():
        for ancestor in taxonomy.ancestors(taxid):clade[ancestor]+=n;active.add(ancestor)
    for taxid in active:
        parent=taxonomy.parent[taxid]
        if parent!=taxid and parent in active:children.setdefault(parent,[]).append(taxid)
    codes={'superkingdom':'D','domain':'D','kingdom':'K','phylum':'P','class':'C','order':'O','family':'F','genus':'G','species':'S'}
    with open(path,'w') as f:
        def emit(taxid,depth):
            code='R' if taxid==1 else codes.get(taxonomy.rank[taxid],'-')
            f.write(f'{100*clade[taxid]/sum(direct.values()):.6f}\t{clade[taxid]}\t{direct.get(taxid,0)}\t{code}\t{taxid}\t'+
                    '  '*depth+taxonomy.names.get(taxid,str(taxid))+'\n')
            for child in sorted(children.get(taxid,[]),key=lambda t:(-clade[t],t)):emit(child,depth+1)
        if direct:emit(1,0)


def link(filter_id,broad_bam,carpdm_bam,broad_counts,carpdm_metrics,kraken,taxonomy_dir,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    eligible={'broad':{r['reference'] for r in read_tsv(broad_counts) if int(r['count'])>0},
              'carpdm':{r['reference'] for r in read_tsv(carpdm_metrics) if r['primary']=='True'}}
    database=sqlite3.connect(out/'associations.sqlite')
    database.execute('PRAGMA journal_mode=OFF');database.execute('PRAGMA synchronous=OFF')
    database.execute('CREATE TABLE associations(catalogue TEXT,reference TEXT,query TEXT,PRIMARY KEY(catalogue,reference,query)) WITHOUT ROWID')
    for catalogue,bam in [('broad',broad_bam),('carpdm',carpdm_bam)]:
        with pysam.AlignmentFile(bam,'rb') as f:
            batch=[]
            for read in f.fetch(until_eof=True):
                if read.flag&2308 or read.mapping_quality<30 or read.reference_name not in eligible[catalogue]:continue
                batch.append((catalogue,read.reference_name,query_id(read.query_name)))
                if len(batch)>=10000:database.executemany('INSERT OR IGNORE INTO associations VALUES(?,?,?)',batch);batch=[]
            database.executemany('INSERT OR IGNORE INTO associations VALUES(?,?,?)',batch)
    database.commit();wanted={r[0] for r in database.execute('SELECT DISTINCT query FROM associations')}
    assignments={};direct=Counter();retained=0;classified=0;total_pairs=0
    with open(kraken) as f:
        for line in f:
            fields=line.rstrip('\n').split('\t')
            if len(fields)<5:raise ValueError('Malformed Kraken output')
            total_pairs+=1
            if fields[0]!='C':continue
            classified+=1;taxid=taxid_from_label(fields[2]);score=support_score(taxid,fields[4])
            if taxid<=0 or (score is not None and score<.2):continue
            retained+=1;direct[taxid]+=1;query=query_id(fields[1])
            if query in wanted:
                if query in assignments and assignments[query]!=taxid:raise ValueError('Conflicting paired-query taxonomy')
                assignments[query]=taxid
    taxonomy=Taxonomy(taxonomy_dir)
    kraken_report(direct,taxonomy,out/'community.report')
    database.execute('CREATE TABLE taxonomy(query TEXT PRIMARY KEY,taxid INTEGER) WITHOUT ROWID')
    database.executemany('INSERT INTO taxonomy VALUES(?,?)',assignments.items());database.commit()
    summary=Counter();annotated=[]
    for catalogue,reference,taxid,count in database.execute('SELECT a.catalogue,a.reference,t.taxid,COUNT(*) FROM associations a LEFT JOIN taxonomy t USING(query) GROUP BY a.catalogue,a.reference,t.taxid'):
        if taxid is None:category='no_retained_taxid'
        elif taxid==1:category='root_only'
        else:
            info=taxonomy.annotate(taxid);category='genus_resolved' if info['genus'] else 'above_genus'
            annotated.append(dict(filter_id=filter_id,catalogue=catalogue,reference=reference,
                                  association_count=count,**info))
        summary[(catalogue,category)]+=count
    fields=['filter_id','catalogue','reference','association_count','taxid','name','rank','domain','lineage',
            'phylum','class','order','family','genus','species','taxon','taxon_id','taxon_rank']
    write_tsv(out/'taxonomy.tsv',annotated,fields)
    completeness=[]
    for catalogue in ['broad','carpdm']:
        row=dict(filter_id=filter_id,catalogue=catalogue)
        for category in ['genus_resolved','above_genus','root_only','no_retained_taxid']:row[category]=summary[(catalogue,category)]
        row['eligible_associations']=sum(row[c] for c in ['genus_resolved','above_genus','root_only','no_retained_taxid'])
        completeness.append(row)
    write_tsv(out/'completeness.tsv',completeness)
    (out/'community.json').write_text(json.dumps(dict(filter_id=filter_id,total_pairs=total_pairs,
        classified_pairs=classified,retained_pairs=retained,custom_score_threshold=.2),indent=2)+'\n')
    database.close();(out/'associations.sqlite').unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['filter-id','broad-bam','carpdm-bam','broad-counts','carpdm-metrics','kraken','taxonomy','out']:p.add_argument('--'+key,required=True)
    a=p.parse_args();link(a.filter_id,a.broad_bam,a.carpdm_bam,a.broad_counts,a.carpdm_metrics,a.kraken,a.taxonomy,a.out)
