"""Shared identifiers, tabular IO and frozen NCBI taxonomy handling."""
import csv
import hashlib
from pathlib import Path


def read_tsv(path):
    with open(path,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def write_tsv(path,rows,fields=None):
    rows=list(rows);path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        if not rows:raise ValueError(f'No schema supplied for empty table {path}')
        fields=list(rows[0])
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        while data:=f.read(8*1024*1024):h.update(data)
    return h.hexdigest()


def fasta_lengths(path):
    result={};name=None
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in result:raise ValueError(f'Duplicate FASTA identifier {name}')
                result[name]=0
            elif line.strip():
                if name is None:raise ValueError('Invalid FASTA')
                result[name]+=len(line.strip())
    if not result or min(result.values())<=0:raise ValueError('Empty FASTA sequence')
    return result


def query_id(value):
    value=value.split()[0]
    return value[:-2] if value.endswith(('/1','/2')) else value


def reference_annotation(reference,catalogue,classes=None):
    fields=reference.split('|')
    if catalogue=='broad':
        if len(fields)<5:raise ValueError(f'Expected MEGARes annotation in {reference}')
        return {'gene':fields[4],'type':fields[1],'class':fields[2],'mechanism':fields[3]}
    gene=fields[-1]
    return {'gene':gene,'type':'','class':(classes or {}).get(gene,'Unassigned'),'mechanism':''}


class Taxonomy:
    """Load one frozen taxdump; resolve ranks by ancestry rather than position."""
    def __init__(self,directory):
        self.parent={};self.rank={};self.names={};self.cache={}
        with (Path(directory)/'nodes.dmp').open() as f:
            for line in f:
                fields=[x.strip() for x in line.split('|')]
                self.parent[int(fields[0])]=int(fields[1]);self.rank[int(fields[0])]=fields[2]
        with (Path(directory)/'names.dmp').open() as f:
            for line in f:
                fields=[x.strip() for x in line.split('|')]
                if fields[3]=='scientific name':self.names[int(fields[0])]=fields[1]

    def ancestors(self,taxid):
        seen=set();result=[]
        while taxid and taxid not in seen:
            if taxid not in self.parent:raise ValueError(f'Taxid {taxid} absent from frozen taxonomy')
            seen.add(taxid);result.append(taxid)
            p=self.parent[taxid]
            if p==taxid:break
            taxid=p
        return result

    def annotate(self,taxid):
        if taxid in self.cache:return self.cache[taxid]
        lineage=self.ancestors(taxid);ranks={self.rank[t]:t for t in lineage}
        roots={2:'Bacteria',2157:'Archaea',2759:'Eukarya',10239:'Viruses'}
        domain=next((roots[t] for t in lineage if t in roots),'Unresolved')
        result={'taxid':taxid,'name':self.names.get(taxid,''),'rank':self.rank.get(taxid,''),
                'domain':domain,'lineage':';'.join(self.names.get(t,str(t)) for t in reversed(lineage))}
        for rank in ['phylum','class','order','family','genus','species']:
            result[rank]=self.names.get(ranks.get(rank),'')
        node=ranks.get('class',taxid)
        result.update(taxon=self.names.get(node,str(node)),taxon_id=node,taxon_rank=self.rank.get(node,''))
        self.cache[taxid]=result
        return result
