#!/usr/bin/env python3
import sys
from pathlib import Path

root=Path(sys.argv[1]);figures=['Figure1B_CARPDM','Figure2','FigureS1_temporal','FigureS2_filter_support',
 'FigureS3_genus_affiliations','FigureS4_gene_heatmap','FigureS4_gene_overlap','FigureS5_depth_diversity',
 'FigureS6_networks','FigureS7_sensitivity','FigureS8_partial_RDA']
missing=[name for name in figures if not (root/'figures'/(name+'.pdf')).is_file()]
missing.extend(f'TableS{i}' for i in range(1,13) if not list((root/'tables').glob(f'TableS{i}_*.tsv')))
if missing:raise SystemExit('Incomplete report outputs: '+', '.join(missing))
