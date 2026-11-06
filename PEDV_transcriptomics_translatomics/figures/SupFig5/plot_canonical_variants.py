#!/usr/bin/env python3
"""Reproduce canonical PEDV sgRNA exact-breakpoint variant figure from Fig2D-G arcs CSV.
Usage: python plot_canonical_variants.py --input /path/to/Fig2D-G_arcs_source_data.csv --output /path/to/output
No raw sequencing files are accessed. Events are exact donor/acceptor pairs, not unique RNA molecules.
"""
# >>> release-candidate path shim (added during repository curation; no analytical change) >>>
import os as _os, sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
while _d != _os.path.dirname(_d) and not _os.path.isdir(_os.path.join(_d, "config")):
    _d = _os.path.dirname(_d)
_sys.path.insert(0, _os.path.join(_d, "config"))
from pedv_paths import PEDV_P, PEDV_K, PEDV_R, PEDV_COV, PEDV_SERVER_ROOT, REPO_ROOT, SCRIPT_DIRS, ORIGIN
_sys.path.extend(_x for _x in SCRIPT_DIRS if _x not in _sys.path)  # sibling modules of other repo modules (own directory keeps priority)
__REPO_FILE__ = globals().get("__REPO_FILE__") or _os.path.abspath(__file__)  # real location of this script inside the repository (sibling-script loading must use this, not the virtual __file__)
_o = ORIGIN.get(_os.path.relpath(_os.path.abspath(__file__), REPO_ROOT).replace(chr(92), "/"))
if _o:
    __file__ = _o  # virtual ORIGINAL location inside the working tree: keeps the scripts' own dirname(__file__)-relative path logic unchanged
# <<< release-candidate path shim <<<
import argparse, csv, math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

ORDER=['S','ORF3','E','M','N']
COLORS={'S':'#9ecae1','ORF3':'#6baed6','E':'#4292c6','M':'#2171b5','N':'#08306b'}


def read_events(path):
    records=[]
    with open(path, newline='', encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            if r['class']!='canonical': continue
            if r['canonical_gene'] not in ORDER: continue
            records.append(dict(gene=r['canonical_gene'],donor=int(r['donor']),acceptor=int(r['acceptor']),reads=int(r['reads']),n_samples=int(r['n_samples'])))
    return records


def save_csv(path, rows, fields):
    with open(path,'w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');wr.writeheader();wr.writerows(rows)


def main(input_file, output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    recs=read_events(input_file)
    assert len(recs)==118, f'Expected 118 canonical events, found {len(recs)}'
    assert sum(r['reads'] for r in recs)==2555133, 'Canonical read sum mismatch'
    groups={g:sorted((r for r in recs if r['gene']==g),key=lambda x:-x['reads']) for g in ORDER}
    summaries=[]
    enriched=[]
    for gene in ORDER:
        gg=groups[gene]; total=sum(r['reads'] for r in gg);mode=gg[0]
        summaries.append({'gene':gene,'n_exact_events':len(gg),'total_junction_observations':total,'top1_reads':mode['reads'],'top1_fraction':mode['reads']/total,'top3_fraction':sum(r['reads'] for r in gg[:3])/total,'events_in_all_20_libraries':sum(r['n_samples']==20 for r in gg),'events_in_ge_10_libraries':sum(r['n_samples']>=10 for r in gg)})
        for rank,r in enumerate(gg,1):
            enriched.append(dict(**r,rank_within_gene=rank,percent_within_gene=100*r['reads']/total,modal_donor=mode['donor'],modal_acceptor=mode['acceptor'],donor_offset=r['donor']-mode['donor'],acceptor_offset=r['acceptor']-mode['acceptor']))
    save_csv(out/'canonical_exact_breakpoint_variants.csv',enriched,list(enriched[0]))
    save_csv(out/'canonical_variant_summary.csv',summaries,list(summaries[0]))

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(8.1,10.5), facecolor='white')
    gs=fig.add_gridspec(3,1,height_ratios=[2.55,1.85,2.15],top=.95,bottom=.075,left=.10,right=.96,hspace=.55)
    fig.text(.095,.978,'Canonical PEDV sgRNA junctions exhibit preferred regions and exact-coordinate variants',weight='bold',size=11)

    # A
    ga=gs[0].subgridspec(1,5,wspace=.28)
    for i,g in enumerate(ORDER):
        ax=fig.add_subplot(ga[0,i]); gg=[r for r in enriched if r['gene']==g]
        maxval=max(r['reads'] for r in gg)
        for r in gg:
            sz=10+210*(r['reads']/maxval)**0.43
            ax.scatter(r['donor_offset'],r['acceptor_offset'],s=sz,color=COLORS[g],edgecolor='#1b2832',lw=.35,alpha=.82,zorder=3)
        ax.axhline(0,color='#bbbbbb',lw=.7);ax.axvline(0,color='#bbbbbb',lw=.7)
        ax.plot([-6,8],[-6,8],ls='--',lw=.6,color='#bbbbbb')
        ax.set(xlim=(-6.5,8),ylim=(-6.5,8),xticks=[-5,0,5],yticks=[-5,0,5],title=g)
        ax.set_aspect('equal',adjustable='box')
        if i==0:ax.set_ylabel("3′ acceptor offset (nt)")
        ax.set_xlabel("5′ donor offset (nt)",fontsize=7)
    fig.text(.064,.949,'A',weight='bold',size=14)
    fig.text(.10,.704,'Offsets are relative to the most-supported exact junction in each sgRNA class; bubble area encodes supporting observations.',fontsize=7,color='#505050')

    # B
    ax=fig.add_subplot(gs[1]); shades=['#08519c','#3182bd','#6baed6','#dddddd']
    for i,g in enumerate(ORDER):
        gg=groups[g];total=sum(r['reads'] for r in gg);left=0
        for rank in range(3):
            r=gg[rank]; frac=r['reads']/total
            ax.barh(i,frac,left=left,color=shades[rank],height=.71,edgecolor='white',linewidth=.8)
            if frac>.17: ax.text(left+frac/2,i,f"{r['donor']}→{r['acceptor']}\n{100*frac:.1f}%",ha='center',va='center',color='white',fontsize=7)
            left+=frac
        ax.barh(i,1-left,left=left,color=shades[3],height=.71,edgecolor='white',linewidth=.8)
        ax.text(1.015,i,f"Top1 {100*gg[0]['reads']/total:.1f}%  |  {len(gg)} variants",va='center',fontsize=7)
    ax.set(xlim=(0,1.34),yticks=list(range(5)),yticklabels=ORDER,xlabel='Fraction of junction observations within each canonical sgRNA class')
    ax.set_xticks([0,.25,.5,.75,1]);ax.invert_yaxis()
    ax.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(shades,['Top 1','Top 2','Top 3','Other'])],ncol=4,frameon=False,loc='upper left',bbox_to_anchor=(.0,1.18),fontsize=7)
    fig.text(.064,.656,'B',weight='bold',size=14)

    # C
    ax=fig.add_subplot(gs[2])
    for g in ORDER:
        gg=groups[g]; m=max(r['reads'] for r in gg)
        ax.scatter([r['n_samples'] for r in gg],[r['reads'] for r in gg],s=[12+145*(r['reads']/m)**.45 for r in gg],color=COLORS[g],edgecolors='#333333',linewidths=.4,alpha=.85,label=g)
    # Exact variant labels are provided in panel B and the source CSV; avoid cluttering the recurrence plot.
    ax.set(yscale='log',xlim=(.5,23.3),ylim=(.75,7e6),xticks=[1,5,10,15,20],xlabel='Libraries with supporting junction observations (out of 20)',ylabel='Total junction observations (log scale)')
    ax.grid(axis='y',color='#ededed',lw=.7,zorder=0)
    ax.legend(ncol=5,frameon=False,loc='upper left',fontsize=7)
    fig.text(.064,.335,'C',weight='bold',size=14)
    fig.text(.10,.025,'Counts are exact breakpoint-coordinate calls; local microhomology may produce alignment-coordinate ambiguity. CHX/HAR are not biological replicates.',size=7,color='#505050')
    fig.savefig(out/'canonical_junction_variants.pdf',bbox_inches='tight')
    fig.savefig(out/'canonical_junction_variants.jpg',dpi=300,bbox_inches='tight',pil_kwargs={'quality':94,'subsampling':0})
    fig.savefig(out/'canonical_junction_variants.png',dpi=300,bbox_inches='tight')
    plt.close(fig)
    print('Saved figure and tables to',out)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='Fig2D-G_arcs_source_data.csv');ap.add_argument('--output',default='.')
    a=ap.parse_args();main(a.input,a.output)
