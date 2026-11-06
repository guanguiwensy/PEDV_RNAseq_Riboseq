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
import sample_registry as SR
import csv
O=PEDV_K+r"Supplementary_Figure_L/"
rows=[
("A","families drawn","76 non-canonical-frame families (80 total)","73 families outside the frame of annotated stop codons (80 total; 7 canonical-stop families not drawn)","3 families (E-x1 = M in-frame extension; S-frame extension through ORF1b; ORF3-frame extension through S) had been labelled internal alternative-frame because classification used the CDS containing the first start","Recompute + replot","Figure 6 featured marker: markers kept as provisional (original Figure 6 A/C: ORF1a-x22, N-x1, N-x2; D: ORF1a-x3, ORF1b-x13, ORF1b-x8); old Panel-A red set (S-x1, S-x2, ORF1a-x4 added) replaced by this record; revised-rule list reported as pending"),
("A","frame labels","genomic frame 0/1/2","genomic frame with canonical/+1 relation of each CDS (frame 1 = ORF1a/ORF3/E/N canonical and ORF1b +1)","terminology clarification only","Annotation/legend","unchanged"),
("A","featured marker colour","red outline on 7 families chosen by me (top-5 HAR + N-x2 + ORF1a-x22)","red outline = original Figure 6 featured (provisional); grey label = examined in Panel C but not originally featured","previous selection did not follow the project record","Replot","changed (provisional); pending revised-rule list in SupL_featured_marker_status.tsv"),
("B","families in region table","77 families / 798 start variants","73 families / 583 start variants (>40 aa)","upstream ORF counted in both E and UTR/gap; mislabelled canonical-stop families; start variants summed from by-family counts","Recompute + replot","n/a"),
("B","E region","8.3 families/kb","0 families / kb (11 start sites located in E belong to families whose first start is elsewhere)","E-x1 is the M-frame extension, not an E-internal family","Recompute","n/a"),
("B","N vs other CDS","3.77 vs 2.72 families/kb (RR 1.38, P 0.42); starts/family median 9 vs 5 (P 0.44)","3.77 vs 2.57 (RR 1.47, exact binomial P 0.40); start sites/kb 33.9 vs 20.8; starts/family median 9 vs 5 (Mann-Whitney P 0.37; permutation P 0.82)","comparison set lost canonical-stop families","Recompute + wording","unchanged; wording 'numerically higher, not significant'"),
("C1","start-class denominators","all >=8 aa start variants of each family","only >40 aa start variants (Figure 6 definition)","made consistent with Figure 6 start definition","Recompute + replot","unchanged"),
("C2","family retained definition","any eligible start (including newly gained) with intact stop","primary: >=1 REFERENCE start variant still eligible with intact stop and no in-frame stop; any-start definition kept as sensitivity column","C2 and C3 now use one definition; legend states exactly what is checked","Recompute + replot","unchanged"),
("C3","families in null","74","71","3 canonical-stop families removed","Recompute + replot","unchanged"),
("D","columns / per-genome denominator","80 families","73 families (canonical-stop families excluded: their retention measures canonical CDS integrity)","same reason as A","Recompute + replot","unchanged"),
("D","lineage labels","GI-b n=30; historical 'G1b' 168","unchanged; historical label kept as separate sensitivity row and named 'historical G1b label' only","naming not to be confused with phylogenetic GI-b","Legend","unchanged"),
("all","retention values of featured families","ORF1a-x3 0.985 ... N-x2 0.995 (all PEDV)","identical to 3 decimals except GII-a for S-x1/ORF1a-x22 (0.004/0.032 vs 0.004/0.033) and N-x1 0.9935 vs 0.9945","definition change C2/C3 consistency","Recompute","unchanged"),
]
with open(O+"SupL_old_vs_revised.tsv","w",encoding="utf-8",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["panel","metric","original","revised","reason_for_change","change_type","Figure6_featured_marker"]); w.writerows(rows)
fe=[("ORF1a-x22","original Figure 6 A/C","provisional red outline",f"dropped in revised-rule list (no resolved HAR response with {SR.lib('RIBO_VERO_CHX_1MOI_16h')})"),
("N-x1","original Figure 6 A/C","provisional red outline","dropped in revised-rule list"),
("N-x2","original Figure 6 A/C","provisional red outline","kept (Strong)"),
("ORF1a-x3","original Figure 6 D","provisional red outline","kept (Suggestive)"),
("ORF1b-x13","original Figure 6 D","provisional red outline","dropped (source not locatable; replaced by rule)"),
("ORF1b-x8","original Figure 6 D","provisional red outline","dropped"),
("ORF1a-x9","-","not marked","added by revised rule (Strong)"),
("ORF1a-x23","-","not marked","added by revised rule (Strong)"),
("S-x1","-","examined in C (grey label)","added by revised rule (Suggestive)"),
("S-x2","-","examined in C (grey label)","added by revised rule (Suggestive)"),
("N-x5","-","not marked","added by revised rule (Suggestive)")]
with open(O+"SupL_featured_marker_status.tsv","w",encoding="utf-8",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["family","original_status","status_in_REVISED_SupL","status_in_Figure6_REVISED_rule_draft_(pending)"]); w.writerows(fe)
print("ok")
