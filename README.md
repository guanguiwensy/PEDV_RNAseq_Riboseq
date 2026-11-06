# PEDV transcriptomics and translatomics

Code for integrated RNA-seq, Ribo-seq, and comparative-genomics analyses of porcine epidemic diarrhea virus (PEDV).

## Repository layout

The analysis package is in `PEDV_transcriptomics_translatomics/`:

- `rnaseq/`: alignment, junction calling, classification, abundance, and junction dynamics.
- `riboseq/`: alignment, P-site offsets, periodicity, PRICE, HARTIS, terminal translation, and PK1 comparisons.
- `comparative_genomics/`: sequence collection, lineage analysis, ORF scans, conservation, and null models.
- `figures/`: figure-generation scripts.
- `config/`: shared paths, sample registry, and fixed analysis inputs.
- `metadata/`: the 30-library sample table, public-data sources, and genome accessions.
- `environment/`: Python and server environment definitions.
- `tests/`: reproduction entry points and the sample-registry audit utility.

## Data configuration

Raw data and large intermediate files are not stored in this repository. Scripts resolve data roots in this order:

1. environment variables;
2. `config/paths.local.yaml`;
3. directories below `PEDV_transcriptomics_translatomics/data_root/`.

Copy `config/paths.example.yaml` to `config/paths.local.yaml` when custom locations are needed. The configured roots must be writable because some analyses place outputs beside their inputs.

## Environments

For workstation analyses:

```bash
python -m pip install -r PEDV_transcriptomics_translatomics/environment/requirements.txt
```

Alignment, BAM processing, and other server workflows use the Conda environment in `environment/environment.yml`. These workflows require Linux tools such as STAR, samtools, minimap2, bowtie2, and pysam. PRICE additionally requires Gedi and OpenJDK 11.

## Samples and public data

`metadata/sample.txt` is the authoritative registry for the 30 in-house libraries. External accessions used by the analyses are listed in `metadata/public_data_sources.tsv`.

The in-house sequencing data are associated with NCBI BioProject `PRJNA1554149`; run accessions and release dates are recorded in `metadata/sample.txt`.

## Reproduction

After configuring populated data roots, run the relevant module scripts or use `tests/repro_commands.sh` as an entry-point list. The commands are not a single end-to-end workflow because server-generated inputs are required by several downstream analyses.

## License

MIT License. See `LICENSE`.
