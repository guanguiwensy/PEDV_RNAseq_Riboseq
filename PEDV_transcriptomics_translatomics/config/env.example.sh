#!/usr/bin/env bash
# Environment variables expected by the server-side shell scripts (copy, edit, then `source` it).
# NEVER place caches/temp files on the root partition ( / , /opt , /usr , /tmp ).
export PEDV_SERVER_ROOT=/data/pedv/server             # data root that holds refdata_GIb_v1.1, analysis_GIb, ...
export PEDV_TMPDIR=/data/pedv/tmp                     # scratch space, must be on a data disk
export STAR_BIN=STAR                                   # STAR 2.7.11b
export STAR_AVX2_BIN=STAR-avx2                         # same version, AVX2 build (alignment of Ribo-seq)
export MINIMAP2_BIN=minimap2                           # minimap2 2.30
export SAMTOOLS_BIN=samtools                           # samtools 1.22.1
export PYTHON_BIN=python3                              # Python 3.10/3.11 with numpy, pandas, pysam
export BOWTIE2_BIN_DIR=/opt/bowtie2/bin               # directory that contains bowtie2 / bowtie2-build
export GEDI_ENV=/opt/conda/envs/gedi                   # PRICE 1.0.5 (Gedi 1.0.6f, JDK 11), built from source
export GEDI_HOME=/opt/gedi                             # Gedi installation
export PEDV_ENV_ORF=/data/pedv/config/env_orf.sh      # environment file sourced by the PRICE scripts
export PEDV_REPO=/opt/pedv_transcriptomics_translatomics
