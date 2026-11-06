#!/usr/bin/env bash
# Environment variables expected by the server-side shell scripts (copy, edit, then `source` it).
# NEVER place caches/temp files on the root partition ( / , /opt , /usr , /tmp ).
export PEDV_SERVER_ROOT=/path/to/PEDV_ribo           # data root that holds refdata_GIb_v1.1, analysis_GIb, ...
export PEDV_TMPDIR=/path/to/data_disk/tmp             # scratch space, must be on a data disk
export STAR_BIN=STAR                                   # STAR 2.7.11b
export STAR_AVX2_BIN=STAR-avx2                         # same version, AVX2 build (alignment of Ribo-seq)
export MINIMAP2_BIN=minimap2                           # minimap2 2.30
export SAMTOOLS_BIN=samtools                           # samtools 1.22.1
export PYTHON_BIN=python3                              # Python 3.10/3.11 with numpy, pandas, pysam
export BOWTIE2_BIN_DIR=/path/to/bowtie2/bin            # directory that contains bowtie2 / bowtie2-build
export GEDI_ENV=/path/to/conda_env_gedi                # PRICE 1.0.5 (Gedi 1.0.6f, JDK 11), built from source
export GEDI_HOME=/path/to/gedi                         # Gedi installation
export PEDV_ENV_ORF=/path/to/env_orf.sh                # environment file sourced by the PRICE scripts (NOT included in the repository, see docs/SERVER_SCRIPTS_REQUIRED.md)
export PEDV_REPO=/path/to/this/repository                  # root of this repository (shell scripts resolve sample numbers through config/sample_registry.py)
