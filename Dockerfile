# Variant-calling pipeline: BWA + samtools + bcftools + Snakemake,
# installed via bioconda for exact, reproducible tool versions.

FROM continuumio/miniconda3:latest

WORKDIR /pipeline
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Activate the environment for all subsequent RUN/CMD instructions
SHELL ["conda", "run", "-n", "variant-calling-pipeline", "/bin/bash", "-c"]

COPY . /pipeline

# Default: download the real data, run the pipeline, then validate the
# output against the published benchmark counts and the real causal mutation
CMD conda run -n variant-calling-pipeline bash -c "\
    bash workflow/get_data.sh && \
    snakemake -s workflow/Snakefile --cores 2 && \
    pytest tests/ -v && \
    python3 tests/validate_results.py results \
    "
