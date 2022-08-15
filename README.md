# Substitution counts for SARS-CoV-2 proteins
This repo contains an analysis by Jesse Bloom Will Hannon that counts the number of occurrences of each amino-acid substitution across the [UShER phylogenetic tree](https://usher-wiki.readthedocs.io/).
These counts provide a way to assess mutational tolerance / flexibility at different sites.

The basic approach is simply to get the counts of all amino-acid mutations across the SARS-CoV-2 global tree from `UShER` and then plot them.

The pipeline is run by the [snakemake](https://snakemake.readthedocs.io/) file [Snakefile](Snakefile), which uses the configuration in [config.yaml](config.yaml).
To use a newer version of the `UShER` matrix-annotated tree, update [config.yaml](config.yaml).

To run the pipeline, build the `conda` environment in [environment.yml](environment.yml), then activate it and run the pipeline:

    conda activate SARS2-prot-sub-counts
    snakemake -j 1

A few results are tracked in this repo:

 - The substitution counts are in [results/sub_counts/sub_counts.csv](results/sub_counts/sub_counts.csv)
 - The proteins in the reference genome are in [results/ref/prots.fa](results/ref/prots.fa)
 - Interactive plots of the substitution coutns are in [results/sub_counts/plots](results/sub_counts/plots)
t

The plots are visualized via GitHub pages served from [./docs/](docs) at [https://jbloomlab.github.io/SARS2-prot-sub-counts/](https://jbloomlab.github.io/SARS2-prot-sub-counts/)
