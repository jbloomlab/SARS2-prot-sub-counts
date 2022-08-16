"""``snakemake`` file that runs analysis."""


import pandas as pd


configfile: "config.yaml"


rule all:
    """Target rule."""
    input:
        "results/sub_counts/sub_counts.csv",
        "results/ref/prots.fa",
        "results/sub_counts_plots",
        "docs/_includes/vega_plots.html"


rule get_mat_tree:
    """Get the pre-built mutation-annotated tree."""
    params:
        url=config["mat_tree"],
    output:
        mat="results/mat/mat_tree.pb.gz"
    shell:
        "curl {params.url} > {output.mat}"


rule get_ref_fasta:
    """Get the reference FASTA."""
    params:
        url=config["ref_fasta"],
    output:
        ref_fasta="results/ref/ref.fa",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_fasta}"


rule get_ref_gtf:
    """Get the reference FASTA."""
    params:
        url=config["ref_gtf"],
    output:
        ref_gtf="results/ref/ref.gtf",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_gtf}"


rule translate_mat:
    """Translated amino-acid mutations on mutation-annotated tree."""
    input:
        mat=rules.get_mat_tree.output.mat,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        ref_gtf=rules.get_ref_gtf.output.ref_gtf,
    output:
        tsv="results/sub_counts/translated_mutations.tsv",
    shell:
        """
        matUtils summary \
            -i {input.mat} \
            -g {input.ref_gtf} \
            -f {input.ref_fasta} \
            -t {output.tsv}
        """


rule get_sub_counts:
    """Counts of amino-acid substitutions."""
    input:
        mut_tsv=rules.translate_mat.output.tsv,
    output:
        sub_counts_csv="results/sub_counts/sub_counts.csv",
    script:
        "scripts/get_sub_counts.py"


rule get_ref_prots:
    """Get sequence of proteins in reference."""
    input:
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        ref_gtf=rules.get_ref_gtf.output.ref_gtf,
    output:
        ref_prots="results/ref/prots.fa",
    script:
        "scripts/get_ref_prots.py"


rule plot_sub_counts:
    """Plot counts of amino-acid substitutions."""
    input:
        sub_counts_csv=rules.get_sub_counts.output.sub_counts_csv,
        ref_prots=rules.get_ref_prots.output.ref_prots,
    output:
        sub_count_plotsdir=directory("results/sub_counts_plots"),
        sub_count_plots_concat="docs/_includes/vega_plots.html"
    log:
        notebook="results/logs/plot_sub_counts.ipynb",
    notebook:
        "notebooks/plot_sub_counts.py.ipynb"
