"""``snakemake`` file that runs analysis."""


import pandas as pd


configfile: "config.yaml"


rule all:
    """Target rule."""
    input:
        sub_counts_csv="results/sub_counts/sub_counts.csv",


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
