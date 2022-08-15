"""Get counts of amino-acid substitutions."""


import pandas as pd


sub_counts = (
    pd.read_csv(snakemake.input.mut_tsv, sep="\t")
    .assign(aa_mutations=lambda x: x["aa_mutations"].str.split(";"))
    [["aa_mutations"]]
    .explode("aa_mutations")
    .groupby("aa_mutations", as_index=False)
    .aggregate(count=pd.NamedAgg("aa_mutations", "count"))
    .assign(
        protein=lambda x: x["aa_mutations"].str.split(":").str[0],
        substitution=lambda x: x["aa_mutations"].str.split(":").str[1],
        site=lambda x: x["substitution"].str[1: -1].astype(int),
    )
    .sort_values(["protein", "site"])
    .drop(columns="aa_mutations")
    [["protein", "site", "substitution", "count"]]
)

assert not len(sub_counts.query("'X' in substitution")), "ambiguous amino acids"

sub_counts.to_csv(snakemake.output.sub_counts_csv, index=False)
