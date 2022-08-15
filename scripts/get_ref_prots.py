"""Get the proteins from the reference."""


import Bio.Seq
import Bio.SeqIO

import pandas as pd


ref_fasta = snakemake.input.ref_fasta
ref_gtf = snakemake.input.ref_gtf

ref = str(Bio.SeqIO.read(ref_fasta, "fasta").seq)

prot_coords = (
    pd.read_csv(
        ref_gtf,
        sep="\t",
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    .query("feature == 'CDS'")
    .assign(
        protein=lambda x: x["attribute"].str.extract(
            'gene_id "(\w+)";',
            expand=False,
        ),
    )
    .sort_values(["protein", "start"])
    .assign(
        nt_sequence=lambda x: x.apply(lambda r: ref[r["start"] - 1: r["end"]], axis=1),
    )
    .groupby(["protein", "start", "end"], as_index=False)
    .aggregate(nt_sequence=pd.NamedAgg("nt_sequence", lambda s: "".join(s.values)))
    .assign(
        prot_sequence=lambda x: x["nt_sequence"].map(
            lambda s: str(Bio.Seq.Seq(s).translate()),
        )
    )
)

with open(snakemake.output.ref_prots, "w") as f:
    for tup in prot_coords.itertuples():
        prot = tup.prot_sequence
        assert len(prot) == (tup.end - tup.start + 1) / 3
        print(f"Writing {tup.protein} of length {len(prot)}")
        f.write(f">{tup.protein}\n{prot}\n")
