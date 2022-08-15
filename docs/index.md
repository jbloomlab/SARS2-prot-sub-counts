---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
---

## Overview
This page allows you to assess the mutational tolerance of SARS-CoV-2 proteins by examining the counts of amino acid substitutions at each site.
These counts are over the [UShER](https://usher-wiki.readthedocs.io/) global phylogenetic tree.
The counts represent how many independent times each substitution is inferred to have occurred, **not** alignment frequencies (where a particular mutation might occur just once but then rise to high frequency).
Note the [UShER](https://usher-wiki.readthedocs.io/) tree contains some sequencing errors / or bad base calls, so some of the mutations will spurious (reversions to reference are especially common due to bad base calls).
Nonetheless, these counts should be a good measure of mutational tolerance.

## Plots for each protein
Here are interactive plots for each SARS-CoV-2 protein.
You can mouse over points for statistics and use the zoom bar to zoom into specific sites.


## Data and code
The code used to obtain the counts and generate this site are [on GitHub](https://github.com/jbloomlab{{ site.baseurl }}).

The counts are tabulated in [this CSV file](https://github.com/jbloomlab{{ site.baseurl }}/blob/main/results/sub_counts/sub_counts.csv).

