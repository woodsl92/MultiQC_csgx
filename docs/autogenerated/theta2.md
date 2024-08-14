---
name: THetA2
urls: ["http://compbio.cs.brown.edu/projects/theta/"]
summary: >
  Estimates tumour purity and clonal / subclonal copy number
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/theta2/theta2.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

The module plots the % germline and % tumour subclone for each sample.
Note that each sample can have multiple maximum likelihood solutions - the MultiQC
module plots proportions for the first one in the results file (`*.BEST.results`).
Also note that if there are more than 5 tumour subclones, their percentages are summed.

### File search patterns

```yaml
theta2:
  fn: "*.BEST.results"
```