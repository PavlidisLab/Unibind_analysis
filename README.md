---
editor_options: 
  markdown: 
    wrap: 80
---

# Unibind_analysis

This projects downloads mouse and human ChIP-seq data from the Unibind database,
organizing this information and gene-scoring each experiment for various
projects in the Pavlab.

<https://unibind.uio.no/>

These are the main deliverables:

1.  A gene score vector is created for each experiment, following the method
    used in Morin et al., 2023 (<https://github.com/PavlidisLab/TR_aggregation>).
    This information is bound into a gene by experiment matrix of binding scores.

    Unibind experiments may be "duplicated" in that the same experiment is matched
    to multiple relevant TF motifs. In all such cases, I average the resulting
    binding scores so that each experiment is only represented once.

    Experiments that have fewer than 100 peaks/regions are discarded. This is a
    fairly lenient filter.

    Finally, I log2 transform and quantile normalize the de-duplicated binding
    matrix. This is what I use for most downstream analysis.

2.  A list of basic metadata is saved out. This only includes the ChIP'd symbol,
    the experiment/file name, the ID, whether the given experiment was
    "duplicated", and the count of peaks in the experiment.

    The "File" column is really the unique identifier as it includes info on which
    motif was used. But given that I de-duplicate/average the experiments in the
    processed data, each experiment is represented by one ID which is used to name
    the columns of the processed matrix.

3.  All experiments are saved out as a list of Genomic Ranges objects. This
    includes "duplicated" experiments (they are not collapsed). This is to
    provide convenience for overlap functions.

4.  An aggregate binding profile for each available TF, which averages the TF's
    set of binding vectors. These aggregate profiles were used for the ChIP-seq
    comparison in my single cell coexpression paper 
    (<https://github.com/PavlidisLab/TR_singlecell/>).

    This list object also contains two further "global" summaries that can be used
    to identify which genes are commonly bound, regardless of the ChIP'd TF. The
    first averages binding scores across all experiments -- this is influenced by
    imbalanced dataset counts across TFs. The second averages within each set of TF
    experiments and then averages these aggregates, such that each TF has equal
    weight towards the global profile.

5.  A linear model is fit per TF, wherein a one-vs-rest contrast is fit to
    compare a TF's binding scores versus all other TF's binding scores. The
    model controls for the number of peaks for each experiment, but does not
    account for imbalanced coverage (i.e., the "rest" comparison will be
    influenced by the most commonly represented TFs.)

The processed bind score matrices and associated metadata can also be found in
this data repository (09_Scored_Unibind_ChIPseq_data.RDS) <https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP3/HJ1B24>
