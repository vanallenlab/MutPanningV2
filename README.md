# MutPanningV2
This repository contains the source code of the revised version of MutPanning. MutPanning is publicly available under the BSD3-Clause open source license.

MutPanning is designed to detect rare cancer driver genes from aggregated whole-exome sequencing data. Most approaches detect cancer genes based on their mutational excess, i.e. they search for genes with an increased number of nonsynonymous mutations above the background mutation rate. MutPanning further accounts for the nucleotide context around mutations and searches for genes with an excess of mutations in unusual sequence contexts that deviate from the characteristic sequence context around passenger mutations.

INTRODUCTION

MutPanning analyzes aggregated DNA sequencing data of tumor patients to identify genes that are likely to be functionally relevant, based on their abundance of nonsynonymous mutations or their increased number of mutations in unusual nucleotide contexts that deviate from the background mutational process. 
 
The name MutPanning is inspired by the words "mutation" and "panning". The goal of the MutPanning algorithm is to discover new tumor genes in aggregated sequencing data, i.e. to "pan" the few tumor-relevant driver mutations from the abundance of functionally neutral passenger mutations in the background. Previous approaches for cancer gene discovery were mostly based on mutational recurrence, i.e. they detected cancer genes based on their excess of nonsynonymous mutation above the local background mutation rate.  Further, they search for mutations that occur in functionally important genomic positions, as predicted by bioinformatical scores). These approaches are highly effective in tumor types, for which the average background mutation rate (i.e., the total mutational burden) is low or moderate.
 
The ability to detect driver genes can be increased by considering the nucleotide context around mutations in the statistical model. MutPanning utilizes the observation that most passenger mutations are surrounded by characteristic nucleotide sequence contexts, reflecting the background mutational process active in a given tumor. In contrast, driver mutations are localized towards functionally important positions, which are not necessarily surrounded by the same nucleotide contexts as passenger mutations. Hence, in addition to mutational excess, MutPanning searches for genes with an excess of mutations in unusual sequence contexts that deviate from the characteristic sequence context around passenger mutations. That way, MutPanning actively suppresses mutations in its test statistics that are likely to be passenger mutations based on their surrounding nucleotide contexts. Considering the nucleotide context is particularly useful in tumor types with high background mutation rates and high nucleotide context specificity (e.g., melanoma, bladder, endometrial, or colorectal cancer).

ALGORITHM
 
Most passenger mutations occur in characteristic nucleotide contexts that reflect the mutational process active in a given tumor. MutPanning searches for mutations in “unusual” nucleotide contexts that deviate from this background mutational process. In these positions, passenger mutations are rare and mutations are thus a strong indicator of the shift of driver mutations towards functionally important positions.
 
The main steps of MutPanning are as follows (adopted from Dietlein et al.): 
(i) Model the mutation probability of each genomic position in the human exome depending on its surrounding nucleotide context and the regional background mutation rate. 
(ii) Given a gene with n nonsynonymous mutations, use a Monte Carlo simulation approach to simulate a large number of random “scenarios” in which n or more nonsynonymous mutations are randomly distributed along the same gene . 
(iii) Compare the number and positions of mutations in each random scenario with the observed mutations in gene . Based on these comparisons, derive a p-value for the gene. 
(iv) Combine this p-value with additional statistical components that account for insertions and deletions, the abundance of deleterious mutations, and mutational clustering.

USAGE
 
You can run the algorithm for multiple cancer types at the same time. All you need to run MutPanning is a mutation file (*.maf) and a sample file (*.txt):
- Mutation File (*.maf): mutations that you would like to include in the analysis 
- Sample File (*.txt): contains sample IDs and assocates them with cancer types 
Unless you are familiar with the MutPanning algorithm, we recommend running MutPanning with standard parameters. MutPanning was only tested to run with standard parameters. A precompiled desktop version of MutPanning is available on cancer-genes.org. In addition, MutPanning can be run as a module on genepattern.org

REFERENCES
 
Dietlein F, Weghorn D, Taylor-Weiner A, Richters A, et al. Identification of cancer driver genes based on nucleotide context. Under review. (preprint available on biorxiv)
