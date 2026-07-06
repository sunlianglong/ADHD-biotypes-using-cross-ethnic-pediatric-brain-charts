# **Population-specific pediatric brain charts reveal replicable ADHD biotypes across ethnicities**

This repository provides the main code, relevant toolboxes, pretrained normative models, and shareable data accompanying the manuscript:

> Bu X#, Sun L#, et al. (2026). *Cross-ethnic pediatric brain charts reveal replicable ADHD biotypes with distinct developmental, transcriptomic, and therapeutic profiles*.

## Overview

This repository includes software dependencies, source code, pretrained pediatric brain chart models, and shareable data used in the manuscript.

For the Chinese cohort, we provide the data necessary to reproduce the main analyses, including regional cortical thickness, pediatric normative model, cortical thickness deviations estimated from normative modeling, and demographic and clinical information.

For the US cohort based on the ABCD Study, individual-level cortical thickness and derived deviation data are not redistributed in this repository because ABCD data access and reuse are governed by the NIH/NBDC Data Use Certification process. Users who wish to reproduce analyses involving ABCD individual-level data should apply for access through the official ABCD/NBDC data access process: [ABCD Study Data Sharing](https://abcdstudy.org/scientists/data-sharing/) and [NBDC Data Access Process](https://www.nbdc-datahub.org/data-access-process). To facilitate reproducibility without redistributing controlled individual-level ABCD-derived data, we provide the pretrained US pediatric normative model in the `Model` folder.

## Code

The analyses were carried out using the following open-source packages and resources:

1. **Multiscale Desikan–Killiany parcellation** files were downloaded using the `netneurotools` toolbox (https://github.com/netneurolab/netneurotools). We used 219 cortical parcels in the analysis [1].

2. **Normative modeling** was implemented using `PCNtoolkit` (https://pcntoolkit.readthedocs.io/en/latest/) [2].

3. **Spectral clustering** was applied to the similarity matrix of cortical thickness deviations to identify ADHD subgroups. This analysis was performed using the R package `SNFtool` (http://compbio.cs.toronto.edu/SNF/SNF/Software.html) [3].

4. **ComBat harmonization** was used to harmonize raw cortical thickness across sites for group-level analyses (https://github.com/Jfortin1/ComBatHarmonization) [4].

5. **Allen Human Brain Atlas** datasets were preprocessed using the Python toolbox `abagen` (https://github.com/rmarkello/abagen) [5].

6. **Gene Ontology enrichment analysis** was conducted using the online tool Metascape (https://metascape.org) [6].

7. **Visualization** was performed using Python `nilearn` and R `ggplot2`.

The main analyses were carried out using the following scripts.

### 1. Normative modeling analysis

* `1.1.NormativeModelTraining.py`
  Training pediatric normative models.

* `1.2.plot_curves.py`
  Plotting normative developmental curves.

### 2. Clustering analysis

* `2.1.running_parameters_iter.R`
  Determining optimal parameters and cluster number.

* `2.2.spectral_clustering_robust.R`
  Running spectral clustering.

* `2.3.clustering_stability.R`
  Evaluating clustering stability using bootstrapping.

### 3. Subgroup analysis

* `3.1.group_analysis.R`
  Univariate subgroup analyses.

* `3.2.pls_symptom_cognition.m`
  Multivariate brain–behavior association analysis.

* `3.3.drug_effect.R`
  Medication response analysis.

* `3.4.AHBA_gene_expression.py` and `3.5.pls_gene_brain.m`
  Imaging transcriptomic analysis.

* Spin test
  Spatial permutation testing was performed using the `spin-test toolbox`: https://github.com/frantisekvasa/rotate_parcellation.

Data analyses were performed using MATLAB R2021b, Python v3.8, and R v4.2.2.

## Data

The `Data` folder contains the following shareable data files of the Chinese cohort.

1. **Demographic and clinical information for ADHD subtypes**

   * `ADHD_subtypes_demographic.csv`
   * `ADHD_Clinic_Info.xlsx`

2. **Follow-up treatment information for ADHD subtypes**

   * `ADHD_subtypes_treatment.csv`

3. **Raw regional cortical thickness**

   * `cortical_thickness_ADHD.csv`
   * `cortical_thickness_TDC.csv`

4. **Cortical thickness deviations estimated from normative modeling**

   * `cortical_thickness_deviation_ADHD.csv`

Please note that individual-level cortical thickness and derived cortical thickness deviation data from the ABCD Study are not included in this repository because of ABCD/NBDC data sharing and controlled-access policies. Researchers who need to reproduce the US cohort analyses should obtain ABCD data access through the official ABCD/NBDC Data Use Certification process.

## Model

The `Model` folder contains the pretrained pediatric normative models used in this study. These models can be applied to estimate individual-level cortical thickness deviations for participants from new sites using `PCNtoolkit` (https://pcntoolkit.readthedocs.io/en/latest/) [2].

1. **Chinese pediatric brain chart**

   This model was trained using the Chinese typically developing cohort and was used to estimate individual cortical thickness deviations for Chinese children with ADHD.

2. **US pediatric brain chart**

   This model was trained using the US typically developing cohort from the ABCD Study and was used to estimate individual cortical thickness deviations for children with ADHD in the US cohort.

## References

1. Cammoun L, Gigandet X, Meskaldji D, et al. Mapping the human connectome at multiple scales with diffusion spectrum MRI. *Journal of Neuroscience Methods* 2012; 203: 386–397.

2. Rutherford S, Kia SM, Wolfers T, et al. The normative modeling framework for computational psychiatry. *Nature Protocols* 2022; 17: 1711–1734.

3. Wang B, Mezlini AM, Demir F, et al. Similarity network fusion for aggregating data types on a genomic scale. *Nature Methods* 2014; 11: 333–337.

4. Fortin JP, Cullen N, Sheline YI, et al. Harmonization of cortical thickness measurements across scanners and sites. *NeuroImage* 2018; 167: 104–120.

5. Markello RD, Arnatkeviciute A, Poline JB, et al. Standardizing workflows in imaging transcriptomics with the abagen toolbox. *eLife* 2021; 10: e72129.

6. Zhou Y, Zhou B, Pache L, et al. Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. *Nature Communications* 2019; 10: 1523.
