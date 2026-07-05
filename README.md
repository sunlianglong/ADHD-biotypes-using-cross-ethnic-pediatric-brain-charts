# **Cross-ethnic pediatric brain charts reveal replicable ADHD biotypes with distinct developmental, transcriptomic, and therapeutic profiles**

This repository provides main code, relevant toolbox, and data for manuscript "Bu X#, Sun L#, et al. (2026). Cross-ethnic pediatric brain charts reveal replicable ADHD biotypes with distinct developmental, transcriptomic, and therapeutic profiles"

## **Overview**
Following contents include software, source code, and demo data. All data necessary to replicate our results of Chinese cohort have been made publicly available, including regional cortical thickness, cortical thickness deviations from normative modeling, and demographic and clinical information.

## **Code**
Our analyses were carried out using following open source packages:
1. **Multiscale Desikan-Kiliany parcellation** files [1] were downloaded using the netneurotools toolbox (https://github.com/netneurolab/netneurotools). We use 219 parcellations for our analysis.

2. **Normative modeling** was implemented with the PCNtoolkit (https://pcntoolkit.readthedocs.io/en/latest/) [2].

3. **Spectral clustering** was applied to the similarity matrix of cortical thickness deviation to cluster the ADHD into subgroups. This analysis was performed using an open R package SNFtool (http://compbio.cs.toronto.edu/SNF/SNF/Software.html) [3].

4. **ComBat** was used to harmonize raw cortical thickness across sites for group-level analyses (https://github.com/Jfortin1/ComBatHarmonization) [4].

5. The **Allen Human Brain Atlas** (AHBA) datasets were preprocessed using the Python toolbox abagen (https://github.com/rmarkello/abagen) [5].

6. The **Gene Ontology enrichment analysis** was conducted using online tool Metascape (https://metascape.org) [6].

7. **Visualization**: Python nilearn and R ggplot.

The main analyses were carried out using following codes:
### **1. Normative modeling analysis:**
   running normative modeling: 1.1.NormativeModelTraining.py
   
   plot normative curves: 1.2.plot_curves.py

### **2. Clustering analysis:**
   determining best parameters and cluster number: 2.1.running_parameters_ireta.R
   
   running spectral clustering: 2.2.spectral_clustering_robust.R
   
   evaluating clutering stability by bootstrapping: 2.3.clustering_stability.R

### **3. Subgroup analysis:**
   univariate analysis: 3.1.group_analysis.R
   
   multivariate analysis (brain-behavior relationship): 3.2.pls_symptom_cognition.m
   
   medication response: 3.3.drug_effect.R
   
   gene expression analysis: 3.4.AHBA_gene_expression.py, 3.5.pls_gene_brain.m

   spin test: https://github.com/frantisekvasa/rotate_parcellation.

Data analysis was performed using MATLAB R2021b, Python v3.8, and R v4.2.2.
   
## **Data**
1. Demographic and clinical information for ADHD subtypes: *ADHD_subtype_demographic.csv* and *ADHD_Clinic_Info.xlsx*

2. Follow-up treatment information for ADHD subtypes: *ADHDsubtypes_treatment.csv*

3. Raw regional cortical thickness for each participant: *cortical_thickness_ADHD.csv* and *cortical_thickness_TDC.csv*

4. Cortical thickness devations estimated form normative modeling: *cortical_thickness_deviation_ADHD.csv*

## **References**
1. Cammoun L, Gigandet X, Meskaldji D, et al., Mapping the human connectome at multiple scales with diffusion spectrum MRI. J Neurosci Meth 2012; 203: 386-97.

2. Rutherford S, Kia SM, Wolfers T, et al., The normative modeling framework for computational psychiatry. Nat Protoc 2022; 17: 1711-34.

3. Wang B, Mezlini AM, Demir F, et al., Similarity network fusion for aggregating data types on a genomic scale. Nat Methods 2014; 11: 333-U19.

4. Fortin JP, Cullen N, Sheline YI, et al., Harmonization of cortical thickness measurements across scanners and sites. Neuroimage 2018; 167: 104-20.

5. Markello RD, Arnatkeviciute A, Poline JB, et al., Standardizing workflows in imaging transcriptomics with the abagen toolbox. Elife 2021; 10.

6. Zhou YY, Zhou B, Pache L, et al., Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nat Commun 2019; 10.
