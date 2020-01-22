# TransportingMultipleTrials
Example R code to implement estimators in *Towards causally interpretable meta-analysis: transporting inferences from multiple randomized trials to a new target population.* Issa J Dahabreh, Lucia C Petito, Sarah E Robertson, Miguel A Hernán, and Jon A Steingrimsson. [arxiv link here] 

Required R package is geex (to get the sandwich variance)

1. 01_run_haltc_analysis_forestplot_table_github.R: Code that replicates analysis used to create table in paper. The data used in the applied example (halt-c) is not available to share, but the code here can be modified to use in a different dataset.

2. Source Code: Three files of source code for unadjusted (no transportability), one-at-a-time analysis transportability, and transportability using the collection of trials. 

