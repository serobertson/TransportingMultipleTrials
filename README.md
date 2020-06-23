# TransportingMultipleTrials
*Towards causally interpretable meta-analysis: transporting inferences from multiple randomized trials to a new target population.* Issa J Dahabreh, Lucia C Petito, Sarah E Robertson, Miguel A Hernán, and Jon A Steingrimsson. Epidemiology (2020). https://journals.lww.com/epidem/Abstract/2020/05000/Toward_Causally_Interpretable_Meta_analysis_.4.aspx

Example R code to implement estimators. 

Required R package is geex (to get the sandwich variance).

1. 01_run_haltc_analysis_forestplot_table_github.R: Code that replicates analysis used to create table in paper. The data used in the applied example (halt-c) is not available to share, but the code here can be modified to use in a different dataset.

2. Source Code: Three files of source code for unadjusted (no transportability), one-at-a-time analysis transportability, and transportability using the collection of trials. 

