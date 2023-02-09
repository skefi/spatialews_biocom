# README

This folder contains additional analyses which can be performed on the data (files biocom-grps.rda and data_biocom.rda in the 'outputs' folder), on the spatial metrics (contained in the files indics-data_Nperm_199.rda and indics-model_Nperm_199.rda in the 'outputs' folder) as well the scripts to generate the figures of the paper. 

The files calculate_indicator_trends.R and calculate_model_indicator_trends.R calculate the trends of the spatial metrics along the stress gradient (aridity in the data, parameter b in the model). They generate the files 'trends_one_group_Nperm_199.rda', 'trends_two_groups_Nperm_199.rda', 'trends_three_groups_Nperm_199.rda', 'trends_model_Nperm_199_BOOTN_2999.rda' 

The plot_fig_xx.R files generate the figures of the paper: 
plot_fig_2groups allows to plot figures 2, 3, S7 and S8
plot_fig_slope-indics.R allows to plot figures 2, S3, S10, S14
plot_fig_raw-indics.R allows to plot figures S2, S9, S13

The file clustering_analyses.R explores different clustering approaches to cluster the sites in different groups. 

The file compare_groups_stats.R contains the statistics performed to compare the different group of sites identified. 

