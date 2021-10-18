# KidneyDiseasesExplorer

Kidney Diseases Eplorer is a shiny app, that makes it easy to visualize some data about proteins enrolled in kidney diseases.

## Functionalities
The application proposes 4 main functionalities : 
- Comparing proteins from 2 datasets : display a scatter-plot containing one point for each protein that is in both datasets. One can choose to display log2 of t-test’s differences (log2(Fold-changes)), or -log2 of t-test’s p-values (-log2(p-values)). Plot can be downloaded as PNG and SVG.
- Comparing statistics for a selected set of proteins and datasets : display a heatmap with datasets on the y-axis and gene names on the y-axis. One can choose to display log2 of t-test’s differences (log2(Fold-changes)), or -log2 of t-test’s p-values (-log2(p-values)). Plot can be downloaded as PNG and SVG.
- Visualizing selected data in a datatable and downloading them as CSV, Microsoft Excel or PDF.
- Downloading the original dataset as CSV : for each dataset there are 3 variables : Gene name, ttest_difference1 which is log2 of t-test’s fold changes and ttest_pvalue1 which is -log2 of t-test’s p-value.

## How to use it?
<strong> Parameters selection : one can select parameters in the sidebar. </strong> <br>
There are 4 parameters that can be set : <br>
- Selection of proteins : here the user can select proteins he desires to study. There are 4 ways of selected proteins (that can be combined)
  -  By Gene names : select or enter name of gene(s) associated to proteins one wants to study (if enter manually, it has to be the human gene name except if it is only known for mouse or rat)
  - By keywords : select uniprot keywords with which proteins must be annotated
  - By gene ontology : select GO term(s) with which proteins must be annotated (GO are divided into the 3 GO axis : biological process, cellular component and molecular function)
- Set a minimum number of datasets in which proteins must be. By default it is 1. <br>
<strong>WARNING: </strong> We recommend not to use this last criterion alone, since for all choices there would be many selected proteins, so computation time would be very high. <br>

For the 3 first criterions, in each case it is possible to do multiple selection. <br>
When several criterions are selected, or when multiple selection is made, only proteins included in their intersection are selected. There is an exception for gene names selection for the scatter plot : in that case, proteins that meet all the other requirements are selected and gene names are only used to be printed on the scatter plot. <br>
<br>
<strong> The 3 others parameters that can be modified in the sidebar are : </strong> <br>
- Confidence level for t-test : One can choose between “no threshold”, “95%” and “99%”, default is “no threshold”. Only proteins whose confidence level is higher or equal to this threshold are selected. <br>
- Show only one copy per gene? : Default is “No”. In some datasets, one protein appears several times (several rows correspond to the same proteins). If one wants to display one row per protein, one can set this parameter to “Yes”. In that case, one row is randomly selected. <br>
- Selection of dataset for the heatmap : One can select  datasets to display on the heatmap. By default all datasets are selected. <br>
In order to apply the selection, one must click on the “Click here to apply changes” button. <br>

<strong> Choose parameters for the scatter plot : </strong> <br>
In this box, one can choose parameters to apply to the scatter plot, in addition to parameters already selected in the sidebar. <br>
- Select datasets to use as axis <br>
- Select whether to display density : default is “No” <br>
- Select whether to display only selected proteins or to display all proteins and annotate the selected ones <br>
In order to apply changes, one must click on the “Display” button.

## Outputs
The outputs of the app consist of 4 boxes in the application : <br>
- A scatter-plot containing one point for each protein that is in both datasets
- A heatmap : gene_names of the heatmap terminate by ”_k”, the “k” means that it is the k th copy of that protein. 
- Datatables of selected data 
- 10 first rows original datasets <br>

<strong> Special cases : </strong> <br>
- Scatter plot : if no proteins fulfill the requirements, a message is printed explaining that there is nothing to display
- Heatmap : 
  - if no proteins fulfill the requirements, a message is printed explaining that there is nothing to display
  - If no selection criterion is selected, a message is printed to ask for selection
  - If only one protein is selected, the heatmap cannot be computed, and a message explains it to the user
  - if more than 1000 proteins are selected, the heatmap is not displayed and a message invites the user to add more filters to its selection

## Download outputs
- Scatter plots and heatmaps cans be download in PNG and SVG format
- Datatable of proteins under study can be downloaded in CSV, Excel and PDF format
- Original datasets can be downloaded as CSV.




