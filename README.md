# KidneyDiseasesExplorer

Kidney Diseases Explorer is a shiny app, that makes it easy to visualize some data about proteins enrolled in kidney diseases.

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

## Datasets
<strong> cPodocyte_G195D_WT_human	hPodocult_Actn4Mutation	The effect of G195D Mutation in cultured human podocytes vs control </strong><br>
This dataset is the first of two derived from PMID: 26740551 (https://pubmed.ncbi.nlm.nih.gov/26740551/). The study describes a case of juvenile familial focal segmental glomerulosclerosis (FSGS) due to a G195D mutation in the ACTN4 gene coding for Alpha-actinin-4. The comparison in the paper and in the app is between stable podocyte cell line transfected with G195D expressing plasmid and wild-type podocyte cell line. Proteins significantly perturbed by G195D presence included the PDLIM proteins PDLIM1, 2, 4 and 7 as well as the LIM-domain protein ZNF185. <br>

<strong> PUC_G195D_WT_con_human	hPUC_ACTN4Mutation	The effect of G195D Mutation in G195D PUCs vs control	 </strong> <br>
This dataset is the second of two derived from PMID: 26740551 (https://pubmed.ncbi.nlm.nih.gov/26740551/). The study describes a case of juvenile familial focal segmental glomerulosclerosis (FSGS) due to a G195D mutation in the ACTN4 gene coding for Alpha-actinin-4. The comparison in the paper and in the app is between primary renal epithelial cells from urine (hPUCs) of the patient carrying the G195D mutation compared with hPUCs of healthy controls. Proteins significantly decreased by G195D presence included LIM-domain proteins ZNF185, PDLIM4, LMO7, LIMA1 and CRSP.<br>

<strong> cPodocyte_PAN_cultured_human	hPodocult_PAN_33	The effect of PAN on cultured human podocytes vs control (33 degree)	</strong> <br>
This dataset is the first of four derived from PMID: 28400537 (https://pubmed.ncbi.nlm.nih.gov/28400537/). The study describes a proteomic analysis before induction of injury of podocytes and at specific time points after induction but before the development of proteinuria to define the cellular effects of puromycin aminonucleoside (PAN). In this particular dataset the effect of puromycin aminonucleoside (PAN) on in-vitro undifferentiated podocytes (cultured at 33 °C) is documented. The comparison in the paper and in the app is between non-differentiated (proliferating) PAN-treated podocytes and vehicle-treated control podocytes. After PAN-treatment, 14 proteins decreased significantly, among these were YAP/TAZ target genes, including ANLN, ANKRD1, CYR61, CTGF, and DIAPH3. This effect depends on concentration and duration of PAN-treatment, which suggest that PAN-induced damage in non-differentiated podocytes inhibits YAP/TAZ function in vitro.<br>

<strong> cPodocyte_PAN_cultured_human	hPodocult_PAN_37	The effect of PAN on cultured human podocytes vs control (37 degree)	</strong> <br>
This dataset is the second of four derived from PMID: 28400537 (https://pubmed.ncbi.nlm.nih.gov/28400537/). The study describes a proteomic analysis before induction of injury of podocytes and at specific time points after induction but before the development of proteinuria to define the cellular effects of puromycin aminonucleoside (PAN). In this particular dataset the effect of puromycin aminonucleoside (PAN) on in-vitro differentiated podocytes (cultured at 37 °C) is documented. The comparison in the paper and in the app is between differentiated (postmitotic) PAN-treated podocytes and vehicle-treated control podocytes. After PAN-treatment, 10 proteins decreased significantly, among these were YAP/TAZ target genes. LATS phosphorylation in response to PAN was not altered. This effect depends on concentration and duration of PAN-treatment, which suggests that PAN-induced damage in postmitotic podocytes inhibits YAP/TAZ function in vitro. <br>

<strong> Glom_PAN_rat	rGlom_PAN_d2	The effect of PAN on rat glomeruli in vivo vs control (day 2)	</strong> <br>
This dataset is third of four derived from PMID: 28400537 (https://pubmed.ncbi.nlm.nih.gov/28400537/). The study describes a proteomic analysis before induction of injury of podocytes and at specific time points after induction but before the development of proteinuria to define the cellular effects of puromycin aminonucleoside (PAN). In this particular dataset the effect on rat glomeruli after two days of in vivo puromycin aminonucleoside (PAN) administration is documented. The comparison in the paper and in the app is between the effect of PAN on rat glomeruli and vehicle-treated rat glomeruli. Increased activity of Yap1 and expression of YAP/TAZ target genes as well as diminished activity of Lats were localized in PAN-treated rat podocytes. <br>

<strong> Glom_PAN_rat	rGlom_PAN_d4	The effect of PAN on rat glomeruli in vivo vs control (day 4)	</strong> <br>
This dataset is the fourth of four derived from PMID: 28400537 (https://pubmed.ncbi.nlm.nih.gov/28400537/). The study describes a proteomic analysis before induction of injury of podocytes and at specific time points after induction but before the development of proteinuria to define the cellular effects of puromycin aminonucleoside (PAN). In this particular dataset the effect on rat glomeruli after four days of in vivo puromycin aminonucleoside (PAN) administration is documented. The comparison in the paper and in the app is between the effect of PAN on rat glomeruli and vehicle-treated rat glomeruli. Increased activity of Yap1 and expression of YAP/TAZ target genes as well as diminished activity of Lats were localized in PAN-treated rat podocytes. Also, the glomerular proteome was dominated by serum proteins (i.e. albumin and complement). <br>

<strong> Human_FSGS_1_FSGS_2_vs_con2	hGlom_FSGS_1	The proteome of FSGS single glomeruli in human </strong> <br>
This dataset is the first of four derived from PMID: 29530281 (https://pubmed.ncbi.nlm.nih.gov/29530281/). The study describes a scalable method for direct extraction of quantitative proteomic information from individual units of the kidney (i.e. individual glomeruli and individual tubule segments). In this particular dataset Wt1-mutation induced focal segmental glomerulosclerosis (FSGS) is analyzed by mass spectroscopy. The comparison in the paper and in the app is between a single glomeruli of a Wt1 heterozygous mice which develop proteinuria and partially sclerotic lesions at the age of 14 weeks resembling FSGS and wild-type mice (3 mice per group, 7 single glomeruli/mouse). Albumin, laminin, collagen, LAMP1 and laminin are significantly increased in single glomeruli from Wt1het (FSGS) mice. ACTN1/4 and CD2AP are decreased. A strong correlation in the comparison occurred between LAMP1 and extracellular matrix proteins and glomerular albumin. <br>

<strong> Human_FSGS_1_FSGS_2_vs_con2	hGlom_FSGS_2	The proteome of FSGS single glomeruli in human	</strong> <br>
This dataset is the second of two derived from PMID: 29530281 (https://pubmed.ncbi.nlm.nih.gov/29530281/). The study describes a scalable method for direct extraction of quantitative proteomic information from individual units of the kidney (i.e. individual glomeruli and individual tubule segments). In this particular dataset Doxorubicin induced focal segmental glomerulosclerosis (FSGS) is analyzed by mass spectroscopy. The comparison in the app is between Doxorubicin induced FSGS mice glomeruli and wild-type mice. Albumin, laminin, collagen, and LAMP1 were significantly increased in single glomeruli Doxorubicin (FSGS) mice. Nephrin, ACTN1/4 and CD2AP were decreased. A strong correlation in the comparison occurred between LAMP1 and extracellular matrix proteins and glomerular albumin. <br>

<strong> Human_NPHS_pat1vscon	hGlom_NephrinMutation_1	The proteome of nephrin mutation in human patient 1	</strong> <br>
This dataset is the third of four derived from PMID: 29530281 (https://pubmed.ncbi.nlm.nih.gov/29530281/). The study describes a scalable method for direct extraction of quantitative proteomic information from individual units of the kidney (i.e. individual glomeruli and individual tubule segments). In this particular dataset nephrin mutation in a human patient (patient 1), is analyzed by mass spectroscopy. The comparison in the paper and in the app is between a single glomeruli from a patient with congenital nephrotic syndrome (NPHS1-patient 1) and minimal change disease caused by nephrin (NPHS1) mutation and glomeruli from a healthy control kidney. In glomeruli obtained from the kidney with NPHS1 mutation (patient 1), Nephrin was strongly reduced and LAMP1 and collagen were increased. <br>

<strong> Human_NPHS_pat2vscon	hGlom_NephrinMutation_2	The proteome of nephrin mutation in human patient 2 </strong> <br>
This dataset is the fourth of four derived from PMID: 29530281 (https://pubmed.ncbi.nlm.nih.gov/29530281/). The study describes a scalable method for direct extraction of quantitative proteomic information from individual units of the kidney (i.e. individual glomeruli and individual tubule segments). In this particular dataset nephrin mutation in a human patient (patient 2), is analyzed by mass spectroscopy. The comparison in the paper and in the app is between a single glomeruli from a patient with congenital nephrotic syndrome (NPHS1-patient 2) and minimal change disease caused by nephrin (NPHS1) mutation and  glomeruli from a healthy control kidney. In glomeruli obtained from the kidney with NPHS1 mutation, Nephrin was strongly reduced and LAMP1 was increased. <br>

<strong> NATIVE_Podocyte_Doxorubicine_vs_control_v2_mouse	mPodo_Doxorubicin	The proteome of Doxorubicin treated native podocytes vs control </strong> <br>
This dataset is the first of two derived from PMID: 32047005 (https://pubmed.ncbi.nlm.nih.gov/32047005/). The study describes a proteomic approach to identify participants in the early and late disease responses of podocytes. In this particular dataset the proteome of Doxorubicin treated native podocytes is analyzed by mass spectroscopy. The comparison in the paper and in the app is between isolated podocytes from R26mTmG mice, which were mated with hNphs2.PodCre mice and treated with Doxorubicin and vehicle-treated controls. Several glomerular disease–associated proteins, such as the transcription factors Stat-1 and Stat-3, were significantly increased. Also, proteins involved in metabolic processes were upregulated and adhesion process proteins were downregulated. <br>

<strong> NATIVE_Podocyte_LPS_vs_Control_mouse	mPodo_LPS	The proteome of LPS treated native podocytes vs control  </strong> <br>
This dataset is the second of two derived from PMID: 32047005 (https://pubmed.ncbi.nlm.nih.gov/32047005/). The study describes a proteomic approach to identify participants in the early and late disease responses of podocytes. In this particular dataset the proteome of LPS treated native podocytes is analyzed by mass spectroscopy. The comparison in the app is between isolated podocytes from Podocin.2A.iCre.2A.mTomato mice, which were injected with LPS and vehicle-treated controls. The transcription factor Stat-1 as well as proteins in essential metabolic pathways were significantly increased. Adhesion process proteins, associated with the actin-cytoskeleton, were downregulated. <br>

### In addition in the organoid version :

<strong> Cultured_kidney_organoids_day25_vs_day21	hOrganoid_d25_d21 The proteome of kidney organoids differentiated to 25 days compared with 21 days (in publication) </strong> <br>
Project description: Kidney organoids are a valuable and innovative model to understand genetic diseases, kidney development and transcriptomic dynamics. However, their proteome has not been analyzed so far. Here, we analyzed the organoid proteome trajectory during differentiation. Genes involved in podocytopathies and cystic kidney diseases were abundantly expressed on protein level, distinguishing organoids from almost every available cell culture model. On their pathway to terminal differentiation, organoids developed increased deposition of extracellular matrix. Single cell transcriptomic analysis suggests that most changes locate to podocytes and early podocyte progenitors. This matrix deposition is different from commonly used animal models of glomerular disease. We grew organoids from two independent batches according to the Freedman protocol, and performed proteomic profiling (Freedman, Brooks et al. 2015, Czerniecki, Cruz et al. 2018). The IPSCs were differentiated for a three-week period until first spheroids formed. From day 21 of the culture they were used in our experiments up until day 29, where off-target differentiation of organoids becomes an issue.This dataset compares day 25 of organoid differentiation with day 21 of organoid differentiation. <br>

<strong> Cultured_kidney_organoids_day29_vs_day21	hOrganoid_d29_d21 The proteome of kidney organoids differentiated to 29 days compared with 21 days (in publication) </strong> <br>
Project description: Kidney organoids are a valuable and innovative model to understand genetic diseases, kidney development and transcriptomic dynamics. However, their proteome has not been analyzed so far. Here, we analyzed the organoid proteome trajectory during differentiation. Genes involved in podocytopathies and cystic kidney diseases were abundantly expressed on protein level, distinguishing organoids from almost every available cell culture model. On their pathway to terminal differentiation, organoids developed increased deposition of extracellular matrix. Single cell transcriptomic analysis suggests that most changes locate to podocytes and early podocyte progenitors. This matrix deposition is different from commonly used animal models of glomerular disease. We grew organoids from two independent batches according to the Freedman protocol, and performed proteomic profiling (Freedman, Brooks et al. 2015, Czerniecki, Cruz et al. 2018). The IPSCs were differentiated for a three-week period until first spheroids formed. From day 21 of the culture they were used in our experiments up until day 29, where off-target differentiation of organoids becomes an issue.This dataset compares day 29 of organoid differentiation with day 21 of organoid differentiation.<br>

<strong> Cultured_kidney_organoids_TNFa_vs_VC_24h	hOrganoid_24h_TNFa The proteome of kidney organoids treated with TNFa for 24 hours (in publication) </strong> <br>
Kidney organoids are a valuable and innovative model to understand genetic diseases, kidney development and transcriptomic dynamics. However, their proteome has not been analyzed so far. Here, we analyzed the organoid proteome after treatment of organoids with 5ng/mL TNFalpha for 24h and 48h compared with vehicle control (VC). Incubation of organoids (day 25 of differentiation) with TNFalpha led to an activation of NFkappaB signaling, and, interestingly, secretion of cytokines and complement components, alongside with extracellular matrix components. Interestingly, this signaling system directly links inflammatory signaling, production of cytokines and complement; and production of extracellular matrix. Thus, we provide a repository of kidney organoid proteins that revealed the potential to model pathophysiological pathways beyond genetic diseases. Organoids were grown according to the Freedman protocol (Freedman, Brooks et al. 2015, Czerniecki, Cruz et al. 2018). The IPSCs were differentiated for a three-week period until first spheroids from. We started TNFa stimulation at day 25, with the 24h stimulation ending on day 26 and the 48h stimulation ending on day 27.  We chose day 25 because it lies centrally in the day 21 to day 29 window, where we observe reproducible spheroids with limited off-target differentiation of organoids, which becomes an issue after day 29. <br>

<strong> Cultured_kidney_organoids_TNFa_vs_VC_48h	hOrganoid_48h_TNFa The proteome of kidney organoids treated with TNFa for 48 hours (in publication) </strong> <br>
Kidney organoids are a valuable and innovative model to understand genetic diseases, kidney development and transcriptomic dynamics. However, their proteome has not been analyzed so far. Here, we analyzed the organoid proteome after treatment of organoids with 5ng/mL TNFalpha for 24h and 48h compared with vehicle control (VC). Incubation of organoids (day 25 of differentiation) with TNFalpha led to an activation of NFkappaB signaling, and, interestingly, secretion of cytokines and complement components, alongside with extracellular matrix components. Interestingly, this signaling system directly links inflammatory signaling, production of cytokines and complement; and production of extracellular matrix. Thus, we provide a repository of kidney organoid proteins that revealed the potential to model pathophysiological pathways beyond genetic diseases. Organoids were grown according to the Freedman protocol (Freedman, Brooks et al. 2015, Czerniecki, Cruz et al. 2018). The IPSCs were differentiated for a three-week period until first spheroids from. We started TNFa stimulation at day 25, with the 24h stimulation ending on day 26 and the 48h stimulation ending on day 27.  We chose day 25 because it lies centrally in the day 21 to day 29 window, where we observe reproducible spheroids with limited off-target differentiation of organoids, which becomes an issue after day 29.

## Add your own dataset
We have made it possible to add a new dataset. However, to do that, one will have to use the app locally, that it to say, one have to download all the files and the database. <br>
The new dataset must be in either <strong> .txt </strong> or <strong> .csv </strong> format. It <strong> MUST </strong> contain : 
- A column named "Gene_names" which contains the gene name
- A column which contains the log2 fold-changes
- A column which contains the -log10 p-values associated to the t-test 

The new dataset <strong> CAN </strong> contain : 
- A column which contains log2 fold-changes of another test
- A column which contains -log10 p-values of the t-test on the second column of fold-changes <br>

The function [add_table_function.R](add_table_function.R) is to be used to add the dataset into the database.<br>
One must enter the following argument in that exact order : 
- Path to the database
- Path to the new dataset
- Name of the dataset to be used in the database
- Decimal character used in the dataset
- Separator character used in the dataset
- Name of the log2 fold-change column
- Name of the -log10 p-value colum
- If there is one, name of the second log2 fold-change column and then name of the second -log10 p-value column

Note : By default, if there is a second test, it will create two tables in the database, one named as passed in the arguments, the second one with the same name following by "2". If one wish to change that, one can follow the procedure under "# Rename a table" in the function file.

