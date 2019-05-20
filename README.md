# GPD
Gene-to-Protein-to-Disease (GPD): a segmentation-based approach for association analysis with whole exome sequencing data

# 1 Preparation
GPD can be downloaded and installed in R as follows. As a prerequisite, devtools must be installed:

```{r}
install.packages("devtools")
```
Next, install GPD:

```{r}
library("devtools")
devtools::install_github("ginnyintifa/GPD")
library(GPD)
```

Files need to be prepared are 

* mutation file
* protein information file 
* clinical information file 

Format of mutaiton file:

```
Hugo_Symbol	Gene	Chromosome	Start_Position	End_Position	Variant_Classification	Variant_Type	HGVSc	HGVSp	barcode
OPN4	ENSG00000122375	10	88419681	88419681	Missense_Mutation	SNP	c.863G>A	p.Gly288Asp	TCGA-OR-A5J1
GOLGA3	ENSG00000090615	12	133360652	133360652	Intron	SNP	c.3267+98G>A	.	TCGA-OR-A5J1
KLRB1	ENSG00000111796	12	9760409	9760409	Missense_Mutation	SNP	c.27G>C	p.Glu9Asp	TCGA-OR-A5J1
SALL2	ENSG00000165821	14	21991730	21991730	Missense_Mutation	SNP	c.2132G>A	p.Arg711Gln	TCGA-OR-A5J1
PLCB2	ENSG00000137841	15	40587141	40587141	Silent	SNP	c.1902G>A	p.%3D	TCGA-OR-A5J1
C15orf27	ENSG00000169758	15	76467904	76467904	Frame_Shift_Del	DEL	c.657delC	p.Tyr219Ter	TCGA-OR-A5J1
KLHDC4	ENSG00000104731	16	87742934	87742934	Missense_Mutation	SNP	c.1384G>A	p.Asp462Asn	TCGA-OR-A5J1
SCARF1	ENSG00000074660	17	1543262	1543262	Silent	SNP	c.1083G>A	p.%3D	TCGA-OR-A5J1
NOL11	ENSG00000130935	17	65734087	65734087	Missense_Mutation	SNP	c.1528A>C	p.Ser510Arg	TCGA-OR-A5J1
LAMA3	ENSG00000053747	18	21461949	21461949	Missense_Mutation	SNP	c.5162G>A	p.Gly1721Asp	TCGA-OR-A5J1
WDR7	ENSG00000091157	18	54348602	54348602	Missense_Mutation	SNP	c.325T>C	p.Cys109Arg	TCGA-OR-A5J1
THOP1	ENSG00000172009	19	2790453	2790453	Silent	SNP	c.51G>A	p.%3D	TCGA-OR-A5J1
ANKRD24	ENSG00000089847	19	4210319	4210319	Silent	SNP	c.1009C>T	p.%3D	TCGA-OR-A5J1
LYPD3	ENSG00000124466	19	43968541	43968541	Missense_Mutation	SNP	c.147G>T	p.Met49Ile	TCGA-OR-A5J1
LYPD3	ENSG00000124466	19	43968550	43968550	Silent	SNP	c.138G>A	p.%3D	TCGA-OR-A5J1
PTCHD2	ENSG00000204624	1	11561526	11561526	Silent	SNP	c.477G>A	p.%3D	TCGA-OR-A5J1
VPS13D	ENSG00000048707	1	12309384	12309384	Silent	SNP	c.552T>G	p.%3D	TCGA-OR-A5J1
PHC2	ENSG00000134686	1	33820015	33820015	Silent	SNP	c.1542G>A	p.%3D	TCGA-OR-A5J1
FOXD2	ENSG00000186564	1	47905517	47905517	3'UTR	SNP	c.*222C>T	.	TCGA-OR-A5J1
RBL1	ENSG00000080839	20	35696518	35696518	Missense_Mutation	SNP	c.362G>A	p.Arg121His	TCGA-OR-A5J1
TTN	ENSG00000155657	2	179501161	179501161	Missense_Mutation	SNP	c.41293A>G	p.Lys13765Glu	TCGA-OR-A5J1
```



Format of protein information file 


```
uniprot_accession	start_position	end_position	center_position	unit_name	gene_name	gene_id	unit_label
Q14D04	211	221	216	py	VEPH1	ENSG00000197415	PTM
Q14D04	375	385	380	ps	VEPH1	ENSG00000197415	PTM
Q14D04	392	402	397	pt	VEPH1	ENSG00000197415	PTM
Q14D04	417	427	422	pt	VEPH1	ENSG00000197415	PTM
Q14D04	425	435	430	ps	VEPH1	ENSG00000197415	PTM
Q14D04	444	454	449	ps	VEPH1	ENSG00000197415	PTM
Q14D04	524	534	529	ps	VEPH1	ENSG00000197415	PTM
Q14D04	717	819	768	PH	VEPH1	ENSG00000197415	Domain
Q14D04	766	776	771	ps	VEPH1	ENSG00000197415	PTM
Q14D04	778	788	783	ps	VEPH1	ENSG00000197415	PTM
Q14D33	50	150	100	zf-3CxxC	RTP5	ENSG00000277949	Domain

```

Please note that users' own files for mutation data and protein information data should follow the above formats with exact columns and column names. 

Format of clinical information file 

```
bcr_patient_barcode	type	age_at_initial_pathologic_diagnosis	gender	race	ajcc_pathologic_tumor_stage	clinical_stage	histological_type	histological_grade	initial_pathologic_dx_year	menopause_status	birth_days_to	vital_status	tumor_status	last_contact_days_to	death_days_to	cause_of_death	new_tumor_event_type	new_tumor_event_site	new_tumor_event_site_other	new_tumor_event_dx_days_to	treatment_outcome_first_course	margin_status	residual_tumor	OS	OS.time	DSS	DSS.time	DFI	DFI.time	PFI	PFI.time	Redaction
TCGA-ACC-1	ACC	65	FEMALE	WHITE	Stage II	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2011	[Not Available]	-24017	Alive	WITH TUMOR	383	NA	[Not Available]	Distant Metastasis	Lung	#N/A	166	[Not Available]	NA	NA	0	383	0	383	NA	NA	1	166	NA
TCGA-ACC-2	ACC	42	FEMALE	WHITE	Stage IV	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	1998	[Not Available]	-15536	Dead	WITH TUMOR	NA	436	[Not Available]	Distant Metastasis	Bone	#N/A	61	Progressive Disease	NA	NA	1	436	1	436	NA	NA	1	61	NA
TCGA-ACC-3	ACC	32	FEMALE	WHITE	Stage III	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2010	[Not Available]	-11970	Dead	WITH TUMOR	NA	994	[Not Available]	Distant Metastasis	Lung	#N/A	97	Progressive Disease	NA	NA	1	994	1	994	NA	NA	1	97	NA
TCGA-ACC-4	ACC	37	FEMALE	[Not Evaluated]	Stage II	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2009	[Not Available]	-13574	Alive	TUMOR FREE	1857	NA	[Not Available]	#N/A	#N/A	#N/A	NA	Complete Remission/Response	NA	NA	0	1857	0	1857	0	1857	0	1857	NA
TCGA-ACC-5	ACC	53	FEMALE	WHITE	Stage IV	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2011	[Not Available]	-19492	Alive	WITH TUMOR	1171	NA	[Not Available]	Distant Metastasis	Lung	#N/A	351	Stable Disease	NA	NA	0	1171	0	1171	NA	NA	1	351	NA


```
Clinical information is supplied for survival analysis. Not all the columns in the above file are necessary, however these columns are required:

```
barcode 
age
gender
race
OS
OS.time
``` 
*`OS` refers to the survival status; `OS.time` refers to survival time.* 


We provide in the R package the mutation file and clinical data for a subset of TCGA ACC cohort and protein information file. Users can access these data after installing the package. 
```
sel_acc_mutation
sel_acc_cdr
ptm_pfam_combine
```

# 2 Mutation Extraction 

GPD is originally designed for cancer specific analysis extracting information from a Pan-cancer file, it processes a designated list of patients (barcodes) in a certain cancer type at a time. A list of barcodes should be provided, it can be from the matching clinical file or simply reading from the mutation file in the following way:

```{r}

mutation_df = read.table("your_path_to_file/user_mutation.tsv",
                header = T)
cancer_barcode = unique(mutation_df$barcode)

```
These objects will be inputs in the `extraction_annotation_pos`funciton. 

An example of running a subset of ACC somatic mutations from TCGA Pan-Cancer dataset is as the following:

```{r}




extraction_annotation_pos(mutation_df = sel_acc_mutation,
                                  cancer_type = "ACC",
                                  cancer_barcode = unique(sel_acc_mutation$barcode),
                                  output_dir = "your_output_dir1/")

```


This process will produce four output files:

* "your_output_dir1/ACC_mutation.tsv"
* "your_output_dir1/ACC_mutation_pc.tsv"
* "your_output_dir1/ACC_mutation_npc.tsv"
* "your_output_dir1/ACC_mutation_pc_pos.tsv"

The last two will be used in the subsequent function. 

# 3 Mutation Mapping 

GPD maps mutations to protein information (piu) data. Users can load their own piu data in the following way:

```{r}
piu_df = read.table("your_path_to_file/user_piu.tsv",
                header = T)

```



Mapping to the provided piu data with the annotated subset of ACC somatic mutations:

```{r}
piu_mapping (piu_df =  ptm_pfam_combine,
             pc_data_name  = "your_output_dir1/ACC_mutation_pc_pos.tsv",
             npc_data_name = "your_output_dir1/ACC_mutation_npc.tsv",
             cancer_barcode = unique(sel_acc_mutation$barcode),
             output_dir = "your_output_dir2/")

```

This function will produce three output files:

* "your_output_dir2/piu_mapping.tsv"
* "your_output_dir2/lu_summarising_count.tsv"
* "your_output_dir2/ncu_summarising_count.tsv"


These files present mutation mapping results which can be used in subsequent statistical analysis. 

# 4 Survival analysis

After mapping, users can perform survival analysis to study the association between mapped mutation counts on PIU, LU and NCU and patient survival.

Users can upload their clinical data in the following way:

```{r}
clinical_df = read.table("your_path_to_file/user_clinical.tsv",
                header = T)

```
With the subset of ACC data mapping results and clinical data, we perform the survival analysis.

```{r}

univariate_cox_model(piu_filename = "your_output_dir2/piu_mapping_count.tsv",
                     lu_filename = "your_output_dir2/lu_summarising_count.tsv",
                     ncu_filename = "your_output_dir2/ncu_summarising_count.tsv",
                     clinical_df = sel_acc_cdr,
                     gender_as_covariate = T,
                     race_group_min = 6,
                     min_surv_days = 90,
                     min_surv_people = 5,
                     patient_sum_min = 3,
                     mutation_type = "somatic",
                     output_dir = "your_output_dir3/")


```


This function will produce these output files:

* "your_output_dir3/somatic_piu_cdr_univariate.tsv"
* "your_output_dir3/somatic_lu_cdr_univariate.tsv"
* "your_output_dir3/somatic_ncu_cdr_univariate.tsv"


