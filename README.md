# GPD
Gene-to-Protein-to-Disease (GPD): a segmentation-based approach for statistical analysis of whole exome sequencing data


# 1 Preparation
GPD can be downloaded and installed in R. 

As a prerequisite, devtools must be installed:

```{r}
install.packages("devtools")
```
Next, install GPD:

```{r}
library("devtools")
devtools::install_github("ginnyintifa/GPD")
library(GPD)
```

Input files need to be prepared are 

* Mutation data file
* Protein information data file 
* Clinical information data file 

Format of mutation data:

```
Hugo_Symbol	Gene	Chromosome	Start_Position	End_Position	Variant_Classification	Variant_Type	HGVSc	HGVSp	Tumor_Sample_Barcode
GFRA1	ENSG00000151892	10	117884831	117884831	Missense_Mutation	SNP	c.671G>A	p.Arg224Gln	TCGA-ACC-35
LRP6	ENSG00000070018	12	12301919	12301919	Missense_Mutation	SNP	c.3163G>A	p.Glu1055Lys	TCGA-ACC-35
KIF5A	ENSG00000155980	12	57957244	57957244	Missense_Mutation	SNP	c.152G>A	p.Arg51His	TCGA-ACC-35
BBS10	ENSG00000179941	12	76741456	76741459	Frame_Shift_Del	DEL	c.306_309delCAGA	p.Asp102GlufsTer6	TCGA-ACC-35
SUGT1	ENSG00000165416	13	53237236	53237236	Nonsense_Mutation	SNP	c.484G>T	p.Glu162Ter	TCGA-ACC-35
HECTD1	ENSG00000092148	14	31598273	31598273	Missense_Mutation	SNP	c.4304T>A	p.Val1435Asp	TCGA-ACC-35
TLN2	ENSG00000171914	15	63017231	63017231	Silent	SNP	c.3183G>A	p.%3D	TCGA-ACC-35
PDE8A	ENSG00000073417	15	85664025	85664025	Intron	SNP	c.1735-3C>A	.	TCGA-ACC-35
ANKS3	ENSG00000168096	16	4747041	4747041	Missense_Mutation	SNP	c.1959G>C	p.Trp653Cys	TCGA-ACC-35
CCL11	ENSG00000172156	17	32614719	32614719	3'UTR	SNP	c.*10A>G	.	TCGA-ACC-35
TP53	ENSG00000141510	17	7578415	7578428	Frame_Shift_Del	DEL	c.502_515delCACATGACGGAGGT	p.His168CysfsTer8	TCGA-ACC-35
LAMA3	ENSG00000053747	18	21464756	21464756	Missense_Mutation	SNP	c.5242G>T	p.Val1748Leu	TCGA-ACC-35
ZNF230	ENSG00000159882	19	44515611	44515611	Missense_Mutation	SNP	c.1420A>T	p.Met474Leu	TCGA-ACC-35
IZUMO1	ENSG00000182264	19	49245556	49245556	Silent	SNP	c.510G>A	p.%3D	TCGA-ACC-35
PTCHD2	ENSG00000204624	1	11561189	11561189	Missense_Mutation	SNP	c.140G>A	p.Arg47Gln	TCGA-ACC-35
IGSF3	ENSG00000143061	1	117150771	117150771	Missense_Mutation	SNP	c.1015G>A	p.Glu339Lys	TCGA-ACC-35
NEK2	ENSG00000117650	1	211846871	211846871	Missense_Mutation	SNP	c.509C>T	p.Thr170Met	TCGA-ACC-35

```



Format of protein information data: 


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

Please note that user's own files of mutation data and protein information data should follow the above formats with exact column names.  Indeed, the format of an MAF file is compatible with the format required for the mutation data file, and columns present MAF file but not listed in the above example are not needed in GPD.  

Format of clinical information data:

```
Tumor_Sample_Barcode	type	age	gender	race	ajcc_pathologic_tumor_stage	clinical_stage	histological_type	histological_grade	initial_pathologic_dx_year	menopause_status	birth_days_to	vital_status	tumor_status	last_contact_days_to	death_days_to	cause_of_death	new_tumor_event_type	new_tumor_event_site	new_tumor_event_site_other	new_tumor_event_dx_days_to	treatment_outcome_first_course	margin_status	residual_tumor	OS	OS.time	DSS	DSS.time	DFI	DFI.time	PFI	PFI.time	Redaction
TCGA-ACC-1	ACC	65	FEMALE	WHITE	Stage II	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2011	[Not Available]	-24017	Alive	WITH TUMOR	383	NA	[Not Available]	Distant Metastasis	Lung	#N/A	166	[Not Available]	NA	NA	0	383	0	383	NA	NA	1	166	NA
TCGA-ACC-2	ACC	42	FEMALE	WHITE	Stage IV	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	1998	[Not Available]	-15536	Dead	WITH TUMOR	NA	436	[Not Available]	Distant Metastasis	Bone	#N/A	61	Progressive Disease	NA	NA	1	436	1	436	NA	NA	1	61	NA
TCGA-ACC-3	ACC	32	FEMALE	WHITE	Stage III	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2010	[Not Available]	-11970	Dead	WITH TUMOR	NA	994	[Not Available]	Distant Metastasis	Lung	#N/A	97	Progressive Disease	NA	NA	1	994	1	994	NA	NA	1	97	NA
TCGA-ACC-4	ACC	37	FEMALE	[Not Evaluated]	Stage II	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2009	[Not Available]	-13574	Alive	TUMOR FREE	1857	NA	[Not Available]	#N/A	#N/A	#N/A	NA	Complete Remission/Response	NA	NA	0	1857	0	1857	0	1857	0	1857	NA
TCGA-ACC-5	ACC	53	FEMALE	WHITE	Stage IV	[Not Applicable]	Adrenocortical carcinoma- Usual Type	[Not Available]	2011	[Not Available]	-19492	Alive	WITH TUMOR	1171	NA	[Not Available]	Distant Metastasis	Lung	#N/A	351	Stable Disease	NA	NA	0	1171	0	1171	NA	NA	1	351	NA


```
Clinical information is an input required for survival analysis. Not all the columns in the above file are required, however these columns are compulsory:

```
Tumor_Sample_Barcode 
age
gender
race
OS
OS.time
``` 
*`OS` refers to the survival status; `OS.time` refers to survival time.* 


We provide in our R package the mutation data file and clinical data file for a subset of TCGA ACC cohort, as well as the protein information data file holding PhosphoSitePlus PTM records and pfam domain records. User can access these data after installing the package. 

```
sel_acc_mutation
sel_acc_cdr
ptm_pfam_combine
```

# 2 Mutation Extraction 

GPD requires a list of interested barcodes to be provided, it can be from the matching clinical data file or simply reading from the mutation data file in the following way (if all samples are of interest):

```{r}

library(data.table)
mutation_df = fread("your_path_to_file/user_mutation.tsv",
                header = T,
                select = c("Hugo_Symbol","Gene","Chromosome","Start_Position","End_position","Variant_Classification","Variant_Type","HGVSc","HGVSp","Tumor_Sample_Barcode"),
                stringsAsFactors = F)
                
cancer_barcode = unique(mutation_df$Tumor_Sample_Barcode)  

```
These objects will be part of the inputs in the `extraction_annotation_pos`function. 

An example of running a subset of ACC somatic mutations is the following:

```{r}




extraction_annotation_pos(mutation_df = sel_acc_mutation,
                                  cancer_type = "ACC",
                                  cancer_barcode = unique(sel_acc_mutation$Tumor_Sample_Barcode),
                                  output_dir = "your_output_dir1/")

```


This will generate two output files that will be used in the subsequent function:


* "your_output_dir1/ACC_mutation_npc.tsv"
* "your_output_dir1/ACC_mutation_pc_pos.tsv"

# 3 Mutation Mapping 

GPD maps mutations to protein information  data. Users can load their own protein information data file in the following way:

```{r}
piu_df = fread("your_path_to_file/user_piu.tsv", stringsAsFactors = F)

```



Mapping to the provided protein information data with the annotated subset of ACC somatic mutations:

```{r}
piu_mapping (piu_df =  ptm_pfam_combine,
             pc_data_name  = "your_output_dir1/ACC_mutation_pc_pos.tsv",
             npc_data_name = "your_output_dir1/ACC_mutation_npc.tsv",
             cancer_barcode = unique(sel_acc_mutation$Tumor_Sample_Barcode),
             output_dir = "your_output_dir2/")

```

This function will produce three output files:

* "your_output_dir2/piu_mapping.tsv"
* "your_output_dir2/lu_summarising_count.tsv"
* "your_output_dir2/ncu_summarising_count.tsv"


These files carry mutation mapping results which can be used in subsequent statistical analysis. 

# 4 Survival analysis

After mapping, users can perform survival analysis to study the association between mapped mutation counts on PIU, LU, NCU and patient survival.

Users can upload their clinical data in the following way:

```{r}
clinical_df = fread("your_path_to_file/user_clinical.tsv", stringsAsFactors = F)

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


