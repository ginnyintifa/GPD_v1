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

Files need to be prepared are 1 mutation file; 2 protein information file 

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
Please note that Users' own mutation file should follow this format with exact columns and column names. 


# 2 Mutation Extraction 

GPD is designed for cancer specific analysis extracting information from a Pan-cancer file, it processes a designated list of patients (barcodes) in a certain cancer type at a time. A list of barcodes should be provided, it can be from the matching clinical file or simple read from the mutation file in the following way:
```{r}

acc_mut = fread("/Users/ginny/Google Drive/R_GPD/GPD_package_0401/package_files/acc_mutation.tsv",
                stringsAsFactors = F)
acc_barcode = acc_mut$barcode
```



An example of running all ACC sample is as the following:
```{r}




extraction_annotation_pos(mutation_file = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/GPD/package_files/acc_mutation.tsv",
                                  cancer_type = "ACC",
                                  cancer_barcode = acc_barcode,
                                  output_dir = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/")

```


This process will produce four output files:
1 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/ACC_mutation.tsv"
2 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/ACC_pc.tsv"
3 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/ACC_npc.tsv"
4 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/ACC_pc_pos.tsv"
The last two will be used in the subsequent function. 

# 3 Mutation Mapping 

```{r}
piu_mapping (piu_filename =  "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/GPD/package_files/ptm_pfam_combine.tsv",
             pc_data_name  = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/ACC_pc_pos.tsv",
             npc_data_name = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/ACC_npc.tsv",
             cancer_barcode = acc_barcode,
             output_dir = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/")

```

This function will produce three output files:

1 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/piu_mapping.tsv"
2 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/bpiu_summarising_count.tsv"
3 "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/npc_summarising_count.tsv"


These files present mutation mapping results which can be used in subsequent statistical analysis. 







