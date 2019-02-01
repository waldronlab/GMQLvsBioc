# GMQLvsBioc
Alternative approaches to https://doi.org/10.1093/bioinformatics/bty688 analysis

[Quick link](https://github.com/waldronlab/GMQLvsBioc/files/2728079/genometric2018.pdf) to the manuscript PDF, which states:

> Also packages of R/
Bioconductor (https://www.bioconductor.org/) have been proposed
for tertiary analysis (Huber et al., 2015); they facilitate typical
specific operations, but require to perform them through
scripts and are not suitable for big data processing.

# Use case (as specified by authors)

The use case and its GMQL query are formulated as follows: “In
TCGA data of BRCA patients, find the DNA somatic mutations
within the first 2000 bp outside of the genes that are both expressed
with FPKM > 3 and have at least a methylation in the same patient
biospecimen, and extract these mutations of the top 5% patients
with the highest number of such mutations.”

EXPRESSED_GENE = SELECT(manually_curated__cases__
disease_type = “Breast Invasive Carcinoma”; region: fpkm >
3.0) GRCh38_TCGA_gene_expression;

METHYLATION = SELECT(manually_curated__cases__disease_
type = “Breast Invasive Carcinoma”) GRCh38_TCGA_
methylation;

MUTATION = SELECT(manually_curated__cases__disease_
type = “Breast Invasive Carcinoma”) GRCh38_TCGA_
somatic_mutation_masked;

GENE_METHYL = JOIN(distance < 0; output: left_distinct;
joinby: biospecimen__bio__bcr_sample_barcode) EXPRESSED
_GENEMETHYLATION;

MUTATION_GENE = JOIN(distance <. 2000, distance >.
0; output: left_distinct; joinby: biospecimen__bio__bcr_sample_
barcode) MUTATION GENE_METHYL;

MUTATION_GENE_count = EXTEND(mutation_count AS
COUNT()) MUTATION_GENE;

MUTATION_GENE_top = ORDER(mutation_count DESC;
meta_topp: 5) MUTATION_GENE_count;
MATERIALIZE MUTATION_GENE_top INTO MUTATION_
GENE_top;

# BigQuery approach

```sql
with tbl as (SELECT 
  m.Reference_Allele,
  m.Tumor_Seq_Allele1,
  m.Tumor_Seq_Allele2,
  m.Start_Position,
  m.End_Position,
  m.Chromosome,
  HTSeq__FPKM,
  CASE WHEN e.strand='+' THEN e.start-2000 ELSE e.end END as region_start,
  CASE WHEN e.strand='-' THEN e.end+2000 ELSE e.start END as region_end,
  e.seq_name,
  e.strand,
  e.gene_name, 
  e.gene_id,
  r.sample_barcode
FROM `isb-cgc.genome_reference.Ensembl_GRCh38_87` e 
JOIN `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression` r on r.Ensembl_gene_id=e.gene_id
JOIN `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation` m on r.sample_barcode=m.sample_barcode_tumor
WHERE r.HTSeq__FPKM>3 
  AND concat('chr',e.seq_name) = m.Chromosome
  AND m.Start_Position > (CASE WHEN e.strand='+' THEN e.start-2000 ELSE e.end END) 
  AND m.Start_Position < (CASE WHEN e.strand='-' THEN e.end+2000 ELSE e.start END))
select tbl.*, c.mutation_count 
from tbl 
join (select sample_barcode, count(*) as mutation_count from tbl group by sample_barcode) c
on c.sample_barcode=tbl.sample_barcode limit 100;
```
