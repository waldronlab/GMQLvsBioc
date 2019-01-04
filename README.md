# GMQLvsBioc
Alternative approaches to https://doi.org/10.1093/bioinformatics/bty688 analysis

[Quick link](https://github.com/waldronlab/GMQLvsBioc/files/2728079/genometric2018.pdf) to the manuscript.

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
