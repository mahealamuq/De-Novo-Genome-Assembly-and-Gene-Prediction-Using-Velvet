# Bacterial Genome Assembly, Gene Prediction, and Phylogenomic Analysis

 ## 1.Project Overview

This repository demonstrates a complete bacterial genome analysis workflow using Illumina paired-end sequencing data (ERR048385). The project begins with **de novo genome assembly** using **Velvet**, followed by assembly evaluation and selection of the optimal assembly based on quality metrics. The best assembly is then used for **ab initio gene prediction** with **Prodigal**, functional annotation using **BLAST**, and phylogenomic analysis based on the **16S rRNA gene** and **AccA protein**.

The workflow covers the major stages of bacterial genome analysis, including sequence quality assessment, genome assembly, assembly evaluation, gene prediction, homolog identification, multiple sequence alignment, and Bayesian phylogenetic inference. The final phylogenetic trees were generated using **MrBayes** and visualised with **MEGA 12** to investigate the evolutionary relationships between the assembled genome and closely related bacterial species.

This repository was developed for  Bioinformatics practical exercises and demonstrates a complete end-to-end bioinformatics pipeline from raw sequencing reads to biological interpretation.

---
## 2. Workflow

```text
Illumina FASTQ
      │
      ▼
FastQC
      │
      ▼
Velvet Assembly
      │
      ▼
Assembly Evaluation
      │
      ▼
Prodigal Gene Prediction
      │
      ▼
BLAST Annotation
      │
      ▼
Extract 16S rRNA & AccA
      │
      ▼
NCBI BLAST
      │
      ▼
MUSCLE Alignment
      │
      ▼
Readseq
      │
      ▼
MrBayes
      │
      ▼
MEGA Tree Visualization
```

## 3.Project Objectives

The objectives of this project are to:

1. Assess the quality of Illumina paired-end sequencing reads.
2. Assemble a bacterial genome de novo using Velvet.
3. Compare genome assemblies generated with different k-mer lengths.
4. Evaluate assembly quality using metrics such as N50, contig length, genome size, coverage, and percentage of assembled reads.
5. Select the optimal genome assembly for downstream analyses.
6. Predict protein-coding genes using Prodigal.
7. Calculate genome statistics, including GC content and the number of predicted genes.
8. Annotate predicted genes using BLAST searches against the NCBI database.
9. Identify homologous 16S rRNA genes and AccA proteins for phylogenetic analysis.
10. Perform multiple sequence alignment using MUSCLE.
11. Construct Bayesian phylogenetic trees using MrBayes.
12. Visualise and interpret phylogenetic trees using MEGA 12.
13. Investigate the evolutionary relationship between the assembled genome and closely related *Staphylococcus* species.
---



## 4.Dataset

The dataset used in this project is:

```text
ERR048385_1.fastq
ERR048385_2.fastq
```

These are paired-end Illumina sequencing reads. The practical notes state that the sequence reads are two FASTQ files and are used for de novo genome assembly with Velvet. 

---

## 5. Software Requirements


| Software | Version |
|----------|---------|
| Ubuntu | 24.04 LTS |
| Velvet | 1.2.10 |
| Prodigal | 2.6.3 |
| BLAST+ | 2.17 |
| MUSCLE | 5 |
| Readseq | 2.1 |
| MrBayes | 3.2 |
| MEGA | 12 |
| SeqKit | 2.x |
| R | 4.6 |

## 6. Genome Assembly

### 6.1 Downloading the FASTQ Files

The accession **ERR048385** can be downloaded from the European Nucleotide Archive.

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR048/ERR048385/ERR048385_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR048/ERR048385/ERR048385_2.fastq.gz
```

Unzip the files:

```bash
gunzip ERR048385_1.fastq.gz
gunzip ERR048385_2.fastq.gz
```

Move files into a raw data folder:

```bash
mkdir raw_data
mv ERR048385_1.fastq ERR048385_2.fastq raw_data/
```

---

### 6.2 Inspect the FASTQ Files

**Count the number of reads**

```bash
grep "^@ERR048385" raw_data/ERR048385_1.fastq | wc -l
grep "^@ERR048385" raw_data/ERR048385_2.fastq | wc -l
```

Result:

```text
ERR048385_1.fastq = 4,363,304 reads
ERR048385_2.fastq = 4,363,304 reads
```

**Check whether reads are paired-end**

```bash
head -n 8 raw_data/ERR048385_1.fastq
head -n 8 raw_data/ERR048385_2.fastq
```

The presence of matching `/1` and `/2` read identifiers confirms that the files contain paired-end reads.

Example:

```text
@ERR048385.1 ... /1
@ERR048385.1 ... /2
```

Therefore, the reads are **paired-end reads**.

---

**Check Read Length**

```bash
sed -n 2p raw_data/ERR048385_1.fastq | wc -L
```

Result:

```text
75
```

Therefore, each read is **75 bp** long.

---

**Identify FASTQ Quality Encoding**

The FASTQ quality string was compared with the standard FASTQ quality encoding chart.

The quality characters matched the **Illumina 1.8+ Phred+33** range.

Example quality characters:

```text
?, D, E, I
```

Using Phred+33 encoding:

```text
? = ASCII 63 → 63 - 33 = Q30
D = ASCII 68 → 68 - 33 = Q35
E = ASCII 69 → 69 - 33 = Q36
I = ASCII 73 → 73 - 33 = Q40
```

These scores fall within the Illumina 1.8+ range of approximately Q0 to Q41.

**Conclusion:**  
The reads were generated using **Illumina 1.8+ chemistry with Phred+33 encoding**.

---

### 6.3 Quality Control Using FastQC

Before assembly, read quality was checked using FastQC.

```bash
mkdir fastqc

fastqc raw_data/ERR048385_1.fastq -o fastqc/
fastqc raw_data/ERR048385_2.fastq -o fastqc/
```

Open the FastQC reports:

```bash
firefox fastqc/ERR048385_1_fastqc.html
firefox fastqc/ERR048385_2_fastqc.html
```

FastQC checks include:

- Per-base sequence quality
- Per-sequence quality scores
- Adapter contamination
- GC content
- Sequence duplication levels

**FastQC Interpretation**

FastQC showed that the reads had high overall quality. Although some tile-level variation was observed, the reads were considered suitable for downstream assembly. Therefore, no trimming was performed.

---

### 6.4 Velvet

**Why Velvet Was Used**

Velvet is a de novo genome assembler designed for short-read sequencing data. It is suitable for assembling bacterial genomes from Illumina reads.

Velvet works in two main steps:

1. **velveth**  
   Reads sequence data, stores k-mers, and creates a hash table.

2. **velvetg**  
   Builds and simplifies the de Bruijn graph and generates assembled contigs.

`velveth` produces `Sequences` and `Roadmaps`, while `velvetg` produces the final assembly outputs including contigs, statistics, and a log file. 

---

**Install Velvet with Larger K-mer Support**

The default Velvet installation may only support k-mers up to 31. Since this project uses k = 59 and k = 69, Velvet must be compiled with a larger `MAXKMERLENGTH`.

```bash
wget https://github.com/dzerbino/velvet/archive/refs/heads/master.zip
unzip master.zip
cd velvet-master
make 'MAXKMERLENGTH=109'
```

Check Velvet:

```bash
./velveth
```

---

**Choosing K-mer Lengths**

The read length is 75 bp, so the k-mer length must be smaller than 75.

Suitable k-mer values:

```text
k = 59
k = 69
```

Unsuitable values:

```text
k = 79
k = 109
```

These are unsuitable because they are longer than the read length. Compare different k-mer lengths and explains that k is specified when running `velveth`. 

---

**Velvet Assembly Using k = 59**

Run `velveth`:

```bash
./velveth ../velvet_assembly/Assembly59 59 \
-shortPaired \
-fastq \
-separate \
../raw_data/ERR048385_1.fastq \
../raw_data/ERR048385_2.fastq
```

Run `velvetg`:

```bash
./velvetg ../velvet_assembly/Assembly59 \
-ins_length 500 \
-exp_cov 20
```

**Explanation of flags**

| Flag / Parameter | Meaning |
|---|---|
| `Assembly59` | Output directory for k = 59 assembly |
| `59` | K-mer length |
| `-shortPaired` | Input reads are short paired-end reads |
| `-fastq` | Input files are in FASTQ format |
| `-separate` | Forward and reverse reads are in separate files |
| `ERR048385_1.fastq` | Forward reads |
| `ERR048385_2.fastq` | Reverse reads |
| `-ins_length 500` | Expected insert size is approximately 500 bp |
| `-exp_cov 20` | Expected sequencing coverage is 20× |

The ERR048385 reads were generated from fragments of approximately 500 bp with expected coverage of 20×.

---

**Velvet Assembly Using k = 69**

Run `velveth`:

```bash
./velveth ../velvet_assembly/Assembly69 69 \
-shortPaired \
-fastq \
-separate \
../raw_data/ERR048385_1.fastq \
../raw_data/ERR048385_2.fastq
```

Run `velvetg`:

```bash
./velvetg ../velvet_assembly/Assembly69 \
-ins_length 500 \
-exp_cov 20
```

---

**Extract Assembly Statistics**

Move into the assembly directory:

```bash
cd velvet_assembly/Assembly69
```

**Count number of contigs**

```bash
grep -c "^>" contigs.fa
```

**Extract N50, maximum contig length, total contig length, and assembled reads**

```bash
grep "n50" Log
grep "max" Log
grep "total" Log
```

**Calculate maximum coverage**

```bash
awk 'NR>1 {print $6}' stats.txt | sort -n | tail -1
```

Or in R:

```r
A <- read.table("stats.txt", header=TRUE)
max(A$short1_cov)
```

**Calculate average coverage**

```bash
awk 'NR>1 {sum+=$6; n++} END {print sum/n}' stats.txt
```

Or in R:

```r
mean(A$short1_cov)
```

Use the `contigs.fa`, `stats.txt`, and `Log` files to calculate N50, maximum contig length, total contig length, assembled reads, number of contigs, and coverage statistics.

---

**Assembly Results**

| Metric | k = 59 | k = 69 |
|---|---:|---:|
| N50 (bp) | 39,640 | 115,016 |
| Maximum Contig Length (bp) | 135,167 | 305,186 |
| Total Contig Length (bp) | 2,904,625 | 2,872,901 |
| Reads Assembled | 7,701,965 | 7,237,537 |
| Total Reads | 8,726,608 | 8,726,608 |
| Percentage of Reads Assembled (%) | 88.26 | 82.94 |
| Number of Contigs | 458 | 93 |
| Maximum Coverage (×) | 1070.31 | 426.77 |
| Average Coverage (×) | 69.60 | 18.03 |

---

**Calculation of Percentage Reads Assembled**

- k = 59

```text
Percentage assembled = (7,701,965 / 8,726,608) × 100
                     = 88.26%
```

- k = 69

```text
Percentage assembled = (7,237,537 / 8,726,608) × 100
                     = 82.94%
```

---

**Coverage Explanation**

Coverage means the average number of sequencing reads supporting each base in a contig.

For example:

```text
Coverage = total sequenced bases / genome length
```

A coverage of 69.60× means that each base in the contigs is supported by approximately 70 reads on average.

**Coverage summary for k = 59**

```r
A <- read.table("stats.txt", header=TRUE)
summary(A$short1_cov)
```

Output:

```text
Min.  1st Qu.  Median  Mean  3rd Qu.  Max.
1.00   39.32    44.71  69.60  89.23   1070.31
```

| Statistic | Value | Meaning |
|---|---:|---|
| Minimum | 1.00× | Lowest coverage contig |
| 1st Quartile | 39.32× | 25% of contigs have coverage below 39.32× |
| Median | 44.71× | Half of contigs have coverage below 44.71× |
| Mean | 69.60× | Average coverage across all contigs |
| 3rd Quartile | 89.23× | 75% of contigs have coverage below 89.23× |
| Maximum | 1070.31× | Highest coverage contig |

The maximum coverage is much higher than the mean. This may occur because repetitive regions or multi-copy genes, such as rRNA operons, can attract many reads and produce unusually high coverage.

---

### 6.5 Assembly Evaluation

**k = 59 Assembly**

Advantages:

- Higher percentage of reads assembled: 88.26%
- Slightly larger total contig length: 2,904,625 bp
- Higher average coverage: 69.60×

Disadvantages:

- Lower N50: 39,640 bp
- Shorter maximum contig: 135,167 bp
- More contigs: 458
- More fragmented assembly

**k = 69 Assembly**

Advantages:

- Higher N50: 115,016 bp
- Longer maximum contig: 305,186 bp
- Fewer contigs: 93
- More continuous assembly

Disadvantages:

- Lower percentage of reads assembled: 82.94%
- Lower average coverage: 18.03×

---

**Which Assembly Is Better?**

The **k = 69 assembly** was selected as the best assembly.

Although the k = 59 assembly assembled more reads, it produced many more contigs and a lower N50. This indicates that the assembly is more fragmented.

The k = 69 assembly produced:

- Higher N50
- Longer maximum contig
- Fewer contigs

These are important indicators of a more continuous and better-quality genome assembly.

The shorter k-mers can allow more overlaps and produce many short contigs, while longer k-mers can reduce assembled reads. Therefore, the best assembly requires a balance between N50 and the percentage of assembled reads.

---

**Best Assembly Summary**

The best assembly is:

```text
Assembly69
```

| Question | Answer |
|---|---:|
| Total number of contigs | 93 |
| Maximum average coverage for a contig | 426.77× |
| Total length of all contigs | 2,872,901 bp |
| Average coverage across all contigs | 18.03× |

Commands:

```bash
cd velvet_assembly/Assembly69

grep -c "^>" contigs.fa
awk 'NR>1 {print $6}' stats.txt | sort -n | tail -1
grep "total" Log
awk 'NR>1 {sum+=$6; n++} END {print sum/n}' stats.txt
```

---

### 6.6 Files Produced by Velvet

Each Velvet assembly directory contains:

```text
Graph2
LastGraph
Log
PreGraph
Roadmaps
Sequences
contigs.fa
stats.txt
```

Important files:

| File | Description |
|---|---|
| `contigs.fa` | Final assembled contigs in FASTA format |
| `stats.txt` | Contig length and coverage statistics |
| `Log` | Assembly summary including N50, max contig length, total length, and reads assembled |

These three files as the key outputs to inspect after assembly.

---

This project successfully assembled paired-end Illumina reads using Velvet and compared the effect of two k-mer lengths, k = 59 and k = 69. The k = 69 assembly was selected as the best assembly because it produced a higher N50, a longer maximum contig, and fewer contigs. These results suggest that k = 69 generated a more continuous and less fragmented draft genome assembly, making it more suitable for downstream genome annotation and phylogenomic analysis..

---


## 7. Gene Prediction from contigs Using Prodigal

### 7.1 Prodigal

The prodigal was used to identify protein-coding genes in assembled genome contigs. Prodigal is specifically designed for bacterial and archaeal genomes and is widely used for accurate prediction of protein-coding genes from assembled contigs. It does not require a pre-trained species model and performs gene prediction  directly from the genomic sequence.


**Create Gene Prediction Directory**

Create a directory for gene prediction analysis and move the assembled contigs into the directory.

```bash
mkdir Gen_prediction

mv Assembly69/contigs.fa Gen_prediction/

cd Gen_prediction
```

Rename the assembly file:

```bash
mv contigs.fa ERR048385_assembled_contigs.fa
```

---

**Install Prodigal**

Create a Conda environment and install Prodigal from Bioconda.

```bash
conda create -n prodigal_env -c bioconda -c conda-forge prodigal

conda activate prodigal_env
```

Verify installation:

```bash
prodigal -v
```

---

**Predict Protein-Coding Genes**

Run Prodigal on the assembled genome contigs.

```bash
prodigal \
-i ERR048385_assembled_contigs.fa \
-o genes.gff \
-a proteins.faa \
-d genes.fna
```

**Output Files**

| File           | Description                         |
| -------------- | ----------------------------------- |
| `genes.gff`    | Gene coordinates and annotations    |
| `proteins.faa` | Predicted protein sequences         |
| `genes.fna`    | Predicted gene nucleotide sequences |

---

### 7.2 Gene Prediction Results

**G+C Content of the Genome Contigs**

The GC content was calculated directly from the assembled contig sequences.

```bash
grep -v "^>" ERR048385_assembled_contigs.fa | \
awk '{gc+=gsub(/[GCgc]/,""); total+=length($0)}
END{print gc/total*100}'
```

- Result

```text
48.5055%
```

- Interpretation

The assembled genome has a GC content of approximately **48.51%**. GC content is a useful genomic characteristic that can provide insights into genome composition and taxonomic identity.

---

**Total Number of Predicted Genes**

Count the number of predicted genes in the nucleotide FASTA file.

```bash
grep -c "^>" genes.fna
```

- Result

```text
2704 genes
```

- Interpretation**

Prodigal predicted **2,704 protein-coding genes** within the assembled genome.

---

**Putative Function of the First Predicted Gene**

The first predicted gene sequence was extracted from `genes.fna` and searched against the NCBI non-redundant nucleotide database using BLASTN.

- Top BLAST Hit

| Feature   | Result                            |
| --------- | --------------------------------- |
| Protein   | Site-specific integrase (partial) |
| Organism  | *Staphylococcus aureus*           |
| Accession | WP_310651144.1                    |


### 7.3 Functional Annotation

The first predicted gene was identified as a **site-specific integrase**.

Site-specific integrases are recombinase enzymes that catalyse:

* DNA integration
* DNA excision
* Site-specific recombination

These proteins are commonly associated with:

* Bacteriophages
* Transposons
* Genomic islands
* Horizontal gene transfer

Integrases play an important role in genome evolution and the movement of mobile genetic elements between bacterial genomes.

---

**Does the First Predicted Gene Support a Staphylococcus Genome?**



Yes, The top BLASTN hit for the first predicted gene was a site-specific integrase from *Staphylococcus aureus*. This strongly supports the hypothesis that the assembled genome is derived from a member of the genus *Staphylococcus*.

Additional evidence supporting this conclusion includes:

* Genome size of approximately **2.87 Mb**
* Presence of a Staphylococcus-associated integrase gene
* Gene content consistent with a bacterial genome

Therefore, the gene prediction and BLAST analysis support the conclusion that the assembled genome is likely a **Staphylococcus species**.

---

**Summary of Results**

| Question                            | Result                  |
| ----------------------------------- | ----------------------- |
| GC content                          | 48.51%                  |
| Total predicted genes               | 2,704                   |
| First predicted gene                | Site-specific integrase |
| Top BLAST organism                  | *Staphylococcus aureus* |
| Supports Staphylococcus hypothesis? | Yes                     |

---

# 8: Phylogenomics


This project performs phylogenomic analysis of a de novo assembled bacterial genome. The main aim is to identify homologous sequences for two molecular markers and infer evolutionary relationships using Bayesian phylogenetic analysis.

The two markers used were:

* **16S rRNA gene**
* **AccA protein**, acetyl-CoA carboxylase subunit alpha






**Prepare Working Directory**

```bash
mkdir Phylogenomics
cd Phylogenomics
```

Copy the files from the gene prediction directory:

```bash
cp ../Gen_prediction/ERR048385_assembled_contigs.fa .
cp ../Gen_prediction/genes.fna .
cp ../Gen_prediction/proteins.faa .
```

The three FASTA files used were:

| File                             | Description                                       |
| -------------------------------- | ------------------------------------------------- |
| `ERR048385_assembled_contigs.fa` | Genome contigs from the best assembly             |
| `genes.fna`                      | Predicted gene nucleotide sequences from Prodigal |
| `proteins.faa`                   | Predicted protein sequences from Prodigal         |

---

## 8.1. Create Local BLAST Databases

```bash
makeblastdb -in ERR048385_assembled_contigs.fa -dbtype nucl
makeblastdb -in genes.fna -dbtype nucl
makeblastdb -in proteins.faa -dbtype prot
```

These commands create searchable BLAST databases from the genome contigs, predicted genes, and predicted proteins.

---

## 8.2 NCBI BLAST Search
**Search for 16S rRNA in Genome Contigs**

```bash
esearch -db nucleotide -query "Staphylococcus aureus 16S ribosomal RNA" | efetch -format fasta > 16S.fa
```

```bash
blastn \
-query 16S.fa \
-db ERR048385_assembled_contigs.fa \
-outfmt 6 \
-out results/16S_vs_contigs.txt
```

The BLAST outputs were saved using `-outfmt 6`, which gives tabular output without headers. The columns name for **16S_vs_contigs.txt** file:

| Column Number | Column Name | Meaning                      |
| ------------: | ----------- | ---------------------------- |
|             1 | `qseqid`    | Query sequence ID            |
|             2 | `sseqid`    | Subject/database sequence ID |
|             3 | `pident`    | Percentage identity          |
|             4 | `length`    | Alignment length             |
|             5 | `mismatch`  | Number of mismatches         |
|             6 | `gapopen`   | Number of gap openings       |
|             7 | `qstart`    | Query start position         |
|             8 | `qend`      | Query end position           |
|             9 | `sstart`    | Subject start position       |
|            10 | `send`      | Subject end position         |
|            11 | `evalue`    | E-value                      |
|            12 | `bitscore`  | Bit score                    |


- Result

Finding Best Significant hits from **16S_vs_contigs.txt**

```R
#Read the table
blast <- read.table("16S_vs_contigs.txt",header = FALSE,sep = "\t",stringsAsFactors = FALSE)

# Add BLAST outfmt 12 column names
colnames(blast) <- c( "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Filter significant 16S hits
sig_hits <- subset(blast, evalue == 0 & pident >= 97 & length >= 1000)
best_hits <- sig_hits[order(-sig_hits$bitscore), ]

# Show best hits
head(best_hits, )
```


Significant 16S rRNA hits were found in the assembled genome contigs.



The most useful 16S rRNA hit was:

| Metric              | Result                              |
| ------------------- | ----------------------------------- |
| Contig              | `NODE_76_length_1557_cov_77.104691` |
| Alignment length    | approximately 1476 bp               |
| Percentage identity | approximately 100%                  |
| E-value             | 0.0                                 |
| Strand              | Reverse strand                      |

The reverse strand was identified because the subject start coordinate was larger than the subject end coordinate.

---


**Search for 16S rRNA in Predicted Genes**

```bash
blastn \
-query 16S.fa \
-db genes.fna \
-outfmt 6 \
-out results/16S_vs_genes.txt
```

**Result**

No reliable significant 16S rRNA hit was found in the predicted gene set.

This result is expected because Prodigal predicts protein-coding genes, while 16S rRNA is a non-coding RNA gene. Therefore, 16S rRNA is present in the genome contigs but is not expected to appear in the predicted protein-coding gene file.

---

**Search for AccA Protein in Predicted Proteins**

- Download Staphylococcus aureus AccA protein from NCBI

```bash
esearch -db protein -query "Staphylococcus aureus acetyl-CoA carboxylase subunit alpha" | efetch -format fasta > AccA.fa
```

```bash
blastp \
-query AccA.fa \
-db proteins.faa \
-outfmt 6 \
-out results/AccA_vs_proteins.txt
```

- Result

A significant AccA hit was found in the predicted protein set.

| Metric            | Result                                            |
| ----------------- | ------------------------------------------------- |
| Predicted protein | `NODE_25_length_256936_cov_16.179901_45`          |
| Function          | Acetyl-CoA carboxylase biotin carboxylase subunit |
| Query coverage    | Full-length protein match                         |
| E-value           | 0.0                                               |
| Interpretation    | Strong AccA homolog                               |

---

** Did I find hits in the three searches?**

Yes, hits were found in the genome contigs and predicted proteins, but no meaningful 16S rRNA hit was found among the predicted genes.

| Search                             | Result          | Explanation                                                                      |
| ---------------------------------- | --------------- | -------------------------------------------------------------------------------- |
| 16S rRNA vs genome contigs         | Hit found       | The genome assembly contains a 16S rRNA region                                   |
| 16S rRNA vs predicted genes        | No reliable hit | 16S rRNA is non-coding and is not predicted by Prodigal as a protein-coding gene |
| AccA protein vs predicted proteins | Hit found       | AccA is a protein-coding housekeeping gene                                       |

A significant 16S rRNA hit was found on `NODE_76_length_1557_cov_77.104691`, indicating that the assembled genome contains a 16S rRNA gene. However, no significant 16S hit was found among the predicted genes because 16S rRNA is not a protein-coding gene. A significant AccA protein hit was found among the predicted proteins, confirming that the assembled genome encodes AccA.

---

**Extract Putative 16S rRNA and AccA Sequences**

- Extract 16S rRNA Gene

The best 16S rRNA hit was on `NODE_76`, approximately from position 77 to 1552 on the reverse strand.

```bash
seqkit grep -p "NODE_76_length_1557_cov_77.104691" ERR048385_assembled_contigs.fa > NODE_76.fa
seqkit subseq -r 77:1552 NODE_76.fa | seqkit seq -t DNA -r -p > my_16S.fa
```

- Check:

```bash
head my_16S.fa
```

---

- Extract AccA Protein

```bash
seqkit grep -p "NODE_25_length_256936_cov_16.179901_45" proteins.faa > my_AccA.faa
```

- Check:

```bash
head my_AccA.faa
```

---

**Download Homologous 16S rRNA Sequences from NCBI**

The selected taxa were:

| Taxon                          | Role                     |
| ------------------------------ | ------------------------ |
| *Staphylococcus aureus*        | Close match              |
| *Staphylococcus epidermidis*   | Related *Staphylococcus* |
| *Staphylococcus haemolyticus*  | Related *Staphylococcus* |
| *Staphylococcus saprophyticus* | Related *Staphylococcus* |
| *Staphylococcus lugdunensis*   | Related *Staphylococcus* |
| *Bacillus subtilis*            | Outgroup                 |

- Download one sequence for each taxon:

```bash
esearch -db nucleotide -query "Staphylococcus aureus[Organism] 16S ribosomal RNA" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Saureus_16S.fa

esearch -db nucleotide -query "Staphylococcus epidermidis[Organism] 16S ribosomal RNA" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Sepidermidis_16S.fa

esearch -db nucleotide -query "Staphylococcus haemolyticus[Organism] 16S ribosomal RNA" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Shaemolyticus_16S.fa

esearch -db nucleotide -query "Staphylococcus saprophyticus[Organism] 16S ribosomal RNA" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Ssaprophyticus_16S.fa

esearch -db nucleotide -query "Staphylococcus lugdunensis[Organism] 16S ribosomal RNA" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Slugdunensis_16S.fa

esearch -db nucleotide -query "Bacillus subtilis[Organism] 16S ribosomal RNA" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Bacillus_16S.fa
```

- Combine with the assembled genome 16S sequence:

```bash
cat my_16S.fa \
Saureus_16S.fa \
Sepidermidis_16S.fa \
Shaemolyticus_16S.fa \
Ssaprophyticus_16S.fa \
Slugdunensis_16S.fa \
Bacillus_16S.fa \
> 16S_gene_set.fa
```

- Check:

```bash
grep -c "^>" 16S_gene_set.fa
```

- Expected result:

```text
7
```

---

**Download Homologous AccA Protein Sequences from NCBI**

```bash
esearch -db protein -query "Staphylococcus aureus[Organism] acetyl-CoA carboxylase biotin carboxylase subunit" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Saureus_AccA.faa

esearch -db protein -query "Staphylococcus epidermidis[Organism] acetyl-CoA carboxylase biotin carboxylase subunit" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Sepidermidis_AccA.faa

esearch -db protein -query "Staphylococcus haemolyticus[Organism] acetyl-CoA carboxylase biotin carboxylase subunit" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Shaemolyticus_AccA.faa

esearch -db protein -query "Staphylococcus saprophyticus[Organism] acetyl-CoA carboxylase biotin carboxylase subunit" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Ssaprophyticus_AccA.faa

esearch -db protein -query "Staphylococcus lugdunensis[Organism] acetyl-CoA carboxylase biotin carboxylase subunit" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Slugdunensis_AccA.faa

esearch -db protein -query "Bacillus subtilis[Organism] acetyl-CoA carboxylase biotin carboxylase subunit" \
| efetch -format fasta | awk '/^>/{if(++n>1)exit} {print}' > Bacillus_AccA.faa
```

- Combine:

```bash
cat my_AccA.faa \
Saureus_AccA.faa \
Sepidermidis_AccA.faa \
Shaemolyticus_AccA.faa \
Ssaprophyticus_AccA.faa \
Slugdunensis_AccA.faa \
Bacillus_AccA.faa \
> AccA_protein_set.faa
```

- Check:

```bash
grep -c "^>" AccA_protein_set.faa
```

- Expected result:

```text
7
```

---

**Which BLAST search yielded greater diversity of source organisms?**

The 16S rRNA BLAST search yielded greater diversity of source organisms.

The 16S rRNA gene is present in almost all bacteria and is highly conserved, so it can identify homologous sequences across a broad range of bacterial taxa. In contrast, the AccA protein search mainly returned close *Staphylococcus* matches because AccA is a conserved housekeeping protein and the query sequence was very similar to *Staphylococcus aureus* AccA.

Therefore, 16S rRNA was better for broad taxonomic diversity, while AccA was more useful for comparing closely related species.

---

**Based on NCBI BLAST Results**

-  How long is each query?

| Query         |                Length |
| ------------- | --------------------: |
| 16S rRNA gene | approximately 1476 bp |
| AccA protein  |       453 amino acids |

- What is the top hit organism?

| Query        | Top hit organism        |
| ------------ | ----------------------- |
| 16S rRNA     | *Staphylococcus aureus* |
| AccA protein | *Staphylococcus aureus* |

- Do the top hits belong to the same species?

Yes. Both the 16S rRNA and AccA top hits belong to *Staphylococcus aureus*.

- Do the top hits belong to the same strain?

Not necessarily. The top hits belong to the same species, *Staphylococcus aureus*, but they do not necessarily correspond to the same strain. The 16S rRNA and AccA results may match different reference strains because they come from different databases and different sequence types.

---

## 8.3 MUSCLE Multiple Sequence Alignment

**Rename Sequence Headers**

Short names were used so that tree labels are readable.

- Example for AccA:

```bash
sed -i 's/>NODE_25_length_256936_cov_16.179901_45.*/>My_AccA/' AccA_protein_set.faa
sed -i 's/>WP_490684228.1.*/>Saureus/' AccA_protein_set.faa
sed -i 's/>YFI15416.1.*/>Sepidermidis/' AccA_protein_set.faa
sed -i 's/>YEU50698.1.*/>Shaemolyticus/' AccA_protein_set.faa
sed -i 's/>WP_135092699.1.*/>Ssaprophyticus/' AccA_protein_set.faa
sed -i 's/>WP_490828678.1.*/>Slugdunensis/' AccA_protein_set.faa
sed -i 's/>WP_327837060.1.*/>Bacillus/' AccA_protein_set.faa
```

- Check:

```bash
grep "^>" AccA_protein_set.faa
```

---

**Multiple Sequence Alignment with MUSCLE**

- Align 16S rRNA

```bash
muscle -align 16S_gene_set.fa -output 16S_gene_set.aln
```

- Align AccA proteins

```bash
muscle -align AccA_protein_set.faa -output AccA_protein_set.aln
```

---
## 8.4 Readseq Format Conversion

**Convert Alignments to NEXUS Format**

```bash
readseq -format=nexus -output=16S_gene_set.nex 16S_gene_set.aln
readseq -format=nexus -output=AccA_protein_set.nex AccA_protein_set.aln
```

---

**How many sequence format options are available in readseq?**


This depends on the installed version of `readseq`. To check, run:

```bash
readseq -h
```

The help menu lists the available input and output sequence formats. Commonly supported formats include FASTA, NEXUS, PHYLIP, GenBank, EMBL, PIR, and others.

---
## 8.5 Bayesian Phylogenetic Analysis Using MrBayes
**Compare MrBayes command blocks for protein and DNA data**

- DNA MrBayes Block

```text
begin mrbayes;
    set autoclose=yes;
    mcmc ngen=1500000 nchains=4 samplefreq=100 burnin=5000;
    lset nucmodel=4by4 rates=gamma ngammacat=4;
    sump burnin=5000;
    sumt burnin=5000 conformat=simple contype=allcompat;
end;
```

Append this block to end of 16S_gene_set.nex file

- Protein MrBayes Block

```text
begin mrbayes;
    set autoclose=yes;
    mcmc ngen=500000 nchains=4 samplefreq=100 burnin=2000;
    prset aamodel=mixed;
    lset rates=gamma ngammacat=4;
    sump burnin=2000;
    sumt burnin=2000 conformat=simple contype=allcompat;
end;
```
Append this block to end of AccA_protein_set.nex file


The obvious difference is the substitution model. The DNA analysis uses a nucleotide model through `lset nucmodel=4by4`, whereas the protein analysis uses an amino acid model through `prset aamodel=mixed`.

- Both analyses use:

* MCMC sampling
* four chains
* gamma-distributed rate variation
* burn-in removal
* posterior probability summary
* consensus tree generation

However, DNA and protein sequences require different evolutionary models because nucleotides and amino acids evolve differently.

---

**Run MrBayes**

- Run MrBayes:

```bash
mb
```

- Inside MrBayes:

```text
execute 16S_gene_set.nex
execute AccA_protein_set.nex
```

---

**MrBayes Output Files**

| Question                                               | File type  | Meaning                         |
| ------------------------------------------------------ | ---------- | ------------------------------- |
| (a) Log likelihood for each sampled tree               | `.p`       | Parameter and likelihood values |
| (b) Phylogenetic tree in each sampled chain            | `.t`       | Sampled trees from the MCMC run |
| (c) All possible trees sorted by posterior probability | `.trprobs` | Tree probabilities              |
| (d) Consensus final tree                               | `.con.tre` | Final consensus tree            |

---

**What parameters need to be changed?**

The practical asks for:

* 1,000,000 MCMC generations
* burn-in after 500,000 generations
* sampling every 100 generations

Calculation:

```text
Total samples = 1,000,000 / 100 = 10,000
Burn-in samples = 500,000 / 100 = 5,000
```

Therefore, the correct burn-in value is:

```text
burnin = 5000
```

The parameters to change are:

```text
ngen=1000000
samplefreq=100
burnin=5000
```

Correct command:

```text
mcmc ngen=1000000 nchains=4 samplefreq=100 burnin=5000;
sump burnin=5000;
sumt burnin=5000 conformat=simple contype=allcompat;
```

---

## 8.6 Tree Visualisation Using MEGA 12

The MrBayes consensus trees were converted into clean Newick format for viewing in MEGA 12:

| Tree              | File                   |
| ----------------- | ---------------------- |
| 16S rRNA tree     | `16S_clean_named.nwk`  |
| AccA protein tree | `AccA_clean_named.nwk` |

The trees were rooted using **Bacillus** as the outgroup.

---

**Do the two trees share the same topology?**


No. The 16S rRNA tree and AccA protein tree do not have exactly the same topology.

Both trees group the *Staphylococcus* species together and separate **Bacillus** as the outgroup, but the internal branching relationships among the *Staphylococcus* species differ.

In the 16S rRNA tree, the assembled 16S sequence clusters closely with *Staphylococcus aureus*. In the AccA protein tree, the predicted AccA sequence clusters more closely with *Staphylococcus epidermidis*.

Therefore, the two trees are similar at the genus level but differ in their detailed topology.

---

**Interpretation of Bayesian posterior probabilities**

Bayesian posterior probability indicates support for each internal node. Values close to 1.0 represent strong support, while lower values indicate weaker confidence.

- 16S rRNA tree

| Posterior probability | Interpretation      |
| --------------------: | ------------------- |
|                 0.992 | Very strong support |
|                 0.776 | Moderate support    |
|                 0.743 | Moderate support    |
|                 0.434 | Weak support        |

The 16S tree strongly supports the close relationship between the assembled sequence and *Staphylococcus aureus*. However, some deeper relationships have weaker support.

- AccA protein tree

| Posterior probability | Interpretation      |
| --------------------: | ------------------- |
|                 0.998 | Very strong support |
|                 0.998 | Very strong support |
|                 0.993 | Very strong support |
|                 0.984 | Very strong support |

The AccA tree has stronger support overall, indicating that the AccA alignment provides clearer phylogenetic signal among the selected taxa.

---

**Which tree has longer branch lengths?**

The AccA protein tree has longer branch lengths than the 16S rRNA tree.

This indicates that AccA has accumulated more substitutions and has evolved faster than 16S rRNA. This is expected because 16S rRNA is highly conserved due to its essential role in ribosome function, while protein-coding genes can evolve more rapidly.

---

**If the tree is rooted using your sequence as the outgroup, would interpretation change?**

Yes. Rooting the tree using the assembled genome sequence would change the interpretation because the root determines the direction of evolutionary relationships.

The assembled sequence is part of the *Staphylococcus* group, so it is not an appropriate outgroup. **Bacillus** is a better outgroup because it is outside the genus *Staphylococcus* and provides a more biologically meaningful root.

---

## 9.Results

The complete bioinformatics workflow successfully assembled the bacterial genome, predicted protein-coding genes, annotated homologous sequences, and inferred phylogenetic relationships.

### Genome Assembly

| Metric                        | Best Assembly (k = 69) |
| ----------------------------- | ---------------------: |
| Number of Contigs             |                     93 |
| N50                           |             115,016 bp |
| Maximum Contig Length         |             305,186 bp |
| Total Assembly Length         |           2,872,901 bp |
| Reads Assembled               |  7,237,537 / 8,726,608 |
| Percentage of Reads Assembled |                 82.94% |
| Maximum Coverage              |                426.77× |
| Average Coverage              |                 18.03× |

### Gene Prediction

| Metric                         | Result                  |
| ------------------------------ | ----------------------- |
| Gene Prediction Tool           | Prodigal                |
| GC Content                     | 48.51%                  |
| Predicted Protein-Coding Genes | 2,704                   |
| First Predicted Gene           | Site-specific integrase |
| Closest BLAST Match            | *Staphylococcus aureus* |

### Phylogenomics

* A putative **16S rRNA gene** and **AccA protein** were successfully identified from the assembled genome.
* Homologous sequences were downloaded from the NCBI RefSeq database.
* Multiple sequence alignments were generated using MUSCLE.
* Bayesian phylogenetic trees were inferred using MrBayes.
* Both phylogenetic analyses placed the assembled genome within the genus *Staphylococcus*.
* The 16S rRNA tree showed the closest relationship to *Staphylococcus aureus*, while the AccA tree produced a slightly different internal topology with stronger posterior support.
* *Bacillus subtilis* was used as the outgroup for both phylogenetic trees.

---

## 10.Discussion

This project demonstrates a complete bacterial genome analysis pipeline, beginning with raw Illumina sequencing reads and ending with phylogenetic interpretation. Genome assembly using Velvet showed that different k-mer lengths substantially influenced assembly quality. Although the k = 59 assembly incorporated a higher percentage of reads, the k = 69 assembly produced fewer contigs, a much larger N50 value, and substantially longer contigs. Consequently, the k = 69 assembly was selected for downstream analyses because it provided a more contiguous draft genome.

Gene prediction using Prodigal identified 2,704 protein-coding genes with a genome GC content of approximately 48.5%, values that are consistent with members of the genus *Staphylococcus*. Functional annotation of the first predicted gene identified a site-specific integrase, a protein commonly associated with bacteriophages and other mobile genetic elements involved in horizontal gene transfer.

Local and remote BLAST analyses successfully identified homologous 16S rRNA genes and AccA proteins. Phylogenetic analyses based on these two molecular markers consistently placed the assembled genome within the genus *Staphylococcus*. Although the 16S rRNA and AccA trees showed slight differences in their internal branching patterns, both strongly supported the taxonomic placement of the assembled genome. The AccA protein tree displayed longer branch lengths and stronger posterior probabilities than the 16S rRNA tree, reflecting the faster evolutionary rate of protein-coding genes compared with highly conserved ribosomal RNA genes.

Overall, this project demonstrates how genome assembly, gene prediction, sequence annotation, and phylogenetic analysis can be integrated into a reproducible workflow for bacterial genome characterisation.

---

## 11.References

1. Zerbino DR, Birney E. **Velvet: Algorithms for de novo short read assembly using de Bruijn graphs.** *Genome Research.* 2008;18(5):821–829.

2. Hyatt D, et al. **Prodigal: Prokaryotic gene recognition and translation initiation site identification.** *BMC Bioinformatics.* 2010;11:119.

3. Camacho C, et al. **BLAST+: Architecture and applications.** *BMC Bioinformatics.* 2009;10:421.

4. Edgar RC. **MUSCLE: Multiple sequence alignment with high accuracy and high throughput.** *Nucleic Acids Research.* 2004;32(5):1792–1797.

5. Huelsenbeck JP, Ronquist F. **MRBAYES: Bayesian inference of phylogenetic trees.** *Bioinformatics.* 2001;17(8):754–755.

6. Tamura K, Stecher G, Kumar S. **MEGA12: Molecular Evolutionary Genetics Analysis.**

7. National Center for Biotechnology Information (NCBI). https://www.ncbi.nlm.nih.gov/

8. European Nucleotide Archive (ENA). https://www.ebi.ac.uk/ena/

9. BINF7001 Bioinformatics Practical Notes.






