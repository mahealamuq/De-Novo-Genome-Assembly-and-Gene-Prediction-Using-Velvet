# De Novo Genome Assembly and Gene Prediction Using Velvet

## Project Overview

This repository documents a de novo genome assembly workflow using paired-end Illumina sequencing reads from accession **ERR048385**. The aim of this project is to assemble an unknown bacterial genome without using a reference genome and evaluate the quality of the assembly using Velvet.

This practical is based on genome assembly and gene prediction, where the main objectives are to:

1. Assemble a genome de novo using high-throughput sequencing data.
2. Assess how different k-mer lengths affect assembly quality.
3. Use the best assembly for downstream gene prediction and genome annotation.

The original practical uses **Velvet**, a de Bruijn graph-based assembler, to assemble Illumina paired-end reads. The practical also explains that the assembled contigs can later be used for ab initio protein-coding gene prediction using GeneMark. :contentReference[oaicite:0]{index=0}

---

## Repository Structure

```text
de-novo-genome-assembly-velvet/
│
├── README.md
│
├── raw_data/
│   ├── ERR048385_1.fastq
│   └── ERR048385_2.fastq
│
├── fastqc/
│   ├── ERR048385_1_fastqc.html
│   └── ERR048385_2_fastqc.html
│
├── velvet_assembly/
│   ├── Assembly59/
│   │   ├── contigs.fa
│   │   ├── stats.txt
│   │   └── Log
│   │
│   └── Assembly69/
│       ├── contigs.fa
│       ├── stats.txt
│       └── Log
│
├── results/
│   ├── assembly_comparison_table.md
│   └── best_assembly_summary.md
│
└── figures/
    └── fastq_quality_encoding.png
```

---

## Dataset

The dataset used in this project is:

```text
ERR048385_1.fastq
ERR048385_2.fastq
```

These are paired-end Illumina sequencing reads. The practical notes state that the sequence reads are two FASTQ files and are used for de novo genome assembly with Velvet. :contentReference[oaicite:1]{index=1}

---

## Downloading the FASTQ Files

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

## Step 1: Inspect the FASTQ Files

### Count the number of reads

```bash
grep "^@ERR048385" raw_data/ERR048385_1.fastq | wc -l
grep "^@ERR048385" raw_data/ERR048385_2.fastq | wc -l
```

Result:

```text
ERR048385_1.fastq = 4,363,304 reads
ERR048385_2.fastq = 4,363,304 reads
```

### Check whether reads are paired-end

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

## Step 2: Check Read Length

```bash
sed -n 2p raw_data/ERR048385_1.fastq | wc -L
```

Result:

```text
75
```

Therefore, each read is **75 bp** long.

---

## Step 3: Identify FASTQ Quality Encoding

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

## Step 4: Quality Control Using FastQC

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

### FastQC Interpretation

FastQC showed that the reads had high overall quality. Although some tile-level variation was observed, the reads were considered suitable for downstream assembly. Therefore, no trimming was performed.

---

## Step 5: Why Velvet Was Used

Velvet is a de novo genome assembler designed for short-read sequencing data. It is suitable for assembling bacterial genomes from Illumina reads.

Velvet works in two main steps:

1. **velveth**  
   Reads sequence data, stores k-mers, and creates a hash table.

2. **velvetg**  
   Builds and simplifies the de Bruijn graph and generates assembled contigs.

The practical notes describe these two Velvet steps and explain that `velveth` produces `Sequences` and `Roadmaps`, while `velvetg` produces the final assembly outputs including contigs, statistics, and a log file. :contentReference[oaicite:2]{index=2}

---

## Step 6: Install Velvet with Larger K-mer Support

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

## Step 7: Choosing K-mer Lengths

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

These are unsuitable because they are longer than the read length. The practical asks students to compare different k-mer lengths and explains that k is specified when running `velveth`. :contentReference[oaicite:3]{index=3}

---

## Step 8: Velvet Assembly Using k = 59

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

### Explanation of flags

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

The practical states that the ERR048385 reads were generated from fragments of approximately 500 bp with expected coverage of 20×. :contentReference[oaicite:4]{index=4}

---

## Step 9: Velvet Assembly Using k = 69

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

## Step 10: Extract Assembly Statistics

Move into the assembly directory:

```bash
cd velvet_assembly/Assembly69
```

### Count number of contigs

```bash
grep -c "^>" contigs.fa
```

### Extract N50, maximum contig length, total contig length, and assembled reads

```bash
grep "n50" Log
grep "max" Log
grep "total" Log
```

### Calculate maximum coverage

```bash
awk 'NR>1 {print $6}' stats.txt | sort -n | tail -1
```

Or in R:

```r
A <- read.table("stats.txt", header=TRUE)
max(A$short1_cov)
```

### Calculate average coverage

```bash
awk 'NR>1 {sum+=$6; n++} END {print sum/n}' stats.txt
```

Or in R:

```r
mean(A$short1_cov)
```

The practical specifically asks students to use the `contigs.fa`, `stats.txt`, and `Log` files to calculate N50, maximum contig length, total contig length, assembled reads, number of contigs, and coverage statistics. :contentReference[oaicite:5]{index=5}

---

## Step 11: Assembly Results

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

## Step 12: Calculation of Percentage Reads Assembled

### k = 59

```text
Percentage assembled = (7,701,965 / 8,726,608) × 100
                     = 88.26%
```

### k = 69

```text
Percentage assembled = (7,237,537 / 8,726,608) × 100
                     = 82.94%
```

---

## Step 13: Coverage Explanation

Coverage means the average number of sequencing reads supporting each base in a contig.

For example:

```text
Coverage = total sequenced bases / genome length
```

A coverage of 69.60× means that each base in the contigs is supported by approximately 70 reads on average.

### Coverage summary for k = 59

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

## Step 14: Assembly Evaluation

### k = 59 Assembly

Advantages:

- Higher percentage of reads assembled: 88.26%
- Slightly larger total contig length: 2,904,625 bp
- Higher average coverage: 69.60×

Disadvantages:

- Lower N50: 39,640 bp
- Shorter maximum contig: 135,167 bp
- More contigs: 458
- More fragmented assembly

### k = 69 Assembly

Advantages:

- Higher N50: 115,016 bp
- Longer maximum contig: 305,186 bp
- Fewer contigs: 93
- More continuous assembly

Disadvantages:

- Lower percentage of reads assembled: 82.94%
- Lower average coverage: 18.03×

---

## Step 15: Which Assembly Is Better?

The **k = 69 assembly** was selected as the best assembly.

Although the k = 59 assembly assembled more reads, it produced many more contigs and a lower N50. This indicates that the assembly is more fragmented.

The k = 69 assembly produced:

- Higher N50
- Longer maximum contig
- Fewer contigs

These are important indicators of a more continuous and better-quality genome assembly.

The practical notes explain that shorter k-mers can allow more overlaps and produce many short contigs, while longer k-mers can reduce assembled reads. Therefore, the best assembly requires a balance between N50 and the percentage of assembled reads. :contentReference[oaicite:6]{index=6}

---

## Step 16: Best Assembly Summary

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

## Step 17: Files Produced by Velvet

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

The practical identifies these three files as the key outputs to inspect after assembly. :contentReference[oaicite:7]{index=7}

---

## Step 18: Future Work

The best assembly can be used for:

- Gene prediction using GeneMark
- Genome annotation
- BLAST searches
- Phylogenomic analysis
- Identification of genes potentially acquired by lateral genetic transfer

The practical states that assembled contigs from the best assembly should be used as the draft genome assembly for gene prediction. :contentReference[oaicite:8]{index=8}

---

## Conclusion

This project successfully assembled paired-end Illumina reads using Velvet and compared the effect of two k-mer lengths, k = 59 and k = 69. The k = 69 assembly was selected as the best assembly because it produced a higher N50, a longer maximum contig, and fewer contigs. These results suggest that k = 69 generated a more continuous and less fragmented draft genome assembly, making it more suitable for downstream genome annotation and phylogenomic analysis..
