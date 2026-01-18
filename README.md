# Tools

Collection of bioinformatics tools for various NGS (Next-Generation Sequencing) data processing tasks.

## Table of Contents

- [PERL Scripts](#perl-scripts)
- [Python Scripts](#python-scripts)
- [R Scripts](#r-scripts)
- [Shell Scripts](#shell-scripts)

---

## PERL Scripts

### Demux_Intensities_And_Link_Basecalls

**Location:** `PERL/Demux_Intensities_And_Link_Basecalls/`

**Script:** `Demultiplex_HITS_V2.pl`

**Description:** Demultiplexes sequencing data from tab-delimited files. Assumes the barcode is in the 4th column. Supports gzipped input files. Used for tail-seq data processing on HiSeq platforms.

**Usage:**
```bash
perl Demultiplex_HITS_V2.pl <Input file> <LIST OF BARCODES SEPARATED BY COMMA> <NUMBER OF MISMATCHES ALLOWED IN BARCODE> <Output Directory>
```

**Features:**
- Handles both compressed (.gz) and uncompressed input files
- Allows user-defined number of mismatches in barcode matching
- Creates separate output files for each barcode plus an "unknown" file for unmatched reads

---

### FASTQ_Format_Conversion

**Location:** `PERL/FASTQ_Format_Conversion/`

**Script:** `FASTQ_Format_Conversion_v4.pl`

**Description:** Multi-purpose FASTQ processing tool that can perform three different tasks:
1. **changeidformat** - Converts FASTQ read header/identifier format between different formats
2. **makeread2umi** - Extracts UMI sequences from index reads into a separate FASTQ file
3. **demultiplex** - Demultiplexes single FASTQ file into multiple files based on expected barcodes

**Usage:**
```bash
FASTQ_Format_Conversion_v4.pl --task <changeidformat|makeread2umi|demultiplex> --read1fastq <path> --outputfolder <path> [options]
```

**Key Options:**
- `--task`: Required. One of "changeidformat", "makeread2umi", or "demultiplex"
- `--read1fastq`: Required. Path to read 1 FASTQ file
- `--read2fastq`: Required for paired-end data
- `--listofbarcodes`: Required for demultiplex task (comma-separated)
- `--mismatches`: Optional. Maximum mismatches allowed (default: 2)
- `--barcodelength`: Optional. Length of barcode (default: 6, required for makeread2umi)

For detailed usage, see the README in the script directory or use `--man` option.

---

## Python Scripts

### ModifyGTF

**Location:** `Python/`

**Script:** `ModifyGTF.py`

**Description:** Modifies GTF (Gene Transfer Format) files to be compatible with 10x Genomics Cell Ranger reference building. Specifically designed for Arabidopsis GTF files, but can be adapted for other organisms. Ensures proper attribute formatting required by Cell Ranger.

**Features:**
- Parses and reformats GTF attributes
- Adds transcript features
- Ensures proper formatting for Cell Ranger compatibility

**Note:** Currently has hardcoded input/output paths that should be modified before use.

---

### star-qualimap-summary

**Location:** `Python/star-qualimap-summary/`

**Script:** `scripts/summarizeoutputs/parse.py`

**Description:** Parses and summarizes output files from STAR aligner and Qualimap RNA-seq quality control. Generates summary Excel files with alignment statistics and quality metrics.

**Usage:**
```bash
python parse.py <full path to directory with output files from STAR and Qualimap>
```

**Features:**
- Automatically locates STAR log files and Qualimap summary files
- Assumes STAR files are in the main input directory
- Qualimap output can be in subdirectories
- Generates consolidated summary reports in Excel format

**Dependencies:**
- pandas
- numpy
- Custom modules: `STAR.py`, `qualimap.py`, `basicio.py`

---

## R Scripts

### CRISPR_CountSummary

**Location:** `R/`

**Script:** `CRISPR_CountSummary_v1.2.R`

**Description:** Merges multiple CRISPR count files into a single summary table. Reads all `.count` files matching the pattern `^[ACGT].*.count` from a directory and combines them into a single matrix with sample names as columns.

**Usage:**
```bash
Rscript CRISPR_CountSummary_v1.2.R -i <input directory> -o <output directory>
```

**Options:**
- `-i, --indir`: Path to input directory containing count files
- `-o, --outdir`: Path to output directory

**Features:**
- Recursively searches for count files matching the pattern
- Merges all count files based on ID column
- Converts ID column to row names in final output
- Handles missing values with `all=TRUE` merge option

---

## Shell Scripts

### Alignment.sh

**Location:** `Shell/`

**Description:** Automated RNA-seq alignment and quantification pipeline using STAR aligner. Designed for LSF cluster environments. Performs alignment, indexing, counting (using both featureCounts and HTSeq), and quality control with Qualimap.

**Usage:**
```bash
./Alignment.sh -1 <read1fastq.gz> -i <STAR_Index> -g <GTF_File> -o <Output_Directory> -s <stranded>
```

**Options:**
- `-1`: Read 1 FASTQ file (gzipped)
- `-i`: STAR genome index directory
- `-g`: GTF annotation file (sorted)
- `-o`: Output directory
- `-s`: Strandedness ("yes", "no", or "reverse")

**Features:**
- Currently supports single-end RNA-seq reads only
- Submits jobs to LSF cluster with proper dependencies
- Generates BAM files sorted by coordinate
- Produces counts using both featureCounts and HTSeq
- Creates Qualimap quality control reports
- Creates log files for all steps

**Dependencies:**
- STAR
- HTSeq Count
- featureCounts
- samtools
- qualimap

---

### CRISPR-Generate-Count.sh

**Location:** `Shell/`

**Description:** Generates count files for CRISPR screening experiments. Aligns reads to a reference guide RNA library and counts perfect matches. Assumes each guide RNA is 20 bases long.

**Usage:**
```bash
./CRISPR-Generate-Count.sh -1 <read1fastq.gz> -i <Reference_Fasta> -o <Output_Directory> -3 <bases_to_trim_3prime> -5 <bases_to_trim_5prime>
```

**Options:**
- `-1`: Read 1 FASTQ file (gzipped)
- `-i`: Reference FASTA file containing guide RNAs
- `-o`: Output directory
- `-3`: Bases to trim from 3' end
- `-5`: Bases to trim from 5' end

**Features:**
- Automatically builds Bowtie index from reference FASTA
- Validates that trimmed read length equals 20 bases (guide RNA length)
- Only counts perfect matches (MD:Z:20)
- Filters for uniquely mapped reads
- Outputs count file with guide RNA ID and count

**Dependencies:**
- Bowtie
- samtools

---

### downsample.sh

**Location:** `Shell/`

**Description:** Downsamples paired-end FASTQ files to a specified number of reads using seqtk. Maintains read pairing by using the same random seed for both forward and reverse reads.

**Usage:**
```bash
./downsample.sh <input_directory> <output_directory> <number_of_reads> <file_pattern>
```

**Arguments:**
1. Input directory containing FASTQ files
2. Output directory for downsampled files
3. Number of reads to sample
4. File pattern prefix (e.g., sample name)

**Features:**
- Automatically finds forward (`*_R1_*.fastq.gz`) and reverse (`*_R2_*.fastq.gz`) files
- Uses same random seed for both files to maintain pairing
- Runs forward file processing in background for efficiency
- Outputs compressed FASTQ files

**Dependencies:**
- seqtk
- gzip

---

### PreProcess-SLURM.sh

**Location:** `Shell/`

**Description:** Preprocessing pipeline for sequencing run folders. Submits FastQC and FastQ Screen jobs to SLURM cluster for quality control of all FASTQ files in a run folder.

**Usage:**
```bash
./PreProcess-SLURM.sh -r <runfolder>
```

**Options:**
- `-r`: Path to sequencing run folder

**Features:**
- Creates organized directory structure (FASTQ/, FASTQC/, FASTQSCREEN/)
- Submits parallel jobs for each FASTQ file
- Uses SLURM job scheduler
- FastQC with 4 threads, no extraction
- FastQ Screen with 1 million read subset

**Dependencies:**
- FastQC
- FastQ Screen
- SLURM cluster environment
- Bowtie2 (for FastQ Screen)

---

### PreProcess-Summary.sh

**Location:** `Shell/`

**Description:** Generates MultiQC summary reports from FastQC and FastQ Screen outputs. Creates both overall summaries and lane-specific interactive reports.

**Usage:**
```bash
./PreProcess-Summary.sh -r <runfolder>
```

**Options:**
- `-r`: Path to sequencing run folder

**Features:**
- Generates overall MultiQC report for all FastQC results
- Creates lane-specific (L001-L008) MultiQC reports for both FastQC and FastQ Screen
- Handles both `L00X` and `LX` lane naming conventions
- Outputs interactive HTML reports
- Configures maximum table rows for large datasets

**Dependencies:**
- MultiQC
- FastQC output files
- FastQ Screen output files

---

## General Notes

- Most scripts include error checking and usage messages
- Cluster-specific scripts (LSF/SLURM) may need modification for different cluster environments
- Some scripts have hardcoded paths that should be updated before use
- Always check script headers and comments for specific requirements and dependencies
- For detailed usage of individual tools, refer to README files in respective subdirectories