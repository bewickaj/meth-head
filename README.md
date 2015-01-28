##Meth-head

A series of programs, commands, and scripts to parse, analyze, and plot bisulfite sequencing (BS-seq) data. Some programs, and scripts are mine, while others are not – why reinvent the wheel, right?

---

###Janitor work: Cleaning, trimming, aligning, and formatting

####FASTX-Toolkit: Quality filter your reads

There are many programs out there to do this, but here is an example using fastq_quality_filter from the FASTX-Toolkit. See [a link](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html) for more information. Run as so:

````bash
fastq_quality_filter -Q 33 -q 20 -p 80 -i read1.fastq -o read1_filtered.fastq
fastq_quality_filter -Q 33 -q 20 -p 80 -i read2.fastq -o read2_filtered.fastq
````

| Flag | Description |
| --- | --- |
| -Q | ASCII quality offset |
| -q | Minimum quality score to keep |
| -p | Minimum percent of bases that must have [-q] quality |
| -i | Input fastq |
| -o | Output fastq |

---

####Trimmomatic: Trim your reads

Quality along the read changes, with the beginnings, and especially ends, being particularly bad. Again, there are lots of programs out there, but here is an example. See [a link](http://www.usadellab.org/cms/?page=trimmomatic) for more information. Run as so:

````bash
java -classpath trimmomatic-0.30.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 -trimlog trimlog.txt read1_filtered.fastq.gz read1_filtered.fastq.gz read1_filtered_paired.fastq.gz read1_filtered_unpaired.fastq.gz read2_filtered_paired.fastq.gz read2_filtered_unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
````

| Flag | Description |
| --- | --- |
| LEADING: | Cut bases off the start of a read, if below a threshold quality |
| TRAILING: | Cut bases off the end of a read, if below a threshold quality |
| SLIDINGWINDOW: | Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold |
| MINLEN: | Drop the read if it is below a specified length |
| LEADING: | Remove leading low quality or N bases (below quality 3) |
| TRAILING: | Remove trailing low quality or N bases (below quality 3) |
| SLIDINGWINDOW: | Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 |
| MINLEN:50 Drop reads below 50 bases long; <60% of read length (75bp) |

---

####Bismark: Align your reads

...now we're getting somewhere. Below are the commands and flags commonly used for Bismark. The first step is to prepare your genome (fasta file of your chromosomes/pseudomolecules/contigs/scaffolds), which generates two files in your genome folder; converted C->T and G->A fasta files. Next, your reads are aligned to the converted genomes. Finally, methylated and un-methylated sequence coverage information at C sites is extracted. See [a link](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) for more information. Run as so (some flags need to be filled in):

````bash
bismark_genome_preparation --bowtie2 --path_to_bowtie /usr/local/bowtie2/ --verbose genome/
bismark -q --path_to_bowtie /usr/local/bowtie2/ -p 4 --bowtie2 -N 1 genome/ -1 read1_filtered_paired.fastq -2 read2_filtered_paired.fastq
bismark_methylation_extractor -p --no_overlap --comprehensive --report --bedGraph --CX_context --buffer_size --cytosine_report --counts --genome_folder genome/ --split_by_chromosome read1_read2_filtered_paired.sam
````

---

####pbinom.pl: Call sites as methylated or un-methylated

A simple script I wrote that calculates the cumulative probability function for a binomial distribution at each C. At any given C (with more than >2 sequence coverage) we can calculate the probability of the C being methylated by using the number of reads methylated, number of reads that are un-methylated, and the error rate (i.e., non-conversion rate of the chloroplast genome, since it is un-methylated, or error rate from lambda DNA sequences). For example, at a C with 7 out of 10 reads methylated, and a non-conversion rate of 0.005, what is the probability that we would get either 7, 8, 9, or 10 methylated bases out of 10 total reads based on errors alone. Below is a generic example of how to run the script, but there is also a help flag (--help):

````bash
./pbinom.pl --bismark_CX_report=CX_report.txt --outfile=out.vcf --probability=0.005 --confidence=0.05
````

The results are printed to an outfile in VCF format with the following columns:

| Column name | Description |
| --- | --- |
| chr | Name of chromosome |
| pos | Base pair position of the C |
| strand | What strand the C was sequenced from |
| mc_class | In plants, the different methylation contexts (or classes): CG, CHG, or CHH (H = A/C/T) |
| mc_count | The number of times the C was sequenced as methylated |
| total | Total number of reads (methylated and un-methylated) |
| methylated | A binary code for the methylation state of the C – 0 = un-methylated, or 1 = methylated – based on the user defined confidence threshold |

---

####bedtools: Find the intersect between genes/TEs/features and C's

bedtools is an extremely handy program with many functions, but one of my favourites is intersect. "bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported." bedtools only works for certain file types – mainly BED/GFF/VCF and BAM. We have already generated a VCF file of C sites, so all we need is a file (BED in this case) of our features. Since, gene-body methylation is a popular topic of interest I have provided an example using exon locations. See [a link](https://bedtools.readthedocs.org/en/latest/content/tools/intersect.html) for more information.

Your BED file should look something like this:

| Column number | Description |
| --- | --- |
| Col0 | Chromosome name – should match the naming format used in the VCF file above |
| Col1 | Start of the feature |
| Col2 | End of the feature |
| Col3 | Gene (model) name |
| Col4..n | Not necessary for the following script, but can be whatever you want – I generally include exon number, and strand |

Run as so:

````bash
bedtools intersect -a pbinom.vcf -b file.bed -wb > intersect.out
````

Your outfile should have the following format:

| Column number | Description |
| --- | --- |
| col0 | Col0 from pbinom.vcf, i.e., chromosome name |
| col1 | Col1 from pbinom.vcf, i.e., base pair position of the C |
| col2 | Col2 from pbinom.vcf, i.e., strand the C was sequenced from |
| col3 | Col3 from pbinom.vcf, i.e., methylation contexts (or classes) |
| col4 | Col4 from pbinom.vcf, i.e., number of times the C was sequenced as methylated |
| col5 | Col5 from pbinom.vcf, i.e., total number of reads |
| col6 | Col6 from pbinom.vcf, i.e., binary code for the methylation state |
| col7 | Col0 from file.bed, i.e., chromosome name |
| col8 | Col1 from file.bed, i.e., start of the feature |
| col9 | Col2 from file.bed, i.e., end of the feature |
| col10 | Col3 from file.bed, i.e, gene (model) name |
| col11 | Col4 from file.bed, i.e., exon number |
| col12 | Col5 from file.bed, i.e., strand |

---

####methtable.pl: Generate a summary table of methylated C's, and probability of being methylated for each class for each gene

A perl script I wrote that tallies the number of methylated and un-methylated sites (with coverage >2), and calculates the cumulative probability function for a binomial distribution for each C class, for each gene. The idea is similar to pbinom.pl; on a per C class basis, this script uses the number of C's that are methylated, and un-methylated in the gene, and the proportion of methylated to un-methylated C's across all genes as the error rate. Also, there is an option to perform a false discovery rate (FDR) correction on p-values, which generates an outfile with file extension "_fdr.txt". Run as so:

````bash
./methtable.pl --intersect_infile=intersect.out --outfile=table.txt --fdr=yes
````

Your outfile should have the following format (if --fdr=yes):

| Column name | Description |
| --- | --- |
| gene_model | Name of the gene (model) |
| total_CG | Total number of CG sites with >2 coverage in this gene |
| meth_CG | Number of methylated CG sites with >2 coverage in this gene |
| pvalue_CG | Probability of this gene being CG gene-body methylated |
| fdr_CG | Adjusted p-value – included if --fdr=yes |
| total_CHG | Total number of CHG sites with >2 coverage in this gene |
| meth_CHG | Number of methylated CHG sites with >2 coverage in this gene |
| pvalue_CHG | Probability of this gene being CHG gene-body methylated |
| fdr_CHG | Adjusted p-value – included if --fdr=yes |
| total_CHH | Total number of CHH sites with >2 coverage in this gene |
| meth_CHH | Number of methylated CHH sites with >2 coverage in this gene |
| pvalue_CHH | Probability of this gene being CHH gene-body methylated |
| fdr_CHH | Adjusted p-value – included if --fdr=yes |

---
