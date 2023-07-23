# Simulated Example Data
This directory contains a small simulated dataset for testing and getting familiar with TieBrush.

## Running TieBrush and TieCov

Output of the commands below is provided along with the simulated data and includes the following files:
1. example/t1/t1.bam and example/t2/t2.bam - collapsed representations of simulated tissue #1 and #2 computed with TieBrush
2. example/t1/t1.coverage.bedgraph and example/t2/t2.coverage.bedgraph - Coverage computed in BED/BedGraph format
3. example/t1/t1.sample.bedgraph and example/t2/t2.sample.bedgraph - Approximated number of samples computed for each mapped position in BED/BedGraph format
4. example/t1/t1.junctions.bed and example/t2/t2.junctions.bed - Splice Junction coordinates and coverage computed in BED format

### Collapsing reads and building representative alignment for each tissue with TieBrush

```tiebrush -o t1/t1.bam t1/t1s0.bam t1/t1s1.bam t1/t1s2.bam t1/t1s3.bam t1/t1s4.bam t1/t1s5.bam t1/t1s6.bam t1/t1s7.bam t1/t1s8.bam t1/t1s9.bam```

```tiebrush -o t2/t2.bam t2/t2s0.bam t2/t2s1.bam t2/t2s2.bam t2/t2s3.bam t2/t2s4.bam t2/t2s5.bam t2/t2s6.bam t2/t2s7.bam t2/t2s8.bam t2/t2s9.bam```

### Extracting BED-formatted summaries for each tissue with TieCov

```tiecov -s t1/t1.sample -c t1/tb.coverage -j t1/tb.junctions t1/t1.bam```

```tiecov -s t2/t2.sample -c t2/tb.coverage -j t2/tb.junctions t2/t2.bam```

### Visualizing results in IGV
Upon completion, tiecov will have generated several summary files for each tissue.
Along with the collapsed alignment, these summaries can be viewed and scrutinized using
genome browsers such as [IGV](http://software.broadinstitute.org/software/igv/)
that support BAM (for collapsed alignment), BED (for the junctions track)
and BEDgraph (for the sample and coverage tracks) formats.

## Data
Example dataset contains mapped simulated reads (BAM) from two loci (NEFL and SLC25A3)
which are known to exhibit differential expression and splicing in several tissues.
Data provided is separated into two tissues (t1 and t2), each comprised of 10 samples.

NEFL (chr8:24,950,955-24,956,721 (GRCh38/hg38)) neurofilament-coding gene is known to be primarily expressed in certain cells of the 
Brain, however low expression are typically observed in several cell types in the heart.

![](https://github.com/alevar/tiebrush/blob/master/example/nefl.sim.png)

**Figure 1.** *Comparison of transcription of the NEFL gene in brain (top) and 
heart (bottom) tissues from the simulated dataset. Each tissue is represented by 
three tracks produced by TieCov: read coverage (top), percent of samples containing the 
reads (middle), and splice junctions (bottom). The plot illustrates the higher 
prevalence and expression of the gene in brain tissue*

SLC25A3 (chr12:98,593,591-98,606,379 (GRCh38/hg38)) can be expressed by several dominant
alternative isoforms, with splicing patterns being known to be tissue specific. In case of 
cells from heart and brain - alternative forms used differ in the third exon [1].

![](https://github.com/alevar/tiebrush/blob/master/example/slc25a3.sim.png)

**Figure 2.** *Comparison of transcription of the SLC25A3 gene in brain (top) 
and heart (bottom) tissues from the simulated dataset. While the gene is expressed 
in both tissues, coverage data clearly indicate an exon switch where the 3rd and 
4th exons are expressed at dramatically different levels in the two tissues*

## Simulation Protocol
1. 10 samples were randomly selected from brain tissue and another 10 samples were randomly selected 
from the heart tissue of the GTEx project [2] (Table 1). 
2. Samples were aligned with HISAT2 [3] against the GRCh.38(patch 12) [4].
3. Transcriptomes for each sample were individually assembled and quantified using StringTie2 [5].
4. Using GffCompare [6] samples were individually compared against the RefSeq annotation [7].
5. Transcript structures (GTF) and coverage data of isoforms overlapping the NEFL and SLC25A3 genes 
   were extracted for each sample.
7. Extracted transcriptome (FASTA sequences obtained with gffread [6]) and coverage data 
   were used to simulate 101bp reads using polyester [8].
8. Fasta files produced by polyester were directly converted into SAM format using sim2sam utility [9].

| Brain Samples | Heart Samples |
|---------------|:--------------|
|SRR1368222     | SRR1345221    |
|SRR1370173  	| SRR1441464    |
|SRR1374221 	| SRR1432366    |
|SRR1363810     | SRR661229     |
|SRR1467882     | SRR1337388    |
|SRR1471773  	| SRR1458484    |
|SRR1473590   	| SRR816358     |
|SRR814989   	| SRR1352213    |
|SRR817658   	| SRR1435584    |
|SRR818146   	| SRR1345736    |

**Table 1.** *GTEx samples used to model parameters for the simulation*

## References
[1] Wang, Eric T., et al. "Alternative isoform regulation in human tissue transcriptomes." Nature 456.7221 (2008): 470-476.

[2] Lonsdale, John, et al. "The genotype-tissue expression (GTEx) project." Nature genetics 45.6 (2013): 580-585.

[3] Kim, Daehwan, et al. "Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype." Nature biotechnology 37.8 (2019): 907-915.

[4] Schneider, Valerie A., et al. "Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring quality of the reference assembly." Genome research 27.5 (2017): 849-864.

[5] Kovaka, Sam, et al. "Transcriptome assembly from long-read RNA-seq alignments with StringTie2." Genome biology 20.1 (2019): 1-13.

[6] Pertea, Geo, and Mihaela Pertea. "GFF utilities: GffRead and GffCompare." F1000Research 9 (2020).

[7] O'Leary, Nuala A., et al. "Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation." Nucleic acids research 44.D1 (2016): D733-D745.

[8] Frazee, Alyssa C., et al. "Polyester: simulating RNA-seq datasets with differential transcript expression." Bioinformatics 31.17 (2015): 2778-2784.

[9] Varabyou, Ales, Steven L. Salzberg, and Mihaela Pertea. "Effects of transcriptional noise on estimates of gene and transcript expression in RNA sequencing experiments." Genome Research 31.2 (2021): 301-308.
