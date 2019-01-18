# Transcriptome analysis

**Prerequisites:**  
[Cutadapt 1.18](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.6.1d](https://github.com/alexdobin/STAR)  
[featureCounts 1.6](http://bioinf.wehi.edu.au/featureCounts/)  

### Preparing genome annotation and index files
Human genomic sequences and annotation files (GRCh38.p12) were downloaded from the [NCBI repository](ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/).

| files             | MD5 check sum (unzipped)         | Description                                               |
| ----------------- |:--------------------------------:| ----------------------------------------------------------|
| GRCh38.p12.fa     | a4cac7d7ac4dd31ac68b384b10cf444d | RNA in fasta format, coding + noncoding                   |
| GRCh38.p12.fna    | 860290186a4ee3e95cd48dc528a45363 | Genome sequence, chromosomes and extrachromosomal contigs |
| GRCh38.p12.gbk    | 3c35b07e638485984479d50dd5cfebca | RNA in gene bank format, coding + noncoding               |
| GRCh38.p12.gff    | 56394751c00a5bdfb74152a7ed146855 | Genome annotation                                         | 

### Customizing genome annotation  
<details><summary><b>Edit chromosome names in GRCh38.p12.fna genome file.</b></summary>
STAR manual recommends not having spaces in contig names 
     
```{perl}
#!/usr/bin/perl
open (INPUT, '<GRCh38.p12.fna') or die "Can't open file";

while ($line = <INPUT>)  {
     @line = split('\s+', $line);
     if(substr($line[0],0,1) eq '>') {
           print $line[0]."\n";
           while ($line = <INPUT>) {
                  if (substr($line,0,1) ne '>') { print $line;   }
                  else {last;}                
           }
           redo; 
     } 
}
close(INPUT);
```
</details>

<details><summary><b>Drop 'Gnomon' (Predicted) records from gff file and only keep 'RefSeq' (manually curated).</b></summary>
  
```perl
#!/usr/bin/perl
# Usage: Script.pl $ARGV0
open (INPUT, "<$ARGV[0]"); 
for($i=0; $i < 8; $i++) {$line = <INPUT>; print $line;}
while ($line = <INPUT>)  {
     @fields = split /\t/, $line;
     if($fields[1] ne 'Gnomon' && $fields[1] ne 'Curated Genomic')  {  print $line;  }
}
close(INPUT); 
```
</details> 

<details><summary><b>Drop non-coding features from annotation.</b></summary>

Remove non-coding RNA genes, leave only coding genes with their mRNA, transcript, exon, and CDS children. Fix the gff annotation from previous script by matching gene coordinates with the childern coordinates (occured due to removal of Gnomon features).  
```{bash}
Discard_noncoding_annotation.R
```



Illumina adapters trimming
```{bash}
cutadapt -j 20 -m 75 -O 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o out.1.fastq -p out.2.fastq read.1.fq.gz read.2.fq.gz
```
Build human genomic index for STAR
```{bash}
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ./Human_index/ --genomeFastaFiles ./GRCh38.p12.STAR.fna --sjdbGTFfile ./GRCh38.p12.Refseq.codingSTAR.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 149
```
Read mappping with STAR
```{bash}
STAR --genomeLoad LoadAndExit --genomeDir ../STAR-2.6.1d/Human_index/ 	# load genome once in the shared memory
STAR --runThreadN 40 --outSAMtype BAM Unsorted --outSAMmultNmax 1 --genomeLoad LoadAndKeep --genomeDir ../STAR-2.6.1d/Human_index/ --readFilesIn out.1.fastq out.2.fastq --outFileNamePrefix ./OUT_folder 
STAR --genomeLoad Remove 	# remove loaded genome from shared memory
# ipcs - check shared memory consumption
# ipcrm - remove object from shared memory
