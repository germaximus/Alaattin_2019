# Transcriptome analysis

**Prerequisites:**  
[cutadapt 1.18](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.6.1d](https://github.com/alexdobin/STAR)  
[featureCounts 1.6](http://bioinf.wehi.edu.au/featureCounts/)  
Transcriptome samples were sequenced in paired-end 150 nt mode on Illumina sequencer.
Raw sequencing files are available from [GEO]().

### Preparing genome annotation and index files
Human genomic sequences and annotation files (GRCh38.p12) were downloaded from the [NCBI repository](http://ftp.ncbi.nih.gov/genomes/H_sapiens/).

| files             | MD5 check sum (unzipped)         | Description                                               |
| ----------------- |:--------------------------------:| ----------------------------------------------------------|
| GRCh38.p12.fa     | a4cac7d7ac4dd31ac68b384b10cf444d | RNA in fasta format, coding + noncoding                   |
| GRCh38.p12.fna    | 860290186a4ee3e95cd48dc528a45363 | Genome sequence, chromosomes and extrachromosomal contigs |
| GRCh38.p12.gbk    | 3c35b07e638485984479d50dd5cfebca | RNA in gene bank format, coding + noncoding               |
| GRCh38.p12.gff    | 56394751c00a5bdfb74152a7ed146855 | Genome annotation                                         | 

### Customizing genome annotation  
<details><summary><b>Edit chromosome names in GRCh38.p12.fna genome file.</b></summary>
STAR manual recommends not having spaces in contig names.   
     
```perl
#!/usr/bin/perl
# write output to GRCh38.p12.headers.fna
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

<details><summary><b>Drop extra-chromosomal contigs.</b></summary>
     
```perl
#!/usr/bin/perl
# write output to GRCh38.p12.chromosomes.fna
open (INPUT, '<GRCh38.p12.headers.fna') or die "Can't open file";

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

<details><summary><b>Add PhiX genome, mitochondrial genome and rRNA precursos</b></summary>
     
 ```bash
dos2unix GRCh38.p12.chromosomes.fna # if the file was created on Windows machine
cat GRCh38.p12.chromosomes.fna PhiX.fna Human_rmRNA >GRCh38.p12.custom.fna
```    
</details>


<details><summary><b>Drop annotation of extrachromosomal contigs</b></summary>

```perl
#!/usr/bin/perl
# write output in GRCh38.p12.chromosomes.gff

my %chromosomes =('NC_000001.11',0,
                  'NC_000002.12',0,
                  'NC_000003.12',0,
                  'NC_000004.12',0,
                  'NC_000005.10',0,
                  'NC_000006.12',0,
                  'NC_000007.14',0,
                  'NC_000008.11',0,
                  'NC_000009.12',0,
                  'NC_000010.11',0,
                  'NC_000011.10',0,
                  'NC_000012.12',0,
                  'NC_000013.11',0,
                  'NC_000014.9',0,
                  'NC_000015.10',0,
                  'NC_000016.10',0,
                  'NC_000017.11',0,
                  'NC_000018.10',0,
                  'NC_000019.10',0,
                  'NC_000020.11',0,
                  'NC_000021.9',0,
                  'NC_000022.11',0,
                  'NC_000023.11',0,
                  'NC_000024.10',0 );

open (INPUT, '<', '../NCBI_GRCh38.p12_originals/GRCh38.p12.gff') or die "Can't open file";

for (my $i=0; $i < 6; $i++){$line = <INPUT>; print $line;}
while($line = <INPUT>)
  {
     if(substr($line,0,3) eq '###') {print '###'; last;}
     else {
            @line = split /\s/,$line;
            if(exists($chromosomes{$line[1]}))
               {
                 print $line;
                 $line = <INPUT>;
                 do {print $line; $line = <INPUT>;} until(substr($line,0,2) eq '##');
                 redo;
               }
             else {do {$line = <INPUT>;} until(substr($line,0,2) eq '##'); redo;}
          }
  }
close(INPUT);
```

</details>


<details><summary><b>Drop 'Gnomon' (predicted) records from gff file and only keep 'RefSeq' (manually curated).</b></summary>
  
```perl
#!/usr/bin/perl
# write output in GRCh38.p12.RefSeq.gff

open (INPUT, "<GRCh38.p12.chromosomes.gff"); 
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
```bash
Discard_noncoding_annotation.R
#save new annotation as GRCh38.p12.Refseq.coding.gff
```
</details>


### Sequencing reads filtering and mapping   
<details><summary><b>Illumina adapters trimming.</b></summary>

```bash
cutadapt -j 20 -m 75 -O 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o out.1.fastq -p out.2.fastq read.1.fq.gz read.2.fq.gz
# -j      - number of threads
# -m      - discard read pair if any of the mates if shorter than 75 nucleotides after adapter trimming
# -O      - minimal length of the adapter to be considered for trimming
```
</details>

<details><summary><b>Build human genomic index for STAR.</b></summary>
     
```bash
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ./Human_index/ --genomeFastaFiles ./GRCh38.p12.STAR.fna --sjdbGTFfile ./GRCh38.p12.Refseq.codingSTAR.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 149
```
</details>

<details><summary><b>Read mappping with STAR.</b></summary>

```bash
STAR --genomeLoad LoadAndExit --genomeDir ../STAR-2.6.1d/Human_index/ 	# load genome once in the shared memory
STAR --runThreadN 40 --outSAMtype BAM Unsorted --outSAMmultNmax 1 --genomeLoad LoadAndKeep --genomeDir ../STAR-2.6.1d/Human_index/ --readFilesIn out.1.fastq out.2.fastq --outFileNamePrefix ./OUT_folder 
STAR --genomeLoad Remove 	# remove loaded genome from shared memory
# ipcs - check shared memory consumption
# ipcrm - remove object from shared memory
```
</details>

<details><summary><b>Count reads per gene.</b></summary>
  
 ```bash
featureCounts -T 20 -p -t exon -g gene -s 0 Aligned.out.bam -a GRCh38.p12.Refseq.codingSTAR.gff -o feature.counts #counting gene expression
# -T      -    number of threads
# -t exon -    only count reads over 'exon' subfeatures (3rd column of the GFF3 annotaiton file)
# -g gene -    integrate reads over 'gene' feature
# -s 0    -    library type (unstranded)
# -a      -    genome annotation GFF3 file
# -p      -    paired-end mode (count fragments instead individual reads)
 ```
</details>

