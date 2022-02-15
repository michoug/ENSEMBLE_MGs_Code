# Pipeline for metagenome analysis of the ENSEMBLE datasets

Grégoire Michoud

November 2021

## Trimming and quality

### Conda environment

`mamba create -n trim_galore -c bioconda -c conda-forge trim-galore`

### Run

All paired `fastq.gz` files should be placed in the `Raw` folder and should follow the convention `sampleName_R1.fastq.gz` and `sampleName_R2.fastq.gz`

``` bash
conda activate trim_galore
for i in Raw/*R1.fastq.gz;
do
    trim_galore -j 4 --fastqc --paired $i ${i%R1.fastq.gz}R2.fastq.gz;
done

mkdir Trim
mv Raw/*val_*.fastq.gz Trim
```

## Assembly

### Conda environment

`mamba create -n megahit -c bioconda -c conda-forge megahit`

### Script

``` bash
conda activate megahit
for i in Trim/*R1_val_1.fastq.gz;
do
    megahit -1 $i -2 ${i%R1_val_1.fastq.gz}R2_val_2.fastq.gz -m 2.4e+11 -t 28 --kmin-1pass --min-contig-len 1000 -o ${i%_mtg_sed_R1_val_1.fastq.gz}
done

mkdir Contigs
for i in */final.contigs.fa
do
    mv $i Contigs/${i%/final.contigs.fa}.fa
done
```

## Annotation

### Conda environment

`mamba create -n annotation -c bioconda -c conda-forge prodigal seqkit eggnog-mapper`

Had some issues with the conda installation of `eggnog-mapper` so downloaded the release (https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/2.1.2.tar.gz)
and used this version inside the conda environment

### Script

``` bash
conda activate annotation

cd Contigs

# Rename contigs to avoid duplicates

for i in *fa
do
    perl -pe "s/>/>${i%.fa}\_/g" $i > ${i%.fa}_rename.fa
done

# Split fasta files into 28 parts as prodigal is single threaded

for i in *_rename.fa
do
    seqkit split2 -p 28 -j 28 $i;
done

# Run in parallel as you wish to greatly speed up the process

for i in *split/*fa
do
    prodigal -a $i\a -d ${i%.fa}.ffn -f gff -i $i -p meta -q -o ${i%.fa}.gff;
done

for i in *split
do
    cat $i\*.faa > ${i%.fa.split}.faa
    cat $i\*.ffn > ${i%.fa.split}.ffn
    cat $i\*.gff > ${i%.fa.split}.gff
done

cd ..

mkdir Annotations

mv Contigs/*faa Annotations
mv Contigs/*ffn Annotations
mv Contigs/*faa Annotations

# The --dbmem parameter should be removed if the memory of the computer used is a little low but its addition strongly reduce the run time

for i in Annotations/*faa
do
    /work/sber/Databases/eggnog-mapper-2.1.2/emapper.py --dbmem --cpu 28 -i $i --itype proteins -m diamond --sensmode very-sensitive -o ${i%.faa}_egg.txt
done
```

## MAGs generation

### Conda environment

`mamba create -n MAGS_prep -c bioconda -c conda-forge coverm concoct metabat2 samtools`

### Script

```bash

conda activate MAGS_prep

mkdir BamFiles

# list reads files in the Trim folder

ls ../Trim > listreads.txt

for i in Contigs/*fa
do
    coverm make --reference -i $i -t 28 -o BamFiles --discard-unmapped -c `< listreads.txt`;
done

for i in BamFiles/*bam
do
    samtools index $i
done

mkdir MAGs

cd Contigs

for i in *fa
do
    jgi_summarize_bam_contig_depths --outputDepth ../MAGs/${i%.fa}_metabat.txt ../BamFiles/${i%.fa}*bam
done

for i in *fa
do
    cut_up_fasta.py $i -c 10000 -o 0 --merge_last -b ${i%.fa}_10K.bed > ${i%.fa}_10K.fasta
    concoct_coverage_table.py ${i%.fa}_10K.bed ../BamFiles/${i%.fa}*bam > ../MAGs/${i%.fa}_concoct.txt
done

cd ../MAGs

for i in *metabat.txt
do
    metabat2 -a $i -o ${i%.t.txt}/${i%.t.txt} -i ../Contigs/${i%_metabat.txt}.fa -t 28 -m 1500
done

for i in *concoct.txt
do
    mkdir ${i%.txt}
    concoct --composition_file ../Contigs/${i%_concoct.txt}_10K.fasta --coverage_file $i -b ${i%.txt}/${i%.txt} -t 28
    merge_cutup_clustering.py ${i%.txt}/${i%.txt}_clustering_gt1000.csv > ${i%.txt}/${i%.txt}_clustering_merged.csv
    # extract_fasta_bins.py ../Contigs/${i%_concoct.txt}.fa ${i%.txt}/${i%.txt}_clustering_merged.csv --output_path ${i%.txt}/
done
```

## MAGs dereplication and quality

### Conda environment

`mamba create -n MAGS_qc -c bioconda -c conda-forge drep Das_Tool checkm gtdbtk`


```bash

cd MAGs


for i in *
do
    DAS_Tool -i $i\_metabat_mags.txt,$i\_concoct_mags.txt -l metabat,concoct -c $i.fa --proteins $i.faa --write_bins 1 -o $i --threads 25 --search_engine diamond
done


checkm lineage_wf -x fa --pplacer_threads 40 -t 50 -f allBins_checkm_2.txt allBins allBins_checkm_2

```