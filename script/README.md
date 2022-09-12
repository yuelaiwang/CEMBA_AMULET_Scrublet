# Scripts Run on TSCC

#### For Ren Lab members only. All paths in this README are TSCC server path, i.e. you can copy, paste, and use any path to a file/folder directly in TSCC/silencer server.

Define abbreviations for regularly used paths:
```
gold=/projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard
bin=/projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/bin
cemba_mop=/projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/cemba_mop
qsub=/projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/qsub
archive=/projects/ren-transposon/home/yangli/projects/CEMBA/00.data/archive
temp=/projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/temp
path1=/projects/ps-renlab/yangli/projects/CEMBA/00.data/szu
path2script=/projects/ps-renlab/yangli/scripts/snATACutils/bin
```

- [Step 1. Generate a pool of singlets](#step-1-generate-a-pool-of-singlets)
- [Step 2. Simulate datasets containing artificial doublets](#step-2-simulate-doublets)
- [Step 3. Remove doublets on simulated datasets](#step-3-remove-doublets-on-simulated-datasets)
- [Step 4. Compare between the two tools by PRC and AUPRC](#step-4-quantify-the-performance-of-the-two-tools-using-prc-and-auprc)

### Step 1: Generate a pool of singlets

#### For replications 1 and 2 of regions 2C, 3C, 4B, 5D (no UQ cutoff)

Identify AMULET doublets.

```
for i in `cat $gold/cemba.mop.sample.lst`;
do j=${i##*_};
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N AMULET_${i:5}
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -j oe

echo $i

# add CB attribute and sort the deduplicated bam file
cat <(samtools view -H $archive/$j/$i/processed/filtered_dedup.bam) <(samtools view $archive/$j/$i/processed/filtered_dedup.bam | awk 'BEGIN{FS = \":\"} {CB = \$1} {printf \"%s\\\t%s%s\\\n\", \$0, \"CB:Z:\", CB}') | samtools sort -o $gold/$i/doublet_removal/AMULET/filtered_dedup.srt.bam

# make singlecell.csv need by AMULET
# the meta table is used here. Theoretically, the meta table consists of singlets only. As we want to ensure we have a pool of singlets, we use this to provide information about which cells are actually singlets.
cat /projects/ps-renlab/yuw044/projects/CEMBA/metatable.tsv | awk '(\$2 == \"$i\"){printf \"%s,%s\\\n\", \$3, \"1\"}' > $gold/$i/doublet_removal/AMULET/singlecell.csv

# run AMULET 
java -jar /projects/ps-renlab/yuw044/apps/AMULET-v1.1_0124/snATACOverlapCounter.jar --forcesorted --iscellidx 1 $gold/$i/doublet_removal/AMULET/filtered_dedup.srt.bam $gold/$i/doublet_removal/AMULET/singlecell.csv ~/projects/CEMBA/practice_doublet_removal/AMULET/testing_10D/input_data/mouse_autosomes.txt $gold/$i/doublet_removal/AMULET/ &> $gold/$i/doublet_removal/log/${i}_find_overlaps.log
python3 /projects/ps-renlab/yuw044/apps/AMULET-v1.1_0124/AMULET.py --rfilter ~/projects/CEMBA/practice_doublet_removal/AMULET/testing_10D/input_data/mm10.blacklist.bed $gold/$i/doublet_removal/AMULET/Overlaps.txt $gold/$i/doublet_removal/AMULET/OverlapSummary.txt $gold/$i/doublet_removal/AMULET/ &> $gold/$i/doublet_removal/log/${i}_multiplet_detection.log

" | sed '1d' > $qsub/AMULET_${i}.qsub

qsub $qsub/AMULET_${i}.qsub -o $gold/$i/doublet_removal/log/AMULET_${i}.log
done
```

Identify Scrublet doublets.

```
for i in `cat $gold/cemba.mop.sample.lst`;
do j=${i##*_};
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N Scr_${i}
#PBS -l nodes=1:ppn=1,walltime=3:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -j oe

echo $i;
mkdir $gold/$i/doublet_removal/Scrublet

# generate qc.filter.RData and qc.filter.meta.txt
Rscript ${path1}/bin/snapATAC.qc.filter.R -i /projects/ps-renlab/yangli/projects/CEMBA/00.data/$j/$i/$i.snap \
                                        --tsse2depth /projects/ren-transposon/home/yangli/projects/CEMBA/00.data/archive/$j/$i/processed/stat.txt \
                                        -o $gold/$i/doublet_removal/Scrublet/${i} \
                                        --fragment_num 1000 --tsse_cutoff 10 &> $gold/$i/doublet_removal/log/${i}.qc.filter.log

Rscript $path2script/snapATAC.rmDoublets.R -i $gold/$i/doublet_removal/Scrublet/${i}.qc.filter.RData \
                                                      --mat gmat --rate 0.08 \
                                                      -o $gold/$i/doublet_removal/Scrublet/${i}.gmat &> $gold/$i/doublet_removal/log/${i}.gmat.rmDoublets.log; 

Rscript $path2script/snapATAC.fitDoublets.R -i $gold/$i/doublet_removal/Scrublet/${i}.gmat.qc.cluster.RData \
                                            -m $gold/$i/doublet_removal/Scrublet/${i}.gmat.rmDoublets.txt \
                                            -o $gold/$i/doublet_removal/Scrublet/${i}.gmat &> $gold/$i/doublet_removal/log/${i}.gmat.fitDoublets.log

" | sed '1d' > $qsub/scrublet_${i}.qsub

qsub $qsub/scrublet_${i}.qsub -o $gold/$i/doublet_removal/log/scrublet_${i}.log 
done
```

Find common singlets identified by both tools.

```
for i in `cat $gold/cemba.mop.sample.lst`;
do 
python $bin/find_common_singlets.py --meta $gold/$i/doublet_removal/Scrublet/$i.qc.filter.meta.txt --scrublet $gold/$i/doublet_removal/Scrublet/$i.gmat.fitDoublets.txt --AMULET $gold/$i/doublet_removal/AMULET/MultipletBarcodes_01.txt --output $gold/$i/doublet_removal/SingletBarcodes_01.txt
done
```

#### For the combined MOP datasets replications 1 and 2

For each of the eight cemba mop bam files, modify query name by adding the rep_region information and add CB attribute

```
# replication 1
for i in `cat $cemba_mop/rep1/sample.lst`;
do j=${i##*_};
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N modify_bam_${i}
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -e $gold/cemba_mop/log/modify_bam_$i.err
#PBS -m a
#PBS -A ren-group
#PBS -j oe

cat <( samtools view $archive/$j/$i/processed/filtered_dedup.bam -H ) <( samtools view $archive/$j/$i/processed/filtered_dedup.bam | sed 's/^/$i./' | awk 'BEGIN{FS = \":\"} {CB = \$1} {printf \"%s\\\t%s%s\\\n\", \$0, \"CB:Z:\", CB}') | samtools view -bS > ${cemba_mop}/rep1/$i.filtered.dedup.bam

" | sed '1d' > $qsub/modify_bam_$i.qsub

qsub $qsub/modify_bam_$i.qsub -o $gold/cemba_mop/log/modify_bam_$i.log
done

# replication 2
for i in `cat $cemba_mop/rep2/sample.lst`;
do j=${i##*_};
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N modify_bam_${i}
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -e $gold/cemba_mop/log/modify_bam_$i.err
#PBS -m a
#PBS -A ren-group
#PBS -j oe

cat <( samtools view $archive/$j/$i/processed/filtered_dedup.bam -H ) <( samtools view $archive/$j/$i/processed/filtered_dedup.bam | sed 's/^/$i./' | awk 'BEGIN{FS = \":\"} {CB = \$1} {printf \"%s\\\t%s%s\\\n\", \$0, \"CB:Z:\", CB}') | samtools view -bS > ${cemba_mop}/rep2/$i.filtered.dedup.bam

" | sed '1d' > $qsub/modify_bam_$i.qsub

qsub $qsub/modify_bam_$i.qsub -o $gold/cemba_mop/log/modify_bam_$i.log
done
```

Merge the bam files
```
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N merge_cemba_mop
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -e $gold/cemba_mop/log/merge_cemba_mop.err
#PBS -m ae
#PBS -A ren-group
#PBS -j oe
samtools merge ${cemba_mop}/rep1/CEMBA_MOP.filtered.dedup.bam ${cemba_mop}/rep1/CEMBA180409_2C.filtered.dedup.bam ${cemba_mop}/rep1/CEMBA171206_3C.filtered.dedup.bam ${cemba_mop}/rep1/CEMBA180104_4B.filtered.dedup.bam ${cemba_mop}/rep1/CEMBA180612_5D.filtered.dedup.bam
samtools merge ${cemba_mop}/rep2/CEMBA_MOP.filtered.dedup.bam ${cemba_mop}/rep2/CEMBA180410_2C.filtered.dedup.bam ${cemba_mop}/rep2/CEMBA171207_3C.filtered.dedup.bam ${cemba_mop}/rep2/CEMBA171213_4B.filtered.dedup.bam ${cemba_mop}/rep2/CEMBA180618_5D.filtered.dedup.bam
" | sed '1d' > $qsub/merge_cemba_mop.qsub
qsub $qsub/merge_cemba_mop.qsub -o $gold/cemba_mop/log/merge_cemba_mop.log
```

Identify AMULET & Scrublet doublets/singlets & find common singlets identified by both tools ()

```
# add the rep_region information to and merge the SingletBarcodes_01.txt files for each cemba rep_region to form the equivalent file for rep1/rep2, respectively
for i in `cat $cemba_mop/rep1/sample.lst`;
do
echo -e "
cat $gold/${i}/doublet_removal/SingletBarcodes_01.txt | sed 's/^/${i}./' >> $cemba_mop/rep1/SingletBarcodes_01.txt "
done
for i in `cat $cemba_mop/rep2/sample.lst`;
do 
echo -e "
cat $gold/${i}/doublet_removal/SingletBarcodes_01.txt | sed 's/^/${i}./' >> $cemba_mop/rep2/SingletBarcodes_01.txt "
done

# add the rep_region information to and merge the OverlapSummary.txt files 
head -n 1 $gold/${i}/doublet_removal/AMULET/OverlapSummary.txt > $cemba_mop/rep1/OverlapSummary.txt
head -n 1 $gold/${i}/doublet_removal/AMULET/OverlapSummary.txt > $cemba_mop/rep2/OverlapSummary.txt
for i in `cat $cemba_mop/rep1/sample.lst`;
do
cat $gold/${i}/doublet_removal/AMULET/OverlapSummary.txt | sed '1d' | awk -e '{printf "'$i'.%s\t%s\t%s\t'$i'.%s\t%s\n", $1, $2, $3, $4, $5}' >> $cemba_mop/rep1/OverlapSummary.txt
done
for i in `cat $cemba_mop/rep2/sample.lst`;
do
cat $gold/${i}/doublet_removal/AMULET/OverlapSummary.txt | sed '1d' | awk -e '{printf "'$i'.%s\t%s\t%s\t'$i'.%s\t%s\n", $1, $2, $3, $4, $5}' >> $cemba_mop/rep2/OverlapSummary.txt
done
```

### Step 2. Simulate doublets

#### For replications 1 and 2 of regions 2C, 3C, 4B, 5D (no UQ cutoff)

```
for i in `cat $gold/cemba.mop.sample.lst`;
do 
mkdir $gold/$i/log
for ((k = 1; k <= 10; k ++));
do
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N sim_$i_$k
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -e $gold/$i/log/cemba_simulate_${i}_${k}.err
#PBS -A ren-group
#PBS -j oe

mkdir $gold/$i/dataset${k}

samtools index $gold/$i/doublet_removal/AMULET/filtered_dedup.srt.bam

python $bin/simulate_cemba_doublets.py $gold/$i/doublet_removal/AMULET/filtered_dedup.srt.bam --poolsize 1100 --region $i /projects/ps-renlab/yuw044/projects/CEMBA/metatable.tsv MajorType $gold/$i/doublet_removal/SingletBarcodes_01.txt $gold/$i/doublet_removal/AMULET/OverlapSummary.txt $gold/$i/dataset${k}

" | sed '1d' > $qsub/cemba_simulate_${i}_${k}.qsub

qsub $qsub/cemba_simulate_${i}_${k}.qsub -o $gold/$i/log/cemba_simulate_${i}_${k}.log
done
done
```

#### For the combined MOP datasets replications 1 and 2

Simulate 10 datasets at each UQ interval (1000-4000, 4000-7000, 7000-10000, 10000-13000)
```
for j in rep1 rep2;
do
for (( i = 11; i <= 14; i++ ));
do
for (( k = 1; k <= 10; k ++ ));
do
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N sim_${j}_${i}_${k}
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -e ${cemba_mop}/log/cemba_simulation_${j}_${i}_${k}.err
#PBS -m a
#PBS -A ren-group
#PBS -j oe

mkdir ${cemba_mop}/simulated_datasets_${j}/dataset${i}_${k}

python $bin/simulate_cemba_doublets.py --fclower $((3000*${i}-32000)) --fchigher $((3000*${i}-29000)) --poolsize 1100 --cellID CellID $cemba_mop/$j/CEMBA_MOP.filtered.dedup.bam /projects/ps-renlab/yuw044/projects/CEMBA/metatable.tsv MajorType $cemba_mop/$j/SingletBarcodes_01.txt $cemba_mop/$j/OverlapSummary.txt $cemba_mop/simulated_datasets_$j/dataset${i}_${k}

" | sed '1d' > $qsub/cemba_simulation_${j}_${i}_${k}.qsub

qsub $qsub/cemba_simulation_${j}_${i}_${k}.qsub -o ${cemba_mop}/log/cemba_simulation_${j}_${i}_${k}.log 
done
done
done
```

### Step 3. Remove doublets on simulated datasets

#### For replications 1 and 2 of regions 2C, 3C, 4B, 5D (no UQ cutoff)

#### For replications 1 and 2 of regions 2C, 3C, 4B, 5D (UQ between 1000 and 4000)

#### For the combined MOP datasets replications 1 and 2

For each of the 40 simulated datasets, prepare AMULET input (the coordinate-sorted bam file) and run AMULET on each of them
```
for j in rep1 rep2;
do
for (( i = 11; i <= 14; i++ ));
do
for (( k = 1; k <= 10; k ++ ));
do
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N AMULET_${j}_${i}_${k}
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -j oe

samtools sort $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.bam -o $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.srt.bam

rm $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.bam

java -jar /projects/ps-renlab/yuw044/apps/AMULET-v1.1_0124/snATACOverlapCounter.jar --forcesorted --bcidx 1 --cellidx 1 --iscellidx 2 $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.srt.bam $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.singlecell.csv ~/projects/CEMBA/practice_doublet_removal/AMULET/testing_10D/input_data/mouse_autosomes.txt $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/ &> ${cemba_mop}/log/AMULET_${j}_${i}_${k}_find_overlaps.log

python3 /projects/ps-renlab/yuw044/apps/AMULET-v1.1_0124/AMULET.py --rfilter ~/projects/CEMBA/practice_doublet_removal/AMULET/testing_10D/input_data/mm10.blacklist.bed $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/Overlaps.txt $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/OverlapSummary.txt $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/ &> ${cemba_mop}/log/AMULET_${j}_${i}_${k}_multiplet_detection.log

" | sed '1d' > $qsub/AMULET_${j}_${i}_${k}.qsub

qsub $qsub/AMULET_${j}_${i}_${k}.qsub -o $cemba_mop/log/AMULET_${j}_${i}_${k}.log -e $cemba_mop/log/AMULET_${j}_${i}_${k}.err
done
done
done
```

Prepare Scrublet pipeline input:

```
for j in rep1 rep2;
do
for (( i = 11; i <= 14; i++ ));
do
for (( k = 1; k <= 10; k ++ ));
do
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N bamtosnap_${i}
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -j oe

# sort the simulated bam file(s) by query name
samtools sort -n -@ 10 -m 1G $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.bam -o $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.nsrt.bam

# convert the nsrt bam file(s) to snap file(s)
snaptools snap-pre --input-file=$cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.nsrt.bam --output-snap=$cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.snap --genome-name=mm10 --genome-size=$gold/CEMBA180104_4B_old/simulated_datasets2/mm10.chrom.sizes --min-mapq=30 --min-flen=50 --max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=False --overwrite=True --max-num=20000 --min-cov=500 --verbose=True

# add cell by bin matrix to the snap object
snaptools snap-add-bmat \
    --snap-file=$cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.snap \
    --bin-size-list 1000 5000 10000 \
    --verbose=True

" | sed '1d' > $qsub/bamtosnap_${j}_${i}_${k}.qsub

qsub $qsub/bamtosnap_${j}_${i}_${k}.qsub -o $cemba_mop/log/bamtosnap_${j}_${i}_${k}.log -e $cemba_mop/log/bamtosnap_${j}_${i}_${k}.err
done
done
done
```

Add cell-by-gene matrix to the snap file

```
for j in rep1 rep2;
do
for (( i = 11; i <= 14; i++ ));
do
for (( k = 1; k <= 10; k ++ ));
do
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N add_GM_$j_$i_$k
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -j oe

/home/yangli1/apps/anaconda2/bin/snaptools snap-add-gmat \
    --snap-file=$cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.snap \
    --gene-file /projects/ps-renlab/yangli/genome/mm10/gencode.vM16.geneUp2k.bed \
    --verbose=True

" | sed '1d' > $qsub/add_GM_${j}_${i}_${k}.qsub

qsub $qsub/add_GM_${j}_${i}_${k}.qsub -o $cemba_mop/log/add_GM_${j}_${i}_${k}.log -e $cemba_mop/log/add_GM_${j}_${i}_${k}.err
done
done
done
```

Generate a proper stat.txt file for each of the two 2*40 simulated datasets (first combine, then append doublets)

```
head -n 1 $archive/4B/CEMBA180104_4B/processed/stat.txt > $cemba_mop/rep1/stat.txt
head -n 1 $archive/4B/CEMBA180104_4B/processed/stat.txt > $cemba_mop/rep2/stat.txt
for i in `cat $cemba_mop/rep1/sample.lst`;
do
j=${i##*_}
cat $archive/$j/$i/processed/stat.txt | sed '1d' | awk -e '{printf "'$i'.%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6}' >> $cemba_mop/rep1/stat.txt
done
for i in `cat $cemba_mop/rep2/sample.lst`;
do
j=${i##*_}
cat $archive/$j/$i/processed/stat.txt | sed '1d' | awk -e '{printf "'$i'.%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6}' >> $cemba_mop/rep2/stat.txt
done
```

Run the Scrublet pipeline on each of the two 2*40 simulated datasets

```
for j in rep1 rep2;
do
for (( i = 11; i <= 14; i++ ));
do
for (( k = 1; k <= 10; k ++ ));
do
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N Scrublet_${j}_${i}_${k}
#PBS -l nodes=1:ppn=1,walltime=5:00:00
#PBS -V
#PBS -M yuw044@ucsd.edu
#PBS -m a
#PBS -A ren-group
#PBS -j oe

python $bin/generate.stat.txt.py $cemba_mop/${j}/stat.txt $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/ground.truth.tsv

Rscript ${path1}/bin/snapATAC.qc.filter.R -i $cemba_mop/simulated_datasets_$j/dataset${i}_${k}/simulated.snap \
                                        --tsse2depth $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/stat.txt \
                                        -o $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/simulated \
                                        --fragment_num 1 --tsse_cutoff 1 &> $cemba_mop/log/qc.filter_${j}_${i}_${k}.log

Rscript $path2script/snapATAC.rmDoublets.R -i $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/simulated.qc.filter.RData \
                                                      --mat gmat --rate 0.08 \
                                                      -o $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/simulated.gmat &> $cemba_mop/log/gmat.rmDoublets_${j}_${i}_${k}.log; 

Rscript $path2script/snapATAC.fitDoublets.R -i $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/simulated.gmat.qc.cluster.RData \
                                            -m $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/simulated.gmat.rmDoublets.txt \
                                            -o $cemba_mop/simulated_datasets_${j}/dataset${i}_${k}/simulated.gmat &> $cemba_mop/log/gmat.fitDoublets_${j}_${i}_${k}.log

" | sed '1d' > $qsub/scrublet_${j}_${i}_${k}.qsub

qsub $qsub/scrublet_${j}_${i}_${k}.qsub -o $cemba_mop/log/scrublet_${j}_${i}_${k}.log 
done
done
done
```

### Step 4. Quantify the performance of the two tools using PRC and AUPRC

#### For replications 1 and 2 of regions 2C, 3C, 4B, 5D (no UQ cutoff)

#### For replications 1 and 2 of regions 2C, 3C, 4B, 5D (UQ between 1000 and 4000)

#### For the combined MOP datasets replications 1 and 2
