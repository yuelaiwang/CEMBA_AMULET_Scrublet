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

### Step 1: Generate a pool of singlets

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

