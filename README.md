# Comparison between AMULET & Scrublet on Doublet Removal for CEMBA2.0 
Comparison of the performance of [AMULET](https://github.com/UcarLab/AMULET) and a [modified](https://www.nature.com/articles/s41586-021-03604-1#Sec15) [Scrublet](https://github.com/swolock/scrublet) pipeline on doublet removal in the CEMBA2.0 project from [Dr.Bing Ren's lab](http://renlab.sdsc.edu/renlab_website/bing/).

**Publication**: the manuscript will be submitted soon 
## About
In the whole project, we performed single-nucleus (sn)ATAC-seq assays using single-cell combinatorial indexing (sci)ATAC-seq for more than 2,300,000 cells from 117 dissected regions in the adult mouse brain to produce comprehensive maps of cCREs in distinct cerebral cell types.

There are limited number of doublet removal tools developed specifically for snATAC-seq. Up to now, there are [SnapATAC](https://github.com/r3fang/SnapATAC), [ArchR](https://github.com/GreenleafLab/ArchR), and [AMULET](https://github.com/UcarLab/AMULET) for snATAC-seq doublet removal to the best of our knowledge. Researchers generally adopt algorithms from snRNA-seq doublet removal tools to achieve the same goal in snATAC-seq data.

We modified [Scrublet](https://github.com/swolock/scrublet), a simulation-based doublet removal tool for snRNA-seq data, and applied it on our datasets. It turned out the its performance was basically satisfactory on our [previous dataset](https://www.nature.com/articles/s41586-021-03604-1#Sec15). 

In 2021, [AMULET](https://github.com/UcarLab/AMULET), a read-count based doublet removal tool for snATAC-seq was published on Genome Biology. We believed that this novel tool could potentially outperform existing simulation-based packages including Scrublet and ArchR. As we hope to remove doublets more accurately, we compared the performance between [AMULET](https://github.com/UcarLab/AMULET) and our [modified](https://www.nature.com/articles/s41586-021-03604-1#Sec15) [Scrublet](https://github.com/swolock/scrublet) pipeline.

**Our results show that [AMULET](https://github.com/UcarLab/AMULET) does not outperform the [modified](https://www.nature.com/articles/s41586-021-03604-1#Sec15) [Scrublet](https://github.com/swolock/scrublet) pipeline on our typical datasets. The reason is that AMULET does a best job on datasets sequenced deeply, like 25K median valid read pairs per nucleus, but that for our datasets is typically 3K~4K.**

## Workflow
In summary, we use our datasets to generate "simulated datasets" in which we know to our best capability if a barcode represents a singlet or a doublet. Then we apply both of the tools to each of the simulated datasets. We quantify the performance of a tool by plotting its precision-recall curve for distinguishing between singlets and doublets. We keep the example dataset in this repo is very small just as an demo, so the results do not really mean something.

- [Step 1. Generate a singlet pool](#step-1-generate-a-singlet-pool)
- [Step 2. Simulatte doublets](#step-2-simulate-doublets)
- [Step 3. Remove doublets on simulated datasets](#step-3-remove-doublets-on-simulated-datasets)
- [Step 4. Compare between the two tools by PRC and AUPRC](#step-4-quantify-the-performance-of-the-two-tools-using-prc-and-auprc)

### Step 1. Generate a singlet pool

We first obtain a pool of singlets applying AMULET and the modified Scrublet pipeline on a raw dataset and keeping the common singlets identified by both tools at their default level (excluding the union set of doublets called by both tools). Here we make an assumption that the common singlets are real singlets without doublets.

#### Identify AMULET doublets

To apply AMULET on a bam file, you need to add an CB attribute (barcode) to each record of the bam file as AMULET identifies which cell a read belongs to by the CB attribute. Usually bam files from 10X Genomics already have the CB attribute, but ours come from a different platforms. AMULET only accepts coordinate-sorted and works on deduplicated bam file. 
```
# add CB attribute and sort the deduplicated bam file
cat <(samtools view -H $temp/data/singlet_pool_generation/AMULET/example.dedup.bam) \
<(samtools view $temp/data/singlet_pool_generation/AMULET/example.dedup.bam | \
awk 'BEGIN{FS = ":"} {CB = $1} {printf "%s\t%s%s\n", $0, "CB:Z:", CB}') \
| samtools sort -o $temp/data/singlet_pool_generation/AMULET/example.dedup.srt.bam
```

AMULET needs a barcode to cell_id map in CSV format (e.g., singlecell.csv from CellRanger). We adopt it from our [nuclei metatable](http://renlab.sdsc.edu/yangli/downloads/mousebrain/Supplementary_Tables/Supplementary%20Table%202%20-%20Metatable%20of%20nuclei.tsv.zip).
```
cat $yourpath1/metatable.tsv | awk '($2 == "CEMBA180104_4B"){printf "%s,%s\n", $3, "1"}' \
> $temp/data/singlet_pool_generation/AMULET/singlecell.csv
```

You also need a text file (included here) documenting all the chromosomes for the species you are analyzing. Then, you are ready to apply AMULET.

Apply the first step of AMULET:
```
java -jar $temp/script/AMULET-v1.1_0124/snATACOverlapCounter.jar --forcesorted --iscellidx 1 \
$temp/data/singlet_pool_generation/AMULET/example.dedup.srt.bam $temp/data/singlet_pool_generation/AMULET/singlecell.csv \
$temp/data/singlet_pool_generation/AMULET/mouse_autosomes.txt $temp/data/singlet_pool_generation/AMULET/ 
```
Apply the second step of AMULET using two resulting files from the first step:
```
python3 $temp/script/AMULET-v1.1_0124/AMULET.py --rfilter $temp/data/singlet_pool_generation/AMULET/mm10.blacklist.bed \
$temp/data/singlet_pool_generation/AMULET/Overlaps.txt $temp/data/singlet_pool_generation/AMULET/OverlapSummary.txt \
$temp/data/singlet_pool_generation/AMULET/
```
Please refer to [AMULET](https://github.com/UcarLab/AMULET) to see how to interpret the results and identify singlets/doublets. Briefly, a ```MultipletProbabilities.txt``` will be generated containing a q-value for each barcode being a singlet. The smaller the q-value is, the more likely the barcode is a doublet.

#### Identify Scrublet doublets

Start with a snap file object (used by Ren lab). See [here for conversion between a bam file and a snap file](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap)

*The following scripts are not disclosed. Please contact us for futher information.*

First convert the snap file into an R data file and perform QC, filtering out cells with less than 1000 uniquely mapped fragments or tssc less than 10.
```
Rscript /projects/ps-renlab/yangli/projects/CEMBA/00.data/szu/bin/snapATAC.qc.filter.R \
-i /projects/ps-renlab/yangli/projects/CEMBA/00.data/4B/CEMBA180104_4B/CEMBA180104_4B.snap \
--tsse2depth /projects/ren-transposon/home/yangli/projects/CEMBA/00.data/archive/4B/CEMBA180104_4B/processed/stat.txt \
-o /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B --fragment_num 1000 --tsse_cutoff 10
```

Then apply the modified Scrublet pipeline.
```
Rscript /projects/ps-renlab/yangli/scripts/snATACutils/bin/snapATAC.rmDoublets.R \
-i /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B.qc.filter.RData \
--mat gmat --rate 0.08 
-o /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B.gmat
Rscript /projects/ps-renlab/yangli/scripts/snATACutils/bin/snapATAC.fitDoublets.R \
-i /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B.gmat.qc.cluster.RData \
-m /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B.gmat.rmDoublets.txt \
-o /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B.gmat
```

This would result in a ".gmat.fitDoublets.txt" file in which each barcode is assigned a doublet score. The higher the score is, the more likely the barcode is to be a doublet.

Last step is to find common singlets identified by both tools, i.e. excluding the union set of doublets found by both tools.

```
python $temp/script/find_common_singlets.py \
--meta $temp/data/singlet_pool_generation/Scrublet/$i.qc.filter.meta.txt \
--scrublet $temp/data/singlet_pool_generation/Scrublet/$i.gmat.fitDoublets.txt \
--AMULET $temp/data/singlet_pool_generation/AMULET/MultipletBarcodes_01.txt \
--output $temp/data/singlet_pool_generation/SingletBarcodes_01.txt
```

The resulting file "SingletBarcodes_01.txt" will serve as the input to step 2.

### Step 2. Simulate doublets

We introduce artificial doublets by randomly picking 2/11 of the nuclei in the dataset and forming nucleus pairs by adding their read count profiles together (repeated 10 times per dataset and thus generating 10 simulated dataset replicates). This will result in a doublet proportion of 10% in each simulated dataset. 

To generate such datasets, use ```$temp/script/simulate_cemba_doublets.py``` and provide the following input arguments:

```BAM``` The path to the parent bam file

```ANNOTATION``` The path to the metatable.tsv file (or other tsv annotation file)

```CELLTYPE``` The column name in the annotation file used to identify cell type; used for doublet type annotation, homotypic/heterotypic

```SINGLETS``` The path to the ```SingletBarcodes_01.txt``` file

```SUMMARY``` The path to the ```OverlapSummary.txt``` file given by running AMULET on the parent bam file

```LOCATION``` The location to store the simulated dataset (bam file) and the ground truth file of artificial doublet barcodes

Example:
```
for ((k = 1; k <= 10; k ++))
do
python $temp/script/simulate_cemba_doublets.py --poolsize 1100 \
$temp/data/singlet_pool_generation/AMULET/example.dedup.srt.bam \
$yourpath1/metatable.tsv MajorType $temp/data/singlet_pool_generation/SingletBarcodes_01.txt \
$temp/data/singlet_pool_generation/AMULET/OverlapSummary.txt \
$temp/data/simulated_datasets/dataset${k}
done
```

Optional arguments:

```--fclower``` An integer number specifying the fragment count lower limit of any singlet going to the singlet pool (default: 0)

```--fchigher``` An integer number specifying the fragment count upper limit of any singlet going to the singlet pool (default: 100000000)

```--poolsize``` An integer specifying the number of singlets to keep in the pool used for downsampling (default: 100000)

```--proportion``` A float specifying the proportion of doublets in a simulated dataset (default: 0.1)

This step would generate a bam file of the simulated dataset containing singlets and artificial doublets, a ground truth tsv file, and ```singlecell.csv``` file for AMULET input. An example of the ground truth file can be seen here ```$temp/data/simulated_datasets/dataset1/ground.truth.tsv```.

### Step 3. Remove doublets on simulated datasets

First run AMULET on each of the simulated datasets

```
# the first step of AMULET
java -jar $temp/script/AMULET-v1.1_0124/snATACOverlapCounter.jar --forcesorted \
--bcidx 1 --cellidx 1 --iscellidx 2 $temp/data/simulated_datasets/dataset1/simulated.bam \
$temp/data/simulated_datasets/dataset1/singlecell.csv \
$temp/data/singlet_pool_generation/AMULET/mouse_autosomes.txt \
$temp/data/simulated_datasets/dataset1/ 

# the second step of AMULET
python3 $temp/script/AMULET-v1.1_0124/AMULET.py --rfilter \
$temp/data/singlet_pool_generation/AMULET/mm10.blacklist.bed \
$temp/data/simulated_datasets/dataset1/Overlaps.txt \
$temp/data/simulated_datasets/dataset1/OverlapSummary.txt \
$temp/data/simulated_datasets/dataset1/
```

This step would generate the ```MultipletProbabilities.txt``` file containing a q-value for each barcode being a singlet. By comparing the data in this file and the ground truth file, we can clearly calculate AMULET's recall and precision at all possible q-value thresholds.

Then run the modified Scrublet pipeline on each of the simulated datasets. To run the pipeline, first convert each bam file containing a simulated dataset to a snap file.

```
# sort the bam file based on query name
samtools sort -n -@ 10 -m 1G $temp/data/simulated_datasets/dataset1/simulated.bam \
-o $temp/data/simulated_datasets/dataset1/simulated.nsrt.bam

# convert the bam file to a snap file
snaptools snap-pre --input-file=$temp/data/simulated_datasets/dataset1/simulated.nsrt.bam \
--output-snap=$temp/data/simulated_datasets/dataset1/simulated.dedup.snap --genome-name=mm10 \
--genome-size=$temp/data/simulated_datasets/mm10.chrom.sizes --min-mapq=30 --min-flen=50 \ 
--max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=False --overwrite=True \
--max-num=20000 --min-cov=500 --verbose=True

# add cell-by-bin matrix to the snap object
snaptools snap-add-bmat --snap-file=$temp/data/simulated_datasets/dataset1/simulated.dedup.snap \
--bin-size-list 1000 5000 10000   --verbose=True

# add cell-by-gene matrix to the snap object
snaptools snap-add-gmat --snap-file=$temp/data/simulated_datasets/dataset1/simulated.dedup.snap \
--gene-file $temp/data/simulated_datasets/gencode.vM16.geneUp2k.bed    --verbose=True
```

The modified Scrublet pipeline needs a ```stat.txt``` file input in which each barcode takes a line. The pipeline collect the fourth column "TSSe" of each barcode and ignores the remaining. Run ```$temp/script/generate.stat.txt.py``` to generate a new ```stat.txt``` file from the parent ```stat.txt``` file. This script assigns a random TSSE between 10~15 to each artificial doublet. Provide the following input arguments:

```STATTXT``` The path to the parent ```stat.txt``` file
```GROUNDTRUTH``` The path to the ```ground.truth.tsv``` file of the simulated dataset

By default, the script will generate a new ```stat.txt``` file under the same directory as ```GROUNDTRUTH```.

Example:
```
python $temp/script/generate.stat.txt.py $yourpath2/stat.txt \
$temp/data/simulated_datasets/dataset1/ground.truth.tsv
```
Optional Arguments:

```--location``` The path to the directory that stores the resulting stat.txt file

Then you are ready to run the modified Scrublet pipeline on each of the simulated datasets.

First convert the snap file into an R data file, this time without doing QC as barcodes in the simulated datasets all have passed QC or are artificial doublets
```
Rscript /projects/ps-renlab/yangli/projects/CEMBA/00.data/szu/bin/snapATAC.qc.filter.R \
-i /projects/ps-renlab/yangli/projects/CEMBA/00.data/4B/CEMBA180104_4B/CEMBA180104_4B.snap \
--tsse2depth /projects/ren-transposon/home/yangli/projects/CEMBA/00.data/archive/4B/CEMBA180104_4B/processed/stat.txt \
-o /projects/ps-renlab/yuw044/projects/CEMBA/practice_doublet_removal/compare/gold_standard/CEMBA180104_4B/doublet_removal/Scrublet/CEMBA180104_4B --fragment_num 0 --tsse_cutoff 1
```

Follow the same workflow to [apply the modified scrublet pipeline](#identify-scrublet-doublets). By comparing the resulting file ```XXX.gmat.fitDoublets.txt``` with ```ground.truth.tsv``` file, we can clearly calculate Scrublet's recall and precision at all possible q-value thresholds.

### Step 4. Quantify the performance of the two tools using PRC and AUPRC

We use the ggplot package in R to make PRCs and plot area under PRCs as bar/box plots. 

#### PRC

Run ```$temp/script/calculate_prc_coordinates.py``` to calculate the coordinates for the precision-recall curves. This will result in a [feather file](https://www.rstudio.com/blog/feather/) storing the coordinate information serving as an input to a following R script. Provide the following arguments:

```DATASETS``` The path to a set of simulated datasets (replicates of each other)

Example: 
```
for i in `cat $temp/data/simulated_datasets/cemba.mop.sample.lst` 
do
echo ${i}
python $temp/script/calculate_prc_coordinates.py --end 10 $temp/data/simulated_datasets/${i}/dataset \
-f $temp/data/simulated_datasets/${i}/${i}_no_cutoff.feather
done
```


