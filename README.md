# Comparison between AMULET & Scrublet on Doublet Removal for CEMBA2.0 
Comparison of the performance of [AMULET](https://github.com/UcarLab/AMULET) and a [modified](https://www.nature.com/articles/s41586-021-03604-1#Sec16) [Scrublet](https://github.com/swolock/scrublet) pipeline on doublet removal in the CEMBA2.0 project from [Dr.Bing Ren's lab](http://renlab.sdsc.edu/renlab_website/bing/).

**Publication**: the manuscript will be submitted soon 
## About
In the whole project, we performed single-nucleus (sn)ATAC-seq assays using single-cell combinatorial indexing (sci)ATAC-seq for more than 2,300,000 cells from 117 dissected regions in the adult mouse brain to produce comprehensive maps of cCREs in distinct cerebral cell types.

There are limited number of doublet removal tools developed specifically for snATAC-seq. Up to now, there are [SnapATAC](https://github.com/r3fang/SnapATAC), [ArchR](https://github.com/GreenleafLab/ArchR), and [AMULET](https://github.com/UcarLab/AMULET) for snATAC-seq doublet removal to the best of our knowledge. Researchers generally have been adopting algorithms from snRNA-seq doublet removal tools to achieve the same goal in snATAC-seq data.

We modified [Scrublet](https://github.com/swolock/scrublet), a simulation-based doublet removal tool for snRNA-seq data, and applied it on our datasets. It turned out that its performance was basically satisfactory on our [previous datasets](https://www.nature.com/articles/s41586-021-03604-1#Sec16). 

In 2021, [AMULET](https://github.com/UcarLab/AMULET), a read-count based doublet removal tool for snATAC-seq was published on Genome Biology. We believed that this novel tool could potentially outperform existing simulation-based packages including Scrublet and ArchR. As we hope to remove doublets more accurately, we compared the performance between [AMULET](https://github.com/UcarLab/AMULET) and our [modified](https://www.nature.com/articles/s41586-021-03604-1#Sec16) [Scrublet](https://github.com/swolock/scrublet) pipeline.

**Our results show that [AMULET](https://github.com/UcarLab/AMULET) does not outperform the [modified](https://www.nature.com/articles/s41586-021-03604-1#Sec16) [Scrublet](https://github.com/swolock/scrublet) pipeline on our typical datasets. The reason is that AMULET does a best job on datasets sequenced deeply, e.g. median valid read pairs per nucleus equal to 25K, but that for our datasets is typically 3K~4K.**

## Workflow
In summary, we use our datasets to generate "simulated datasets" in which we know to our best capability if a barcode represents a singlet or a doublet (ground truth). Then we apply each tool respectively on each of the simulated datasets. Then the doublet removal results are compared with the ground truth. We quantify the performance of a tool by plotting its precision-recall curve for distinguishing between singlets and doublets. We keep the example dataset in this repo very small just as a demo.

- [Step 1. Generate a singlet pool](#step-1-generate-a-singlet-pool)
- [Step 2. Simulate datasets containing artificial doublets](#step-2-simulate-doublets)
- [Step 3. Remove doublets on simulated datasets](#step-3-remove-doublets-on-simulated-datasets)
- [Step 4. Compare between the two tools by PRC and AUPRC](#step-4-quantify-the-performance-of-the-two-tools-using-prc-and-auprc)

### Step 1. Generate a singlet pool

We first obtain a pool of singlets by applying AMULET and the modified Scrublet pipeline on a raw dataset and keeping common singlets identified by both tools at their default levels (excluding the union set of doublets called by both tools). Here we make an assumption that the common singlets are ground truth singlets, i.e. the union set of doublets identified by both tools contains ground truth doublets.

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
-o l --fragment_num 1000 --tsse_cutoff 10
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

We introduce artificial doublets by randomly picking 2/11 of the nuclei in the singlet pool, randomly pairing up 2 nuclei from the selected nuclei, and adding their read count profiles together for each pair (repeated 10 times per dataset and thus generating 10 simulated dataset replicates). This will result in a final doublet proportion of 10% in each simulated dataset. 

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

PRCs are in ```$temp/result```. Run ```$temp/script/calculate_prc_coordinates.py``` to calculate the coordinates for the precision-recall curves. This will result in a [feather file](https://www.rstudio.com/blog/feather/) storing the coordinate information serving as an input to a following R script. Basically a feather file stores a dataframe that can be read by both Python and R. Provide the following arguments:

```DATASETS``` The path to a set of simulated datasets (replicates of each other)
```FEATHERPATH``` The path of the resulting feather file

Example: 
```
for i in `cat $temp/data/simulated_datasets/cemba.mop.sample.lst` 
do
echo ${i}
python $temp/script/calculate_prc_coordinates.py --end 10 $temp/data/simulated_datasets/${i}/dataset \
$temp/data/simulated_datasets/${i}/${i}_no_cutoff.feather
done
```

Optional arguments:
```--end``` The number of replicates simulated based on a parent bam file (default: 9)

```--AMULETorder``` A negative integer representing the starting order of magnitude for q-value of AMULET results (default: -20). Sometimes AMULET finds barcodes extremely unlikely to be a singlet, and then this value needs to be adjusted down to probably -40.

```--scrhigh``` A float between zero and one representing the highest doublet score for Scrublet results (default: 0.8)

Then, make PRC plots by running ```$temp/script/make_prc_plot.R```. The following arguments are required:

```--dir``` the directory that contains the feather files

```--samples``` the path to the sample.lst file containing the sample name of each parent dataset

```--suffix``` the suffix of each feather file 

```-o``` the path to the output PRC figure

Example: 
```
# directly comparing the performance of AMULET and Scrublet on 8 sets of simulated datasets
Rscript $temp/script/make_prc_plot.R --labels $temp/data/simulated_datasets/simple.cemba.mop.sample.lst \
--dir $temp/data/simulated_datasets/feather_collection/ \
--samples $temp/data/simulated_datasets/cemba.mop.sample.lst \
--suffix _no_cutoff.feather -o $temp/result/performance/cemba.mop.prc.nocutoff.pdf

# illustrating why AMULET does not outperform Scrublet
Rscript $temp/script/make_prc_plot.R --label $temp/data/simulated_datasets/cemba_mop_combined/uq.lst \
--row 1 --dir $temp/data/simulated_datasets/cemba_mop_combined/simulated_datasets_rep1 \
--samples $temp/data/simulated_datasets/cemba_mop_combined/simulated_datasets_rep1/sample.lst \
--suffix .feather -o $temp/result/why_AMULET_not_good/cemba.mop.combined.rep1.prc.pdf
```

Optional Arugments:

```--labels``` the path to the label.lst file. If specified, will be used to label each PRC; otherwise, will label each PRC using ```--samples```

```-n``` the number of simulated datasets (replicates) in the feather file (default: 10)

```--row``` the number of rows in the resulting figure panel (default: 2)

```--column``` the number of columns in the resulting figure panel (default: 4)

#### AUPRC

AUPRC bar plots are in ```$temp/result```. Run ```$temp/script/calculate_auprc.py``` to calculate and store the area under the precision-recall curves in a text file. Provide the following arguments:

```DATASETS``` The path to a set of simulated datasets (replicates of each other)
```PREFIX``` the prefix of where to store the resulting text file (a .txt will be appended automatically to the end of the prefix)

Example:
```
# calculate auPRC for simulated datasets based on each CEMBA mop region 
for i in `cat $temp/data/simulated_datasets/cemba.mop.sample.lst`;
do 
echo ${i}
python $script/calculate_auprc.py $temp/data/simulated_datasets/${i}/dataset \
$temp/data/simulated_datasets/${i}/prc
done
```

Optional Arguments:
```--end``` The number of replicates simulated based on a parent bam file (default: 9)

```--AMULETorder``` A negative integer representing the starting order of magnitude for q-value of AMULET results (default: -20). Sometimes AMULET finds barcodes extremely unlikely to be a singlet, and then this value needs to be adjusted down to probably -40.

```--scrhigh``` A float between zero and one representing the highest doublet score for Scrublet results (default: 0.8)

Then combine the prc.auc.txt for simulated datasets from each CEMBA mop region. 

```
for i in `cat $temp/data/simulated_datasets/cemba.mop.sample.lst`;
do 
echo $i
cat $temp/data/simulated_datasets/${i}/prc.auc.txt >> $temp/data/simulated_datasets/cemba.mop.prc.auc.txt
done
```

Then run ```$temp/script/auprc_to_tsv.py``` to convert the combined text file to a tsv file for a following R script. Provide the following arguments:

```AUCTXT``` the path to the combined ```prc.auc.txt``` file

```SAMPLES``` the path to the ```sample.lst``` file

```TOOLS``` the path to the ```tools.lst``` file

```OUTPUT``` the path of the output file

```--ggplot``` a "store truth" option that asks the script to produce a tsv file for R ggplot instead of for Excel

Example:

```
python $script/auprc_to_tsv.py --ggplot $temp/data/simulated_datasets/cemba.mop.prc.auc.txt \
$temp/data/simulated_datasets/cemba.mop.sample.lst $temp/data/simulated_datasets/tools.lst \
$temp/data/simulated_datasets/cemba.mop.prc.auc.ggplot.tsv
```

The output tsv would serve as an input to the following R ggplot script. Run ```$temp/plot_auprc.R``` to plot AUPRC as bar/box plots. The following arguments are required:

```-i```/```--auprc``` the path to the ggplot.tsv file where the first column is sample, the second column is tool, and the third column is auprc

```-o```/```--output``` output file prefix (a ".barplot.png" or ".boxplot.png" will be automatically appended)

Example:

```
Rscript $script/plot_auprc.R -i $temp/data/simulated_datasets/cemba.mop.prc.auc.ggplot.tsv \
--labels $temp/data/simulated_datasets/simple.sample.lst --sig --box --bar \
-o $temp/result/performance/cemba.mop.auprc.nocutoff
```

Optional Arguments:

```--labels``` the path to the label.lst file. If specified, will be used to label each PRC; otherwise, will label each PRC using ```--samples```

```--sig``` a "store true" option, if specified, add the significance label to the box/bar plot

```--box``` a "store true" option, if specified, generate a box plot

```--bar``` a "store true" option, if specified, generate a bar plot

```-f/--xtickfont``` the font size for xticks (default: 8)

The figures in the ```result``` directory have been converted pdf format.
