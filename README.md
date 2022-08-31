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
- [Step 2. Simulate doublets](#step-2-simulate-doublets)
- [Step 3. Remove doublets on simulated datasets](#step-3-remove-doublets-on-simulated-datasets)
- [Step 4. Compare between the two tools by PRC and AUPRC](#step-4-compare-between-the-two-tools-by-prc-and-auprc)

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
Please refer to [AMULET](https://github.com/UcarLab/AMULET) to see how to interpret the results and identify singlets/doublets.

#### Identify Scrublet doublets

Start with a snap file object (used by Ren lab). See [here for conversion between a bam file and a snap file](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap)



### Step 2. Simulate doublets

and introducing artificial doublets by randomly selecting 2/11 of nuclei in the dataset and forming nuclei pairs by adding their read count profiles together (repeated 10 times per dataset). 

### Step 3. Remove doublets on simulated datasets



### Step 4. Compare between the two tools by PRC and AUPRC
