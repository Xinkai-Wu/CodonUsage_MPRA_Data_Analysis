# Preprocessing of Massive Parallel Report Assays Sequencing Data

## Introduction
The purpose of this script is to calculate the distribution of each sequence variant in each component from the original sequencing data. Please refer to our article for specific experimental design and description.

## Processing Steps

### Step 1. Localization of variant sequences from raw sequencing data.

``` bash
for gzfile in $(ls *gz)
do
name=$(echo $gzfile|awk -F '.' '{print $1}')
gunzip -c $gzfile | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $name.fa
done

# Region 1
for sample in $(cat ./Data/samples.txt|grep 'NA1|hp1')
do
for file in $(ls $sample*R1*)
do
name=$(echo $file|awk -F 'R1' '{print $1}')
tail=$(echo $file|awk -F 'R1' '{print $2}')
R2File=$(echo -e ${name}R2$tail|awk '{print $1}')
python ./Script/CodonUsageBias_01dataSplit.py $file $R2File $file 1 > $file.log
done
done

# Region 2
for sample in $(cat ./Data/samples.txt|grep 'NA2|hp2')
do
for file in $(ls $sample*R1*)
do
name=$(echo $file|awk -F 'R1' '{print $1}')
tail=$(echo $file|awk -F 'R1' '{print $2}')
R2File=$(echo -e ${name}R2$tail|awk '{print $1}')
python ./Script/CodonUsageBias_01dataSplit_region2.py $file $R2File $file 1 > $file.log
done
done
```
\# Note: When processing sequencing data, it is recommended to divide the fastq file into several small files and use multiple threads for processing.

### Step 2. Using barcode to distinguish libraries from different components.

``` bash
for sample in $(cat ./Data/samples.txt|grep 'NA1|hp1')
do
python ./Script/rawSeqR1Process.py ${sample} region1 > $sample.umi.log
done
for sample in $(cat ../samples.txt|grep 'NA2|hp2')
do
python ./Script/rawSeqR1Process.py ${sample} region2 > $sample.umi.log
done
```

### Step 3. The counting results of experimental libraries were integrated.

Please see the **Step3_DataCombination.ipynb** for details.

### Step 4. Count normalization and calculation of RiboScore and GFPScore.
The sequencing depths of all libraries were normalized using DESeq2. 
The number of ribosomes bound and the efficiency of GFP synthesis were used to quantify the transltaion efficiency.

Please see the **Step4.1_RiboScore_Calculate.ipynb** and **Step4.2_GFPScore_Calculate.ipynb** for details.

## License
Each file included in this repository is licensed under the MIT License.
