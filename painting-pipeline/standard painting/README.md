# Pipeline to paint a target dataset with a reference panel without reference correction

Here we describe the pipeline for painting a bio-bank scale target dataset using a reference dataset on a single chromosome through [SparsePainter](https://github.com/YaolingYang/SparsePainter). Below are the files we have before painting:  

-   ``chr1_ref.vcf.gz``: The phased, non-missing reference dataset.

-   ``chr1_target.vcf.gz``: The phased, non-missing target dataset. Assume this dataset has 500,000 individuals and 50,000 SNPs.

-   ``chr1_map.txt``: The genetic map file as required by SparsePainter.

-   ``popnames.txt``: The population file of reference individuals as required by SparsePainter.

-   ``names.txt``: The names of all the target individuals.

## Step1: Split the input file into subfiles

To paint huge target samples, it is suggested to split the huge target file into small target subfiles, then we can submit multiple small jobs to run efficiently with HPC. It is usually more efficient to split phase files than vcf files, so we first convert vcf to phase files. All the commands below, unless specifically stated, are run on **Linux bash shell**.

```
pbwt -readVcfGT chr1_ref.vcf.gz -writePhase chr1_ref.phase  
pbwt -readVcfGT chr1_target.vcf.gz -writePhase chr1_target.phase
```

Then we create a folder ‘chr1split’ to store the subfiles.

```mkdir chr1split```

Now we run below commands in **Python** to split target files into 500 subfiles, each subfile contains 1000 individuals. Note that this is an example code, more efficient codes for splitting files are possible.

```
import os
chr = '1'
input_file = f"chr{chr}_target.phase"
total_files = 500
lines_per_file = 2000

def split_file(input_file, total_files, lines_per_file, chr):
    with open(input_file, 'r') as file:
       next(file)  # Skip the first line
       second_line = next(file)
       third_line = next(file)

    with open(input_file, 'r') as file:
        for _ in range(3):
            next(file)  

        file_number = 1
        line_count = 0
        output_file = gzip.open(f"chr{chr}split/chr{chr}_target{file_number}.phase", 'wt')
        output_file.write(f"{lines_per_file}\n")  
        output_file.write(second_line)           
        output_file.write(third_line)             

        for line in file:
            if line_count == lines_per_file:
                output_file.close()
                file_number += 1
                line_count = 0
                output_file = gzip.open(f"chr{chr}split/chr{chr}_target{file_number}.phase", 'wt')
                output_file.write(f"{lines_per_file}\n")  
                output_file.write(second_line)
                output_file.write(third_line)

            output_file.write(line)
            line_count += 1

        output_file.close()

split_file(input_file, total_files, lines_per_file, chr)
```

To save storage space, we could compress the reference file, and remove ``chr1_target.phase`` which has already been split into subfiles.

```
gzip chr1_ref.phase  
rm chr1_target.phase
```

Also, we should update the name files for the individuals in each subfile. These files are stored in the folder ``namefile``.  

```
input_file="names.txt"
output_dir="namefile"
mkdir -p "$output_dir"

lines_per_file=1000
file_count=1
line_count=0

while IFS= read -r line; do
    if [ $line_count -eq 0 ]; then
        output_file="${output_dir}/target${file_count}.txt"
        > "$output_file"
    fi

    echo "$line" >> "$output_file"

    ((line_count++))

    if [ $line_count -eq $lines_per_file ]; then
        ((file_count++))
        line_count=0
        echo $file_count
    fi
done < "$input_file"
```

Now we have finished data preprocessing. ``chr1_target.vcf.gz`` has been split into 1000 subfiles in phase format in folder ``chr1split``; ``names.txt`` has been split into 1000 subfiles in folder ``namefile``. ``chr1_ref.vcf.gz`` has been converted to ``chr1_ref.phase.gz``. Then the painting starts.

## Step2: Perform painting for all the target individuals

Before we paint all the target individuals, we need to determine the recombination scaling constant lambda, and use the fixed lambda to paint all the individuals. To estimate lambda, we only need to paint one target subset:

```
mkdir chr1  
 ./SparsePainter -reffile chr1_ref.phase.gz -targetfile chr1split/chr1_target1.phase -popfile popnames.txt -mapfile chr1_map.txt -namefile namefile/target1.txt -indfrac 1 -prob -chunklength -chunkcount -probstore linear -out chr1/chr1_target1
```

This generates ``chr1/chr1_target1_fixlambda.txt`` and other files (described below). Assume the estimated recombination scaling constant is 100 from ``chr1/chr1_target1_fixlambda.txt``, then we use this fixed lambda to paint all the other subfiles (and other chromosomes if it applies). We usually submit the remaining 499 array jobs on HPC. Let ``SLURM_ARRAY_TASK_ID`` denote the array task IDs, then we run the below command to paint all the subfiles:

```
 ./SparsePainter -reffile chr1_ref.phase.gz -targetfile chr1split/chr1_target${SLURM_ARRAY_TASK_ID}.phase -popfile popnames.txt -mapfile chr1_map.txt -namefile namefile/target${SLURM_ARRAY_TASK_ID}.txt -fixlambda 100 -prob -chunklength -chunkcount -probstore linear -out chr1/chr1_target${SLURM_ARRAY_TASK_ID}
```

It is optional to paint with ``-LDAS -AAS -aveSNP -aveind``, etc.

Now we finish the painting, and get the below output files for each target subfile indexed ``SLURM_ARRAY_TASK_ID``.

``chr1/chr1_target${SLURM_ARRAY_TASK_ID}_prob.txt.gz``: The local ancestry probabilities for each target sample at each SNP stored in linear form (see SparsePainter manual).
``chr1/chr1_target${SLURM_ARRAY_TASK_ID}_chunklength.txt.gz``: The expected length (in centiMorgan) of copied chunks of each local ancestry for each target sample.
``chr1/chr1_target${SLURM_ARRAY_TASK_ID}_chunkcount.txt.gz``: The expected number of copied chunks of each local ancestry for each target sample.

## Step3: Post-processing output files

Please follow the [instructions](https://github.com/YaolingYang/SparsePainter) to merge these subfiles, extract paintings for analysis, etc. Below we provide an example in **Python** to merge the probability (``chr1_target_prob.txt.gz``), chunk length (``chr1_target_chunklength.txt.gz``) and chunk count (``chr1_target_chunkcount.txt.gz``) files.

```
import gzip
import os
import sys

chr = '1'
total_files = 500

#Function to merge files
def merge_probfiles(file_type, output_file_name):
    with gzip.open(output_file_name, 'wt') as output_file:
        for i in range(1, total_files+1):
            print(f"Processing file {i} for {file_type}")
            full_path = f"chr{chr}/chr{chr}_target{i}_{file_type}.txt.gz"

            if os.path.isfile(full_path):
                with gzip.open(full_path, 'rt') as input_file:
                    if i != 1:
                        next(input_file, None)
                        next(input_file, None)
                    for line in input_file:
                        output_file.write(line)
            else:
                print(f"File not found: {full_path}")

def merge_chunkfiles(file_type, output_file_name):
    with gzip.open(output_file_name, 'wt') as output_file:
        for i in range(1, total_files+1):
            print(f"Processing file {i} for {file_type}")
            full_path = f"chr{chr}/chr{chr}_target{i}_{file_type}.txt.gz"

            if os.path.isfile(full_path):
                with gzip.open(full_path, 'rt') as input_file:
                    if i != 1:
                        next(input_file, None)

                    for line in input_file:
                        output_file.write(line)
            else:
                print(f"File not found: {full_path}")

#Merge 'prob' files
merge_probfiles('prob', f"chr{chr}_target_prob.txt.gz")

#Merge 'chunklength' files
merge_chunkfiles('chunklength', f"chr{chr}_target_chunklength.txt.gz")

#Merge 'chunkcount' files
merge_chunkfiles('chunkcount', f"chr{chr}_target_chunkcount.txt.gz")
```
