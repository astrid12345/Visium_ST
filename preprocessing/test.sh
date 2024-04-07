#! /bin/bash
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
# Author    : Astrid Alsema
# Date      : March 2021
# Datasets  : Fastq files from libraries generated with 10x Visium kit.
# Purpose   : Automated SpaceRanger script for the preproccesing of multiple samples in a dataset. 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

# Set variables
dir=$(pwd)
SpaceDir=/data/bcn/Pipelines/SpaceRanger/

# Prompt user for sample IDs
printf "\nWhat are the IDs of the sample names? (e.g.ST31 ST32 ST33 ST34 for the visium MS study) \n"
read samples
# example input 
# samples=ST31 ST32 ST33 ST34 ST35 ST36 ST37 ST38

# Prompt user for slide IDs
printf "\nWhat are the slide numbers? input every slide number once (e.g. V19T12-044 V19S23-115) \n"
read -a slide_nr
# example input 
# slide_nr=V19T12-044 V19S23-115

# Expand the array of slides. 4 capture areas per slides are assumed to be filled.
slides=($(for k in "${slide_nr[@]}"; do printf "$k%.0s " {1..4}; done))

printf "\nWhat are the capture areas (in the SAME ORDER as sample names)? (e.g. with 2 slides you need to type A1 B1 C1 D1 A1 B1 C1 D1 for the MS study) \n"
read -a areas
# example input 
# areas=A1 B1 C1 D1 A1 B1 C1 D1

# Choose reference genome based on species
printf "\nWhat is the species of you samples? (human/mouse)\n"
read species

if [ $species = "human" ]; then
    ref=${SpaceDir}refdata-cellranger-GRCh38-3.0.0
    elif [ $species = "mouse" ]; then
    ref=${SpaceDir}refdata-cellranger-mm10-3.0.0 
    fi 

# Prompt user for input directories
printf "\nName of the folder containing the fastqs to be analysed? \n"
read dir1
fastq_dir=${dir}/${dir1}

printf "\nName of the folder containing the images? \n"
read dir2
image_dir=${dir}/${dir2}

# Display input information to user
i=0
for sample in $samples; do
	echo sample: $sample 
	echo on slide: ${slides[i]}
	echo on capture area: ${areas[i]}
	echo with image: ${image_dir}/${sample}.jpg 
	echo with fastq: ${fastq_dir}/${sample}
	((i=i+1))
	done

# Confirm input with user
printf "\nOK, almost there. Is the info ABOVE correct? (yes/no) \n"

read decision	
if [  $decision = "yes" ]; then
		echo Great, you are ready to go!
	else
		echo Please check the path to your images and fastq files
	fi

# Prompt user to start SpaceRanger
printf "\nDo you want to start SpaceRanger? (yes/no) \n"
read decision2
mkdir output

if [ $decision2 = "yes" ]; then
  # Run SpaceRanger for each sample
	i=0
	for sample in $samples; do
		echo running pipeline for $sample 
		spaceranger count --id=${sample}-out \
		--transcriptome=${ref} \
		--fastqs=${fastq_dir}/${sample} \
		--sample=${sample} \
		--image=${image_dir}/${sample}.jpg \
		--slide=${slides[i]} \
		--localcores=70 \
		--localmem=155 \
		--area=${areas[i]} \
		--loupe-alignment=${image_dir}/${sample}.json

		# Move output to output directory
		mv -n ${sample}-out output
		((i=i+1))
		done
	elif [ "$decision2" = "no" ]; then
    printf "\nPipeline stopped.\n"
	fi

# Concatenate metrics summary files of each sample and append together in one file total_metrics_summary
find ${dir}/output/*-out/outs -name metrics_summary.csv -exec cat {} + > ${dir}/output/total_metrics_summary.txt
