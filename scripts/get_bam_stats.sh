#!/bin/bash
# This script will run samtools to extract BAM header information
# Uses manifest file from Portal as input to create directory structure and to access mounted files

MANIFEST=""

while getopts "m:" OPTION;
do
   case "$OPTION" in
      m)
         MANIFEST=${OPTARG?}
      ;;
      :)
         echo "Option -$OPTION requires the name of the manifest file" >&2
         exit 1
         ;;
   esac
done

OUT="completed_ids.tsv"
{
read
#while IFS=$'\t' read -r -a data
sed 's/\t/;/g'  | while IFS=';' read -r -a data
do
   file_id=${data[1]}
   echo ${file_id}
   object_id=${data[2]}
   donor_id=${data[8]}
   project_id=${data[9]}
   bundle_id=${data[11]}
   filename=${data[4]}
   # create directory structure to store bam header information
   if [[ ! -e "bam_data/$project_id" ]]; then
      mkdir "bam_data/$project_id"
   fi
   if [[ ! -e "bam_data/$project_id/$donor_id" ]];  then
      mkdir "bam_data/$project_id/$donor_id"
   fi
   if [[ ! -e "bam_data/$project_id/$donor_id/$file_id" ]]; then
      mkdir "bam_data/$project_id/$donor_id/$file_id"
   fi
   output="bam_data/${project_id}/${donor_id}/${file_id}/${object_id}.header"
   echo "donor_id = ${donor_id}" 
   echo "project_id = ${project_id}"
   echo "bundle_id = ${bundle_id}"
   echo "object_id = ${object_id}"
   echo `samtools view -H /mnt/${bundle_id}/${filename} 1>${output} 2>>stderr.log`
   echo "samtools view -H /mnt/${bundle_id}/${filename} 1>${output} 2>>stderr.log"
   echo "${project_id}	${donor_id}	${object_id}	${bundle_id}	${filename}	${file_id}"
   echo "${project_id}	${donor_id}	${object_id}	${bundle_id}	${filename}	${file_id}" >> ${OUT}
done } < ${MANIFEST}
