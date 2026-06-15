#!/bin/bash
# Usage: bash level_analysis_adapt.sh <call_level_output_directory>
LEVEL_DIR=$1
if [ -z "$LEVEL_DIR" ]; then
    echo "Usage: $0 <directory containing subdirs with .all.levels.txt>"
    exit 1
fi

ALU_BED="/path/to/hg19_Alu.bed"
NONALU_REP_BED="/path/to/human_nonAlu.bed"

SAMPLES=($(find ${LEVEL_DIR} -maxdepth 2 -name "*.all.levels.txt" -type f | sort))
mkdir -p ${LEVEL_DIR}/tmp_dir

for file in "${SAMPLES[@]}"; do
    sample=$(basename ${file} .all.levels.txt)
    echo "Processing ${sample} ..."
    sed '1d' ${file} | awk -v OFS='\t' '{print $1,$2,$2+1,$3,$4,$5,$6}' | \
        bedtools intersect -a - -b ${ALU_BED} -wa | \
        awk -v OFS='\t' '{print $1"_"$4"_"$2,$5,$6,$7,"ALU"}' > ${LEVEL_DIR}/tmp_dir/${sample}.all.levels.txt
    sed '1d' ${file} | awk -v OFS='\t' '{print $1,$2,$2+1,$3,$4,$5,$6}' | \
        bedtools intersect -a - -b ${NONALU_REP_BED} -wa | \
        awk -v OFS='\t' '{print $1"_"$4"_"$2,$5,$6,$7,"nonALU_rep"}' >> ${LEVEL_DIR}/tmp_dir/${sample}.all.levels.txt
    sed '1d' ${file} | awk -v OFS='\t' '{print $1,$2,$2+1,$3,$4,$5,$6}' | \
        bedtools intersect -a - -b ${ALU_BED} -b ${NONALU_REP_BED} -v -wa | \
        awk -v OFS='\t' '{print $1"_"$4"_"$2,$5,$6,$7,"other"}' >> ${LEVEL_DIR}/tmp_dir/${sample}.all.levels.txt
done
echo "Classification results saved in ${LEVEL_DIR}/tmp_dir/"
