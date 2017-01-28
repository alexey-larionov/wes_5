#!/bin/bash

# s01_export_txt.sh
# data export
# Alexey Larionov, 06Sep2016

# Not used options:
# -M 1000 output the first 1000 variants only (may be used for debugging)
# -SMA split multi-allelic variants
# -raw keep filtered (all variants have passed filters)

# Stop at errors
set -e

# Read parameters
job_file="${1}"
scripts_folder="${2}"

# Update pipeline log
echo "Started s01_export_txt: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Set parameters
source "${scripts_folder}/a01_read_config.sh"
echo "Read settings"
echo ""

# Make tmp folder
mkdir -p "${tmp_folder}"

# Go to working folder
init_dir="$(pwd)"
cd "${export_folder}"

# --- Copy source vcf to cluster --- #

# Progress report
echo "Started copying source data"

# Source files and folders (on source server)
source_vcf_folder="${dataset}"
source_vcf="${dataset}.vcf"

# Suspend stopping on errors
set +e

# Copy source vcf
rsync -thrqe "ssh -x" "${data_server}:${project_location}/${project}/${source_vcf_folder}/${source_vcf}" "${tmp_folder}/"
exit_code="${?}"

# Stop if copying failed
if [ "${exit_code}" != "0" ] 
then
  echo ""
  echo "Failed getting source data from NAS"
  echo "Script terminated"
  echo ""
  exit
fi

# Restore stopping on erors
set -e

source_vcf="${tmp_folder}/${source_vcf}"

# Progress report
echo "Completed copying source data: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export raw VCF-VEP table --- #

# Progress report
echo "Started exporting raw VCF-VEP table (non-split VEP field)"

# File names
VV_raw_txt="${tmp_folder}/${dataset}_VV_raw.txt"
VV_raw_log="${logs_folder}/${dataset}_VV_raw.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID -F LocID -F TYPE -F ID \
  -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F AS_VQSLOD -F FILTER -F AC -F AF -F AN \
  -F NEGATIVE_TRAIN_SITE -F POSITIVE_TRAIN_SITE \
  -F ANN \
  -o "${VV_raw_txt}" \
  -AMD &>  "${VV_raw_log}"

# -AMD : allow missed data
# -nda : number of discovered alt alleles

# Progress report
echo "Completed exporting raw VCF-VEP table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Split VEP field in VCF-VEP table --- #

# Progress report
echo "Started splitting VEP field in VCF-VEP table"

# Get number of the ANN column, which contains VEP annotations

found_ANN="no"
colnum=0

# Get names of columns in raw VV file
colnames=$(head -n 1 ${VV_raw_txt})

# For each column
for col in $colnames
do

  # Increment column number
  colnum=$(( $colnum + 1 ))
  
  # Check if the column name is ANN
  if [ $col == "ANN" ]
  then
      found_ANN="yes"
      break
  fi
  
done

if [ $found_ANN == "no" ]
then
  echo ""
  echo "Can not find ANN column with VEP annotations"
  echo "Script terminated"
  echo ""
  exit
fi

# File names
VV_txt="${export_folder}/${dataset}_VV.txt"

# Update header line
sed -i "1 s/""ANN""$/"${vep_fields}"/" "${VV_raw_txt}"
echo "Updated header"

# Split VEP records in the table (note the $colnum use)
awk 'BEGIN {OFS="\t"}{gsub(/\|/,"\t",$'"$colnum"'); print}' "${VV_raw_txt}" > "${VV_txt}"
echo "Updated VEP fields"

# Progress report
echo "Completed splitting VEP field: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export 1k table --- #

# Progress report
echo "Started exporting 1k table"

# File names
kgen_txt="${export_folder}/${dataset}_kgen.txt"
kgen_log="${logs_folder}/${dataset}_kgen.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID \
  -F kgen.AC \
  -F kgen.AN \
  -F kgen.AF \
  -F kgen.AFR_AF \
  -F kgen.AMR_AF \
  -F kgen.EAS_AF \
  -F kgen.EUR_AF \
  -F kgen.SAS_AF \
  -o "${kgen_txt}" \
  -AMD &>  "${kgen_log}"

# Progress report
echo "Completed exporting 1k table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export Exac table --- #

# Progress report
echo "Started exporting Exac table"

# File names
exac_txt="${export_folder}/${dataset}_exac.txt"
exac_log="${logs_folder}/${dataset}_exac.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID \
  -F exac_non_TCGA.AF \
  -F exac_non_TCGA.AC \
  -F exac_non_TCGA.AN \
  -F exac_non_TCGA.AC_FEMALE \
  -F exac_non_TCGA.AN_FEMALE \
  -F exac_non_TCGA.AC_MALE \
  -F exac_non_TCGA.AN_MALE \
  -F exac_non_TCGA.AC_Adj \
  -F exac_non_TCGA.AN_Adj \
  -F exac_non_TCGA.AC_Hom \
  -F exac_non_TCGA.AC_Het \
  -F exac_non_TCGA.AC_Hemi \
  -F exac_non_TCGA.AC_AFR \
  -F exac_non_TCGA.AN_AFR \
  -F exac_non_TCGA.Hom_AFR \
  -F exac_non_TCGA.Het_AFR \
  -F exac_non_TCGA.Hemi_AFR \
  -F exac_non_TCGA.AC_AMR \
  -F exac_non_TCGA.AN_AMR \
  -F exac_non_TCGA.Hom_AMR \
  -F exac_non_TCGA.Het_AMR \
  -F exac_non_TCGA.Hemi_AMR \
  -F exac_non_TCGA.AC_EAS \
  -F exac_non_TCGA.AN_EAS \
  -F exac_non_TCGA.Hom_EAS \
  -F exac_non_TCGA.Het_EAS \
  -F exac_non_TCGA.Hemi_EAS \
  -F exac_non_TCGA.AC_FIN \
  -F exac_non_TCGA.AN_FIN \
  -F exac_non_TCGA.Hom_FIN \
  -F exac_non_TCGA.Het_FIN \
  -F exac_non_TCGA.Hemi_FIN \
  -F exac_non_TCGA.AC_NFE \
  -F exac_non_TCGA.AN_NFE \
  -F exac_non_TCGA.Hom_NFE \
  -F exac_non_TCGA.Het_NFE \
  -F exac_non_TCGA.Hemi_NFE \
  -F exac_non_TCGA.AC_SAS \
  -F exac_non_TCGA.AN_SAS \
  -F exac_non_TCGA.Hom_SAS \
  -F exac_non_TCGA.Het_SAS \
  -F exac_non_TCGA.Hemi_SAS \
  -F exac_non_TCGA.AC_OTH \
  -F exac_non_TCGA.AN_OTH \
  -F exac_non_TCGA.Hom_OTH \
  -F exac_non_TCGA.Het_OTH \
  -F exac_non_TCGA.Hemi_OTH \
  -o "${exac_txt}" \
  -AMD &>  "${exac_log}"

# Progress report
echo "Completed exporting exac table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export GT table --- #

# Progress report
echo "Started exporting GT table"

# File names
GT_txt="${export_folder}/${dataset}_GT.txt"
GT_log="${logs_folder}/${dataset}_GT.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID -GF GT \
  -o "${GT_txt}" \
  -AMD &>  "${GT_log}"  

# Progress report
echo "Completed exporting GT table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Translate GTs from alphabetic to numeric notations --- #

# Progress report
echo "Translate GTs from alphabetic to numeric notations"

# File names
translate_gt_html="${logs_folder}/${dataset}_gt2num.html"
translate_gt_log="${logs_folder}/${dataset}_gt2num.log"

# Paremeters for R script
VV=$(basename "${VV_txt}") 
GT=$(basename "${GT_txt}")

# Prepare R script for translation with html report
r_script_translate_gt="library('rmarkdown', lib='"${r_lib_folder}"'); render('"${scripts_folder}"/r01_translate_gt_html.Rmd', params=list(dataset='"${dataset}"', working_folder='"${export_folder}"', vv_file='"${VV}"', gt_file='"${GT}"', file_out_base='"${GT%.txt}"'), output_file='"${translate_gt_html}"')"

# Execute R script for html report
echo "-------------- Preparing html report -------------- " > "${translate_gt_log}"
echo "" >> "${translate_gt_log}"
"${r_bin_folder}/R" -e "${r_script_translate_gt}" &>> "${translate_gt_log}"
echo "" >> "${translate_gt_log}"

# Names of created files 
#(hardwired within r01_translate_gt_html.Rmd script used for translation above)
GT_add="${GT%.txt}_add.txt"
GT_dom="${GT%.txt}_dom.txt"
GT_rec="${GT%.txt}_rec.txt"

# Progress report
echo "Completed translating GTs to numeric notations: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export DP table --- #

# Progress report
echo "Started exporting DP table"

# File names
DP_txt="${export_folder}/${dataset}_DP.txt"
DP_log="${logs_folder}/${dataset}_DP.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID -GF DP \
  -o "${DP_txt}" \
  -AMD &>  "${DP_log}"  

# Progress report
echo "Completed exporting DP table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export AD table --- #

# Progress report
echo "Started exporting AD table"

# File names
AD_txt="${export_folder}/${dataset}_AD.txt"
AD_log="${logs_folder}/${dataset}_AD.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID -GF AD \
  -o "${AD_txt}" \
  -AMD &>  "${AD_log}"  

# Progress report
echo "Completed exporting AD table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export GQ table --- #

# Progress report
echo "Started exporting GQ table"

# File names
GQ_txt="${export_folder}/${dataset}_GQ.txt"
GQ_log="${logs_folder}/${dataset}_GQ.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID -GF GQ \
  -o "${GQ_txt}" \
  -AMD &>  "${GQ_log}"  

# Progress report
echo "Completed exporting GQ table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export PL table --- #

# Progress report
echo "Started exporting PL table"

# File names
PL_txt="${export_folder}/${dataset}_PL.txt"
PL_log="${logs_folder}/${dataset}_PL.log"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -F SplitVarID -GF PL \
  -o "${PL_txt}" \
  -AMD &>  "${PL_log}"  

# Progress report
echo "Completed exporting PL table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Check biallelic tables --- #

# Progress report
echo "Preparing report to check exported tables"

# File names
check_html="${logs_folder}/${dataset}_check_txt.html"
check_log="${logs_folder}/${dataset}_check_txt.log"

# Paremeters for Rmarkdown scripts
VV=$(basename "${VV_txt}")
kgen=$(basename "${kgen_txt}") 
exac=$(basename "${exac_txt}")
GT=$(basename "${GT_txt}")
GT_add=$(basename "${GT_add}")
GT_dom=$(basename "${GT_dom}")
GT_rec=$(basename "${GT_rec}")
DP=$(basename "${DP_txt}")
AD=$(basename "${AD_txt}")
GQ=$(basename "${GQ_txt}")
PL=$(basename "${PL_txt}")

# Prepare script for html report
r_script_to_check_tables="library('rmarkdown', lib='"${r_lib_folder}/"'); render('"${scripts_folder}"/r02_check_update_tables_html.Rmd', params=list(dataset='"${dataset}"', working_folder='"${export_folder}"', vv_file='"${VV}"', kgen_file='"${kgen}"', exac_file='"${exac}"', gt_file='"${GT}"', gt_add_file='"${GT_add}"', gt_dom_file='"${GT_dom}"', gt_rec_file='"${GT_rec}"', dp_file='"${DP}"', ad_file='"${AD}"', gq_file='"${GQ}"', pl_file='"${PL}"'), output_file='"${check_html}"')"

# Execute R script for html report
echo "-------------- Preparing html report -------------- " > "${check_log}"
echo "" >> "${check_log}"
"${r_bin_folder}/R" -e "${r_script_to_check_tables}" &>> "${check_log}"
echo "" >> "${check_log}"

# Progress report
echo "Completed check for exported tables: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Make md5 sum for all tables --- #

# Progress report
echo "Started making md5 sums for exported text tables"

# File names
txt_md5="${dataset}_txt.md5"
 
md5sum \
  "${VV}" \
  "${kgen}" \
  "${exac}" \
  "${GT}" \
  "${GT_add}" \
  "${GT_dom}" \
  "${GT_rec}" \
  "${DP}" \
  "${AD}" \
  "${GQ}" \
  "${PL}" \
  > "${txt_md5}"
  
# Progress report
echo "Completed making md5 sums: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Copy results to NAS --- #

# Progress report
echo "Started copying results to NAS"

# Remove temporary data
rm -fr "${tmp_folder}"

# Suspend stopping on errors
set +e

# Copy files to NAS
rsync -thrqe "ssh -x" "${export_folder}" "${data_server}:${project_location}/${project}/" 
exit_code="${?}"

# Stop if copying failed
if [ "${exit_code}" != "0" ]  
then
  echo ""
  echo "Failed copying results to NAS"
  echo "Script terminated"
  echo ""
  exit
fi

# Restore stopping on errors
set -e

# Progress report
echo "Completed copying results to NAS: ${timestamp}"
echo ""
  
# Progress report to log on nas
log_on_nas="${project_location}/${project}/${dataset}_txt/logs/${dataset}_export_txt.log"
timestamp="$(date +%d%b%Y_%H:%M:%S)"
ssh -x "${data_server}" "echo \"Completed copying results to NAS: ${timestamp}\" >> ${log_on_nas}"
ssh -x "${data_server}" "echo \"\" >> ${log_on_nas}"

# Remove results from cluster
#rm -fr "${logs_folder}"
#rm -f "${VV_txt}"
#rm -f "${kgen_txt}"
#rm -f "${exac_txt}"
#rm -f "${GT_txt}"
#rm -f "${GT_add}"
#rm -f "${GT_dom}"
#rm -f "${GT_rec}"
#rm -f "${DP_txt}"
#rm -f "${AD_txt}"
#rm -f "${GQ_txt}"
#rm -f "${PL_txt}"
#rm -f "${txt_md5}"

#echo $(ssh -x "${data_server}" "echo \"Removed results from cluster\" >> ${log_on_nas}")
#ssh -x "${data_server}" "echo \"\" >> ${log_on_nas}"

# Return to the initial folder
cd "${init_dir}"
