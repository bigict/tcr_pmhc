#!/bin/bash
#
CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

export PYTHONPATH=${CWD}/profold2

help() {
  echo "usage: `basename $0` [OPTIONS] [CSV_FILE]..."
  echo "Options:"
  echo "    -o OUTPUT_DIR, --output_dir OUTPUT_DIR"
  echo "               output dir. (default: .)"
  echo "    -m MODEL, --model MODEL {fold[0-4]}"
  echo "               model. (default: {fold[0-4]})"
  echo "    -t MHC_ALIGN_RATIO_THRESHOLD --mhc_align_ratio_threshold MHC_ALIGN_RATIO_THRESHOLD"
  echo "               filter MHC with threshold=t. (default: 0.85)"
  echo "    -h, --help show this help message and exit"
  exit $1
}

output_dir="."
output_params="?attr_idx=attr.idx_all&mapping_idx=mapping.idx_all"
mhc_align_ratio_threshold=0.85
model_args=""

ARGS=$(getopt -o "o:m:h" -l "output_dir:,model:,help" -- "$@") || help 1
eval "set -- ${ARGS}"
while true; do
  case "$1" in
    (-o | --output_dir) output_dir="$2"; shift 2;;
    (-m | --model) model_args="${model_args} --model $2"; shift 2;;
    (-t | --mhc_align_ratio_threshold) mhc_align_ratio_threshold="$2"; shift 2;;
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

if [ $# -eq 0 ]; then
  help 1
fi

csv_file=$*

# convert csv to fasta files
python ${CWD}/main.py csv_to_fasta \
    --target_uri "${output_dir}${output_params}" \
    --pid_prefix tcr_pmhc_test_ \
    --default_y=1.0 \
    --verbose \
    ${csv_file}

# make chain.idx
cat ${output_dir}/mapping.idx_all | \
  cut -f2 | \
  awk -F _ '{printf("%s",$1);for (i=2;i<NF;++i) printf("_%s", $i); printf(" %s\n", $NF);}' | \
  sort -T . | \
  awk -f ${CWD}/scripts/collapse.awk  > ${output_dir}/chain.idx_all

# filter out ones that has only one chain
#   1. load dict a (in test dataset) from attr.idx_all
#   2. filter out those that:
#      i.  has no peptide
#      ii. only have peptide & MHC and in dict a
#      iii.has only one chain
#
cat ${output_dir}/chain.idx_all | \
  awk -v attr_idx=${output_dir}/attr.idx_all 'BEGIN{
      while(getline<attr_idx) {
        if ($1~/^tcr_pmhc_test_[0-9]+/) {  // We DO NOT test pMHCs
          a[$1]=1;
        }
      }
    }{
      has_P = 0;
      has_M = 0;
      for (i=2; i<=NF; ++i) {
        if ($i == "P")
          has_P +=1;
        else if ($i=="M")
          has_M +=1;
      }
      if (has_P==0 || (NF==3 && has_P>0 && has_M>0 && ($1 in a)) || NF<=2)
        print $0;
    }' > ${output_dir}/chain.idx_all_blacklist

# make attr.idx for test fold_i
cat ${output_dir}/attr.idx_all | \
  awk -v blacklist=${output_dir}/chain.idx_all_blacklist 'BEGIN{
      a["xxxxxxxx"] = 1;
      while(getline<blacklist)
        a[$1] = 1;
    }{
      if (!($1 in a)) {
        if ($1~/^tcr_pmhc_test_[0-9]+/) {
          print $0;
        }
      }
    }' > ${output_dir}/attr.idx

python ${CWD}/main.py attr_update_weight_and_task \
    --weight 1.0 \
    data/tcr_pmhc_db/attr.idx  >> ${output_dir}/attr.idx

# build the dataset (test data included) mapping.idx and chain.idx
cat ${CWD}/data/tcr_pmhc_db/mapping.idx ${output_dir}/mapping.idx_all > ${output_dir}/mapping.idx
cat ${output_dir}/mapping.idx | \
  cut -f2 | \
  awk -F _ '{printf("%s",$1);for (i=2;i<NF;++i) printf("_%s", $i); printf(" %s\n", $NF);}' | \
  sort -T . | \
  awk -f ${CWD}/scripts/collapse.awk  > ${output_dir}/chain.idx


# build fasta for each chain
for c in "A" "B" "P" "M"; do
  python ${CWD}/main.py fasta_extract \
      --target_uri ${output_dir} \
      --chain ${c} > ${output_dir}/tcr_pmhc_${c}.fa
  #find ${output_dir}/fasta -name "*_${c}.fasta" -exec awk '{print $0}' {} \; > ${output_dir}/tcr_pmhc_${c}.fa
  n=$(cat ${output_dir}/tcr_pmhc_${c}.fa | wc -l)
  if [ ${n} -eq 0 ]; then
    echo ">fake1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" > ${output_dir}/tcr_pmhc_${c}.fa
  fi
done

# align A B M with jackhmmer
for c in "A" "B" "M"; do
  find ${CWD}/data/tcr_pmhc_db/fasta -name "*_${c}.fasta" > ${output_dir}/tcr_pmhc_db_${c}
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "uniref90_db=${output_dir}/tcr_pmhc_${c}.fa mgnify_db=${CWD}/data/tcr_pmhc_db_${c}.fa sh ${CWD}/scripts/run_jackhmmer.sh -o ${output_dir}/a3m" -c 10
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "PIPELINE_UNIREF_MAX_HITS=1000000 PIPELINE_MGNIFY_MAX_HITS=1000000 PIPELINE_DEDUPLICATE=0 sh ${CWD}/scripts/run_pipeline.sh -o ${output_dir}/a3m" -c 10
done

# align P with equal length
python ${CWD}/main.py peptide_align \
  --output_dir ${output_dir}/a3m \
  --target_db ${output_dir}/tcr_pmhc_P.fa \
  --target_db ${CWD}/data/tcr_pmhc_db_P.fa \
  --verbose \
  ${CWD}/data/tcr_pmhc_db/fasta/*_P.fasta \

for c in "P"; do
  find ${CWD}/data/tcr_pmhc_db/fasta -name "*_${c}.fasta" > ${output_dir}/tcr_pmhc_db_${c}
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "PIPELINE_UNIREF_MAX_HITS=1000000 PIPELINE_MGNIFY_MAX_HITS=1000000 PIPELINE_DEDUPLICATE=0 sh ${CWD}/scripts/run_pipeline.sh -o ${output_dir}/a3m" -c 10
done

# filter a3m with threshold=t
if [ -d ${output_dir}/var ]; then
  rm -rf ${output_dir}/var
fi

echo "filter a3m (MHC): align_ratio>=${mhc_align_ratio_threshold}"
cp -r ${output_dir}/a3m ${output_dir}/var
python ${CWD}/main.py a3m_filter \
    --output_dir ${output_dir}/var \
    --aligned_ratio_threshold ${mhc_align_ratio_threshold} \
    --trim_gap \
    ${CWD}/data/tcr_pmhc_db/fasta/*_M.fasta

echo "predict ${csv_file}"
python main.py predict \
    ${model_args} \
    --output_dir ${output_dir}/pred \
    --data_dir ${output_dir}
