#!/bin/bash
#

set -e

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

CWD=`dirname ${CWD}`

help() {
  echo "usage: `basename $0` -x <conda_env> -o <output_dir>"
  exit $1
}

CONDA_ENV="pf2"
OUTPUT_DIR="."
while getopts 'x:o:h' OPT; do
  case $OPT in
    x) CONDA_ENV="$OPTARG";;
    o) OUTPUT_DIR="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))

. ${HOME}/.bashrc

cat ${HOME}/.bashrc

echo ${PATH}

USER_HOME=${HOME}
CONDA_HOME=${CONDA_HOME:-${HOME}/.local/anaconda3}

. ${CONDA_HOME}/bin/activate ${CONDA_ENV} || exit 1

export PATH=${PATH}:${USER_HOME}/bin:${CONDA_HOME}/bin

uniref90_db=${uniref90_db:-db/uniref90/uniref90.fasta}
mgnify_db=${mgnify_db:-db/mgnify/mgy_clusters.fa}

for f in $(cat); do
  echo "msa_build: ${f}"
  python ${CWD}/profold2/profold2/data/tools/jackhmmer.py \
      --output_dir=${OUTPUT_DIR} \
      --uniref90_database_path=${uniref90_db} \
      --bfd_database_path=db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
      --uniclust30_database_path=db/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
      --mgnify_database_path=${mgnify_db} \
      --pdb70_database_path=db/pdb70/pdb70 \
      --fasta_paths ${f} \
      $*
  pid=$(echo ${f}|awk -F '/' '{print $NF}' | awk -F '.fasta' '{print $1}')
  # find ${OUTPUT_DIR}/${pid}/msas -name "*.sto" -exec gzip -f {} \;
done
