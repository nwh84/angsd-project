#! /bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="run qorts"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=40G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=naw4005@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

qc_dir=$1
arg_count=$#

#if [ $arg_count -ne 1 ]; then
#  echo "Not enough command line arguments. Exiting ..."
#  echo "Usage: run_qorts.sh <qc_directory>"
#  exit
#fi

while read line
do
  if [ ! -d /athena/angsd/scratch/naw4005/qort_out/$line/ ]; then
    mkdir /athena/angsd/scratch/naw4005/qort_out/$line/
  fi
  echo /athena/angsd/scratch/naw4005/alignments/$line.Aligned.sortedByCoord.out.bam
  qorts QC \
  -Xmx5G \
  --stranded \
  /athena/angsd/scratch/naw4005/alignments/$line.Aligned.sortedByCoord.out.bam \
  /athena/angsd/scratch/naw4005/hg38.ensGene.gtf \
  /athena/angsd/scratch/naw4005/qort_out/$line/
done < "uniqueID_list.txt"

