#! /bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="align"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=40G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=naw4005@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

mamba activate angsd

# read in comman line arguments
gen_dir=$1
fastq_files=$2
out_dir_prefix=$3
arg_count=$#
out_dir="$(dirname "${out_dir_prefix}")"

# check we have the right number of command line arguments
if [ $arg_count -ne 3 ]; then
  echo "Not enough command line arguments. Exiting ..."
  echo "Usage: align.sh <genome_directory> <fastq_files> <out_directory>"
  exit
fi

# make directory if it does not exist
if [ ! -d ${out_dir} ]; then
  mkdir ${out_dir}
fi

STAR --runMode alignReads \
 --runThreadN 8 \
 --genomeDir ${gen_dir} \
 --readFilesIn ${fastq_files} \
 --readFilesCommand zcat \
 --outFileNamePrefix ${out_dir_prefix} \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMattributes All \