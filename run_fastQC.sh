#! /bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --job-name="run fastqc"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=8G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=naw4005@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

arg_count=$#

if [ $arg_count -ne 1 ]; then
  echo "Not enough command line arguments. Exiting ..."
  echo "Usage: run_fastQC.sh <out_fastqc_directory>"
  exit
fi

mamba activate angsd

out_dir=$1
if [ -d "$out_dir" ]; then
    echo "Directory already exists. Exiting script."
    exit
else
    echo "Creating out directory..."
    mkdir $out_dir
fi

# run fastqc on all files
fastqc -t 6 $SCRATCH_DIR/*_1.fastq.gz -o $out_dir
fastqc -t 6 $SCRATCH_DIR/*_2.fastq.gz -o $out_dir

# index files
samtools index $SCRATCH_DIR/alignments/*.Aligned.sortedByCoord.out.bam

# run samtools flagstat on all our samples
for bam_file in $SCRATCH_DIR/alignments/*.Aligned.sortedByCoord.out.bam
do
    directory="$(dirname "$bam_file")"
    filename="$(basename "$bam_file")"
    output_file=$out_dir/$filename.flagstats
    if [ ! -e $output_file ]
    then
        samtools flagstat $bam_file > $output_file
    fi
done

# run multiqc
mamba activate multiqc
multiqc -o $out_dir -n multiqc.ouput $out_dir