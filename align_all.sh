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

for file in $SCRATCH_DIR/*_1.fastq.gz
do
    echo $file
    basename=$(basename "$file")  # Extracts file name with extension
    basename="${basename%%_*}"  # Removes extension
    echo "$basename"
    # perform alignment if alignment file doesn't exist
    if [ ! -e ${SCRATCH_DIR}/"alignments"/${basename}".Aligned.sortedByCoord.out.bam" ]; then
        echo "does not exist"
        echo "$file /athena/angsd/scratch/naw4005/"${basename}"_2.fastq.gz"
        echo "/athena/angsd/scratch/naw4005/alignments/"${basename}"."
        sbatch /home/naw4005/align.sh "/athena/angsd/scratch/naw4005/hg38_ind" "$file /athena/angsd/scratch/naw4005/"${basename}"_2.fastq.gz" "/athena/angsd/scratch/naw4005/alignments/"${basename}"."
    else
        echo "exists"
    fi
done

