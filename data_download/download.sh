#! /bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="download_gen"
#SBATCH --time=24:0:00 # HH/MM/SS
#SBATCH --mem=4G # memory requested, units available: K,M,G,T
#SBATCH --mail-user=naw4005@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

# process csv file and only take needed files 
# cut -d',' -f2-3 data_url.csv | sed '1d;2d;8d' | sed 's/,/ /g' > to_download.txt

cd /athena/angsd/scratch/naw4005
while read line
do
    wget $line
done < /home/naw4005/project/to_download.txt
