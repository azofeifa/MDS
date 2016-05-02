#PBS -S /bin/bash


#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/
#PBS -o /Users/azofeifa/qsub_stdo/

#PBS -l walltime=48:00:00
#PBS -l nodes=20:ppn=64
#PBS -l mem=1000gb
#PBS -N simulations_scanner
hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

#================================================================
#paths to config and src
src=/Users/azofeifa/Lab/gTFIv2/CPP_src/SE
DB=/Users/azofeifa/Lab/gTFIv2/PSSM_DB/HOCOMOCOv10_HUMAN_mono_meme_format_filtered_duplicates.meme
#bed=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/before_1_1/Allen2014_Nutlin2_3-3_divergent_classifications.bed
bed=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/before_1_1/Allen2014_DMSO2_3-19_divergent_classifications.bed
fasta=/Users/azofeifa/Lab/genome_files/hg19.fa
out=/Users/azofeifa/motif_distances_new/PSSM_DBS/generated_1000000_null_simulations.txt
br=100
pv=0.000001
log_out=/Users/azofeifa/
ID=DMSO_1000_high
sim_N=1000000
#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src -DB $DB -fasta $fasta -bed $bed -o $out -br $br -pv $pv -log_out $log_out -ID $ID -sim_N $sim_N
#================================================================