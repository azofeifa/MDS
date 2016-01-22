#PBS -S /bin/bash


#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/
#PBS -o /Users/azofeifa/qsub_stdo/

#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=64
#PBS -l mem=20gb

hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

#================================================================
#paths to config and src
src=/Users/azofeifa/Lab/gTFIv2/CPP_src/SE
DB=/Users/azofeifa/Lab/gTFIv2/PSSM_DB/HOCOMOCOv10_HUMAN_mono_meme_format.meme
bed=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/before_1_1/Allen2014_Nutlin2_3-3_divergent_classifications.bed
fasta=/Users/azofeifa/Lab/genome_files/hg19.fa
out=/Users/azofeifa/
br=500
pad=1500
log_out=/Users/azofeifa/
#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src -DB $DB -fasta $fasta -bed $bed -o $out -br $br -pad $pad -log_out $log_out
#================================================================