#PBS -S /bin/bash

#PBS -N gTFIv2
#PBS -m ae

#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=64
#PBS -l mem=10gb
hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4



cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"


#===========================================================
install_location=/Users/azofeifa/Lab/gTFIv2/ #wherever you installed gTFI, change this accordingly
#===========================================================


src=${install_location}/CPP_src/MDS

ID=unit_test

#===========================================================
#necessary input files
input_interval_file=${install_location}examples/test_intervals_example.bed
fasta_file=${install_location}examples/test_hg19.fa
PSSM_DB=${install_location}examples/test_motif_db.txt
TSS=${install_location}examples/hg19_TSS.bed
log_out=${install_location}/examples/
#===========================================================
#files that will output
out_dir=${install_location}/examples/
#===========================================================
NP=3
sim_N=100
pv=0.00001
bsn=150
H=1500
h=150



#mpirun -np $NP $src DB -bed $input_interval_file -TSS $TSS -fasta $fasta_file -DB $PSSM_DB -o $out_dir -log_out $log_out -ID $ID -sim_N $sim_N -pv $pv -H $H ;


mpirun -np $NP $src EVAL -bed $input_interval_file -fasta $fasta_file -DB ${install_location}/examples/${ID}.db -o $out_dir -log_out $log_out -ID $ID -bsn $bsn -TSS $TSS

