
#===========================================================
install_location=/Users/joazofeifa/Lab/gTFIv2 #whereever you installed gTFI
#change this accordingly
#===========================================================
#input parameters and necessary paths

src=${install_location}/CPP_src/SE #path to src

input_interval_file=${install_location}/examples/test_intervals_example.bed
fasta_file=${install_location}/examples/test_hg19.fa
PSSM_DB=${install_location}/examples/test_motif_db.txt
out_file_DB=${install_location}/examples/test_null_database.txt
out_file_stats=${install_location}/examples/test_interval_stats.txt
log_out=${install_location}/examples/
ID=unit_test_run_DB
ID2=unit_test_run_EVAL
NP=3
sim_N=100
pv=0.00001
bsn=150

echo '\n\n-------------------------Running Unit Tests-------------------------\n\n\n'

echo 'Installation Location: ' $install_location ' is this correct?' 


echo '\n\n-------------------------Running DB Module-------------------------\n'

mpirun -np $NP $src DB -bed $input_interval_file -fasta $fasta_file -DB $PSSM_DB -o $out_file_DB -log_out $log_out -ID $ID -sim_N $sim_N -pv $pv ; 

echo '\n\n-----------------------Running EVAL Module---------------------------\n'

mpirun -np $NP $src EVAL -bed $input_interval_file -fasta $fasta_file -DB $out_file_DB -o $out_file_stats -log_out $log_out -ID $ID2 -bsn $bsn










