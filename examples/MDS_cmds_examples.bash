
#===========================================================
install_location=/Users/joazofeifa/Lab/ #whereever you installed MDS
#change this accordingly
#===========================================================
#input parameters and necessary paths

src=${install_location}/MDS/src/MDS #path to src

input_interval_file=${install_location}/MDS/examples/test_intervals_example.bed
fasta_file=${install_location}/MDS/examples/test_hg19.fa
PSSM_DB=${install_location}/MDS/examples/HOCOMOCO_human.txt
out_file_DB_dir=${install_location}/MDS/examples/
out_file_stats_dir=${install_location}/MDS/examples/
log_out=${install_location}/MDS/examples/

ID=unit_test_run_DB
ID2=unit_test_run_EVAL
NP=4
sim_N=100
pv=0.000001
bsn=150
H=1500


echo '-------------------------Running Unit Tests-------------------------'

echo 'Installation Location: ' $install_location ' is this correct?' 


echo '-------------------------Running DB Module--------------------------'

mpirun -np $NP $src DB -bed $input_interval_file -fasta $fasta_file -DB $PSSM_DB -o $out_file_DB_dir -log_out $log_out -ID $ID -sim_N $sim_N -pv $pv  -H $H  -t 10
wait
echo '-------------------------Running EVAL Module------------------------'

#mpirun -np $NP $src EVAL -bed $input_interval_file -fasta $fasta_file -DB $PSSM_DB -o $out_file_stats_dir -log_out $log_out -ID $ID2 -bsn $bsn  -H $H -pv $pv -t 3










