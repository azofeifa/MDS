DB=/Users/joazofeifa/Lab/gTFIv2/PSSM_DB/HOCOMOCO_DB_reformatted_human.txt
DB2=/Users/joazofeifa/Desktop/generated_DB.txt
fasta=/Users/joazofeifa/Lab/genome_files/hg19.fa
log_out=/Users/joazofeifa/Lab/log_out/
bed=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/test_Allen2014_DMSO2_3-19_divergent_classifications.bed
out=/Users/joazofeifa/Desktop/generated_DB.txt
out2=/Users/joazofeifa/Desktop/stats.txt
out3=/Users/joazofeifa/Desktop/test_stuff/
out_bed=/Users/joazofeifa/Desktop/hits.bed
TSS=/Users/joazofeifa/Lab/genome_files/chr1_TSS.bed


TSS=/Users/joazofeifa/Lab/TF_predictions/peak_files/ENCFF001UEO.bed
bed=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/HCT116_EMG.bed
bed=/Users/joazofeifa/Lab/genome_files/hg19_chromosome_sizes.bed
#DB=/Users/joazofeifa/Lab/gTFIv2/PSSM_DB/TEAD4_motif.txt





br=100
pv=0.00001
ID=test
simN=100
bed_out=/Users/joazofeifa/Lab/motif_distances/test_bed_motif_files/
NP=2
site_br=10
boot_strap_number=1000
test=1
PSSM_test=7
order=1
MD_window=250
#mpirun -np $NP ./SE DB -DB $DB -fasta $fasta -bed $bed -o $out -sim_N $simN -br $br -pv $pv -log_out $log_out -t $PSSM_test -order $order -TSS $TSS

#mpirun -np $NP ./SE EVAL -DB $DB2 -fasta $fasta -bed $bed -o $out2 -bsn $boot_strap_number -log_out $log_out -ID $ID -window $MD_window -TSS $TSS -o2 $out_bed 

mpirun -np $NP ./SE GENOME -DB $DB -fasta $fasta -o $out3 -log_out $log_out -t $PSSM_test -bed $bed
