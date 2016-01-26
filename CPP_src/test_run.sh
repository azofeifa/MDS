DB=/Users/joazofeifa/Lab/gTFIv2/PSSM_DB/HOCOMOCOv10_HUMAN_mono_meme_format.meme
fasta=~/Lab/genome_files/hg19.fa
bed=~/Lab/gro_seq_files/Allen2014/EMG_out_files/test_Allen2014_DMSO2_3-19_divergent_classifications.bed
out=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/
br=5
pv=0.000001
log_out=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/
ID=test
simN=0

NP=2


mpirun -np $NP ./SE -DB $DB -fasta $fasta -bed $bed -o $out -br $br  -log_out $log_out -pv $pv -ID $ID -simN $simN