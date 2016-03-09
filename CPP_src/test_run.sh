DB=/Users/joazofeifa/Lab/gTFIv2/PSSM_DB/HOCOMOCOv10_HUMAN_mono_meme_format.meme
fasta=~/Lab/genome_files/hg19.fa
bed=~/Lab/gro_seq_files/Allen2014/EMG_out_files/test_Allen2014_DMSO2_3-19_divergent_classifications.bed
out=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/
br=100
pv=0.00001
log_out=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/
ID=test
simN=10
bed_out=/Users/joazofeifa/Lab/motif_distances/test_bed_motif_files/
NP=5
site_br=1

mpirun -np $NP ./SE -DB $DB -fasta $fasta -bed $bed -o $out -br $br -bed_out $bed_out  -log_out $log_out -pv $pv -ID $ID -simN $simN -site_br $site_br