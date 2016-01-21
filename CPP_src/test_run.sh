DB=/Users/joazofeifa/Lab/gTFI/files/HOCOMOCOv10_HUMAN_mono_meme_format.meme
fasta=~/Lab/genome_files/hg19.fa
bed=~/Lab/gro_seq_files/Allen2014/EMG_out_files/test_Allen2014_DMSO2_3-19_divergent_classifications.bed
out=/Users/joazofeifa/Lab/gro_seq_files/Allen2014/
br=500
pad=1500
log_out=XX

./SE -DB $DB -fasta $fasta -bed $bed -o $out -br $br -pad $pad -log_out $log_out