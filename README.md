#Motif Displacement Calculator

#Modules

##DB

| Flag | Type | Description |
|------|------|-------------| 
|-ID| some string |An identifier, all output files will begin with this prefix
|-bed  |/path/to/.bed|A bed file over which GC content will be average (Tfit, MACs) 
|-o|/path/to/ |Will output <-ID>.db; this will be required for the EVAL module
|-fasta|/path/to/genome.fasta|fasta file of the same genome build as <-bed> 
|-TSS  |/path/to/.bed|A bed file of promoter start sites  
|||

##EVAL

