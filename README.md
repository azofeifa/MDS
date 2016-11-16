#Motif Displacement Calculator
This package provides the necessary algorithms to scan for significant sites of TF-binding motifs at locations of regulatory DNA; i.e. enhancers and promoters. To compute a measure of co-occurrence between motifs and regulatory DNA, this package implements the so called motif displacement (MD) score which computes the proportion of motifs falling within some radius (-h) of all regulatory DNA centers against a larger local background (-H). 

This package consists of two modules (DB/EVAL) and is invoked as below.

```
mpirun -np <.> SE  DB <paramater flags and values>

mpirun -np <.> SE  EVAL <paramater flags and values> 
```

Evident from the run command, this c++ package requires three dependencies:

```
1. c++11
2. openmp (include <omp.h>)
3. mpi (include <mpi.h>)
```

To build 

```
$ cd /CPP_src/
$ make clean
$ make
```

The make file requires the path to mpic++ (install and config openMPI) to be in your PATH. Installing and configuring your gcc compilers will likely be cause for a headache. I found these sites useful

1. https://www.open-mpi.org
2. https://www.cyberciti.biz/faq/howto-apple-mac-os-x-install-gcc-compiler/



#Modules

##DB
A fair warning, running this module will likely take upwards of a week on a single node machine. In short, this module requires access to a large compute cluster. DB files that can be used for the eval module, in both human and mouse, are located within the PSSM_DB/. Fortunately, this file type is genome-build-free and so does not to be regenerated even if a new fasta file is used during the EVAL module. However, if you would like to re-estimate the GC distribution at your regulatory element of interest than go for it! 

| Flag | Type | Description |
|------|------|-------------| 
|-ID| some string |An identifier, all output files will begin with this prefix
|-bed  |/path/to/.bed|A bed file over which GC content will be average (Tfit, MACs) 
|-DB|/path/to/file/from/<PSSM_DB>|A specific file from PSSM_DB/; gives the motif models
|-o|/path/to/ |Will output <-ID>.db; this will be required for the EVAL module
|-fasta|/path/to/genome.fasta|fasta file of the same genome build as <-bed> 
|-TSS  |/path/to/.bed|A bed file of promoter start sites  
|-log_out|/path/to |Where temporary and final log files will be generate <-ID>.log
|-H|numerical|distance around which sequence will collected (default = 1500bp)
|-pv|numerical|pvalue threshold under which a motif will be considered significant
|-sim_N|numerical|number of random sequence generations; (default=10,000,000)

###Output file type (-o)

![alt tag](https://github.com/azofeifa/gTFIv2/blob/master/images/DB_FILE_OUT.png)

Above: A screen shot of a small porition of the db file that is outputed from running DB module. The file type is broken up into blocks according to the PSSM model (641 in human).  Each block (delimited by the ~ symbol) contains the probability distribution matrix of each PSSM model. The final two lines is the empiracle distribution of motif displacement estimated from the non-stationary GC content surrounding the regulatory element.  

##EVAL
The EVAL module computes the so called motif displacement (MD) score. The module follows this generic pipeline:


1. Takes in as input a bed file corresponding to the regulatory DNA
2. Extracts the underlying nucleotide sequence from a provide fasta file (of the same build as the bed file) 
3. Scans for significant motif sites (gathered from the file.db either from the PSSM_DB/ directory or generated from the DB module) 
4. assess signficance of the MD score under binomial model (stationary background model) and (non-stationary background model estimated from the GC bias from the DB module).

| Flag | Type | Description |
|------|------|-------------|
|-ID| some string |An identifier, all output files will begin with this prefix
|-bed  |/path/to/.bed|A bed file over which motif displacement will be calculated
|-DB|/path/to/file.db|output file from the db module; gives the simulation,parameter info
|-o|/path/to/ |Will output <-ID>.tsv; this provides information on motif displacements and scores
|-fasta|/path/to/genome.fasta|fasta file of the same genome build as <-bed>
|-TSS  |/path/to/.bed|A bed file of promoter start sites
|-log_out|/path/to |Where temporary and final log files will be generate <-ID>.log
|-h|numerical|distance around which the MD score will be computed (default = 150bp)
|-bsn|numerical|number of random draws from the empiracle distribution estimated from DB module; (default=10,000)

![alt tag](https://github.com/azofeifa/gTFIv2/blob/master/images/Enrichment_FILE_OUT.png)

Above: A screen shot of output from the EVAL module (-o). To explain the file type briefly, lines starting with a \# indicate descriptive headings and should be ignored in downstream analysis.  The file is broken up into three parts. 

(1) The first portion of the MD score file provides for each motif model the computed MD score from the set of promoter associated, non-promoter associated and combined bed file intervals. Along with the MD score, each line provides the number of motif occurrences in the larger (-H) window. Finally, p-values computed under two different background null models (stationary and non-stationary; binomial and empirical simulations; are provided for each motif model and each promoter association type. 

(2) The second portion of the MD score file-following the header: #Binned Observation statistics range-provides the raw motif displacement frequency (histogram) according to each motif model. 

(3) The third and final portion of the MD score file-following the header: #Empiracle Bootstrapped Distribution-provides the empirical samples of the MD-score and used to calculate the non-stationary p-value.


#Advanced HPC usage
Those using compute cluster may find this below qsub script useful!


```bash
#PBS -S /bin/bash

#PBS -N gTFIv2


#PBS -l walltime=72:00:00
#PBS -l nodes=10:ppn=64
#PBS -l mem=10gb
hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )

# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
src=/path/to/gTFIv2/CPP_src/repo/SE


$cmd $src <module> <parameter flags and values>

```

#Questions/Comments
Please email me (Joey) at joseph[.]azofeifa[@]colorado[.]edu if you have any questions on usage or bug reports. Or open an Issue.

