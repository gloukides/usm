This is the code for the paper 
     Utility-Oriented String Mining, to appear in the SIAM International Conference on Data Mining (SDM), 2024. 
Authors:
Giulia Bernardini (Univ. of Trieste), Huiping Chen (Univ. of Birmingham), Alessio Conte (Univ. of Pisa), 
Roberto Grossi (Univ. of Pisa), Veronica Guerrini (Univ. of Pisa), Grigorios Loukides (King's College London), 
Nadia Pisanti (Univ. of Pisa), Solon Pissis (CWI).
    
To complile the code run each of the following commands:
    make -f Makefile.mba.64.gcc 
    make -f Makefile.mia.64.gcc 
    make -f Makefile.msa.64.gcc 

To execute the different methods (this may need ulimit -v and sudo)

    with weights provided in a "weights_file" file: 
    [./mba | ./mia | ./run_msa.sh]  <text_file> <weights_file> <U> <L> 

    with random weights chosen uniformly at random from {0.7, 0.75, . . . , 1}: 
    [./mba | ./mia | ./run_msa.sh] <text_file> random <U> <L> 

    with unit weights: 
    [./mba | ./mia | ./run_msa.sh] <text_file> unit <U> <L>

    with weights read from a FASTQ file: 
    [./mba | ./mia | ./run_msa.sh] <fastq_file> fastq <U> <L>


CAVEAT: This program is a work-in-progress prototype; do not use it on very small inputs (maybe sudo is needed in some cases).
