Publication: 

    Giulia Bernardini (University of Trieste), Huiping Chen (King's College London), Alessio Conte (University of Pisa), Roberto Grossi (University of Pisa), Veronica Guerrini (University of Pisa), Grigorios Loukides (King's College London), Nadia Pisanti (University of Pisa), Solon Pissis (CWI).
    Utility-Oriented String Mining. 
    SIAM International Conference on Data Mining, 2024 (SDM2024).

To complile:
    make -f Makefile.mba.64.gcc 
    make -f Makefile.mia.64.gcc 
    make -f Makefile.msa.64.gcc 

To run (it might need ulimit -v and sudo):

    For getting weights from file: 
    [./mba | ./mia | ./run_msa.sh]  <text_file> <weights_file> <U> <L> 

    For getting random weights: 
    [./mba | ./mia | ./run_msa.sh] <text_file> random <U> <L> 

    For getting unit weights: 
    [./mba | ./mia | ./run_msa.sh] <text_file> unit <U> <L>

    For getting weights from a FASTQ file: 
    [./mba | ./mia | ./run_msa.sh] <fastq_file> fastq <U> <L>


CAVEAT: This program is a work-in-progress prototype, do not use on very small inputs (maybe sudo is needed in some cases)