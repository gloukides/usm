/usr/bin/time -v ../mba dna.txt random 10 2 

/usr/bin/time -v ../mia dna.txt random 10 2 

cp ../msa ../msa_create .
/usr/bin/time -v ../run_msa.sh dna.txt random 10 2 
rm -f msa msa_create 