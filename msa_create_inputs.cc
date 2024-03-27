/* 
  g++ -O3 -o msa_create msa_create_inputs.cpp -lz 
*/
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <chrono>
#include <utility>
#include <random>
#include <iterator>
#include <algorithm>
#include <functional>
#include <malloc.h>
#include <float.h>
#include <time.h>

// mmap 
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>


using namespace std;

using INT = int64_t;

#include <zlib.h>
#include "headers/kseq.h"
KSEQ_INIT(gzFile, gzread) //We use the library kseq.h to read a gzipped FASTQ dataset and to split its components

#define SEP_CHAR 126 //char to delimit 


/* Computes the size of the input alphabet in linear time */
INT alphabet_size(string s)
{
    unordered_map<unsigned char, INT> h;
    for (INT i = 0; i < s.length(); i++)	h[s[i]]++;
    return h.size();
}


double convertQS (char Q){
	return 1.0- pow(10,-((int)Q-33)/10.0);	
}

int main(int argc, char **argv)
{
	if( argc != 3 && argc != 5)
 	{
        	cout<<"Wrong arguments!\n\n";
        	cout<<" For getting weights from file: \n";
 		cout<<  argv[0] << " <text_file> <weights_file> \n";
 		cout<<  argv[0] << " <text_file> <weights_file> <new_text_file> <new_weights_file> \n";
 		cout<<"\n For getting random weights: \n";
 		cout<< argv[0] << " <text_file> random \n";
 		cout<< argv[0] << " <text_file> random <new_text_file> <new_weights_file>\n";
 		cout<<"\n For getting unit weights: \n";
 		cout<< argv[0] << " <text_file> unit \n";
		cout<< argv[0] << " <text_file> unit <new_text_file> <new_weights_file>\n";
		cout<<"\n For getting weights from a FASTQ file: \n";
 		cout<< argv[0] << " <fastq_file> fastq \n";
 		cout<< argv[0] << " <fastq_file> fastq <new_text_file> <new_weights_file>\n";
 		exit(-1);
 	}
	
	string second_arg(argv[2]); // check the 2nd argument "fastq" or weights file or "random"

	string text_string; //text
	INT n; //size of text_string
	vector<double> WEIGHTS; //weights
	vector<string> all_patterns_string; //patterns
			
	if (second_arg.compare("fastq")!=0)
	{
		
		ifstream is;
		is.open (argv[1], ios::in | ios::binary);
		
		ifstream in_file(argv[1], ios::binary);
		in_file.seekg(0, ios::end);
		INT file_size = in_file.tellg();

		vector<unsigned char> text;
		char c = 0;
		for (INT i = 1; i < file_size; i++)
		{
			is.read(reinterpret_cast<char*>(&c), 1);
			text.push_back( (unsigned char) c );
		}
		is.close();
			
		text_string.assign(text.begin(), text.end());

		text.clear();

		INT sigma = alphabet_size(text_string);
		cout<<"The input alphabet is of size "<<sigma<<endl;
		
		cout<<"The text is of length "<< text_string.size()<<endl;
		n = text_string.size();
		  
		if(second_arg.compare("random")==0)
		{
			   cout<<"Random weight generation\n";
				vector<double> weight_domain{0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
				//mt19937 random_engine(time(nullptr));
			mt19937 random_engine(12345);
			uniform_int_distribution<INT> dist(0, weight_domain.size()-1);

			//RANDOM weight generation	
			for(INT i=0;i<n;++i)
			{
				double d=weight_domain.at(dist(random_engine));
				WEIGHTS.push_back(d);			
			}		
		}
		else if(second_arg.compare("unit")==0)
		{
			cout<<"Unit weight assignment\n";
			for(INT i=0;i<n;++i)
			{
				double d = 1.0;
				WEIGHTS.push_back(d);			
			}		

		}
		else
		{
			is.open(argv[2], ios::binary);
			std::istream_iterator<double> start(is), end;
			//WEIGHTS(start, end);
			for(auto it=start; it!=end; ++it)
			{
        WEIGHTS.push_back(*it);

      }
			is.close();
		}
	}	
	else{
		
		//Read a (possibly gzipped) FASTQ file by using kseq.h
		gzFile fastq_file; //input
		fastq_file = gzopen(argv[1], "r");
		INT numReads=1;
		//to read FASTQ 
		kseq_t *read_ele = kseq_init(fastq_file);
		
		//--------First read
		int kseq_res = kseq_read(read_ele);
		//read_ele contains attributes seq.s (DNA seq) and qual.s (per-base qualities)
		text_string.assign(read_ele->seq.s);
		for(size_t i=0; i<read_ele->qual.l; i++)
			WEIGHTS.push_back(convertQS(read_ele->qual.s[i]));
		//--------
		
		//--------Iterate for all reads
		while (( kseq_res = kseq_read(read_ele)) >= 0) {
			text_string.append(1,SEP_CHAR); //append delimiter
			text_string.append(read_ele->seq.s); //append text
			WEIGHTS.push_back(0.0); //append delimiter
			for(size_t i=0; i<read_ele->qual.l; i++)
				WEIGHTS.push_back(convertQS(read_ele->qual.s[i])); //append qs
			numReads++;
		}
		//--------
		kseq_destroy(read_ele);
		gzclose(fastq_file);
		cout << "Number of reads: " << numReads << endl;
		
		INT sigma = alphabet_size(text_string);
		cout<<"The input alphabet is of size "<<sigma<<endl;
		
		cout<<"The text is of length "<< text_string.size()<<endl;
		n = text_string.size();
	}
	
	cout<<"The weights are "<<WEIGHTS.size()<<endl;
  
	if(WEIGHTS.size()!=n)
	{
		cout<<"WEIGHTS must be n, one for each symbol!\n";
		exit(-1);
	}


  unsigned char * seq = ( unsigned char * ) text_string.c_str(); //the string

  cout << "you need free " << 33 * n / 1000000 << "Mb on your disk" << endl;

  string nameS, nameC;
  if (argc == 3){
	srand(time(NULL));
	nameC = to_string(rand());
	nameS = "tempfileText_" + nameC;
	nameC = "tempfileWeights_" + nameC;
  } else {
	nameS = string(argv[3]);
	nameC = string(argv[4]);
  }
  auto fdS = open(nameS.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
  lseek(fdS, n-1, SEEK_SET);
  auto dummy = write(fdS, "", 1);
  unsigned char *S = (unsigned char *)mmap(NULL, n, PROT_READ | PROT_WRITE, MAP_SHARED, fdS, 0);
  
  for (int64_t i = 0; i < n; i++) S[i] = seq[i];

  munmap(S, n);
  close(fdS);

  cout << " created " << nameS << endl;

  size_t sizedoublearr = sizeof(double) * n;
  auto fdC = open(nameC.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
  lseek(fdC, sizedoublearr-1, SEEK_SET);
  dummy = write(fdC, "", 1);
  double *C = (double *)mmap(NULL, sizedoublearr, PROT_READ | PROT_WRITE, MAP_SHARED, fdC, 0);


  for (int64_t i = 0; i < n; i++) C[i] = WEIGHTS[i];

  munmap(C, sizedoublearr);
  close(fdC);

  cout << " created " << nameC << endl;


  return 0;
}
