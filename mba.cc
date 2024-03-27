/**
    CAVEAT: This program is a work-in-progress prototype

    Publication: 

    Giulia Bernardini (University of Trieste), Huiping Chen (King's College London), Alessio Conte (University of Pisa), 
    Roberto Grossi (University of Pisa), Veronica Guerrini (University of Pisa), Grigorios Loukides (King's College London), 
    Nadia Pisanti (University of Pisa), Solon Pissis (CWI).
    Utility-Oriented String Mining. 
    SIAM International Conference on Data Mining, 2024 (SDM2024).
    
    Code for MBA (MUSM Baseline Algorithm for Maximal Useful String Mining), S. P. Pissis, 2024

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#define DEBUG false
#include <cstdint>
#include <vector>
#include <utility> 
#include <bits/stdc++.h>
#include <random>

#include <unordered_set>
#include <unordered_map>

#include <zlib.h>
#include "headers/kseq.h"
KSEQ_INIT(gzFile, gzread) //We use the library kseq.h to read a gzipped FASTQ dataset and to split its components

#define SEP_CHAR '#' //char to delimit reads

using namespace std;

#define PRINT_STRINGS // added by RG
//#define PRINT_EXAMPLE // added by RG

#ifdef _USE_64
typedef int64_t INT;
#include <divsufsort64.h>                                         // include header for suffix sort
#endif
#ifdef _USE_32
typedef int32_t INT;
#include <divsufsort.h>                                           // include header for suffix sort
#endif
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
using namespace sdsl;

/* Comparator of tuples */
bool sort_cluster_id(const tuple<INT, INT, INT>& a, const tuple<INT, INT, INT>& b)
{
    return (get<2>(a) < get<2>(b));
}

/* Merging intervals to make them inclusion-maximal */
void mergeIntervals(vector<tuple<INT,INT,INT> > &arr, vector<tuple<INT,INT,INT> > &res) 
{
    // Test if the given set has at least one interval
    if (arr.size() <= 0)
        return;
 
    // Create an empty stack of intervals
    stack<tuple<INT, INT, INT> > s;
 
    // sort the intervals in increasing order of start time
    sort(arr.begin(), arr.end());
    if(DEBUG)
    {
    	cout<<"\n Sorted: "<<endl;
    	for (INT i = 0; i < arr.size(); i++) {

  		cout<<get<0>(arr[i])<<" "<<get<1>(arr[i])<<" "<<get<2>(arr[i]);
    		cout<<endl;
    	}
    }
    // push the first interval to stack
    s.push(arr[0]);
 
    // Start from the next interval 
    for (INT i = 1; i < arr.size(); i++) {
        // get interval from stack top
        tuple<INT,INT,INT> top = s.top();
 
	tuple<INT,INT,INT> cur=arr[i];
 	if(get<1>(top)>=get<1>(cur)) //if one interval is contained in the current top, continue 
 		continue;
        else
              s.push(arr[i]);
 
    }
   
    // Print contents of stack
    if(DEBUG)cout << "\n The Merged Intervals are: ";
    while (!s.empty()) {
        tuple<INT,INT,INT> t = s.top();
        if(DEBUG)cout << "[" << get<0>(t) << "," << get<1>(t) << "] id: "<<get<2>(t)<<endl;
        res.push_back(t);
        s.pop();
    }      
    
    return;
}

/* Kasai et al algorithm for O(n)-time LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }
        return ( 1 );
}


/*constexpr bool isvalid(double x){
  return std::isfinite(x) && (std::numeric_limits<double>::min() < x) && (x <= 1);
}*/


/* Compute the SCORE array for some kpower using WEIGHTS -- it expects kpower <= n */
void score_array(vector<double> &WEIGHTS, INT n, INT kpower, double * SCORE)//, bool &prec_loss)
{
	/* compute the score for the first k-mer */
	SCORE[0] = WEIGHTS[0];
	for(INT i=1; i < kpower; ++i) 		SCORE[0] *= WEIGHTS[i];

	for(INT i=1; i+kpower-1 < n;++i)
	{	
		/* if what is on the left is zero, I cannot divide by it, so I have to re-set */
		if ( WEIGHTS[i-1] == 0 )
		{
			INT j = i; 					// this is the starting position I have to re-set
			SCORE[j] = WEIGHTS[j];				// initialization
			for(INT j = i + 1; j < i + kpower; ++j )
			{
				if ( j == n || WEIGHTS[j] == 0 )	// if the end of the vector is reached or the first zero is read then the score is set to 0 and we break
				{
					SCORE[i] = 0;
					break;
				}
				SCORE[i] *= WEIGHTS[j];			// otherwise we multiply the j-th value
			}
		}	
		else
		{	
			SCORE[i] = (SCORE[i-1]/WEIGHTS[i-1])*WEIGHTS[i+kpower-1]; // standard 2-finger computation

			/*auto tmp = (SCORE[i-1]/WEIGHTS[i-1])*WEIGHTS[i+kpower-1]; // standard 2-finger computation
			//SCORE[i] = isvalid(tmp) ? tmp : 1.0;  // it's probably better than 0.0 --RG
			
			if(isvalid(tmp))
			{
				SCORE[i]=tmp;
			}
			else	
			{	
				prec_loss=true;
				SCORE[i]=tmp;
			}*/

		}
	}
}

bool compare_float(double x, double y, double epsilon = 0.0000000001f)
{
   if(fabs(x - y) < epsilon)
      return true; //they are same
      return false; //they are not same
}
   

/* The baseline algorithm: it prints the useful maximal patterns of seq with respect to WEIGHTS and the U threshold -- time complexity is O(n^2) in the worst case */
INT useful_pattern_mining_baseline(unsigned char * seq, vector<double> &WEIGHTS, INT n, double U, INT L, char * file)
{
 	INT * SA;
  	INT * LCP;
  	INT * invSA;

	auto start = std::chrono::high_resolution_clock::now();
	 
  	/* Compute the suffix array */
  	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
  	if( ( SA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        	return ( 0 );
  	}

	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

  	/*Compute the inverse SA array */
  	invSA = ( INT * ) calloc( n , sizeof( INT ) );
 	if( ( invSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        	return ( 0 );
  	}
  	for ( INT i = 0; i < n; i ++ )
  	{
  		invSA [SA[i]] = i;
  	}

  	/* Compute the LCP array */
  	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
  	if( ( LCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
  	}
  	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
  	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#ifndef PRINT_EXAMPLE
	free(invSA);
	#endif
	auto end = std::chrono::high_resolution_clock::now();    
	auto duration = (end - start);    
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration); // Milliseconds (as int)
	    
  	cout<<"SA and LCP constructed in "<<ms.count()<<" milliseconds"<<endl;
  	
  	start = std::chrono::high_resolution_clock::now();
  	
	double * SCORE=( double* ) calloc(n , sizeof(double));
	for ( INT i = 0; i < n; i++ ) SCORE[i] = WEIGHTS[i];
	 
  	// PARTIAL LEN, e.g. PARTIAL_LEN[i] = 2^k tells us that the pattern seq[SA[i]..SA[i]+2^k-1] is useful and THIS occurrence has individual score PARTIAL_SCR[i]
  	INT * PARTIAL_LEN;	
  	PARTIAL_LEN = ( INT * ) calloc  ( n, sizeof( INT ) );
  	double * PARTIAL_SCR;
  	PARTIAL_SCR = ( double * ) calloc  ( n, sizeof( double ) );
  	
	#ifdef PRINT_EXAMPLE
	cout << "string  \t";
	for (INT i = 0; i < n; i++){
			cout << seq[i] << "\t";
	}
	cout << endl << "weights \t";
	for( INT i = 0; i < n; i ++) {
			cout << WEIGHTS[i] << "\t";
	}
	cout << endl << endl;
	#endif

	//bool prec_loss=false;
	/* First pass */
        for ( INT len = 1; len <= n; ++len )             		// search over all lengths: 1, 2, ,..., n
        {

                bool succ = 0;                          // if a length is no good for any cluster of occurrences then we will terminate
                INT first = 0;                          // this is the index of the lex-first element (representative) of a cluster of occurrences with an LCP value of 2^k
                double score = 0;                       // this is the cumulative score (TOTAL sum) for a cluster of occurrences
		
                for( INT i = 0; i < n; i++ )            // for all O(n) starting positions we will make a linear (top to bottom) scan of the SA/LCP array
                {
                        if( LCP[i] >= len )               // if there is an LCP value of this length
			{
			          score += SCORE[SA[i]];           // add to the cumulative score for the current cluster of occurrences
			}
                        else
                        {                                  // no such an LCP value and so we need to process the current cluster and start a new one
                        	if ( ( compare_float(score, U) || score > U ) && i > 0 )       // check if the score of the current cluster is sufficient     
                               	{
                                        succ = 1;               // turn succ on because at least one cluster was good among all positions in [0,n-1]
                                        for (INT j = first; j < i; j++ )           // go through all occurrences of the current cluster
                                        {
						if(SA[j]+len <= n)
                                                {
							PARTIAL_LEN[j] = len;           	// store the current length that is good for this occurence of the current cluster
	                                                PARTIAL_SCR[j] = SCORE[SA[j]];      	// store the score of this specific occurrence
						}
                                        }
                               	}			     
                               	score = SCORE[SA[i]];                        // set the new cumulative score to the score of the first occurrence of the new cluster
                        	first = i;                                   // this is the first occurrence of the new cluster                        	
                        }
                        
                        if ( i == n - 1 && ( compare_float(score, U) || score > U ) )	// special case: we need to process separately the last cluster & check if the score of the last cluster is sufficient     
			{               
				succ = 1;               			// turn succ on because at least one cluster was good
				for (INT j = first; j <= n-1; j++ )           	// go through all occurrences of this cluster
				{
					if(SA[j]+len<=n)
					{
						PARTIAL_LEN[j] = len;           	// store the current length that is good
						PARTIAL_SCR[j] = SCORE[SA[j]];      // store the score of this specific occurrence
					}
				}
                        }
		}

		if ( succ == 0 )	break;         	// if len is no good for any cluster, terminate because no bigger len can be good by monotonicity
		else					// we need to examine the next length
		{
			if ( len + 1 <= n )		score_array(WEIGHTS, n, len + 1, SCORE);//, prec_loss);

			
			
        	}

		#ifdef PRINT_EXAMPLE
		cout << ":step: " << len << endl << "lengths \t";
		for (INT i = 0; i < n; i++)
		{
						cout << PARTIAL_LEN[invSA[i]] << "\t";
		}
		cout << endl << "loc.util\t";
		for( INT i = 0; i < n; i ++) {
			cout << PARTIAL_SCR[invSA[i]] << "\t";
		}
		cout << endl << endl;
		#endif

	}
	
	WEIGHTS.clear();
	free(SCORE);
	
	/* We next ensure that the output patterns are inclusion-maximal by sorting and merging the underlying occurrences (intervals) */
	/* In the following code, we split the SA entries into clusters depending on the occurrences of useful patterns: one cluster corresponds to one pattern with many occurrences */
	INT cluster_id = 0;
	vector<tuple<INT,INT,INT> > all_intervals; // occurrences are encoded as intervals and each interval has a cluster id: interval start, interval end, and then cluster_id

	vector<pair<INT,double> > clusters;
	
	for(INT i=0; i<n; ++i)  //iterate over the SA in O(n) time and form the clusters according to the above variables
	{
		INT first = i;
		vector<INT> is_of_cluster;
		double score = PARTIAL_SCR[i];
		is_of_cluster.push_back(i); //adds the index of the first element of cluster
					
		while((i<n-1) && PARTIAL_LEN[i]>0 && PARTIAL_LEN[i+1]==PARTIAL_LEN[i] && LCP[i+1]>=PARTIAL_LEN[i])
		{
			score+=PARTIAL_SCR[i+1];	
			is_of_cluster.push_back(i+1);			
			i++;
		}
		INT last=max(i,first);
		
		
		if(( compare_float(score, U) || score > U ) )
		{
			
			bool no_cluster=false;

			/* We next ensure that the current cluster (pattern) is a candidate for being inclusion-maximal */
			if(i > 0)	//here we check its predecessor in the LCP array
			{
				if(LCP[first]>=PARTIAL_LEN[first] && PARTIAL_LEN[first-1]>PARTIAL_LEN[first]) // this means that the current cluster (pattern) is not inclusion maximal due to its predecessor
				{
					no_cluster=true;
				}
			}
			if(i<n-1)	//here we check its succesor in the LCP array
			{
				if(LCP[last+1]>=PARTIAL_LEN[last] && PARTIAL_LEN[last+1]>PARTIAL_LEN[last]) // this means that the current cluster (pattern) is not inclusion maximal due to its sucessor
					no_cluster=true;
			}
			
			if( !no_cluster ) // the current cluster is indeed a candidate
			{
			 	for(auto &it: is_of_cluster) // we will add all occurrences (intervals) into a vector that is later going to be merged to infer inclusion-maximal patterns
				{
					tuple<INT,INT,INT> t{SA[it], SA[it]+PARTIAL_LEN[it]-1, cluster_id}; //makes the interval plus cluster_id tuples
					all_intervals.push_back(t);  					//adds them to vector
				}
  				//cluster_sizes[cluster_id]=is_of_cluster.size();  			//cluster_id -> how many intervals it has
				//cluster_scores[cluster_id]=score; 					//cluster_id -> score
				pair<INT,double> p(is_of_cluster.size(),score);
				clusters.push_back(p);
				cluster_id++;
			}
		}				
	}
	
	/* free the memory */
	free(SA);
	free(LCP);
	free(PARTIAL_LEN);
	free(PARTIAL_SCR);

	vector<tuple<INT,INT,INT> > res;
	mergeIntervals(all_intervals, res); // merge the intervals using sorting to kill patterns that are not inclusion-maximal -- this takes O(n log n) time in the worst case (unless we use radix sort)
	sort(res.begin(),res.end(),sort_cluster_id); //merged intervals sort w.r.t. cluster id
	
	INT how_many=0; //number of useful patterns	
	end = std::chrono::high_resolution_clock::now();    
	auto duration2 = (end - start);    
	
INT length_longest_useful=0;
#ifdef PRINT_STRINGS
	cout<<"Useful maximal patterns: \n";
#endif
	for(auto it = res.begin();it!=res.end();)
	{
		INT cnt=0;
		vector<tuple<INT,INT,INT>>::iterator it2=it;
		while(it2 != res.end())
		{
			if(get<2>(*it2)==get<2>(*it)) // next cluster_id is same
			    cnt++;
			else	break;  
			it2++;	  					
		}
		
		if(cnt == clusters[get<2>(*it)].first && get<1>(*it) - get<0>(*it) + 1 >= L ) // if for a cluster no interval (occurrence) was killed then it belongs to the output because it is inclusion-maximal && the L constraint is satisfied
		{
			    if(length_longest_useful < get<1>(*it) - get<0>(*it) + 1)
                                length_longest_useful= get<1>(*it) - get<0>(*it) + 1;

#ifdef PRINT_STRINGS
			for(INT j=get<0>(*it); j<=get<1>(*it);++j) // print the useful pattern
				cout<<seq[j];
			cout<<"\t score: "<<clusters[get<2>(*it)].second<<endl;
#endif
			how_many++;
		}
		it =it + cnt;
	}
	cout<<"The number of useful patterns is:"<<how_many<<endl;
	auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(duration2); // Milliseconds (as int)
	
	cout<<"The algorithm took: "<<ms2.count()<<" milliseconds."<<endl;
	cout<<"Total time is: "<<ms.count()+ms2.count()<<" milliseconds."<<endl;
	//cout<<"Prec loss: "<<prec_loss<<endl;
	cout<<"The length of the longest maximal useful:"<<length_longest_useful<<endl;	
	
	std::ofstream out_file;
        string filename = "outBASE_"+string(file); 
        out_file.open(filename,std::ios_base::app);
        out_file<<file<<"	"<<U<<"	"<<L<<"	"<<how_many<<"	"<<(ms.count()+ms2.count())<<"	"<<length_longest_useful <<endl;//"	"<<prec_loss<<endl;
        out_file.close();


	return 0;
}

/* Computes the size of the input alphabet in linear time */
INT alphabet_size(string s)
{
    unordered_map<unsigned char, INT> h;
    for (INT i = 0; i < s.length(); i++)	h[s[i]]++;
    return h.size();
}

double convertQS (char Q){
	return 1.0- pow(10,-((INT)Q-33)/10.0);	
}

int main(int argc, char **argv)
{
	if( argc != 5 )
 	{
        	cout<<"Wrong arguments!\n\n";
        	cout<<" For getting weights from file: \n";
 		cout<<" ./useful <text_file> <weights_file> <U> <L> \n";
 		cout<<"\n For getting random weights: \n";
 		cout<<" ./useful <text_file> random <U> <L> \n";
 		cout<<"\n For getting unit weights: \n";
 		cout<<" ./useful <text_file> unit <U> <L>\n";
		cout<<"\n For getting weights from a FASTQ file: \n";
 		cout<<" ./useful <fastq_file> fastq <U> <L>\n";
 		exit(-1);
 	}
	
	string second_arg(argv[2]); // check the 2nd argument "fastq" or weights file or "random"

	string text_string; //text
	INT n; //size of text_string
	vector<double> WEIGHTS; //weights
	vector<string> all_patterns_string; //patterns
	
	std::string str2(argv[3]);
	double U;			//the utility threshold
	std::stringstream(str2)>>U;

	INT L = stoi(argv[4]);			//the length threshold

	cout<<"The parameter U is set to "<<U<<endl;
	cout<<"The parameter L is set to "<<L<<endl;
	
	if( U <= 0)
	{
		cout<<"U must be >= 0"<<endl;
		return -1;
	}
		
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
		INT kseq_res = kseq_read(read_ele);
		//read_ele contains attributes seq.s (DNA seq) and qual.s (per-base qualities)
		text_string.assign(read_ele->seq.s);
		for(INT i=0; i<read_ele->qual.l; i++)
			WEIGHTS.push_back(convertQS(read_ele->qual.s[i]));
		//--------
		
		//--------Iterate for all reads
		while (( kseq_res = kseq_read(read_ele)) >= 0) {
			text_string.append(1,SEP_CHAR); //append delimiter
			text_string.append(read_ele->seq.s); //append text
			WEIGHTS.push_back(0.0); //append delimiter
			for(INT i=0; i<read_ele->qual.l; i++)
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
	
	
	cout<<"\nBASELINE ALGORITHM\n";
	/* Calling the baseline algorithm */
  	useful_pattern_mining_baseline( seq, WEIGHTS, n, U, L, argv[1]);

	cout<<"\n";
	
	return 0;
}
