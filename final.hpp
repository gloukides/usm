#define DEBUG false
#include <iostream>
#include <cstdint>
#include <vector>
#include <utility> 
#include <string.h>
#include <unordered_set>
#include <unordered_set>

#define SEP_CHAR 126 //char to delimit 
#define MY_USEFUL_MAX_LEN 8192 // a power of 2

using namespace std;

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
typedef int64_t INT;
#endif
#ifdef _USE_32
#include <divsufsort.h>                                           // include header for suffix sort
typedef int32_t INT;
#endif

#include <stack>
#include <math.h> 

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

bool compare_float(double x, double y, double epsilon = 0.0000000001f)
{
   if(fabs(x - y) < epsilon)
      return true; //they are same
      return false; //they are not same
}

INT process(unsigned char *seq, INT n, vector<double> &util, INT L)
{
	
	INT * SA;
  	INT * LCP;
  	INT * invSA;

	vector<INT>  TEMP_PARTIAL_LEN;
	TEMP_PARTIAL_LEN.reserve(n);

 	INT cnt=0;
    for(INT x=n-1;x>=0;--x)
    {
    	if(seq[x]==SEP_CHAR)
    		cnt=0;
    	TEMP_PARTIAL_LEN[x]=cnt++;

    }
    
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
	free(invSA);

    vector<INT> PARTIAL_LEN;
    PARTIAL_LEN.reserve(n);
    for(INT x=0;x<n;++x)
    {
    	PARTIAL_LEN[x]=TEMP_PARTIAL_LEN[SA[x]];
    }

    /* We next ensure that the output patterns are inclusion-maximal by sorting and merging the underlying occurrences (intervals) */
	/* In the following code, we split the SA entries into clusters depending on the occurrences of useful patterns: one cluster corresponds to one pattern with many occurrences */
	INT cluster_id = 0;
	vector<tuple<INT,INT,INT> > all_intervals; // occurrences are encoded as intervals and each interval has a cluster id: interval start, interval end, and then cluster_id

	vector<INT > clusters;
	
	for(INT i=0; i<n; ++i)  //iterate over the SA in O(n) time and form the clusters according to the above variables
	{
		INT first = i;
		vector<INT> is_of_cluster;
		is_of_cluster.push_back(i); //adds the index of the first element of cluster
					
		while((i<n-1) && PARTIAL_LEN[i]>0 && PARTIAL_LEN[i+1]==PARTIAL_LEN[i] && LCP[i+1]>=PARTIAL_LEN[i])
		{
			is_of_cluster.push_back(i+1);			
			i++;
		}
		INT last=max(i,first);
		
			
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
				
				clusters.push_back(is_of_cluster.size());
				cluster_id++;
			}
						
	}
	vector<tuple<INT,INT,INT> > res;
	mergeIntervals(all_intervals, res); // merge the intervals using sorting to kill patterns that are not inclusion-maximal -- this takes O(n log n) time in the worst case (unless we use radix sort)
	sort(res.begin(),res.end(),sort_cluster_id); //merged intervals sort w.r.t. cluster id
	
	INT how_many=0; //number of useful patterns
	INT longest = 0; // longest among them

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
		
		if(cnt == clusters[get<2>(*it)] && get<1>(*it) - get<0>(*it) + 1 >= L ) // if for a cluster no interval (occurrence) was killed then it belongs to the output because it is inclusion-maximal && the L constraint is satisfied
		{
#ifdef PRINT_STRINGS
			for (INT j = get<0>(*it); j <= get<1>(*it); ++j)
			{ // print the useful pattern
				cout << seq[j];
			}
			#if 0
			auto pos = 0;
			for (INT j = 0; j < get<0>(*it); j++)
				if (seq[j] == SEP_CHAR) pos++;
			cout << "\t score: " << util[pos];
			#endif 
			cout<<endl;
#endif
			how_many++;
			//total_length_of_output+=(get<1>(*it)-get<0>(*it)+1);
			auto local = (get<1>(*it) - get<0>(*it) + 1);
			if (longest < local) longest = local;
		}
		it =it + cnt;
	}

	cout << endl << "Found " << how_many << " useful strings" << endl;
	cout << endl << "Longest has length " << longest << endl;
	if (longest >= MY_USEFUL_MAX_LEN)
		cout << "Please increase MY_USEFUL_MAX_LEN and USEFUL_MAX_LEN as the outcome is not reliable" << endl;

	/* free the memory */
	free(SA);
	free(LCP);

	return 0;

}

#if 0
int main(int argc, char **argv)
{
	string text_string; //text
	INT n; //size of text_string

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

  	unsigned char * seq = ( unsigned char * ) text_string.c_str(); //the string
  	cout<<"The text is of length "<< text_string.size()<<endl;
	n = text_string.size();
	cout<<text_string<<"|"<<n<<endl;

	INT L = stoi(argv[2]);			//the length threshold
	process(seq,n, L);

	return 0;
}
#endif