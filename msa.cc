/**
    CAVEAT: This program is a work-in-progress prototype

    Publication: 

    Giulia Bernardini (University of Trieste), Huiping Chen (King's College London), Alessio Conte (University of Pisa), 
    Roberto Grossi (University of Pisa), Veronica Guerrini (University of Pisa), Grigorios Loukides (King's College London), 
    Nadia Pisanti (University of Pisa), Solon Pissis (CWI).
    Utility-Oriented String Mining. 
    SIAM International Conference on Data Mining, 2024 (SDM2024).
    
    Code for MSA (MUSM Scan-based Algorithm for Maximal Useful String Mining), R. Grossi, 2024

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

#undef NDEBUG

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

// hash tables
#include <unordered_map>
#include <unordered_set>

#include "headers/unordered_dense.h"  
// rolling hash
#include "headers/rabinkarphash.h"

// mmap 
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#define ERR(msg)  {fprintf(stderr, "%s\n", msg); exit(2);}

//#define VERBOSE
#define PRINT_STRINGS

using namespace std;

using INT = int64_t;

constexpr auto USEFUL_MAX_LEN = 64; // a power of 2
constexpr auto BLOCKSIZE = 64 * 4096; // 1048576; // 1048576; // a power of 2 > 2 * USEFUL_MAX_LEN


/* ---------------- AUX CLASSES AND FUNCTIONS ------------------- */

constexpr auto isvalid(double x){
  return std::isfinite(x) && (std::numeric_limits<double>::min() < x) && (x <= 1);
}
class RollingWeights final { 
  // read weights on the fly like RollingHash
  // maintain the products of the last m ones
  // the sum/count of these products
public:
  double currentWeight;
  int64_t currentPosition;
  double totalNonZeroWeights;
  int64_t countNonZeroWeights;

  RollingWeights(int_fast16_t _m) : 
    m(_m) 
  { }
  //~RollingWeights() { }

  void reset() {
    currentWeight = 0.0;
    currentNonZeroWeight = 1.0;  
    currentNonZeroWeightGhost = 1.0;
    countGhost = m;
    currentPosition = -1;
    lastZeroCostPos = -1;
    totalNonZeroWeights = 0.0;
    countNonZeroWeights = 0;
  }

  void eat(double w){ 
    if (++currentPosition < m) {
      if (w == 0) {
        lastZeroCostPos = currentPosition;
        currentWeight = 0.0;
      } else {
        currentNonZeroWeight *= w;
        currentWeight = isvalid(currentNonZeroWeight) ? currentNonZeroWeight : 0.0;
        }
      if (currentPosition == m-1 && currentWeight != 0.0){
        totalNonZeroWeights = currentWeight;
        countNonZeroWeights = 1;
      }
      if (currentPosition == m-1){
        currentNonZeroWeightGhost = 1.0;
        countGhost = m;
      }
    } else {
      throw("exceeded interval length " + to_string(m));
    }
  }

  void update(double outw, double inw){ 
    currentPosition++;
    if (outw != 0){
      currentNonZeroWeight /= outw;
    }
    if (inw != 0){
      currentNonZeroWeight *= inw;
      currentNonZeroWeightGhost *= inw;
    } else {
      lastZeroCostPos = currentPosition;
      currentNonZeroWeightGhost = 1.0;
      countGhost = m+1;  // as it is decremented anyway in the next if () statement
    }
    if (!--countGhost) {
      currentNonZeroWeight = currentNonZeroWeightGhost;
      currentNonZeroWeightGhost = 1.0;
      countGhost = m;
    } 
    if (currentPosition >= lastZeroCostPos + m){
      currentWeight = isvalid(currentNonZeroWeight) ? currentNonZeroWeight : 0.0;
    } else
      currentWeight = 0.0;
    if (currentWeight != 0.0){
      totalNonZeroWeights += currentNonZeroWeight;
      countNonZeroWeights++;
    }
  }

  private:
    long m;  // intervals' common length
    double *last_m;
    double currentNonZeroWeight;
    double currentNonZeroWeightGhost;  // better precision: every m steps it is copied into currentNonZeroWeight, instead of rolling the latter (Alessio's idea)
    long countGhost;  //  count the m steps above
    int64_t lastZeroCostPos;

}; // END class declaration

/*  ------------ */

#include "buffermmap.hpp"

/*  ------------ */

double computeUtility(string z, int32_t m, char *S, double *C, int32_t *LB, int32_t *UB, int64_t n){
  const auto ell = z.size();
  double res = 0.0;
  for (int64_t i = 0; i <= n-ell; i++){
    if (LB[i] == 0 || (LB[i]+m) >= UB[i]) continue;  
    if (z.compare(string(S+i, ell)) == 0){
      double prod = 1;
      for (int j = 0; j < ell; j++)
        prod *= C[i + j];
      res += prod;
    }
  }
  return res;
}

void score_array(vector<double> &WEIGHTS, INT n, INT kpower, double * SCORE); 
void intervalWeights(int32_t m, double *IW, int64_t n, vector<double> &C) {
  score_array(C, n, m, IW);
  for (int64_t i = n-m+1; i < n; i++)
    IW[i] = 0;
}

void printStuff(int32_t *LB, int32_t *UB, double *UT, int64_t n){
    cout << "LB: ";
    for (int64_t i = 0; i < n; i++)
      cout << LB[i] << " ";
    cout << endl;
    cout << "UB: ";
    for (int64_t i = 0; i < n; i++)
      cout << UB[i] << " ";
    cout << endl;
      cout << "UT: ";
    for (int64_t i = 0; i < n; i++)
      cout << UT[i] << " ";
    cout << endl;
}

void printDebug(int32_t *LB, int32_t *UB, double *UT, char *S, int64_t n, int64_t j, int64_t offset){
    assert(j < n);
    if (j+offset > n)
      offset = n-j;
    cout << "positions from " << j << " to " << j+offset-1 << endl;
    cout << "LB: ";
    for (int64_t i = j; i < j+offset; i++)
      cout << LB[i] << " ";
    cout << endl;
    cout << "UB: ";
    for (int64_t i = j; i < j+offset; i++)
      cout << UB[i] << " ";
    cout << endl;
    cout << "UT: ";
    for (int64_t i = j; i < j+offset; i++)
      cout << UT[i] << " ";
    cout << endl;
    cout << "S : ";
    for (int64_t i = j; i < j+offset; i++)
      cout << S[i] << " ";
    cout << endl;
}

double searchDebug(string z, double *C, int32_t *LB, int32_t *UB, double *UT, char *S, int64_t n){
  const auto ell = z.size();
  double res = 0.0;
  for (int64_t i = 0; i <= n-ell; i++){
    if (LB[i] == 0) continue; 
    if (z.compare(string(S+i, ell)) == 0){
      double prod = 1;
      for (int j = 0; j < ell; j++)
        prod *= C[i + j];
      res += prod;
      auto ii = i - 1;
      if (ii < 0) ii = 0;
      printDebug(LB, UB, UT, S, n, ii, LB[i] + 5);
    }
  }
  return res;
}


/*void display_mallinfo2(void){
  struct mallinfo2 mi;

  mi = mallinfo2();

  printf("Total non-mmapped kbytes (arena):       %zu\n", mi.arena/1000);
  printf("# of free chunks (ordblks):            %zu\n", mi.ordblks);
  printf("# of free fastbin blocks (smblks):     %zu\n", mi.smblks);
  printf("# of mapped regions (hblks):           %zu\n", mi.hblks);
  printf("kbytes in mapped regions (hblkhd):      %zu\n", mi.hblkhd/1000);
  printf("Max. total allocated space (usmblks):  %zu\n", mi.usmblks/1000);
  printf("Free kbytes held in fastbins (fsmblks): %zu\n", mi.fsmblks/1000);
  printf("Total allocated space (uordblks):      %zu\n", mi.uordblks/1000);
  printf("Total free space (fordblks):           %zu\n", mi.fordblks/1000);
  printf("Topmost releasable block (keepcost):   %zu\n", mi.keepcost/1000);
}
*/

/*-------------------------------------------------------------------*/
/*                        SCAN-BASED METHOD                          */
/*-------------------------------------------------------------------*/


/* ---------------------- INITIALIZATION SCANS --------------------- */

int32_t initializationScans(int32_t m, int64_t n, double v){
  assert(m >= 1 && v >= 1);

  buffer.begin(m);

  constexpr auto N_ROWS = 11;  // error probability delta = e^(-ROWS)

  //----------------------------//
  //cout << "/* SCAN #1-#2 */" << endl;
  //----------------------------//

  // Built Count-MinSketch (CMS) 
  constexpr auto universalHash = [](const uint64_t x, const uint64_t l, const uint64_t a) 
  {     // hashes x universally into l<=64 bits using random ODD seed a > 0
    return static_cast<uint64_t>(a*x) >> (64-l);    // see https://en.wikipedia.org/wiki/Universal_hashing
  };

  auto logLargestPowerTwo = [&](int offset)
  {
    constexpr double EULERC = 2.718281828459; // Euler's constant e, the base of natural logarithm
    constexpr auto MIN_NCOLS = 128;
    // arg is e / epsilon = e * ||counters|| / v <=  e * n / v, where epsilon n = v is the additive error; we force to be a function of min-utility v
    uint64_t value = EULERC * n / v;
    //cout << "default #columns: " << value << endl;
    assert(value >= 1);
    int l = 0;
    uint64_t p = 1;
    while (2*p < value)
    {
      p *= 2;
      l++;
    }
    if (offset >= 0 && l > offset)
      l -= offset;
    if (n <= sizeof(double) * N_ROWS * (1 << l)){ // CMS table will be larger than the text size, N_COLS == (1 << l)
      //cout << "WARNING: reduced #columns: ";
      while (l > 0 && n <= sizeof(double) * N_ROWS * (1 << l))
        l--;
      cout << (1 << l) << endl;
    }
    return l;
  };
  const uint64_t LOG_N_COLS = logLargestPowerTwo(1); 

  // create CMS
  const int64_t N_COLS = (1 << LOG_N_COLS);  //  N_COLS == 2^LOG_N_COLS
  cout << "Actual number of columns in CMS: " << N_COLS << endl; 
  vector<double> counters[N_ROWS];
  uint64_t a[N_ROWS]; //, b[N_ROWS];
  uint64_t currentHash[N_ROWS];

  // CMS init 
  mt19937_64 random_gen64(131); // change random seed here
  for (auto j = 0; j < N_ROWS; j++){
    do {
      a[j] = random_gen64();
    } while(a[j] % 2 == 0); // (a[j] == 0) || (a[j] % 2 == 0)) ;
    assert(a[j] > 0 && a[j] % 2); // a[j] = random ODD number for Universal Hash
    counters[j].reserve(N_COLS);
    for (auto l = 0; l < N_COLS; l++)
      counters[j][l] = 0.0;
  }
  
  // CMS query and update lambdas (use __builtin_prefetch args)
  constexpr auto READ_ACCESS = 0;
  constexpr auto WRITE_ACCESS = 1;
  constexpr auto NO_TEMPORAL_LOCALITY = 0;
  constexpr auto HIGH_TEMPORAL_LOCALITY = 3;
  
  auto queryCMS = [&]() 
  {
    __builtin_prefetch(&(counters[0][currentHash[0]]), READ_ACCESS, NO_TEMPORAL_LOCALITY);
    auto minval = counters[0][currentHash[0]];
    for (auto j = 1; j < N_ROWS; j++){
      __builtin_prefetch(&(counters[j][currentHash[j]]), READ_ACCESS, NO_TEMPORAL_LOCALITY);
      auto current = counters[j][currentHash[j]];
      if (minval > current) minval = current;
    }
    return minval;
  };

  auto updateCMS = [&](double w)
  {
    if (w == 0 || queryCMS() >= v) return;  // modding the original CMS
    for (auto j = 0; j < N_ROWS; j++){
      __builtin_prefetch(&(counters[j][currentHash[j]]), WRITE_ACCESS, HIGH_TEMPORAL_LOCALITY);
      counters[j][currentHash[j]] += w;
    }
  };

  // Store interval hash values and interval weights in RH and UT, actual scan occurs here <-------
  buffer.reset_S_C_RH_UT_IW_IH();
  int64_t imod = 0;
  for (int64_t j = 0; j <= n - m; j++, imod++){
    // scan using buffer.RH, buffer.UT, buffer.IH, buffer.IW, buffer.S, buffer.C
    if (j == buffer.offset){
      imod = 0;
      buffer.remapAllRefillIHIW();
    }
    buffer.UT[imod] = buffer.IW[imod];
    auto item = buffer.IH[imod];
    buffer.RH[imod] = item;
    for (auto k = 0; k < N_ROWS; k++)
      currentHash[k] = universalHash(item, LOG_N_COLS, a[k]);
    updateCMS(buffer.UT[imod]);
  }
  buffer.unmapAll();

  //----------------------------//
  //cout << "/* SCAN #3 */" << endl;
  //----------------------------//

  // Compute stats
  vector<double> V;
  V.clear();
  for (auto r = 0; r < N_ROWS; r++)
    for (auto c = 0; c < N_COLS; c++){
      V.push_back(counters[r][c]);
    }

  sort(V.begin(), V.end());
  //cout << "CM sketch size r x c : " << N_ROWS << " x " << N_COLS << " with bytes to text size ratio " << double(N_ROWS * N_COLS * sizeof(counters[0][0]))/n <<  endl;
  //cout << "minimum CM sketch val: " << V[0] << endl;
  //cout << "median  CM sketch val: " << V[V.size()/2] << endl;
  //cout << "maximum CM sketch val: " << V[V.size()-1] << endl;

  // Store m-strings and their utilities in a dictionary D only if queryCMS() >= v for their hash values, set LB[]
  auto D = ankerl::unordered_dense::map<uint64_t, double>();
  D.clear();

  buffer.reset_S_LB_RH_UT();
  imod = 0;
  for (int64_t j = 0; j <= n - m; j++, imod++){
    // scan using buffer.LB, buffer.RH, buffer.UT, buffer.S
    if (j == buffer.offset){
      imod = 0;
      buffer.remapAll();
    }
    auto item = buffer.RH[imod];
    for (auto k = 0; k < N_ROWS; k++)
      currentHash[k] = universalHash(item, LOG_N_COLS, a[k]);
    if (queryCMS() >= v){
      buffer.LB[imod] = m;
      //D[string(buffer.S + imod, m)] += buffer.UT[imod];
      D[item] += buffer.UT[imod];
    } else {
      buffer.LB[imod] = 0;
      buffer.UT[imod] = 0;
    }
  }
  buffer.unmapAll();
  
  //----------------------------//
  //cout << "/* SCAN #4 */" << endl;
  //----------------------------//

  // Only useful m-string survives (buffer.LB[i] = m)
  // - if buffer.S[i,i+m-1] in D gives utility >= v, update buffer.LB[i] and buffer.UT[i]

  buffer.reset_S_LB_UB_RH_UT();
  imod = 0;
  for (int64_t j = 0; j <= n - m; j++, imod++){
    // scan using buffer.LB, buffer.UB, buffer.RH, buffer.UT, buffer.S
    if (j == buffer.offset){
      imod = 0;
      buffer.remapAll();
    }
    //if (buffer.LB[imod] == 0 || D[string(buffer.S + imod, m)] < v) { // second term evaluated only if buffer.LB[i] != 0
    if (buffer.LB[imod] == 0 || D[buffer.RH[imod]] < v) { // second term evaluated only if buffer.LB[i] != 0
      buffer.RH[imod] = 0;
      buffer.LB[imod] = 0;
      buffer.UB[imod] = 0;
      buffer.UT[imod] = 0;
    }
  }
  buffer.unmapAll();
  
  int64_t count = 0;
  for (const auto &x : D){
     if (x.second >= v) count++;
  }
  //cout << "distinct useful: " << count << endl;
  //cout << "ratio of useful " << m << "-strings to dictionary size: " << double(count) / D.size() << endl;
  D.clear();

  //----------------------------//
  //cout << "/* SCAN #5 */" << endl;
  //----------------------------//

  // Compute buffer.UB[i] as tightly as possible for buffer.LB[i] == m
  int32_t maxrunlen = 0;
  buffer.reset_LB_UB();
  imod = 0;
  int64_t j = 0;
  while(j <= n - m){
    if (j >= buffer.offset){
      imod = j % BLOCKSIZE;
      buffer.remapAll();
    }
    // scan using buffer.LB, buffer.UB with independent indices j <= k
    if (buffer.LB[imod] == m) { // beginning of island
      int64_t k = j; 
      while (j + 1 <= n - m){
        if (buffer.LB[imod+1] != m) break;
        j++, imod++;
      }
      // k..j is an island (i.e. j = n-m or buffer.LB[imod+1] == 0)
      if (maxrunlen < j - k + 1) maxrunlen = j - k + 1;
      int64_t kmod = k % BLOCKSIZE;
      while (k <= j){
        int32_t x = j + m - k;
        if (x > USEFUL_MAX_LEN) x = USEFUL_MAX_LEN;
        buffer.UB[kmod++] = x + 1;  
        k++;
      }
    }
    j++, imod++;
  }
  buffer.unmapAll();

  assert(maxrunlen > 0);
  cout << "max run length: " << maxrunlen << endl;
  if (maxrunlen > USEFUL_MAX_LEN)
  {
    maxrunlen = USEFUL_MAX_LEN;
    //cout << "Warning. Shortened max run length to " << USEFUL_MAX_LEN << endl;
  }

  // largest power of 2 for the binary search
  int32_t p = 1;
  while (2*p <= maxrunlen) p *= 2;
  buffer.end();
  
  return p;
}



/* ---------------------- INTERMEDIATE SCANS -------------------------- */

int64_t intermediateScans(int32_t m, int64_t n, double v){

  buffer.begin(m);

  //----------------------------//
  //cout << "/* SCAN #1 */" << endl;
  //----------------------------//

  // Compute the utility of prefixes of length buffer.LB[i] + m when buffer.LB[i] >= L  using a dictionary D for their hash value and utility

  auto D = ankerl::unordered_dense::segmented_map<uint64_t, double>();
  D.clear();

  buffer.resetAll();
  int64_t imod = 0;
  for (int64_t j = 0; j <= n - m; j++, imod++){
    // scan using buffer.LB, buffer.UB, buffer.RH, buffer.UT, buffer.IH, buffer.IW (and so buffer.S, buffer.C)
    if (j == buffer.offset){
      imod = 0;
      buffer.remapAllRefillIHIW();
    }
    if (buffer.LB[imod] == 0 || (buffer.LB[imod] + m) >= buffer.UB[imod]) continue;   // ***
    uint64_t h = buffer.hashconcat(buffer.RH[imod], buffer.IH[imod + buffer.LB[imod]]);         // imod + LB[i] <= BLOCKSIZE + USEFUL_MAX_LEN < 2 * BLOCKSIZE
    /*if (D[h] < v)*/ D[h] += (buffer.UT[imod] * buffer.IW[imod + buffer.LB[imod]]);   // modded: increase only if D[h] < v
  }
  buffer.unmapAll();
  
  //----------------------------//
  //cout << "/* SCAN #2 */" << endl;
  //----------------------------//

  // Extended prefixes fromn buffer.LB[i] to buffer.LB[i] + m for those useful
  // Shorten buffer.UB[i] to buffer.LB[i] + m for those _not_ useful

  int64_t farthest = -1;
  int64_t farthestcount = 0;
 
  buffer.resetAll();
  imod = 0;
  for (int64_t j = 0; j <= n - m; j++, imod++){
    // scan using buffer.LB, buffer.UB, buffer.RH, buffer.UT, buffer.IH, buffer.IW (and so buffer.S, buffer.C)
    if (j == buffer.offset){
      imod = 0;
      buffer.remapAllRefillIHIW();
    }
    if (buffer.LB[imod] == 0) continue; 
    if ((buffer.LB[imod] + m) < buffer.UB[imod]) {   // ***
      uint64_t h = buffer.hashconcat(buffer.RH[imod], buffer.IH[imod + buffer.LB[imod]]);
      if (D[h] >= v){ // useful: extend the lower bound for position i
        buffer.RH[imod] = h;
        buffer.UT[imod] *= buffer.IW[imod + buffer.LB[imod]];  // do not invert!!
        buffer.LB[imod] += m;
      }
      else {  // not useful: shorten the upper bound for position i as cannot extend by m
        //assert(buffer.UB[i] > buffer.LB[i]+m);  
        buffer.UB[imod] = buffer.LB[imod] + m;
      }
    }
    // begin checkEliminateInconsistent 
    if (farthest >= j + buffer.UB[imod] - 1){  // it should be == but we have inconsistencies
      buffer.LB[imod] = buffer.UB[imod] = 0; // remove position i as it cannot be maximal
      farthestcount++;
    }
    if (farthest < j + buffer.LB[imod])
      farthest = j + buffer.LB[imod];
    // end checkEliminateInconsistent
  }
  buffer.unmapAll();
  
  //cout << "eliminated " << farthestcount << endl; // roughly eliminates half of dictionary entries

  //display_mallinfo2();

  auto dictsize = D.size();
  //cout << "#dictionary_items / text_size = " << (static_cast<double>(dictsize) / n) << endl;
  D.clear();

  buffer.end();

  return dictsize;
}



/* ---------------------- FINALIZATION SCAN ------------------------ */

#include "final.hpp"

void finalizationScan(int32_t m, int64_t n){
  auto T = ankerl::unordered_dense::set<std::string>(); 
  T.clear();
  buffer.reset_S_LB();
  int64_t imod = 0;
  for (int64_t j = 0; j <= n-m; j++, imod++){
    // scan using buffer.S and buffer.LB
    if (j == buffer.offset){
      imod = 0;
      buffer.remapAll();
    }
    if (buffer.LB[imod] == 0 || (j + buffer.LB[imod]) > n) continue;
    auto z = string(buffer.S + imod, buffer.LB[imod]);
    if (!T.contains(z)){ // first time found useful
      T.insert(z);
    }
  }
  buffer.unmapAll();

  string X = "";
  vector<double> util;
  for (const auto &s : T) {
    X += s;
    X += SEP_CHAR;
    //util.push_back(s.second);
  }
  T.clear();
  
  process(reinterpret_cast<unsigned char *>(const_cast<char *>(X.c_str())), X.size(), util, m);
}


/* ---------------------- PUTTING ALL TOGETHER -----------------------------*/

// TO-DO: use template T to have "T* seq"
void useful_pattern_mining(string nameS, string nameC, double v, int L) {
  assert(v >= 1);
  assert(L >= 1);
  assert(BLOCKSIZE > 2 * USEFUL_MAX_LEN);
  assert(BLOCKSIZE >= sysconf(_SC_PAGE_SIZE) && BLOCKSIZE % sysconf(_SC_PAGE_SIZE) == 0); // sysconf(_SC_PAGE_SIZE) is the system page size, e.g. 4096

  /* shared variables */
  auto nameRH = "tempfileRH_" + to_string(rand());
  auto nameLB = "tempfileLB_" + to_string(rand());
  auto nameUB = "tempfileUB_" + to_string(rand());
  auto nameUT = "tempfileUT_" + to_string(rand());
  // auto nameIH = "tempfileIH_" + to_string(rand());  // to-do: circular buffer of size 2*m
  // auto nameIW = "tempfileIW_" + to_string(rand());  // to-do: circular buffer of size 2*m
  //int fdIH, fdIW;

  size_t sizeuint64arr, sizeint32arr, sizedoublearr;
  
  /* lambdas */
  auto openfiles = [&](){
    // open the text and weights input files
    buffer.fdS = open(nameS.c_str(), O_RDONLY, (mode_t)0600);
    if ( buffer.fdS == -1 ) ERR("error opening the text file");
    buffer.fdC = open(nameC.c_str(), O_RDONLY, (mode_t)0600);
    if ( buffer.fdC == -1 ) ERR("error opening the weights file");
    /* file size in byte */
    struct stat stat_buffer;
    if (fstat( buffer.fdS, &stat_buffer) == -1) ERR("fstat error");
    buffer.n = stat_buffer.st_size;
	  cout<<"The text is of length "<< buffer.n <<endl;
    //assert( buffer.n > 2 * BLOCKSIZE);
    if (fstat(buffer.fdC, &stat_buffer) == -1)
      ERR("fstat error");
    assert( buffer.n * sizeof(double) == stat_buffer.st_size);

    // open the aux files (to be removed at the end)  TO DO error handling
    sizeuint64arr = sizeof(uint64_t) * buffer.n;
    sizeint32arr = sizeof(int32_t) * buffer.n;
    sizedoublearr = sizeof(double) * buffer.n;
    buffer.fdRH = open(nameRH.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
    lseek(buffer.fdRH, sizeuint64arr-1, SEEK_SET);
    auto dummy = write(buffer.fdRH, "", 1);
    buffer.fdLB = open(nameLB.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
    lseek(buffer.fdLB, sizeint32arr-1, SEEK_SET);
    dummy = write(buffer.fdLB, "", 1);
    buffer.fdUB = open(nameUB.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
    lseek(buffer.fdUB, sizeint32arr-1, SEEK_SET);
    dummy = write(buffer.fdUB, "", 1);
    buffer.fdUT = open(nameUT.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
    lseek(buffer.fdUT, sizedoublearr-1, SEEK_SET);
    dummy = write(buffer.fdUT, "", 1);
   };

  auto closefiles = [&](){
    close(buffer.fdS);
    close(buffer.fdC);
    close(buffer.fdRH);
    close(buffer.fdLB);
    close(buffer.fdUB);
    close(buffer.fdUT);
    // remove temp fi;es except the two input files for S and C
    std::remove(nameRH.c_str());
    std::remove(nameLB.c_str());
    std::remove(nameUB.c_str());
    std::remove(nameUT.c_str());
  };

  // ---- Open files

  openfiles();

  // ---- Initialization scans

  int64_t totalTime = 0;
  //cout << endl << "*** INITIAL SCAN ***" << endl << endl;
  auto start = std::chrono::high_resolution_clock::now();

  auto maxpower2 = initializationScans(L, buffer.n, v);  // no need for checkShorten()

  auto end = std::chrono::high_resolution_clock::now();    
  auto duration = (end - start);    
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration); // Milliseconds (as int)
  //cout << endl << "### Inital scan: "<< ms.count() << " milliseconds" << endl;
  totalTime += ms.count();

  // ---- Intermediate scans for binary search

  start = std::chrono::high_resolution_clock::now();
  
  int64_t maxDsize = 0;
  while (maxpower2 > 0){
    //cout << endl << "*** STEP " << maxpower2 << " ***" << endl << endl;

    auto x = intermediateScans(maxpower2, buffer.n, v);
    
    if (maxDsize < x) maxDsize = x;
    maxpower2 /= 2;
  }

  end = std::chrono::high_resolution_clock::now();    
  duration = (end - start);    
  ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration); // Milliseconds (as int)
  //cout << endl << "### Binary search: "<< ms.count() << " milliseconds" << endl;
  totalTime += ms.count();

  // --- Finalization scan: discard the few useful strings that survived and are not maximal, and print the maximal ones

  start = std::chrono::high_resolution_clock::now();
  finalizationScan(L, buffer.n);

  end = std::chrono::high_resolution_clock::now();    
  duration = (end - start);    
  ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration); // Milliseconds (as int)
  //cout << endl << "### Maximal useful: "<< ms.count() << " milliseconds" << endl;
  totalTime += ms.count();

  //cout << endl << "max #dictionary_items size / text_size = " << (static_cast<double>(maxDsize) / buffer.n) << endl;

  cout << endl << "### Total time: "<< totalTime << " milliseconds" << endl;


  // --- Close files
  closefiles(); 

}
/* ---------------------- END OF SCAN-BASED METHOD --------------------------*/



/*------------------------------------------------------------------*/
/*                BELOW IS VERY SIMILAR TO useful.cc                */
/*------------------------------------------------------------------*/

#include <zlib.h>
#include "headers/kseq.h"
KSEQ_INIT(gzFile, gzread) //We use the library kseq.h to read a gzipped FASTQ dataset and to split its components

/* Compute the SCORE array for some kpower using WEIGHTS -- it expects kpower <= n */
void score_array(vector<double> &WEIGHTS, INT n, INT kpower, double * SCORE)
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
		else	SCORE[i] = (SCORE[i-1]/WEIGHTS[i-1])*WEIGHTS[i+kpower-1]; // standard 2-finger computation
	}
}

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
	if( argc != 5 )
 	{
    cout<<"Wrong arguments!\n\n";
    cout<<" For getting weights from file: \n";
 		cout<< argv[0] << " <text_file> <weights_file> <U> <L> \n";
 		cout<<"\n For getting random weights: \n";
 		cout<< argv[0] << " <text_file> random <U> <L> \n";
 		cout<<"\n For getting unit weights: \n";
 		cout<< argv[0] << " <text_file> unit <U> <L>\n";
		cout<<"\n For getting weights from a FASTQ file: \n";
 		cout<< argv[0] << " <fastq_file> fastq <U> <L>\n";
 		exit(-1);
 	}
		
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
	
  string nameS(argv[1]);
  string nameC(argv[2]);
  
  useful_pattern_mining( nameS, nameC, U, L);

  return 0;
}
