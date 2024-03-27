/* BUFFERED MMAP-ED FILES */

// Important struct/class for handling access to files in consecutive blocks of size 2 * BLOCKSIZE, to save main memory 
// for each X in {S,C,LB,UB,RH,UT}, all of n entries:
// At any time only X[offset - BLOCKSIZE ... offset + BLOCKSIZE] is accessible, and index j ranging over X
// is actually indexing the first BLOCKSIZE elements, and the last BLOCKSIZE elements are needed as lookaheas
// as the maximal allowed string is USEFUL_MAX_LEN < BLOCKSIZE/2
// Hence global index j = offset - BLOCKSIZE corresponds to local index imod = 0
// Typical usage
//
// buffer.begin(m);  /* m = size of the sliding window */
// /* it has allocaded RollingHash, RollingWeight */
// ...
// buffer.reset_X();
// int64_t imod = 0;
// for (int64_t j = 0; j <= n - m; j++, imod++){
//     /* scan using buffer.X */
//     if (j == buffer.offset){
//       imod = 0;
//       buffer.remapAll();
//     }
//     auto item = buffer.X[imod];  
//     /* it would be X[j] if the entire X would be mmap'ed */
//     ...
//   }
// buffer.unmapAll();
// ...
// buffer.reset_Y();
// int64_t imod = 0;
// for (int64_t j = 0; j <= n - m; j++, imod++){
//     /* scan using buffer.Y */
//     if (j == buffer.offset){
//       imod = 0;
//       buffer.remapAll();
//     }
//     auto item = buffer.Y[imod];  
//     /* it would be Y[j] if the entire Y would be mmap'ed */
//     ...
//   }
//   buffer.unmapAll();
//   ...
//   buffer.end;
//
  struct
  {
    uint64_t IH[2 * BLOCKSIZE];
    double IW[2 * BLOCKSIZE];
    // the ones below are mmaped into blocks of 2 * BLOCKSIZE entries
    char *S;
    double *C;
    int32_t *LB;
    int32_t *UB;
    uint64_t *RH;
    double *UT;
    // their corresponding file numbers
    int fdS, fdC, fdLB, fdUB, fdRH, fdUT;
    // Invariants
    // - S is the input string
    // - C is the input weights (double floats)
    // Invariants:
    // - L <= LB[i] <= (length of the longest useful string occuring at position i) < UB[i], where the right inequality is strict
    // - S[i..i+UB[i]-1] _not_ useful if a _maximal_ useful occurs at position i
    // - S[i..i+LB[i]-1] is _surely_ useful (and so its substrings)
    // - RH[i] = rolling hash of S[i..i+LB[i]-1]
    // - UT[i] = prod_k=i..i+LB[i]-1 C[k] is the utility of the occurence of S[i, i+LB[i]-1] at position i
    // Note that LB[i] == UB[i] == 0 indicates that there is no useful prefix of length >= L at position i, so UT[i] == 0.0
    // Invariants for step m=maxpower2=2^q where q is decreasing:
    // - IH[i] = rolling hash for S[i..i+m-1]
    // - IW[i] = prod_k=i..i+m-1 C[k] is the utility of the occurrence of S[i..i+m-1] at position i
    int64_t offset;
    size_t mapsize, oldmapsize;
    int64_t n;
    int32_t m;
    uint64_t Btom;

    // Rolling weights
    std::unique_ptr<RollingWeights> rollingWeight;

    // Rolling hash implementation by Lemine, see https://github.com/lemire/rollinghashcpp
    std::unique_ptr<KarpRabinHash<uint64>> rollingHash;

    //--- hacked from Lemire's hash to concatenate two hash values ---//
    uint64_t hashconcat(uint64_t x, uint64_t y) // y = hash of m-long string
    {
      auto res = x * Btom;
      res &= rollingHash->HASHMASK;
      res += y;
      return (res & rollingHash->HASHMASK);
    }
    //--- end ----------------------------------------------------------//

    void begin(int32_t _m)
    {
      S = NULL;
      C = NULL;
      RH = NULL;
      LB = NULL;
      UB = NULL;
      UT = NULL;
      m = _m;
      oldmapsize = mapsize = 0;
      offset = 0;

      rollingWeight.reset(new RollingWeights(m));
      rollingHash.reset(new KarpRabinHash<uint64>(m, 64));
      Btom = rollingHash->B;
      for (int32_t i = 1; i < m; i++)
      {
      Btom *= rollingHash->B;
      Btom &= rollingHash->HASHMASK;
    }
  }

  void end(){
    offset = 0;
    oldmapsize = mapsize = 0;
    rollingHash.reset();
    rollingWeight.reset();
    unmapAll();
  }

  void unmapAll(){
    if (S != NULL)
      munmap(S,  mapsize * sizeof(char));
    if (C != NULL)
      munmap(C,  mapsize * sizeof(double));
    if (RH != NULL)
      munmap(RH, mapsize * sizeof(uint64_t));
    if (LB != NULL)
      munmap(LB, mapsize * sizeof(int32_t));
    if (UB != NULL)
      munmap(UB, mapsize * sizeof(int32_t));
    if (UT != NULL)
      munmap(UT, mapsize * sizeof(double));
    S = NULL;
    C = NULL;
    RH = NULL;
    LB = NULL;
    UB = NULL;
    UT = NULL;
  }

  // All reset functions
  void __init(){
    offset = BLOCKSIZE;
    mapsize = 2 * BLOCKSIZE;
    if (mapsize > n) mapsize = n;
    oldmapsize = mapsize;
  }
  void reset_S(){
    __init();
    S = (char *)      mmap(NULL, mapsize * sizeof(char),     PROT_READ,              MAP_SHARED, fdS, 0);
    assert(S != MAP_FAILED && S != NULL);
  }
  void reset_C(){
    __init();
    C =  (double *)   mmap(NULL, mapsize * sizeof(double),   PROT_READ,              MAP_SHARED, fdC, 0); 
    assert(C != MAP_FAILED && C != NULL); 
  }
  void reset_LB(){
    __init();
    LB = (int32_t *)  mmap(NULL, mapsize * sizeof(int32_t),  PROT_READ | PROT_WRITE, MAP_SHARED, fdLB, 0);
    assert(LB != MAP_FAILED && LB != NULL);
  }
  void reset_UB(){
    __init();
    UB = (int32_t *)  mmap(NULL, mapsize * sizeof(int32_t),  PROT_READ | PROT_WRITE, MAP_SHARED, fdUB, 0);
    assert(UB != MAP_FAILED && UB != NULL);
  }
  void reset_RH_UT(){
    __init();
    RH = (uint64_t *) mmap(NULL, mapsize * sizeof(uint64_t), PROT_READ | PROT_WRITE, MAP_SHARED, fdRH, 0);
    cout << strerror(errno) << endl << flush;
    assert(RH != MAP_FAILED && RH != NULL);
    UT = (double *)   mmap(NULL, mapsize * sizeof(double),   PROT_READ | PROT_WRITE, MAP_SHARED, fdUT, 0);
    assert(UT != MAP_FAILED && UT != NULL);
  }
  void reset_S_C_IW_IH(){     
    reset_S();
    reset_C();
    for (int64_t j = 0; j < 2 * BLOCKSIZE; j++){
        IH[j] = 0;
        IW[j] = 0.0;
    }
    int64_t end = 2 * BLOCKSIZE;
    if (end > n) end = n;
    if (m <= end){ // there is at least one m-string inside
      rollingHash->reset();
      rollingWeight->reset();
      for (int64_t j = 0; j < m; j++){
        rollingHash->eat(S[j]);
        rollingWeight->eat(C[j]);
      }
      IH[0] = rollingHash->hashvalue;
      IW[0] = rollingWeight->currentWeight;
      for (int64_t j = 1; j <= end - m; j++){
        rollingHash->update(S[j - 1], S[j + m - 1]);
        IH[j] = rollingHash->hashvalue;
        rollingWeight->update(C[j - 1], C[j + m - 1]);
        IW[j] = rollingWeight->currentWeight;
      }
    }
    assert(offset == BLOCKSIZE);
  }
  
  // the following are combinations of their previous ones

  void reset_S_LB(){
    reset_S();
    reset_LB();
  }
  void reset_S_C_RH_UT_IW_IH(){  // All minus LB UB
    reset_S_C_IW_IH();
    reset_RH_UT();
  }
  void reset_S_LB_RH_UT(){
    reset_S_LB();
    reset_RH_UT();
  }
  void reset_S_LB_UB_RH_UT(){
    reset_S_LB_RH_UT();
    reset_UB();
  }
  void reset_LB_UB(){
    reset_LB();
    reset_UB();
  }
  void resetAll(){
    reset_S_C_RH_UT_IW_IH();
    reset_LB_UB();
 }

  // All remap functions

  void advance(){
    offset += BLOCKSIZE;
    oldmapsize = mapsize;
    mapsize = 2 * BLOCKSIZE;
    if (offset >= n)    
      mapsize = (n % BLOCKSIZE);
    else if (offset + BLOCKSIZE > n) mapsize = BLOCKSIZE + (n % BLOCKSIZE);
  }

  void remapAll(){
    advance();
    if (S != NULL){
      munmap(S,  oldmapsize * sizeof(char));
      S = (char *)mmap(NULL, mapsize * sizeof(char),            PROT_READ,              MAP_SHARED, fdS, (offset - BLOCKSIZE) * sizeof(char));
      assert(S != MAP_FAILED && S != NULL);
    }
    if (C != NULL){
      munmap(C,  oldmapsize * sizeof(double));
      C =  (double *)   mmap(NULL, mapsize * sizeof(double),   PROT_READ,              MAP_SHARED, fdC, (offset - BLOCKSIZE)  * sizeof(double));  
      assert(C != MAP_FAILED && C != NULL);
    }
    if (RH != NULL){
      munmap(RH, oldmapsize * sizeof(uint64_t));
      RH = (uint64_t *) mmap(NULL, mapsize * sizeof(uint64_t), PROT_READ | PROT_WRITE, MAP_SHARED, fdRH, (offset - BLOCKSIZE)  * sizeof(uint64_t));
      assert(RH != MAP_FAILED && RH != NULL);
    }
    if (LB != NULL){
      munmap(LB, oldmapsize * sizeof(int32_t));
      LB = (int32_t *)  mmap(NULL, mapsize * sizeof(int32_t),  PROT_READ | PROT_WRITE, MAP_SHARED, fdLB, (offset - BLOCKSIZE)  * sizeof(int32_t));
      assert(LB != MAP_FAILED && LB != NULL);
    }
    if (UB != NULL){
      munmap(UB, oldmapsize * sizeof(int32_t));
      UB = (int32_t *)  mmap(NULL, mapsize * sizeof(int32_t),  PROT_READ | PROT_WRITE, MAP_SHARED, fdUB, (offset - BLOCKSIZE)  * sizeof(int32_t));
      assert(UB != MAP_FAILED && UB != NULL);
    }
    if (UT != NULL){
      munmap(UT, oldmapsize * sizeof(double));
      UT = (double *)   mmap(NULL, mapsize * sizeof(double),   PROT_READ | PROT_WRITE, MAP_SHARED, fdUT, (offset - BLOCKSIZE)  * sizeof(double));
      assert(UT != MAP_FAILED && UT != NULL);
    }
  }

  // It stores interval hash values and interval weights in IH and IW
  void remapAllRefillIHIW(){
    remapAll();

    // Copy second half in first half
    int64_t k = BLOCKSIZE;
    for (int64_t j = 0; j < BLOCKSIZE; j++, k++){
      IH[j] = IH[k];
      IW[j] = IW[k];
    }
    // Fill the second half IH[BLOCKSIZE...] and IW[BLOCKSIZE...] using S[start ... end-1] and C[start ... end-1]
    int64_t start = offset;
    int64_t end = start + BLOCKSIZE;
    if (end > n) end = n;
    if (start <= end - m){ // there is at least one m-string inside
      rollingHash->reset();
      rollingWeight->reset();
      int64_t imod = BLOCKSIZE;
      for (int64_t j = start; j < start+ m; j++, imod++){
        rollingHash->eat(S[imod]);
        rollingWeight->eat(C[imod]);
      }
      IH[BLOCKSIZE] = rollingHash->hashvalue;
      IW[BLOCKSIZE] = rollingWeight->currentWeight;
      imod = BLOCKSIZE + 1;
      for (int64_t j = start + 1; j <= end - m; j++, imod++){
        rollingHash->update(S[imod - 1], S[imod + m - 1]);
        IH[imod] = rollingHash->hashvalue;
        rollingWeight->update(C[imod - 1], C[imod + m - 1]);
        IW[imod] = rollingWeight->currentWeight;
      }
    }
  }

  } buffer; // END struct