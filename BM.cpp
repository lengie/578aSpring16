/*
MATH578a Spring 2016 Unit2 HW1 - Liana Engie
Purpose: Implementation of Boyer-Moore algorithm for pattern matching.
Input:
    P, pattern
    T, fasta file containing text string to be matched against
Output:
    All starting positions of substrings in T that exactly match P
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cassert>

clock_t begintime = clock();

using namespace std;

static void create_Rlist(const string &P, struct &Rlist){
    Rlist.a; //wait need to revisit this
    Rlist.c;
    Rlist.g;
    Rlist.t;
}

static size_t
match(const string &s, size_t q, const size_t n) {
  for (size_t i = n; max(q, i) < s.length() &&
         (s[i] == s[q]); ++i, ++q);
  return q;
}

static void z_box(const string &Z){
  size_t l = 0, r = 0;
  for (size_t k = 1; k < s.length(); ++k) {
    if (k >= r) { // Case 1: full comparison
      Z[k] = match(s, 0, k);
      if (Z[k] > 0) {
        r = k + Z[k];
        l = k;
      }
    }
    else { // Case 2: (we are inside a Z-box)
      const size_t k_prime = k - l;
      const size_t beta_len = r - k;
      if (Z[k_prime] < beta_len) { // Case 2a: stay inside Z-box
        Z[k] = Z[k_prime];
      }
      else {  // Case 2b: need to match outside the Z-box
        const size_t q = match(s, r, beta_len);
        Z[k] = q - k;
        r = q;
        l = k;
      }
    }
  }
}

static void
read_fasta_file_single_sequence(const string &filename, string &T) {

    ifstream in(filename.c_str());

    string line;
    in >> line;
    while (in >> line)
        T += line;
}

int main(int argc, const char * const argv[]) {
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <PATTERN> <FASTA-FILE>" << endl;
        return EXIT_FAILURE;
    }
    const string P(argv[1]), filename(argv[2]);
    const size_t n = P.length();

    string T;
    read_fasta_file_single_sequence(filename, T);
    const size_t m = T.length();
    assert(n <= m);

    struct Rs {
        vector int a;
        vector int c;
        vector int g;
        vector int t;
    } Rlist;

    create_Rlist(P,Rlist)

    vector<size_t> Z(n); //from the pattern?
    z_box(Z);

    for(int i=0;i<m-n;++i){
        for(int j=n-1;j=0,--j){
            if(P[j]=T[j+i]){
                k=k+1;
            }
        }
        if(k==n){
        cout << "There's a matching occurence starting at " << j+1 << endl;
    }

    clock_t end = clock();
    printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);
}
