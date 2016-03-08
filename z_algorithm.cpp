/*
MATH578a Spring 2016 Unit2 HW1 - Liana Engie
Purpose: Implementation of Z algorithm for pattern matching.
Input:
    P, pattern
    T, fasta file containing text string to be matched against
Output:
    All starting positions of substrings in T that exactly match P
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using std::vector;
using std::string;
using std::cout;
using std::endl;

//clock_t begintime = clock();

static size_t
match(const string &s, size_t q, const size_t n, size_t &comparisons) {
  for (size_t i = n; std::max(q, i) < s.length() &&
         (s[i] == s[q]); ++i, ++q, ++comparisons);
  return q;
}


static void
read_fasta_file_single_sequence(const string &filename, string &T) {

  std::ifstream in(filename.c_str());

  string line;
  in >> line;
  while (in >> line)
    T += line;
}

int main(int argc, const char * const argv[]) {

  if (argc != 3) {
    std::cerr << "Usage: <PATTERN> <FASTA-FILE>" << endl;
    return EXIT_FAILURE;
  }

  const string P(argv[1]), filename(argv[2]);
  const size_t n = P.size();

  string T;
  read_fasta_file_single_sequence(filename, T);
  const string s = P+T;

  vector<size_t> Z(s.length());

  size_t matches=0;
  size_t comparisons = 0;

//cout << "There are occurrences of the pattern at " << endl;

  size_t l = 0, r = 0;
  for (size_t k = 1; k < s.length(); ++k) {
    if (k >= r) { // Case 1: full comparison
      Z[k] = match(s, 0, k, comparisons);
      if (Z[k] > 0) {
        r = k + Z[k];
        l = k;
        if(Z[k]>=n){
            //cout << k+1-n << ", ";
            ++matches;
        }
      }
    }
    else { // Case 2: (we are inside a Z-box)
      const size_t k_prime = k - l;
      const size_t beta_len = r - k;
      if (Z[k_prime] < beta_len) { // Case 2a: stay inside Z-box
        Z[k] = Z[k_prime];
        if(Z[k]>=n){
            //cout << k+1-n << ", ";
            ++matches;
        }
      }
      else {  // Case 2b: need to match outside the Z-box
        const size_t q = match(s, r, beta_len,comparisons);
        Z[k] = q - k;
        if(Z[k]>=n){
            //cout << k+1-n << ", ";
            ++matches;
        }
        r = q;
        l = k;
      }
    }
  }
cout << "Number of matches: "  << matches << endl;
cout << "Number of comparisons: " << comparisons << endl;
//    clock_t end = clock();
//   printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
}
