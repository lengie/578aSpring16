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
#include <string>
#include <vector>
#include <cstdlib>

using std::vector;
using std::string;
using std::cout;
using std::endl;

static size_t
match(const string &s, size_t q, const size_t n) {
  for (size_t i = n; std::max(q, i) < s.length() &&
         (s[i] == s[q]); ++i, ++q);
  return q;
}

int main(int argc, const char * const argv[]) {

  if (argc != 3) {
    std::cerr << "Usage: <PATTERN> <FASTA-FILE>" << endl;
    return EXIT_FAILURE;
  }

  const string P(argv[1]);
  const string T(argv[2]);
  const string s = P+T;

  vector<size_t> Z(s.length());

  string the_case; // A label used to print the case encountered (as
                   // described in Gusfield's book)

  // !!!! The code here is based on indexing strings starting from
  // !!!! position 0, which will make the values used here different
  // !!!! from those in the book and my pseudocode.

  cout << "k" << "\t" << "l" << "\t" << "r" << "\t" << "Z[k]" << endl;

  size_t l = 0, r = 0;
  for (size_t k = 1; k < s.length(); ++k) {
    if (k >= r) { // Case 1: full comparison
      the_case = "1";
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
        the_case = "2a";
        Z[k] = Z[k_prime];
      }
      else {  // Case 2b: need to match outside the Z-box
        the_case = "2b";
        const size_t q = match(s, r, beta_len);
        Z[k] = q - k;
        r = q;
        l = k;
      }
    }
    cout << k + 1 << "\t" << l + 1 << "\t" << r << "\t"
         << Z[k] << "\t" << the_case << endl;
  }

  cout << endl
       << s << endl
       << "i\tZ[i]" << endl
       << "==\t====" << endl;

  for (size_t i = 0; i < Z.size(); ++i)
    cout << i << "\t" << Z[i] << endl;

  return EXIT_SUCCESS;
}
