/*    Purpose: Pattern matching via the Boyer-Moore algorithm.
 *    Created by Liana Engie, February 2016
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

using std::vector;
using std::string;
using std::cout;
using std::endl;

void

struct Rs {
    char nuc;
    string pos;
} Rlist;

static void create_Rlist(const string &P, struct Rlist){
    Rlist.nuc='a'; //wait need to revisit this
    Rlist.nuc='c';
    Rlist.nuc='t';
    Rlist.nuc='g';
}

static size_t
match(const string &s, size_t q, const size_t n) {
  for (size_t i = n; std::max(q, i) < s.length() &&
         (s[i] == s[q]); ++i, ++q);
  return q;
}

void z_box(string s){

  vector<size_t> Z(s.length());

  string the_case;
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
    std::cerr << "usage: " << argv[0] << " <PATTERN> <FASTA-FILE>" << endl;
    return EXIT_FAILURE;
  }



}
