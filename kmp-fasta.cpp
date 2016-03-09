/*
MATH578a Spring 2016 Unit2 HW1 - Liana Engie
Purpose: Implementation of KMP algorithm for pattern matching.
Input:
    P, pattern
    T, fasta file containing text string to be matched against
Output:
    All starting positions of substrings in T that exactly match P

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cassert>

using std::vector;
using std::string;
using std::cout;
using std::endl;


int comparisons = 0;

static void
compute_prefix_function(const string &P, vector<size_t> &sp, int &comparisons) {
  const size_t n = P.length();
  sp.resize(n, 0);

  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {

    size_t inner_loops = 0;
    while (k > 0 && P[k] != P[i])
        k = sp[k - 1];
        ++comparisons;
    if (P[k] == P[i]) ++k;
    sp[i] = k;
  }
}

static void
Knuth_Morris_Pratt(const string &T, const string &P, const vector<size_t> &sp,
                   vector<int> &matches, int &comparisons) {

  const size_t m = P.length();
  const size_t n = T.length();

  size_t j = 0;
  for (size_t i = 0; i < n; ++i) {

    // look for the longest prefix of P that is the same as a suffix
    // of P[1..j - 1] AND has a different next character
    while (j > 0 && P[j] != toupper(T[i]))
      j = sp[j - 1];
      ++comparisons;
    // check if the character matches
    if (P[j] == toupper(T[i])) ++j, ++comparisons;

    // if we have already successfully compared all positions in P,
    // then we have found a match
    if (j == m) {
      matches.push_back(i - m + 1);
      j = sp[j - 1]; // shift by the length of the longest suffix of P
                     // that matches a prefix of P
    }
  }
}

static void
read_fasta_file_single_sequence(const string &filename, string &T) {

  std::ifstream in(filename.c_str());

    string line;
    line.reserve(3e8);
    in >> line;
    while (in >> line)
        T.append(line);
}

int
main(int argc, const char * const argv[]) {

  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " <PATTERN> <FASTA-FILE>" << endl;
    return EXIT_FAILURE;
  }

  const string P(argv[1]), filename(argv[2]);

  string T;
  read_fasta_file_single_sequence(filename, T);

  // make sure pattern not bigger than text
  assert(P.length() <= T.length());

  // preprocess the pattern "P"
  vector<size_t> sp;
  compute_prefix_function(P, sp, comparisons);

  vector<int> matches;
  Knuth_Morris_Pratt(T, P, sp, matches, comparisons);

  // output the results
  cout << endl << "MATCH COUNT:\t" << matches.size() << endl;
  cout << "Number of comparisons: " << comparisons << endl;


  return EXIT_SUCCESS;
}
