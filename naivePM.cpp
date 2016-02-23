/*    Purpose: Naive implementation of pattern matching.
 *    Created by Liana Engie, February 2016
 *
 *    Input: P, pattern and T, fasta file containing text string to be matched against
 *    Output: All starting positions of substrings in T that exactly match P
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

using std::vector;
using std::string;
using std::cout;
using std::endl;

static void
read_fasta_file_single_sequence(const string &filename, string &T) {

  std::ifstream in(filename.c_str());

  string line;
  in >> line;
  while (in >> line)
    T += line;
}

int main(int argc, const char * const argv[]) {
    const string P(argv[1]), filename(argv[2]);
    const size_t n = P.length();

    const string P(argv[1]), filename(argv[2]);

    string T;
    read_fasta_file_single_sequence(filename, T);
    assert(P.length() <= T.length());

    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " <PATTERN> <TEXT>" << endl;
        return EXIT_FAILURE;
    }

    for(int j=0;j<m-n;++j){ //Could(should?) also use size_t for count here
        k = 0;
        for(int i=0;i<n;++i){
            if P[i]=t[j+i-1]{
                k=k+1;
            }
        }
        if(k==n){
            cout << "There's a matching occurence starting at " << j << endl;
        }
    }
}
