/*
MATH578a Spring 2016 Unit2 HW1 - Liana Engie
Purpose: Naive implementation of pattern matching.
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

static void
read_fasta_file_single_sequence(const string &filename, string &T) {

  std::ifstream in(filename.c_str());

  string line;
  in >> line;
  while (in >> line)
    T += line;
}

int main(int argc, const char * const argv[]) {
    /*
    //Testing on a string T, not fasta file
    const string P(argv[1]), T(argv[2]);
    const size_t n = P.length();
    const size_t m = T.length();*/
    cout << "There's a matching occurrence starting at ";

    const string P(argv[1]), filename(argv[2]);
    const size_t n = P.length();

    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " <PATTERN> <FILE-NAME>" << endl;
        return EXIT_FAILURE;
    }

    string T;
    read_fasta_file_single_sequence(filename, T);
    const size_t m = T.length();
    assert(n <= m);

    for(int j=0;j<m-n;++j){ //Could(should?) also use size_t for count here
        size_t k = 0;
        for(int i=0;i<n;++i){
            if (P[i]==T[j+i]){
                k=k+1;
            }
        }
        if(k==n){
            cout << j+1 << ", ";
        }
    }

    clock_t end = clock();
    printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);
}
