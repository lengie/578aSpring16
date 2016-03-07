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
    for(size_t i=P.size()-1; i=0; --i){
        if(s[i]=='A'){
            Rlist.a.push_back(i);
        }
        else if(s[i]=='C'){
            Rlist.c.push_back(i);
        }
        else if(s[i]=='G'){
            Rlist.g.push_back(i);
        }
        else if(s[i]=='T'){
            Rlist.t.push_back(i);
        }
    }
}
static int lookupR(struct &Rlist, size_t &pos, char bp){
 //Is it possible to do a binary search on a struct?
    size_t i = 0;
    if(bp=='A'){
        while(Rlist.a[i]<=pos){
            ++i;
        }
    } else if(bp=='C'){
        while(Rlist.c[i]<=pos){
                ++i;
            }
    } else if(bp=='G'){
        while(Rlist.g[i]<=pos){
            ++i;
        }
    }else{
        while(Rlist.a[i]<=pos){
            ++i;
        }
    }
    return i;
}

static size_t
match(const string &s, size_t q, const size_t n) {
  for (size_t i = n; max(q, i) < s.length() &&
         (s[i] == s[q]); ++i, ++q);
  return q;
}

static void reverseZ(vector<size_t> &N, const string &P){
    n = P.size();
    for(size_t k=0;k<n;++k){
        Pr[n-k-1] = P[i];
    }

    size_t l = 0, r = 0;
    for (size_t k = 1; k < n; ++k) {
        if (k >= r) { // Case 1: full comparison
            N[k] = match(Pr, 0, k);
            if (N[k] > 0) {
                r = k + N[k];
                l = k;
            }
        }
        else { // Case 2: (we are inside a Z-box)
            const size_t k_prime = k - l;
            const size_t beta_len = r - k;
            if (N[k_prime] < beta_len) { // Case 2a: stay inside Z-box
                N[k] = N[k_prime];
            }
            else {  // Case 2b: need to match outside the Z-box
                const size_t q = match(Pr, r, beta_len);
                N[k] = q - k;
                r = q;
                l = k;
            }
        }
    }
    P = Pr; //Difference between &P = &Pr? test this later
}

static void good_suffix(vector<size_t> &Lprime,const string &N, const size_t n){
    for(size_t j = 0; i<n; ++i){ //n-1?
        size_t i = n-N[j];
        Lprime[i] = j;
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

    struct Rs{
        vector<int> a;
        vector<int> c;
        vector<int> g;
        vector<int> t;
    } Rlist;

    //Setting up extended bad character rule
    create_Rlist(P,Rlist)

    //Setting up good suffix rule
    vector<size_t> N(n);
    reverseZ(N,P);

    vector<size_t> Lprime(n,0); //should be a constant or not?
    good_suffix(Lprime,N,n);

    vector<size_t> l(n,0);
    //largest suffix of P[i...n] that is also a prefix of P

    /*for(int i=0;i<m-n;){
        size_t k=0;
        for(int j=n-1;j=0,--j){
            if(P[j]=T[j+i]){
                k=k+1;
            }
            else{ //mismatch
                //Extended bad character rule
                i=lookupR(&Rlist,j+i,T[j+i]);
            }
        }
        if(k==n){
            cout << "There's a matching occurrence starting at " << j+1 << endl;
            i=i+n;
        }
    }*/

    size_t k =0;
    while(k<=m){
        size_t i = n;
        size_t h = k;
        while(i>0 && P[i] = T[h]){
            i = i-1;
            h = h-1;
        }
        if(i=0){
            cout << "There's a matching occurrence starting at " << k+1 << endl;
            k=k_n-l(1);
        }else{
            //shift P (increase k) by the max amunt determined by extended bad char rule and good suffix rule
        }
    }

    clock_t end = clock();
    printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);
}
