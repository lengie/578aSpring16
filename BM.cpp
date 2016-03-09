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
#include <algorithm>

using namespace std;

struct Rs{
    vector<size_t> a;
    vector<size_t> c;
    vector<size_t> g;
    vector<size_t> t;
} Rlist;

static void create_Rlist(const string &P, Rs Rlist){
    for(size_t i=P.size()-1; i>=0; --i){
        if(P[i]=='A'){
            Rlist.a.push_back(i);
        }
        else if(P[i]=='C'){
            Rlist.c.push_back(i);
        }
        else if(P[i]=='G'){
            Rlist.g.push_back(i);
        }
        else if(P[i]=='T'){
            Rlist.t.push_back(i);
        }
    }
    cout << "Rs list was fine" << endl;
}
static int lookupR(struct Rs &Rlist, int pos, char bp){
    size_t j = 0;
    int newshift = pos;
    if(bp=='A'){
        while(Rlist.a[j]>=pos && (j<Rlist.a.size())){
            ++j;
        }newshift = Rlist.a[j];
    } else if(bp=='C'){
        while(Rlist.c[j]>=pos && j<Rlist.c.size()){
            ++j;
            }newshift = Rlist.c[j];
    } else if(bp=='G'){
        while(Rlist.g[j]>=pos && j<Rlist.g.size()){
            ++j;
        }newshift = Rlist.g[j];
    }else{
        while(Rlist.a[j]>=pos && j<Rlist.t.size()){
            ++j;
        }newshift = Rlist.g[j];
    }
    return newshift;
}

static size_t
match(const string &s, size_t q, const size_t n, size_t &comparisons) {
  for (size_t i = n; max(q, i) < s.length() &&
         (s[i] == s[q]); ++i, ++q, ++comparisons);
  return q;
}

static void reverseZ(vector<size_t> &N, const string &P, size_t &comparisons){
    const size_t n = P.size();
    string Pr(n, ' ');
    for(size_t k=0;k<n;++k){
        Pr[n-k-1] = P[k];
    }

    size_t l = 0, r = 0;
    for (size_t k = 1; k < n; ++k) {
        if (k >= r) {
            N[k] = match(Pr, 0, k, comparisons);
            if (N[k] > 0) {
                r = k + N[k];
                l = k;
            }
        }
        else {
            const size_t k_prime = k - l;
            const size_t beta_len = r - k;
            if (N[k_prime] < beta_len) {
                N[k] = N[k_prime];
            }
            else {
                const size_t q = match(Pr, r, beta_len, comparisons);
                N[k] = q - k;
                r = q;
                l = k;
            }
        }
    }
}

static void good_suffix(vector<size_t> &Lprime,const vector<size_t> &N, const size_t n){
    for(size_t j = 0; j<n; ++j){
        size_t i = n-N[j];
        Lprime[i] = j;
    }
}

static void sufpref(vector<int> &lprime,const vector<size_t> &N, const size_t n){
    //largest suffix of P[i...n] that is also a prefix of P
    for(size_t i = 0; i<n;++i){
        vector<int> temp;
        for(size_t j = i; j<n; ++j){
            if(j <= n-1+i && N[j]==j){
                temp.push_back(j);
            }
        }
    lprime[i]=*max_element(temp.begin(),temp.end());
    }
}

static void
read_fasta_file_single_sequence(const string &filename, string &T) {

    ifstream in(filename.c_str());

    string line;
    in >> line;
    while (in >> line)
        T.append(line); //reserve
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

    size_t matches=0;
    size_t comparisons =0;

    //Setting up extended bad character rule
    create_Rlist(P,Rlist);

    //Setting up good suffix rule
    vector<size_t> N(n);
    reverseZ(N,P,comparisons);
    cout << "Got the Ns" << endl;

    vector<size_t> Lprime(n,0);
    good_suffix(Lprime,N,n);
    cout << "Got the L's" << endl;

    vector<int> lprime(n,0);
    //largest suffix of P[i...n] that is also a prefix of P
    sufpref(lprime,N,n);
    cout << "Got the l's" << endl;

    size_t k = n;
    while(k<=m){
        size_t i = n; //is it bad to be defining these each loop?
        size_t h = k;
        size_t goodsuf;
        while(i>0 && P[i] == T[h]){
            i = i-1;
            h = h-1;
            //++comparisons;
            cout << "Matching..." << endl;
        }
        if(i==0){
            //cout << "There's a matching occurrence starting at " << k+1 << endl;
            ++matches;
            k=k+n-lprime[1];
            cout << "Found a match" << endl;
        }else{
            //shift P (increase k, position of end of pattern) by the max amount
            //Good suffix rule shift
            if(Lprime[i] > 0){
                goodsuf = n-Lprime[i];
            }else{
                goodsuf = n-lprime[i];
            }
            //Bad character rule
            int badchar = i-lookupR(Rlist,i,T[h]);
            k = max(k+goodsuf,k+badchar); //or k=k+max(goodsuf,badchar)
        }
    }
    cout << "Number of matches is: " << matches << endl;
    //cout << "Number of comparisons was: " << comparisons << endl;
}
