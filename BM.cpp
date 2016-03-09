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

static void create_Rlist(const string &P, Rs &Rlist){
    for(int i=P.size()-1; i>=0; i--){
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
}
static int lookupR(Rs &Rlist, int pos, char bp){ //This is a messy program
    size_t j = 0;
    int newshift = pos;
    cout << "The basepair is " << bp << endl;
    cout << "position in pattern is " << pos << endl;
    if(bp=='A' && !Rlist.a.empty()){
        while(Rlist.a[j]>=pos && (j<Rlist.a.size())){
            ++j;
        }newshift = Rlist.a[j];
    } else if(bp=='C' && !Rlist.c.empty()){
        while(Rlist.c[j]>=pos && (j<Rlist.c.size())){
            ++j;
            }newshift = Rlist.c[j];
    } else if(bp=='G' && !Rlist.g.empty()){
        while(Rlist.g[j]>=pos && (j<Rlist.g.size())){
            ++j;
        }newshift = Rlist.g[j];
    }else if(bp=='T' && !Rlist.t.empty()) {
        while(Rlist.t[j]>=pos && (j<Rlist.t.size())){
            ++j;
        }newshift = Rlist.t[j];
    }else{
        newshift = pos-1;
    }
    cout << "The next position is " << newshift << endl;
    return newshift;
}

static size_t
match(const string &s, size_t q, const size_t n, size_t &comparisons) {
  for (size_t i = n; max(q, i) < s.length() &&
         (toupper(s[i])==toupper(s[q])); ++i, ++q, ++comparisons);
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
    for(int i = 0; i<n;++i){
        vector<int> temp;
        for(int j = i; j<n; ++j){
            if(j <= n-i && N[j]==j){
                temp.push_back(j);
            }
        }
        if(temp.size()>0){
            lprime[i]=*(max_element(temp.begin(),temp.end()));
        }
    }
}

static void
read_fasta_file_single_sequence(const string &filename, string &T) {

    ifstream in(filename.c_str());

    string line;
    //streampos filesize = filename.tellg();
    line.reserve(3e8);
    in >> line;
    while (in >> line)
        T.append(line);
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

    size_t k = n-1;
    while(k < m){
        size_t i = n-1; //is it bad to be defining these each loop?
        size_t h = k;
        size_t goodsuf;
        cout << "P is " << P[i] << " and T is " << T[h] << endl;
        while(i>0 && P[i]==toupper(T[h])){
            i = i-1;
            h = h-1;
            //++comparisons;
        }
        if(i==0){
            ++matches;
            k=k+n-lprime[1];
            cout << "After match, the k is " << k << endl;
        }else{
            //shift P (increase k, position of end of pattern) by the max amount
            //Good suffix rule shift
            cout << "Have to shift" << endl;
            if(Lprime[i] > 0){
                goodsuf = n-Lprime[i]-1;
            }else{
                goodsuf = n-lprime[i]-1;
            }
            cout << goodsuf << " is the good suffix rule shift" << endl;
            //Bad character rule
            int badchar = i-lookupR(Rlist,i,toupper(T[h]));
            cout << badchar << " bad character shifts" << endl;
            cout << "Comparing..." << endl;
            k = max(k+goodsuf,k+badchar); //or k=k+max(goodsuf,badchar)
            cout << "k is now " << k << endl;
        }
    }
    cout << "Number of matches is: " << matches << endl;
    //cout << "Number of comparisons was: " << comparisons << endl;
}
