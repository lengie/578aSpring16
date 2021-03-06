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
            //cout << "Added an A at " << i << endl;
        }
        else if(P[i]=='C'){
            Rlist.c.push_back(i);
            //cout << "Added a C at " << i << endl;
        }
        else if(P[i]=='G'){
            Rlist.g.push_back(i);
            //cout << "Added a G at " << i << endl;
        }
        else if(P[i]=='T'){
            Rlist.t.push_back(i);
            //cout << "Added a T at " << i << endl;
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
            newshift = Rlist.a[j], ++j;
        }
    } else if(bp=='C' && !Rlist.c.empty()){
        while(Rlist.c[j]>=pos && (j<Rlist.c.size())){
            newshift = Rlist.c[j], ++j;
            }
    } else if(bp=='G' && !Rlist.g.empty()){
        while(Rlist.g[j]>=pos && (j<Rlist.g.size())){
            newshift = Rlist.g[j], ++j;
        }
    }else if(bp=='T' && !Rlist.t.empty()) {
        while(Rlist.t[j]>=pos && (j<Rlist.t.size())){
            newshift = Rlist.t[j], ++j;
        }
    }else{
        newshift = pos-1;
    }
    cout << "The next position is " << newshift << endl;
    return newshift;
}

static size_t
match(const string &s, size_t q, const size_t n, size_t &comparisons) {
  for (size_t i = n; max(q, i) < s.length() &&
         (s[i]==s[q]); ++i, ++q, ++comparisons);
  return q;
}

static void reverseZ(vector<size_t> &N, const string &P, size_t &comparisons){
    const size_t n = P.size();
    string Pr(n, ' ');
    for(size_t k=0;k<n;++k){
        Pr[n-k-1] = P[k];
    }

    vector<size_t> Z(n,0);
    size_t l = 0, r = 0;
    for (size_t k = 1; k < n; ++k) {
        if (k >= r) {
            Z[k] = match(Pr, 0, k, comparisons);
            if (Z[k] > 0) {
                r = k + Z[k];
                l = k;
            }
        }
        else {
            const size_t k_prime = k - l;
            const size_t beta_len = r - k;
            if (Z[k_prime] < beta_len) {
                Z[k] = Z[k_prime];
            }
            else {
                const size_t q = match(Pr, r, beta_len, comparisons);
                Z[k] = q - k;
                r = q;
                l = k;
            }
        }
    }
    for(int i =0;i<n;++i){
        N[i] = Z[n-i-1];
        cout << "N at " << i << " is " << N[i] << endl;
    }
}

static void good_suffix(vector<int> &Lprime,const vector<size_t> &N, const size_t n){
    for(int j = 0; j<n; ++j){
        int i = n - N[j]-1; //debated a lot about subtracting one from this
        Lprime[i] = j;
        cout << "Lprime at " << i << " is " << Lprime[i] << ", ";
    }
    cout << endl;
}

static void sufpref(vector<int> &lprime,const vector<size_t> &N, const size_t n){
    //largest suffix of P[i...n] that is also a prefix of P
    /*for(int i = 0; i<n;++i){
        vector<int> temp;
        for(int j = i; j<n; ++j){
            if(N[j]==j+1){
                temp.push_back(j);
            }
        }
        if(temp.size()>0){
            lprime[i]=*(max_element(temp.begin(),temp.end()));
            cout << "lprime at " << i << " is "<<lprime[i] << ", ";
        }
    }
    cout << endl;*/
    int i = 0,j = n-1;
    while(j >= 0 && i < n){
        if(N[j]==j+1){
            while(j <= n-i+1 && i < n){
                i++;
                lprime[i] = j+1;
                cout << "lprime at " << i << " is "<<lprime[i] << ", ";
            }
        }
    j--;
    }
    cout << endl;
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
        cerr << "usage: " << argv[0] << " <PATTERN> <TEXT>" << endl;
        return EXIT_FAILURE;
    }
    const string P(argv[1]), filename(argv[2]); //, T(argv[2]);
    const size_t n = P.length();
    string T;
    read_fasta_file_single_sequence(filename, T);
    const size_t m = T.length();

    size_t matches=0;
    size_t comparisons =0;

    //Setting up extended bad character rule
    create_Rlist(P,Rlist);

    //Setting up good suffix rule
    vector<size_t> N(n,0);
    reverseZ(N,P,comparisons);

    vector<int> Lprime(n,0);
    good_suffix(Lprime,N,n);

    vector<int> lprime(n,0);
    sufpref(lprime,N,n);

    int k = n-1;
    while(k < m){
        int i = n-1;
        int h = k;
        int goodsuf;
        cout << "P is " << P[i] << " and T is " << T[h] << endl;
        while(i > -1 && (P[i]==toupper(T[h]))){
            i=i-1;
            cout << i << endl;
            h--;
            //++comparisons;
        }
        if(i==-1){
            ++matches;
            k+=n-lprime[1];
            cout << "After match, the k is " << k << endl;
        }else{
            //shift P (increase k, position of end of pattern) by the max amount
            //Good suffix rule shift
            //cout << "Have to shift" << endl;
            if(Lprime[i] > 0){
                cout << "using L' for shift" << endl;
                goodsuf = n-Lprime[i]-1;
            }else{
                cout << "using l' for shift" << endl;
                goodsuf = n-lprime[i]-1;
            }
            cout << goodsuf << " is the good suffix rule shift" << endl;
            //Bad character rule
            int badchar = i-lookupR(Rlist,i,toupper(T[h]))-1;
            cout << badchar << " bad character shifts" << endl;
            //cout << "Comparing..." << endl;
            k = k+max(1,max(goodsuf,badchar));
            cout << "k is now " << k << endl;
        }
    }
    cout << "Number of matches is: " << matches << endl;
    //cout << "Number of comparisons was: " << comparisons << endl;
}
