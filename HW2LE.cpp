/*
MATH578a Spring 2016 HW2 - Liana Engie
This program finds the banded alignment of pairs of sequences from a file with several input sequences
Input:
    File with FASTA sequences,
    w(match),
    w(mismatch),
    w(indel)
Output:
    d, distance between sequences
    T, optimal alignment in a .tex file
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <ctime>
#include <map>
#include <math.h>
#include <stdio.h>

using namespace std;

//clock_t begintime = clock();

int dmatch =  0;
int dmism  =  1;
int dindel = 1;

map <pair<int,int>,pair<int,int> > optPos;
//map <pair<int,int>,int> dist;
string aligned1;
string aligned2;

float k;
int dist;

int distance(int **score,string seq1,string seq2,int n,int m,int c){

    optPos.clear();

    for(int i=0;i<=n;++i)
    {
        score[i][0]= 100; //boundary on left as later rows won't initialize
        int leftlim = max(0,i-c);
        int rightlim = min(m,m-n+i+c);
        for(int j=leftlim; j <= rightlim; ++j){
            int currgap = j-i+c;
            if (i==0){
                score[i][j+c+1] = j*dindel; //replace with argv[4]
                //cout << i << "," << j+c+1 << "a score is: " << score[i][j+c+1] << endl;
            }else if (j==0){
                score[i][c-i+1] = i*dindel;
            }else if (currgap+1>=m-n+2*c+1){
                score[i][currgap+1] = 110;
                //cout << i << "," << j+c << "aa score is: " << score[i][j+c] << endl;
            }else {
                int t1 = score[i-1][currgap+1] + ((seq1[i-1]==seq2[j-1]) ? dmatch:dmism);
                int t2 = score[i][currgap] + dindel;
                int t3 = score[i-1][currgap+2] + dindel;
                score[i][currgap+1] = min(t1,min(t2,t3));
                //cout << i << "," << j-(i-c)+1 << "aaa score is: " << score[i][j-(i-c)+1] << endl;

                //keeping track of the position of optimal alignment and distance
                if(score[i][currgap+1] == t1){
                    optPos[make_pair(i-1,j-1)] = make_pair(i-2,j-2);
                }else if(score[i][currgap+1] == t2){
                    optPos[make_pair(i-1,j-1)] = make_pair(i-1,j-2);
                }else{
                    optPos[make_pair(i-1,j-1)] = make_pair(i-2,j-1);
                }
            }
        }
    }
    /*//Printing the distance matrix to test
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= m-n+2*c+1; j++) {
            printf("%3d ", score[i][j]);
        }
        printf("\n");
    }*/

    return score[n][c+m-n+1];
}

int max_seq(string seq1,string seq2,float k){
    int n = seq1.size(); //should make a sequence class with names and sizes public instead of doing this
    int m = seq2.size();
    int c = ceil(k/2);

    int **score;
    score = new int*[n+1]; //initializations => score has 1 more row than seq size
    int temp = m-n+ 2*c +2;

    for(int i=0;i<n+1;++i)
    {
        score[i]=new int[temp]; //band is m-n+2*k wide
    }
    return distance(score,seq1,seq2,n,m,c);
}

void trace_back(string seq1,string seq2);

int align(string name1, string name2, string seq1, string seq2);

//argv[1]:string file, argv[2] int wmatch, argv[3] int wmism, argv[4] int windel
int main(int argc, char **argv)
{
    vector <vector<string> > species;

    vector<string> specname;
    vector<string> specseq;

    //Each file has multiple sequences, first line is species name and next lines are sequence
    ifstream openfile;
    //openfile.open(argv[1]);
    openfile.open("HW1-File2.txt");

    string read_lines;
    string appending;

    openfile >> read_lines;
    while (!openfile.eof()) // or, if(!openread.fail())?
    {
        /* Note to self: >> is enough to acquire line up to first space.
        Don't use getline because it will make a new line if there's nothing
        on the current line, or may mess with spaces and such
        */

        if(read_lines[0] == '>')
        {
            //Push back the prior species' sequence
            if(appending.size()>0){
                specseq.push_back(appending);
                appending = "";
            }
            specname.push_back(read_lines);
            getline(openfile,read_lines); //just to get the rest of the line out of the way
        }
        else if(read_lines.size() > 0)
        {
            appending.append(read_lines);
        }
        openfile >> read_lines;
    }
    //must pushback or will miss the sequence of the final species
    specseq.push_back(appending);
    openfile.close();

    ofstream resfile;
    resfile.open("HW2-2-output.txt");

    for(int i=0;i<specname.size();++i){
        for(int j=i+1;j<specname.size();++j){
            //The first dimension should be smaller than the second
            if(specseq[i].size()>specseq[j].size()){
                dist = align(specname[j],specname[i],specseq[j],specseq[i]);
                cout << "The distance between sequences is " << dist << endl;

                trace_back(specseq[j],specseq[i]);

            }else{
                dist = align(specname[i],specname[j],specseq[i],specseq[j]);
                cout << "Distance is " << dist << endl;
                //int temp = specseq[j].size()-specseq[i].size()+ceil(k/2);
                //cout << "distance between sequences is " << dist[make_pair(specseq[i].size()-1,specseq[j].size()-1)] << endl;
                trace_back(specseq[i],specseq[j]);
            }
            resfile << "sequence " << i << " and " << j << '\n';
            resfile << aligned1 << '\n';
            resfile << "second sequence of the pair" << '\n';
            resfile << aligned2 << '\n';
        }
    }
    resfile.close();

    //clock_t end = clock();
    //printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);
    return 0;
}

int align(string name1, string name2, string seq1, string seq2){
    cout << name1 << " " << name2 << endl;
    cout << "What is the k value for these sequences? Input an integer or -1 if unknown: " << endl;
    cin >> k;

    if(k == -1){
        k = 1;
        dist = max_seq(seq1,seq2,k);
        while(k<dist){
            k=k*2;
            dist = max_seq(seq1,seq2,k);
        }
        cout << "k was " << k << endl;
        return dist;
    }else{
        return max_seq(seq1,seq2,k);
    }
}


void trace_back(string seq1,string seq2){
    int n = seq1.size(); //just so I don't have to type this a bunch of times
    int m = seq2.size();

    aligned1 = "";
    aligned2 = "";

    //Starting the traceback from the last element
    vector<pair<int,int> > relevant;
    relevant.push_back(make_pair(n-1,m-1)); //corresponds to sequence index, not score index
    //cout << "The first position is " << relevant[0].first << "," << relevant[0].second<<endl;

    while(relevant[0].first + relevant[0].second != 0){
        relevant.insert(relevant.begin(),
                        optPos[relevant[0]]);
        //cout << "The next position is " << optPos[relevant[0]].first << "," << optPos[relevant[0]].second << endl;
        }

    /*//Printing the distance matrix to test
    cout << relevant.size() <<endl;
    for(int i = 0; i < relevant.size(); i++) {
        printf("%3d %3d", relevant[i].first, relevant[i].second);
        printf("\n");
    }*/

    for(int i=0;i<relevant.size();++i){
        if (relevant[i].first == 0 && relevant[i].second == 0){
            aligned1.push_back(seq1[relevant[i].first]);
            aligned2.push_back(seq2[relevant[i].second]);
        }else if (relevant[i].first == relevant[i-1].first){
            aligned1.push_back('-');
            aligned2.push_back(seq2[relevant[i].second]);
        }else if(relevant[i].second == relevant[i-1].second){
            aligned1.push_back(seq1[relevant[i].first]);
            aligned2.push_back('-');
        }else{
            aligned1.push_back(seq1[relevant[i].first]);
            aligned2.push_back(seq2[relevant[i].second]);
        }
    }
    cout << "First sequence aligned: " << endl;
    cout << aligned1 << endl;
    cout << "Second sequence aligned: " << endl;
    cout << aligned2 << endl;
}
//I should probably practice using classes instead of having everything be public
