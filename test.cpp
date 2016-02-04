#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <ctime>

using namespace std;

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

using namespace std;

clock_t begintime = clock();

int wmatch =  3;
int wmism  =  -1;
int windel = -3;

int dmatch =  0;
int dmism  =  1;
int dindel = 1;

map <pair<int,int>,pair<int,int> > optPos;
map <pair<int,int>,int> dist;
string aligned1;
string aligned2;

int alignment_score(int **score,string seq1,string seq2,int n,int m,int c){

    optPos.clear();
    dist.clear();
    for(int i=0;i<=n;++i)
    {
        score[i][0]=-100; //boundary on left as later rows won't initialize
        for(int j=max(0,i-c);j<=min(m,m-n+i+c);++j){
            cout << j << endl;
            if (i==0){
                score[i][j+c+1] = j*windel; //replace with argv[4]
            }else if (j==0){
                score[i][c-i+1] = i*windel;
            }else if (j==m-n+i+c){
                score[i][j+c] = -100;
            }else {
                int t1 = score[i-1][j-(i-c)+1]+ ((seq1[i-1]==seq2[j-1]) ? wmatch:wmism); //argv[2]:argv[3]
                int t2 = score[i][j-(i-c)] + windel;
                int t3 = score[i-1][j-(i-c)+2] + windel;
                score[i][j-(i-c)+1] = max(t1,max(t2,t3));
                cout << i << "," << j << " score is: " << score[i][j-(i-c)+1] << endl;

                //keeping track of the position of optimal alignment and distance
                if(t1==score[i][j-(i-c)+2]){
                    optPos[make_pair(i-1,j-1)] = make_pair(i-2,j-2);
                    dist[make_pair(i-1,j-1)] = dist[make_pair(i-2,j-2)]+((seq1[i-1]==seq2[j-1]) ? dmatch:dmism);
                }else if(t2==score[i][j-(i-c)]){
                    optPos[make_pair(i-1,j-1)] = make_pair(i-1,j-2);
                    dist[make_pair(i-1,j-1)] = dist[make_pair(i-1,j-2)]+dindel;
                }else{
                    optPos[make_pair(i-1,j-1)] = make_pair(i-2,j-1);
                    dist[make_pair(i-1,j-1)] = dist[make_pair(i-2,j-1)]+dindel;
                }
            }
        }
    }

    cout << "getting score at " << n << " , " << c+m-n+1 << endl;
    return score[n][c+m-n+1];
}

int max_seq(string seq1,string seq2,float k)
{
    int n = seq1.size(); //should make a sequence class with names and sizes public instead of doing this
    int m = seq2.size();
    int c = ceil(k/2);

    int **score;
    score = new int*[n+1]; //initializations => score has 1 more row than seq size
    int temp = m-n+ 2*c +2;
    cout << "matrix size is " << temp << endl;

    for(int i=0;i<n+1;++i)
    {
        score[i]=new int[temp]; //band is m-n+2*k wide
    }
    return alignment_score(score,seq1,seq2,n,m,c);
}

void trace_back(string seq1,string seq2)
{
    int n = seq1.size(); //just so I don't have to type this a bunch of times
    int m = seq2.size();

    aligned1 = "";
    aligned2 = "";

    //Starting the traceback from the last element
    vector<pair<int,int> > relevant;
    relevant.push_back(make_pair(n-1,m-1)); //corresponds to sequence index, not score index

    while(relevant[0].first+relevant[0].second != 0){
        relevant.insert(relevant.begin(),
                        optPos[relevant[0]]);
    }
    for(int i=0;i<relevant.size();++i)
    {
        if (relevant[i].first == relevant[i-1].first){
            aligned1.append("-");
            aligned2.push_back(seq2[i]);
        }else if(relevant[i].second == relevant[i-1].second){
            aligned1.push_back(seq1[i]);
            aligned2.append("-");
        }else{
            aligned1.push_back(seq1[i]);
            aligned2.push_back(seq2[i]);
        }
    }
    cout << aligned1.size() << endl;
    cout << aligned2.size() << endl;
}

//argv[1]:string file, argv[2] int wmatch, argv[3] int wmism, argv[4] int windel
int main(int argc, char **argv)
{
    vector <vector<string> > species;

    vector<string> specname;
    vector<string> specseq;

    //Each file has multiple sequences, first line is species name and next lines are sequence
    ifstream openfile;
    //openfile.open(argv[1]);
    openfile.open("HW1-File1.txt");

    string read_lines;
    string appending;

    while (!openfile.eof()) // or, if(!openread.fail())?
    {
        openfile >> read_lines;
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
    }
    //must pushback or will miss the sequence of the final species
    specseq.push_back(appending);
    openfile.close();

    //ofstream resfile;
    //resfile.open("HW2-output-LEngie.txt");

    float k;
    int alignmentscore;
    int distance;

    for(int i=0;i<1;i++){//specname.size();++i){
        for(int j=0;j<1;j++){//specname.size();++j){
                cout << specname[i] << " " << specname[j] << endl;
                cout << "What is the k value for these sequences? Input an integer or -1 if unknown: " << endl;
                cin >> k;
                    if(k == -1){
                        k = 1;
                        int distance = 0;
                        while(k>distance){
                            alignmentscore = max_seq(specseq[i],specseq[j],k);
                            distance = dist[make_pair(specseq[j].size()-1,specseq[j].size()-specseq[i].size()+ceil(k/2)-1)];
                            k=k*2;
                        }
                        cout << "k was " << k/2 << endl;
                    }else{
                        alignmentscore = max_seq(specseq[i],specseq[j],k);
                    }
                cout << "alignment score is " << alignmentscore << endl;
                int temp = specseq[j].size()-specseq[i].size()+ceil(k/2)-1;
                cout << "distance between sequences is " << dist[make_pair(specseq[i].size()-1,temp)] << endl;
                trace_back(specseq[i],specseq[j]);
                //resfile << "first sequence" << endl;
                //resfile << aligned1 << endl;
                //resfile << "second sequence" << endl;
                //resfile << aligned2 << endl;

            }
        }

    //resfile.close();

    clock_t end = clock();
    printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);
    return 0;
}


