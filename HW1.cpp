/*
MATH578a Spring 2016 HW1 - Liana Engie
This program finds the optimal global alignment of two inputted sequences
Input:
    File with FASTA sequences,
    w(match),
    w(mismatch),
    w(indel)
Output:
    r, alignment score
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

using namespace std;

clock_t begintime = clock();

int wmatch =  3;
int wmism  =  -1;
int windel = -3;

map <pair<int,int>,pair<int,int> > optPos;
string aligned1;
string aligned2;

int alignment_score(int **score,string seq1,string seq2){

    optPos.clear();
    for(int i=0;i<=seq1.size();++i)
    {
        for(int j=0;j<=seq2.size();++j){
            if (i==0){
                score[i][j] = j*windel; //replace with argv[4]
            }else if (j==0){
                score[i][j] = i*windel;
                }
            //else if (score[i][j]>0)
                //return score[i][j];
            else {
                //sequence char is off from corresponding score by 1
                int t1 = score[i-1][j-1]+ ((seq1[i-1]==seq2[j-1]) ? wmatch:wmism); //argv[2]:argv[3]
                int t2 = score[i][j-1] + windel;
                int t3 = score[i-1][j] + windel;
                //cout << t1 << " " << t2 << " " << t3 << endl;
                score[i][j] = max(t1,max(t2,t3));

                //keeping track of the position of optimal alignment
                if(t1==score[i][j]){
                    optPos[make_pair(i,j)] = make_pair(i-1,j-1);
                }else if(t2==score[i][j]){
                    optPos[make_pair(i,j)] = make_pair(i,j-1);
                }else{
                    optPos[make_pair(i,j)] = make_pair(i-1,j);}
            }
        }
    }
    return score[seq1.size()][seq2.size()];
}

int max_seq(string seq1,string seq2)
{
    int **score;
    score = new int*[seq1.size()+1]; //initializations mean score has 1 more row than seq size
    for(int i=0;i<seq1.size()+1;++i)
    {
        score[i]=new int[seq2.size()+1];
    }
    return alignment_score(score,seq1,seq2);
}

void trace_back(string seq1,string seq2)
{
    aligned1 = "";
    aligned2 = "";

    //Starting the traceback from the last element
    vector<pair<int,int> > relevant;
    relevant.push_back(make_pair(seq1.size(),seq2.size()));

    while(relevant[0].first + relevant[0].second != 0){
        relevant.insert(relevant.begin(),
                        optPos[relevant[0]]);
        //cout << "The next position is " << optPos[relevant[0]].first << "," << optPos[relevant[0]].second << endl;
        }

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

//argv[1]:string file, argv[2] int wmatch, argv[3] int wmism, argv[4] int windel
int main(int argc, char **argv)
{
    vector <vector<string> > species;

    struct paired{
        string species1;
        string species2;
        int alignmentscore;
        //string optGA;
    };
    paired align;

    vector<string> specname;
    vector<string> specseq;

    //Each file has multiple sequences, first line is species name and next lines are sequence
    ifstream openfile;
    //openfile.open(argv[1]);
    openfile.open("HW1-File2.txt");

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

    /* Testing to see if it worked. Must use a loop; can't cout whole thing
    for(int i=0; i<specname.size();i++){
        cout << specseq[i]<<endl;}*/

    ofstream resfile;
    resfile.open("HW1-output2-LEngie.txt");

    for(int i=0;i<specname.size();++i){
            for(int j=i+1;j<specname.size();++j)
            {
                if(i!=j)
                {
                //Don't actually need to store these, just output
                align.species1 = specname[i];
                align.species2 = specname[j];
                cout << specname[i] << " " << specname[j] << endl;
                align.alignmentscore = max_seq(specseq[i],specseq[j]);
                cout << align.alignmentscore << endl;
                trace_back(specseq[i],specseq[j]);
                resfile << "first sequence" << endl;
                resfile << aligned1 << endl;
                resfile << "second sequence" << endl;
                resfile << aligned2 << endl;
                }
            }
        }

    resfile.close();

    clock_t end = clock();
    printf("Time taken: %.2fs\n", (double)(clock() - begintime)/CLOCKS_PER_SEC);
    return 0;
}
