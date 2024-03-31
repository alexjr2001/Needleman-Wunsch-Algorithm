//Student: Alexander Gómez Del Carpio
//Subject: Computational Molecular Biology
//Description: In this piece of code we can see a use of dynamic programming in order to get the
//alignment of sequences in DNA or RNA with Needleman-Wunsch algorithm, whether you need the best
//result or all possible results.
//Date: 21/03/2024

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>


using namespace std;

//Little class to define a cell 
class Cell{
    public:
        int score;                      //Score depending of previous
        vector<int> xy;                 //Store indexes of previous cells with the greater score
        Cell():score(0),xy(0,0){}
};


//Two functions to manage the txt file
void eraseSpacesNumbers(string& s){         //Erase digit and letter to get the sequence
    for(auto it=s.begin();it!=s.end();){
        if(isdigit(*it) || isspace(*it)){
            s.erase(it);
        }
        else ++it;
    }
    return;
}


vector<string> sequenceFileManager(){           //Returns clean sequences in a vector
    ifstream inputFile("Sequencias.txt");
    if (!inputFile.is_open()) {
        cout << "Can't open" << endl;
    }
    vector<string> sequences;
    string cur_seq="";
    string line; 
    while (getline(inputFile, line)) {
        if(!isalpha(line[0])){
            cur_seq+=line;
        }
        else{
            if(!cur_seq.empty()){
                eraseSpacesNumbers(cur_seq);
                sequences.push_back(cur_seq);
            } 
            cur_seq="";
        }
    }
    eraseSpacesNumbers(cur_seq);
    sequences.push_back(cur_seq);

    return sequences;
}


void printMatrix(vector<vector<Cell>>& M){      //Prints the matrix of DP
    for(auto i:M){
        for(auto j:i){
            if(j.score>=0) cout<<" "; 
            cout<< j.score<<" ";
        }
        cout<<endl;
    }
}

void initialFill(vector<vector<Cell>>& M){      //Initialize the matrix first column and row
    for(int i=1;i<M.size();i++){
        M[i][0].score=i*-2;
        M[i][0].xy={i-1,0};
    }
    for(int i=1;i<M[0].size();i++){
        M[0][i].score=i*-2;
        M[0][i].xy={0,i-1};
    }
}


//Calculate score and previous cell coordinates with best result
pair<int,vector<int>> calculateScore(const vector<vector<Cell>>& M, int i, int j, int diff_score){
    int scores[3];
    scores[0]=M[i-1][j-1].score+diff_score;     //0 -> Diagonal
    scores[1]=M[i-1][j].score-2;                //1 -> Up
    scores[2]=M[i][j-1].score-2;                //2 -> Left

    vector<int> paths;             //Stores 0,1,2 according best score
    vector<int> paths_idx;         //Stores coordinates previous cells with best score
    int max=INT_MIN;

    for(int k=0; k<3; k++){        //Get cell with biggest score
        if(scores[k]>max){
            max=scores[k];
            paths.clear();
            paths.push_back(k);
        }
        else if(scores[k]==max){
            paths.push_back(k);
        }
    }

    for(auto k:paths){              //Give the coordinates of the best score cell
        if(k==0){
            paths_idx.push_back(i-1);
            paths_idx.push_back(j-1);
        }
        else if(k==1){
            paths_idx.push_back(i-1);
            paths_idx.push_back(j);
        }
        else{
            paths_idx.push_back(i);
            paths_idx.push_back(j-1);
        }
    }
    
    return make_pair(max,paths_idx);

}

//Global variables in order to not use too many arguments
int minScore= INT_MAX;                  //Counts how many breaks for gaps occur, less better
pair<string,string> optimas={"",""};    //Store optimal solution of two sequences    
unordered_map<string, int> visited;     //In case the path is already visited with its previous score
string state;                           //Which coordinate we visit to store in visited
pair<bool,bool> gapAtTheEnd={true,true};       //In case the sequences have gaps at the end in order to not count them as ruptures

//Function that finds the optimal solution
void rebuiltPath(const vector<vector<Cell>>& M, const string& a,const string& b, int i, int j, char ins1,char ins2, vector<string> ans, int scoreGaps){
    ans[0].insert(ans[0].begin(),ins1);     //Insert solution at the beginning
    ans[1].insert(ans[1].begin(),ins2);

    if(ins1=='-' && !gapAtTheEnd.first) scoreGaps++;    //In case if there's not gap at the end anymore
    else if(ins1!='-') gapAtTheEnd.first=false;
    if(ins2=='-' && !gapAtTheEnd.second) scoreGaps++;
    else if(ins2!='-') gapAtTheEnd.second=false;

    //Check if is visited and store the new path in visited
    string state=to_string(i)+","+to_string(j);

    if(visited.find(state) != visited.end()){       //We ignore the visited path only if has a lower score
        if(scoreGaps>=visited[state]) return;       
    }

    visited[state]=scoreGaps;
    

    if(i==0 && j==0){       //In case the path arrives to the first cell (the end), verify if it is the best
        if(ans[0][0]=='-') scoreGaps--;
        if(ans[1][0]=='-') scoreGaps--;
        if(scoreGaps<=minScore){
            minScore=scoreGaps;
            optimas.first=ans[0];
            optimas.second=ans[1];
            cout<<optimas.first<<endl<<optimas.second<<endl<<endl<<endl;
        }

        return;
    }

    for(int k=0;k<M[i][j].xy.size();k+=2){      //See what continues in the path
        ins1=a[i];  //Ins1 and 2 next insertions
        ins2=b[j];
        if(M[i][j].xy[k]==i) ins1='-';      //In case is vertical or horizontal there's a gap
        else if(M[i][j].xy[k+1]==j) ins2='-';
        rebuiltPath(M,a,b,M[i][j].xy[k],M[i][j].xy[k+1],ins1,ins2,ans,scoreGaps);
    }

}


//It rebuilds all possible paths not just the optimal, we ignore the gaps
void rebuiltAllPaths(const vector<vector<Cell>>& M, const string& a,const string& b, int i, int j, char ins1,char ins2, vector<string> ans){
    ans[0].insert(ans[0].begin(),ins1);
    ans[1].insert(ans[1].begin(),ins2);
    
    if(i==0 && j==0){
        cout<<ans[0]<<endl<<ans[1]<<endl;
        return;
    }

    for(int k=0;k<M[i][j].xy.size();k+=2){
        ins1=a[i];
        ins2=b[j];
        if(M[i][j].xy[k]==i) ins1='-';
        else if(M[i][j].xy[k+1]==j) ins2='-';
        rebuiltAllPaths(M,a,b,M[i][j].xy[k],M[i][j].xy[k+1],ins1,ins2,ans);
    }

}

void writeResultTxt(){                  //Writes in a new txt file the output
    ofstream outputFile("output.txt");
    if (!outputFile) {
        cout << "Doesn't open" <<endl;
    }
    outputFile << optimas.first.substr(0,150) << std::endl;
    outputFile << optimas.second.substr(0,80) << std::endl;
    outputFile.close();
    return;
}


int main(){
    vector<string> sequences=sequenceFileManager();

    //ALWAYS ADD A GAP FIRST TO HAVE THE RIGHT MATRIX
    //Example sequences
    string a="-AAAC";
    string b="-AGC";

    //In case we want sequences from the txt file
    a="-"+sequences[0];
    b="-"+sequences[1];

    //Matrix init
    vector<vector<Cell>> matrix(a.size(), vector<Cell>(b.size(),Cell()));
    initialFill(matrix);
    
    int diff_score;     //If we come from diagonal cell, how much is the penalty for being different
    pair<int,vector<int>> receiveScore;     //Variable receives the result of every cell.

    //We go from second column and row and see what is the previous cell with the best score
    for(int i=1;i<a.size();i++){
        for(int j=1;j<b.size();j++){
            a[i]==b[j]?diff_score=1:diff_score=-1;
            receiveScore=calculateScore(matrix,i,j,diff_score);
            matrix[i][j].score=receiveScore.first;
            
            matrix[i][j].xy=receiveScore.second;
            
        }
    }

    vector<string> initial_ans(2,"");   //Store results during recursion

    //Uncomment the function that you want, the first one gives you all possible paths and the second
    //one gives you the best alignment

    //rebuiltAllPaths(matrix,a,b,matrix.size()-1,matrix[0].size()-1,' ',' ',initial_ans);
    rebuiltPath(matrix,a,b,matrix.size()-1,matrix[0].size()-1,' ',' ',initial_ans,0);

    writeResultTxt();       //Writing result in txt to plot the point matrix


    return 0;
}

