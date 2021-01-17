#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

#define MAX 50

using namespace std;

//defining a structure for represent each edge from a particular vertex.
typedef struct edge{
    int dest;
    int weight;
    struct edge* next;
}edge;

//defining a type for the adjacency list for the graph.
typedef edge* adjList_t[MAX];


void calcRank(vector<pair<int,long double>>& rank,vector<bool>& visited,adjList_t adjList, int i){
    int n = rank.size();

    long double maxVal = 0;
    long double value;
    edge* child = adjList[i];
    while(child){
        if(!visited[child->dest]){
            calcRank(rank,visited,adjList,child->dest);
        }
        
        value = rank[child->dest].second + child->weight;
        maxVal = max(maxVal,value);

        child = child->next;
    }

    visited[i] = true;
    rank[i].second += maxVal;
}

bool compare(pair<int,long double> a, pair<int,long double> b){
    return a.second>b.second;
}

long double CalcEST(int currTask,int currProc,vector<int> Upr,vector<long double> avail,vector<long double> AFT, adjList_t predList){

    if(predList[currTask]==NULL)return 0;

    long double EST = avail[currProc];

    edge* pred = predList[currTask];
    while(pred){
        if(Upr[pred->dest]==currProc){
            EST = max(EST,AFT[pred->dest]);
        }
        else{
            EST = max(EST,AFT[pred->dest] + pred->weight);
        }
        pred = pred->next;
    }

    return EST;
}

long double ISAECCAlgo(vector<pair<int,long double>>& rank,vector<vector<int>>& timeVec,adjList_t adjList,long double EgGiven,vector<long double> Pkind,vector<long double> Ckef,vector<long double> mk,vector<long double> fklow,vector<long double> fkmax,long double& SLG,adjList_t predList){
    int n = rank.size();
    int u = timeVec.size();

    vector<bool> visited(n,false);    
    calcRank(rank,visited,adjList,0);
    sort(rank.begin(),rank.end(),compare);

    vector<long double> Emin(n,LDBL_MAX);
    vector<long double> Emax(n,LDBL_MIN);
    vector<long double> Eave(n);

    long double EminG = 0;
    long double EmaxG = 0;
    long double EaveG = 0;

    for(int i = 0;i<n;i++){
        for(int j = 0;j<u;j++){
            Emin[i] = min(Emin[i],(Pkind[j]+(Ckef[j]*pow(fklow[j],mk[j])))*timeVec[j][i]*(fkmax[j]/fklow[j]));
        }
        for(int j = 0;j<u;j++){
            Emax[i] = max(Emax[i],(Pkind[j]+(Ckef[j]*pow(fkmax[j],mk[j])))*timeVec[j][i]*(fkmax[j]/fkmax[j]));
        }
        Eave[i] = (Emax[i]+Emin[i])/2;
        EminG += Emin[i];
        EmaxG += Emax[i];
    }

    EaveG = (EminG+EmaxG)/2;

    vector<long double> Epre(n);
    vector<long double> E(n);
    long double EieG = EgGiven - EminG;

    long double sum = 0;

    for(int i = 0;i<n;i++){
        Epre[i] = min(Emax[i],(EieG*(Eave[i]/EaveG))+Emin[i]);
        sum+=Epre[i];
    }

    // cout<<"Sum: "<<sum<<endl;

    vector<long double> AFT(n);
    vector<int> Upr(n);
    vector<long double> fprhz(n);
    vector<long double> avail(u,0);

    cout<<"\n";
    cout<<"EGiven(ni):\n";
    for(int i = 0;i<n;i++){
        int index = rank[i].first;

        AFT[index] = LDBL_MAX;

        long double EgGivencurrent = EgGiven;
        for(int j = 0;j<i;j++){
            EgGivencurrent-=E[rank[j].first];
        }
        for(int j = i+1;j<n;j++){
            EgGivencurrent -= Epre[rank[j].first];
        }

        EgGivencurrent = min(EgGivencurrent,Emax[index]);

        cout<<EgGivencurrent<<endl;

        long double Ecurrent = 0;
        long double ESTcurrent = 0;
        long double EFTcurrent = 0;
        for(int j = 0;j<u;j++){
            for(long double f = fklow[j];f<fkmax[j]+0.01;f+=0.01){
                Ecurrent = (Pkind[j]+(Ckef[j]*pow(f,mk[j])))*timeVec[j][index]*(fkmax[j]/f);

                if(Ecurrent>EgGivencurrent)continue;

                ESTcurrent = CalcEST(index,j,Upr,avail,AFT,predList);
                EFTcurrent = ESTcurrent + (timeVec[j][index]*(fkmax[j]/f));

                if(EFTcurrent<AFT[index]){
                    Upr[index] = j;
                    fprhz[index] = f;
                    E[index] = Ecurrent;

                    AFT[index] = EFTcurrent;
                }
            }
        }

        avail[Upr[index]] = AFT[index];
        //This makes sense. After we pass this task, we cannot change any of its allocations.
    }

    long double totEnergy = 0;
    for(int i = 0;i<n;i++){
        totEnergy += E[i];
    }

    SLG = AFT[rank[n-1].first];

    cout<<"\n";
    cout<<"Processors:\n";
    for(int i = 0;i<n;i++){
        cout<<Upr[rank[i].first]+1<<endl;
    }

    cout<<"\n";
    cout<<"Frequency:\n";
    for(int i = 0;i<n;i++){
        cout<<fprhz[rank[i].first]<<endl;
    }

    return totEnergy;
}

int main(){
    
    ifstream in;
    string A;
    
    in.open("Graph.txt");

    getline(in,A);
    istringstream is1(A);

    int n;
    is1>>n;

    int u;
    is1>>u;

    if(!(0<n && n<=MAX)){
        cout<<"Error: The number of tasks is not valid (minimum 1).\n";
        return 1;
    }

    if(!(0<u)){
        cout<<"Error: The number of processors is not valid (minimum 1).\n";
        return 1;
    }

    cout<<"Assumptions:\n";
    cout<<"1. The graph on input is a DAG.\n";
    cout<<"2. The vertices that are present in each row in the file Graph.txt are not redundant.\n";
    cout<<"3. The vertices that are present in each row in the file Graph.txt are within the range 1-"<<n<<" (both inclusive).\n";

    adjList_t adjList;
    adjList_t predList;
    adjList_t currPtr;
    adjList_t currPtrPred;
    for(int i = 0;i<n;i++){
        currPtr[i] = NULL;
        adjList[i] = NULL;
        predList[i] = NULL;
        currPtrPred[i] = NULL;
    }
    
    vector<vector<int> > timeVec(u,vector<int> (n));

    for(int j = 0;j<u;j++){
        getline(in,A);
        istringstream is2(A);
        for(int i = 0;i<n;i++){
            is2>>timeVec[j][i];
        }
    }

    int index;
    int wt;
    int dest;
    int count;
    while(true){
        getline(in,A);
        istringstream is3(A);
        is3>>index;
        index--;

        if(index==-1)break;
        
        //count for current Index.
        is3>>count;

        for(int i = 0;i<count;i++){
            is3>>dest>>wt;

            edge* nextEdge = (edge*)malloc(sizeof(edge));
            //offset from 0.
            nextEdge->dest = dest-1;
            //making edge weights negative.
            nextEdge->weight = wt;
            nextEdge->next = NULL;

            edge* predEdge = (edge*)malloc(sizeof(edge));
            predEdge->dest = index;
            predEdge->weight = wt;
            predEdge->next = NULL;

            if(currPtr[index]){
                currPtr[index]->next = nextEdge;
            }
            else{
                adjList[index] = nextEdge;
            }

            if(currPtrPred[dest-1]){
                currPtrPred[dest-1]->next = predEdge;
            }
            else{
                predList[dest-1] = predEdge;
            }
            currPtr[index] = nextEdge;
            currPtrPred[dest-1] = predEdge;
        }

    }

    cout<<"\n";
    cout<<"Adjacency List:\n";
    for(int i = 0;i<n;i++){
        cout<<i+1<<": ";
        edge* temp = adjList[i];
        while(temp){
            cout<<" ";
            cout<<"("<<(temp->dest)+1<<","<<temp->weight<<")";
            temp = temp->next;
        }

        cout<<"\n";
    }

    vector<long double> Pkind(u);
    vector<long double> Ckef(u);
    vector<long double> mk(u);
    vector<long double> fklow(u);
    vector<long double> fkmax(u);

    for(int i = 0;i<u;i++){
        getline(in,A);
        istringstream is4(A);
        is4>>Pkind[i]>>Ckef[i]>>mk[i]>>fklow[i]>>fkmax[i];
        // cout<<Pkind[i]<<" "<<Ckef[i]<<" "<<mk[i]<<" "<<fklow[i]<<" "<<fkmax[i]<<endl;
    }

    vector<pair<int,long double>> rank(n);
    

    for(int i = 0;i<n;i++){
        rank[i].first = i;
        rank[i].second = 0;

        long double val = 0;

        for(int j = 0;j<u;j++){
            val+=timeVec[j][i];
        }
        val/=u;

        rank[i].second = val;
    }

    cout<<"\n";
    cout<<"Ranks:\n";
    for(int i = 0;i<n;i++){
        cout<<i+1<<" n"<<rank[i].first+1;
        cout<<"\n";
    }

    getline(in,A);
    istringstream is5(A);

    long double EgGiven;
    is5>>EgGiven;

    long double SLG = 0;
    long double energyConsumption = ISAECCAlgo(rank,timeVec,adjList,EgGiven,Pkind,Ckef,mk,fklow,fkmax,SLG,predList);

    // Code for printing the generated adjList.
    // cout<<"Time in processors:\n";
    // for(int j = 0;j<u;j++){
    //     cout<<"Processor "<<j+1<<": ";
    //     for(int i = 0;i<n;i++){
    //         cout<<timeVec[j][i]<<" ";
    //     }
    //     cout<<"\n";
    // }

    cout<<"\n";

    cout<<"Energy Consumed: "<<energyConsumption<<endl;
    cout<<"Schedule Length: "<<SLG<<endl;
    
    in.close();

    return 0;
}