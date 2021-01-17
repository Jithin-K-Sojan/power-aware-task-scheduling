#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std;

int shiftToBack(vector<int>& Q1,int a,int b){
    int q = Q1.size();
    if((a+b) == q){
        return -1;
    }

    int start = q-b;

    for(int j = 0;j<b;j++){
        int temp = Q1[start];
        Q1[start] = Q1[a];
        Q1[a] = temp;

        a++;
        start++;
    }
    return 1;
}

void getMinNotZero(vector<int>& P,int& min,int& max){
    int p = P.size();
    min = INT_MAX;
    max = INT_MIN;
    for(int i = 0;i<p;i++){
        if(P[i]>max){
            max = P[i];
        }
        if(P[i]==0)continue;
        if(P[i]<min){
            min = P[i];
        }
    }
}

bool backtrack(vector<int>& P,int& minDeadline,int &deadline,int& currE, int E,vector<int>& Q,int alpha,int& minEnergy){

    // This algorithm works because of the way the speeds are distributed.
    int q = Q.size();
    int p = P.size();
    int prevE = currE;
    int max,min;
    vector<int> A(p,0);
    getMinNotZero(P,min,max);
    bool res = false;
    bool tag = false;

    for(int i = 0;i<q;i++){
        if(min<Q[i])continue;
        if(min>=Q[i]){
            for(int j = 0;j<p;j++){
                A[j] = P[j];

                if(P[j]!=0){
                    currE += alpha*(int)pow(Q[i],3)*1;
                    P[j] = P[j]-Q[i];
                }
            }
            deadline++;
            if(min==Q[i] && max==min){
                if(currE<=E){
                    if(deadline<minDeadline){
                        minDeadline = deadline;
                        minEnergy = currE;
                        res = true;
                    }
                    else{
                        minEnergy = currE;
                        tag = true;
                    }
                }
            }
            else{
                res = backtrack(P,minDeadline,deadline,currE, E,Q,alpha,minEnergy);
                if(!res && minEnergy!=INT_MAX){
                    tag = true;
                }
            }

            currE = prevE;
            for(int j = 0;j<p;j++){
                P[j] = A[j];
            }
            deadline--;
            
            if(res){
                return true;
            }

            if(tag){
                return false;
            }
        }

    }

    return false;

}

int main(){

    ifstream in;
    string A;
    
    in.open("inputs.txt");

    getline(in,A);
    istringstream is1(A);

    int p,q,t,E;
    is1>>p>>q>>t>>E;

    int alpha;
    is1>>alpha;

    string A1;

    getline(in,A1);
    istringstream is2(A1);

    int m,n;
    is2>>m;
    is2>>n;

    string B;
    getline(in,B);
    istringstream is3(B);

    vector<int> C(t);
    for(int i = 0;i<t;i++){
        is3>>C[i];
    }

    vector<int> Q(q);
    // cout<<q<<endl;
    int val = q;
    for(int i = 0;i<q;i++){
        Q[i] = val;
        val--;
    }

    struct timeval ts;
    gettimeofday(&ts,NULL);

    vector<int> F(q,0);
    // for(int i = 0;i<q;i++){
    //     F[i] = 0;
    // }

    int sigmaC = 0;
    for(int i = 0;i<t;i++)sigmaC+=C[i];
    
    vector<vector<int>> H(p,vector<int>(q,0));
    vector<int> Q1(q,0);
    for(int i = 0;i<q;i++){
        Q1[i] = Q[i];
        // cout<<Q1[i]<<endl;
    }
    vector<int> Q11(q,0);
    vector<int> F1(q,0);

    int minResidual = INT_MAX;

    int minDeadline = INT_MAX;
    // int totMinEnergy = INT_MAX;
    
    for(int i = 0;i<(int)(pow(p,t));i++){

        // cout<<"i: "<<i<<endl;

        vector<int> P(p,0);

        int temp = i;

        for(int j = 0;j<t;j++){
            int ind = temp%p;
            temp/=p;
            P[ind]+=C[j];
        }

        int deadline = 0;
        int currE = 0;
        int minEnergy = INT_MAX;
        bool res = backtrack(P,minDeadline,deadline,currE,E,Q,alpha,minEnergy);
        
    }

    struct timeval te;
    gettimeofday(&te,NULL);

    double time = (te.tv_sec-ts.tv_sec)*1000 + (te.tv_usec-ts.tv_usec)*1.0/1000;

    cout<<"Deadline: "<<minDeadline<<endl;
    cout<<"Time: "<<time<<endl;


}