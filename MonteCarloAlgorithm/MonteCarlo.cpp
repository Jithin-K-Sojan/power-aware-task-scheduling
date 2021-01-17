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

int main(){

    ifstream in;
    string A;
    
    in.open("inputs.txt");

    getline(in,A);
    istringstream is1(A);

    int p,q,t,Energy;
    is1>>p>>q>>t>>Energy;

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

    vector<int> QArr(q);
    // cout<<q<<endl;
    int val = q;
    for(int i = 0;i<q;i++){
        QArr[i] = val;
        val--;
    }

    vector<int> deadlineArr(10,0);
    vector<int> energies(10,0);
    vector<double> times(10,0);
    int numSols = 0;

    for(int y = 0;y<10;y++){

        struct timeval ts;
        gettimeofday(&ts,NULL);

        int E = Energy;
        vector<int> Q(q);
        for(int i = 0;i<q;i++)Q[i] = QArr[i];

        srand(y);

        int noSol = 0;

        vector<int> F(q,0);

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

        for(int l = 0;l<m;l++){
            for(int i = 0;i<q;i++)F1[i] = 0;
            int totalRemaining = E;
            int residual = 0;
            int M = 0;

            for(int j = 0;j<q;j++){
                // cout<<Q1[j]<<endl;
                F1[j] = (int)(floor(totalRemaining/(p*alpha*(int)pow(Q1[j],3))));
                totalRemaining = totalRemaining - F1[j]*(int)pow(Q1[j],3)*p;
            }
            
            for(int i = 0;i<q;i++){
                M += F1[i]*Q1[i];
            }
            residual = totalRemaining;

            if(sigmaC<=M*p){
                if(residual<minResidual){
                    for(int i = 0;i<q;i++){
                        F[i] = F1[i];
                        Q11[i] = Q1[i];
                        minResidual = residual;
                    }
                }
            }

            while(1){
                int a = rand()%q;
                int b = rand()%(q-a) + 1;
                for(int i = 0;i<q;i++)Q1[i] = Q[i];
                if(shiftToBack(Q1,a,b)!=-1){
                    break;
                }
            }
        }

        if(minResidual==INT_MAX){
            printf("NO SOLUTION AFTER PART 1.\n");
            noSol = 1;
            // exit(0);
        }

        for(int j = 0;j<q;j++)Q[j] = Q11[j];
        sort(C.begin(),C.end(),greater<int>());

        bool flag = false;
        for(int r = 1;r<=n;r++){

            // Index
            int a = rand()%p;
            // Index
            int b;
            while(1){
                b = rand()%q;
                if((int)floor(F[b]/p)>0){
                    break;
                }
            }
            int c = rand()%((int)floor(F[b]/p)) + 1;
            c *= p;

            cout<<F[b]<<" "<<Q[b]<<" "<<c<<" "<<b<<endl;

            F[b] = F[b] + c;
            for(int j = 0;j<p;j++){
                if(j!=a){
                    H[j][b] = H[j][b]+c+c/p;
                }
            }

            vector<int> Comp(p,0);
            for(int k = 0;k<p;k++){
                for(int z = 0;z<q;z++){
                    Comp[k] += Q[z]*(F[z]-H[k][z]);
                }

            }
            cout<<endl;
            sort(Comp.begin(),Comp.end(),greater<int>());

            bool check = true;

            // Assigning cores
            for(int i = 0;i<t;i++){
                if(Comp[0]<C[i]){
                    check = false;
                    break;
                }
                else{
                    Comp[0] = Comp[0]-C[i];
                    cout<<Comp[0]<<endl;
                    for(int k = 1;k<p;k++){
                        if(Comp[k-1]<Comp[k]){
                            int temp = Comp[k-1];
                            Comp[k-1] = Comp[k];
                            Comp[k] = temp;
                        }
                        else{
                            break;
                        }
                    }
                }
            }
            flag = flag | check;
            if(flag){
                break;
            }
        }
        int deadline = 0;
        if(flag){
            deadline = 0;
            for(int i = 0;i<q;i++){
                deadline += F[i];
            }
            cout<<deadline<<endl;
        }
        else{
            printf("NO SOLUTION PART2\n");
            noSol = 1;
            // exit(0);
        }

        struct timeval te;
        gettimeofday(&te,NULL);

        if(!noSol){
            deadlineArr[numSols] = deadline;
            energies[numSols] += (E-minResidual);
            times[numSols] += (te.tv_sec-ts.tv_sec)*1000 + (te.tv_usec-ts.tv_usec)*(1.0)/1000;
            numSols+=1;
        }
    }

    int deadlineMean = 0;
    double timesMean = 0;
    for(int i = 0;i<numSols;i++){
        deadlineMean+=deadlineArr[i];
        timesMean+=times[i];
    }
    deadlineMean/=numSols;
    timesMean/=numSols;

    double deadlineStd = 0;
    for(int i = 0;i<numSols;i++){
        deadlineStd+=(deadlineArr[i]-deadlineMean)*(deadlineArr[i]-deadlineMean);
    }
    deadlineStd/=(numSols-1)*1.0;
    deadlineStd = pow(deadlineStd,0.5);

    cout<<"Means of the dealine: "<<deadlineMean<<endl;
    cout<<"Std. Dev. of the deadline: "<<deadlineStd<<endl;
    cout<<"Means of the time: "<<timesMean<<endl;
    cout<<"Number of SOLUTIONS: "<<numSols<<endl;

}