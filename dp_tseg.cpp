#include "mex.h"
//#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <iostream>
#define Inf 9999999999.0

using namespace std;
typedef struct C_type{
    double e;
    vector <int> S;
}Cache;
Cache *C;
int * V;
int M,T;

inline int a(int i, int j, int k, int h, int w)
{
    return k*h*w + j*h + i;
}
inline int a2(int i, int j, int h)
{
    return j*h + i;
}

vector<int> dpseg(double & e,
                  const vector<vector<double> > & D, int K, int off)
{
    //printf("size of current D: %d\n", D.size());
    //printf("current K : %d\n",K);
    vector<int> S;
    S.clear();
    if(K>D.size())
    {
        e = Inf;
        return S;
    }
    if(K==1)
    {
        vector<double> m;
        
        for(int i = 0; i<M; i++)
        {
            double mean=0.0;
            for(int j = 0; j<D.size();j++)
            {
                mean+=D[j][i];
            }
            mean/=(double)D.size();
            m.push_back(mean);
        }
        e = 0.0;
        for(int i = 0; i<M; i++)
        {
            double error = 0.0;
            for(int j = 0; j<D.size();j++)
            {
                error+=(D[j][i]-m[i])*(D[j][i]-m[i]);
            }
            error = sqrt(error);
            e+=error;
        }
        return S;
    }
    double e1,e2;
    vector<int>S1,S2;
    // let's assume m is empty
    double min_err = Inf;
    for(int j = 0; j<D.size()-1; j++) {
        int check = a(off,off+j,K-2,T,T);
        if (V[check]==1) {
            e1 = C[check].e;
            S1 = C[check].S;
        } else {
            vector<vector<double> >::const_iterator first;
            vector<vector<double> >::const_iterator last;
            first = D.begin();
            last = D.begin()+j+1;
            vector<vector<double> > D1(first,last);
            S1 = dpseg(e1,D1,K-1,off);
	    //cout<<K-1<<endl;

            V[check]=1;
            C[check].e = e1;
	    //cout<<"<->"<<endl;
	    //cout<<S1.size()<<endl;
	    //cout<<C[check].S.size()<<endl;
	    //cout<<"->"<<endl;
            C[check].S = vector<int>(S1);
	    //cout<<"<-"<<endl;
        }
        check = a(off+j+1,off+D.size()-1,0,T,T);
        if(V[check]==1) {
            e2 = C[check].e;
            S2 = C[check].S;
        } else {
            vector<vector<double> >::const_iterator first;
            vector<vector<double> >::const_iterator last;
            first = D.begin()+j+1;
            last = D.end();
            vector<vector<double> > D2(first,last);
            S2 = dpseg(e2,D2,1,off+j+1);
            V[check]=1;
            C[check].e = e2;
            C[check].S = S2;
        }
        if(e1+e2<min_err) {
            min_err = e1+e2;
            // todo: merge S1 and S2
            S.clear();
            for(int i = 0; i<S1.size(); i++) {
                S.push_back(S1[i]);
            }
            S.push_back(j+1);
            for(int i = 0; i<S2.size(); i++) {
                S.push_back(S2[i]);
            }
        }
    }
    e = min_err;
    return S;
}
void mexFunction(int nlhs, mxArray*plhs[],int nrhs, const mxArray*prhs[]) {
    double *A = mxGetPr(prhs[0]);
    double *dK = mxGetPr(prhs[1]);
    M = mxGetM(prhs[0]);
    T = mxGetN(prhs[0]);
    int K = (int)(dK[0]);
    V = (int*)malloc(sizeof(int)*K*T*T);
    for(int i = 0; i<K*T*T; i++) {
        V[i]=0;
    }
    C = new Cache [K*T*T];
    //C = (Cache *)malloc(sizeof(Cache)*K*T*T);
    for(int i = 0; i<K*T*T;i++)
	    C[i].S = vector<int>();
    vector<vector<double> > D;
    // fill D using A
    for(int i = 0; i<T; i++) {
        vector<double> sample;
        for(int j = 0; j<M; j++) {
            sample.push_back(A[a2(j,i,M)]);
        }
        D.push_back(sample);
    }
    double e;
    //printf("here");
    vector<int> S = dpseg(e, D, K, 0);
    //printf("after");
    
    plhs[0] = mxCreateDoubleMatrix(K+1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *ret = mxGetPr(plhs[0]);
    if(S.size()!=K-1) {
        printf("%d %d\n",S.size(),K-1);
        mexErrMsgTxt("not valid length");
    }

    for(int i = 0; i<S.size();i++) {
        ret[i+1] = (double)S[i];
    }
    ret[0] = 1;
    ret[K] = T;
    ret = mxGetPr(plhs[1]);
    ret[0] = e;
    //printf("%lf\n",e);
    delete [] C;
    free(V);
	
    //free(C);
    D.clear();
    return;
}
