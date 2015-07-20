#ifndef HierarchyNode_h
#define HierarchyNode_h
#include <iostream>
#include <map>
#include <string.h>
#include "ComplexMatrix.h"
using namespace std;

const double PI = 3.14159265358979323846264338327950288419716939937510;

struct MultiIndex{
    int *I;
    int N;
    int sum;
    int t; //processing thread;
    int Weight;
    MultiIndex(){
	I = 0;
	N = 0;
	sum = 0;
        Weight = 0;
        t = -1;
    };
    MultiIndex(int m){
	N = m;
	I = new int[m];
	for(int i = 0; i < N; i++)
	    I[i] = 0;
	sum = 0;
        Weight = 0;
        t = -1;
    };
    void Create(int m){
    	N = m;
	if (I == NULL)
	    delete I;
	I = new int[m];
	for(int i = 0; i < N; i++)
	    I[i] = 0;
	sum = 0;
        Weight = 0;
        t = -1;
    };
    void Erase(){
	if (N!=0 || I == NULL)
    	    delete []I;
	I = new int[N];
	for(int i = 0; i < N; i++)
	    I[i] = 0;
	sum = 0;
        Weight = 0;
    };
    int operator[](int n){
	return I[n];
	};

    MultiIndex &operator=(const MultiIndex &rhs){
	if (this == &rhs)
	    return *this;
	if (N!=rhs.N){
	    if(N!=0)    delete []I;
	    else if( I==NULL) delete I;
	    I = new int [N];
    	    N = rhs.N;
	}

        for(int i = 0; i < rhs.N; i++) 
            I[i] = rhs.I[i]; 
        sum = rhs.sum; 
        Weight = rhs.Weight;
        t = rhs.t;
        return *this; 
    }; 
    void Print() const { 
        for(int i = 0; i < N-1; i++) 
            cout << I[i] << ",";
        cout << I[N-1] << endl; 
    }; 
    void inc(int n) { 
        if(n < 0){
            I[-n+1] -= 1; 
            sum -= 1;
        } 
        else if(n > 0){ 
            I[n-1] += 1; 
            sum += 1;
        } 
    }; 
    void dec(int n){ 
        if(n < 0){
            I[-n+1] += 1; 
            sum += 1; 
        } 
        else if (n > 0){ 
            I[n-1] -= 1; 
            sum -= 1; 
        } 
    };
    void updateWeight(){
        Weight = 0;
        for (int i = 0; i < N; ++i) {
            if (I[i] > Weight){
                Weight = I[i];
            }
        }
    }
};

struct CompareMultiIndex{
  bool operator()(const MultiIndex &I1, const MultiIndex &I2) const
  {
    bool res = false;
    if(I1.N > 0 || I2.N > 0){
    	if (I1.N > I2.N)
    	    res = false;
    	else if (I1.N < I2.N)
    	    res = true;
    	else if (I1.sum < I2.sum)
    	    res = true;
    	else if (I1.sum > I2.sum)
    	    res = false;
    	else if (I1.I != NULL && I2.I != NULL){
	    for(int i =0; i < I1.N; i++){
	    	if (I1.I[i] != I2.I[i]){
	    	    res = (I1.I[i] < I2.I[i]);
	    	    return res;
	    	}
	    }
	}
    }
    return res;
  }	
  
};



struct HierarchyNode {
    int *I;  //Index
    int L;   //Heirarchy
    int Ns;  //Number of sites
    int M;   //Number of coupling terms
    int K;   //Number of temperature terms
    int Weight; //Maximum entry in I;
    int id; 

    bool active;

    Complex NZI_prefactor;
    Complex *Same_prefactor;
    Complex *Prev_prefactor_row;
    Complex *Prev_prefactor_col;
    Complex *Next_prefactor; 
    int *SameInd;

    HierarchyNode **Next;
    HierarchyNode **Prev;
    int *PrevInd;       //Entries are 0...M*Kt
    int *PrevIndMatrix; //Entries are 0...M
    int *NextInd;       //Entries are 0...M*Kt
    int *NextIndMatrix; //Entries are 0...0
    int Nprev;
    int Nnext;
    int Ns2;

    bool TL_truncation;	
 
    int size1, size2;
    Complex *rho;
    HierarchyNode();
    void create(int, int, int, int *, Complex *);
    void createVec(int, int, int, int *, Complex *);
    void printI();
    Float max();
    Float maxDiag();
    //HierarchyNode &operator=(const HierarchyNode &);

    //Integration Matrices
    Complex *k0, *k1, *k2, *k3, *k4, *k5;
    Complex *H, *Hdag;

    //Steady-state Liouville-space Matrices:
    Complex **NextLiouville;
    int **NextIndices;
    Complex **PrevLiouville;
    int **PrevIndices;
    Complex *SameLiouville;
    void createSameLiouvilleOperator(const Complex *, const int *, const Index2D *, const int &);
    void createNextLiouvilleOperator(const int &);
    void createPrevLiouvilleOperator(const int &);

    //Steady-state bicgstab
    Complex *r;
    Complex *r0;
    Complex *v;
    Complex *p;
    Complex *t;
    Complex *SteadyStateRho;
    void bicgstabInit();
    void bicgstablInit(int);
    void freeRho();

    Complex dt;
    int step;
};

bool CompAll(int *, int *, int);

Float iter_factorial(int);


#endif
