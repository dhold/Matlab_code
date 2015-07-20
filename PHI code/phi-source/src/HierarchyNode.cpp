#include "HierarchyNode.h"

bool CompAll(int * I1, int * I2, int N)
{
    bool ret = true;
    for(int i = 0; i < N; i++)
    	ret = ret && (I1[i] == I2[i]);
    return ret;
};


Float iter_factorial(int n) {
    Float accu = 1;
    int i;
    for(i = 1; i <= n; i++) {
        accu *= i;
    }
    return accu;
}




HierarchyNode::HierarchyNode(){
    Ns = 0;
    M = 0;
    L = 0;
    Next = NULL;
    Prev = NULL;
    I = NULL;
    Nprev = 0;
    Nnext = 0;
    PrevInd = NULL;
    active = true;
    NZI_prefactor=0;
    id = 0;
    Weight = 0;
    k0 = k1 = k2 = k3 = k4 = k5 = NULL;
    rho = NULL;
};


void HierarchyNode::create(int N, int m, int k, int *In, Complex * s){
    Ns = N;
    K = k;
    M = m;
    L = 0;
    NZI_prefactor=0;
    active = true;
    if (I != NULL){
	delete []I;
    }
    I = new int[M*K];
    for (int i = 0; i < M*K; i++){
	I[i] = In[i];
	L += I[i];
    }
    //rho = s;
    SameInd = new int[M];
    Weight = 0;
    k0 = k1 = k2 = k3 = k4 = k5 = NULL;

};


void HierarchyNode::createVec(int N, int m, int k, int *In, Complex * s){
    Ns = N;
    K = k;
    M = m;
    L = 0;
    if (I != NULL){
	delete []I;
    }
    I = new int[M*K];
    for (int i = 0; i < M*K; i++){
	I[i] = In[i];
	L += I[i];
    }
    //rho = s;
    SameInd = new int[M];
    Weight = 0;
    k0 = k1 = k2 = k3 = k4 = k5 = NULL;
};

void HierarchyNode::printI(){
    for(int i = 0; i < M*K; i++)
	cout << I[i] << " ";
};

Float HierarchyNode::max(){
    Float maxval = 0;
    for(int i = 0; i < Ns*Ns; i++){
        if (norm(rho[i]) > maxval)
            maxval = norm(rho[i]);
    }
    return sqrt(maxval);
};

/*
HierarchyNode & HierarchyNode::operator=(const HierarchyNode &rhs){
    if (this == &rhs)
	return *this;
    L = rhs.L;
    M = rhs.M;
    K = rhs.K;
    Ns = rhs.Ns;
    Nprev = rhs.Nprev;
    Nnext = rhs.Nnext;
    Weight = rhs.Weight;
    I = new int[M*K];
    if (Nprev > 0){ 
	Prev = new HierarchyNode * [Nprev];
	PrevInd = new int[Nprev];
    }
    if (Nnext > 0) Next = new HierarchyNode * [Nnext];
    rho = new Complex[Ns*Ns];
    for (int i = 0; i < M*K; i++){
	I[i] = rhs.I[i];
	if(i < Nnext) Next[i] = rhs.Next[i];
	if(i < Nprev){
	    Prev[i] = rhs.Prev[i];
	    PrevInd[i] = rhs.PrevInd[i];
	}
    }
    for( int j =0; j < Ns*Ns; j++)
        rho[j] = rhs.rho[j];
    id = rhs.id;
    return *this;
};*/

void HierarchyNode::createSameLiouvilleOperator(const Complex * H_l, const int * DiagNum, const Index2D * DiagBlockInd, const int &Nentries) {

    SameLiouville = new Complex[Nentries];
    int n;
    memcpy(SameLiouville,H_l,Nentries*sizeof(Complex));
    for (int i = 0; i < Ns*Ns; ++i){
        n = DiagNum[i];
        SameLiouville[n] -= NZI_prefactor;
        for (int m = 0; m < M; ++m){ 
            if ( DiagBlockInd[n].i != SameInd[m] && DiagBlockInd[n].j == SameInd[m] )
                SameLiouville[n] -= Same_prefactor[m];
            if ( DiagBlockInd[n].j != SameInd[m] && DiagBlockInd[n].i == SameInd[m] )
                SameLiouville[n] -= Same_prefactor[m];
        }
    }
}

void HierarchyNode::createNextLiouvilleOperator(const int & NumEntries){
    NextLiouville = new Complex * [Nnext];
    int m = 0;
    int n = 0;
    for (int nj = 0; nj < Nnext; ++nj) {
        n = 0;
        m = SameInd[NextIndMatrix[nj]];
        NextLiouville[nj] = new Complex[NumEntries];
        for (int i = 0; i < m; ++i) {
             NextLiouville[nj][n++] = +Next_prefactor[nj];
        }
        for (int i = 0 ; i < m; ++i){
            NextLiouville[nj][n++] = -Next_prefactor[nj];
        }
        for (int i = m+1 ; i < Ns; ++i){
            NextLiouville[nj][n++] = -Next_prefactor[nj];
        }
        for (int i = m+1; i < Ns; ++i) {
            NextLiouville[nj][n++] = +Next_prefactor[nj];
        }
    }
    
};

void HierarchyNode::createPrevLiouvilleOperator(const int & NumEntries){
    PrevLiouville = new Complex * [Nprev];
    int m = 0;
    int n = 0;
    for (int nj = 0; nj < Nprev; ++nj) {
        n = 0;
        m = SameInd[PrevIndMatrix[nj]];
        PrevLiouville[nj] = new Complex[NumEntries];

        for (int i = 0; i < m; ++i) {
           PrevLiouville[nj][n++] = Prev_prefactor_col[nj];
        }
        for (int i = 0; i < m; ++i) {
           PrevLiouville[nj][n++] = -Prev_prefactor_row[nj];
        }
        PrevLiouville[nj][n++] = Prev_prefactor_col[nj]-Prev_prefactor_row[nj];
        for (int i = m+1; i < Ns; ++i) {
           PrevLiouville[nj][n++] = -Prev_prefactor_row[nj];
        }
        for (int i = m+1; i < Ns; ++i) {
           PrevLiouville[nj][n++] = Prev_prefactor_col[nj];
        }
    }
};

void HierarchyNode::freeRho(){
        rho = NULL;
};

void HierarchyNode::bicgstabInit(){
    r = new Complex[Ns*Ns];
    r0 = new Complex[Ns*Ns];
    v = new Complex[Ns*Ns];
    p = new Complex[Ns*Ns];
    t = new Complex[Ns*Ns];
    SteadyStateRho = new Complex [Ns*Ns];

    for(int i = 0; i < Ns*Ns; ++i){
        r[i] = 0;
        r0[i] =0;
        v[i] = 0;
        p[i] = 0;
        t[i] = 0;
        SteadyStateRho[i] = 0;
    }
}

void HierarchyNode::bicgstablInit(int l){
    r = new Complex[(l+1)*Ns*Ns];
    v = new Complex[(l+1)*Ns*Ns];
    r0 = new Complex[Ns*Ns];
    SteadyStateRho = new Complex [Ns*Ns];

    for(int i = 0; i < Ns*Ns; ++i){
        for(int k = 0; k <= l; ++k){
            r[i*(l+1)+k] = 0;
            v[i*(l+1)+k] = 0;
        }
        r0[i] =0;
        SteadyStateRho[i] = 0;
    }
}
