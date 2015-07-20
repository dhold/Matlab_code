#include "ComplexMatrix.h"

// Place a copy of A multiplied by each element in B.
// A = NxN complex matrix, B = NxN complex matrix
// returns N^2xN^2 complex matrix
Complex * outerProduct( const Complex * A, const Complex * B, const int N){
    Complex * Result;
    Result = new Complex[N*N*N*N];
    int iiN = 0;
    int jjN = 0;
    for (int ii = 0; ii < N; ++ii, iiN+=N) {
        jjN = 0;
        for (int jj = 0; jj < N; ++jj, jjN+=N) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    Result[(iiN+i)*N*N+jjN+j] = A[i*N+j]*B[iiN+jj];
                }
            }
        }
    } 
    return Result;
}


/////////////////////////////////////////////////////////////////////
//         REPLACEMENT FUNCTIONS IN CASE BLAS NOT AVAILABLE        //
//      THESE ARE SLOW - NO ATTEMPTS AT OPTIMISATION WERE MADE     //
//           THESE ARE *NOT* THE GENERAL BLAS FUNCTIONS:           //
//        THEY **ONLY** APPLY TO SQUARE ROW-MAJOR MATRICES!!!      //
/////////////////////////////////////////////////////////////////////
//void phi_gemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
//                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
//                 const int K, const void *alpha, const void *A,
//                 const int lda, const void *B, const int ldb,
//                 const void *beta, void *C, const int ldc) {
void phi_gemm(const int N, const Complex *alpha, const Complex *A,
              const Complex *B, const Complex *beta, Complex *C) {
#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,alpha,A,N,B,N,beta,C,N);
    #else
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,alpha,A,N,B,N,beta,C,N);
    #endif
#else

    int i,j,n;
       
    //multiply beta:
    //for(i = 0; i < N*N; ++i)
    //    C[i] *= (*beta);
    //cout << "beta= " << *beta << " alpha=" << *alpha << "\n";

    for(i = 0; i < N; ++i) {
        for(j = 0; j < N; ++j) {
            C[i+j*N] *= (*beta);
            for (n = 0; n < N; ++n)
                C[i+j*N] += (*alpha)*(A[i+n*N]*B[j*N+n]);
        }
    } 
#endif

}


void phi_hemm(const enum CBLAS_SIDE Side,
              const int N,
              const Complex *alpha, const Complex *A,
              const Complex *B, const Complex *beta,
              Complex *C) {

#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_chemm(CblasColMajor,Side,CblasUpper,N,N,alpha,A,N,B,N,beta,C,N);
    #else
    cblas_zhemm(CblasColMajor,Side,CblasUpper,N,N,alpha,A,N,B,N,beta,C,N);
    #endif
#else
    
    if (Side == CblasLeft) {
        phi_gemm(N,alpha,A,B,beta,C);
    }
    else if (Side == CblasRight) {
        phi_gemm(N,alpha,B,A,beta,C);
    }

#endif

}

void phi_symm(const enum CBLAS_SIDE Side,
              const int N,
                 const Complex *alpha, const Complex *A,
                 const Complex *B, const Complex *beta,
                 Complex *C){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_csymm(CblasColMajor,Side,CblasUpper,N,N,alpha,A,N,B,N,beta,C,N);
    #else
    cblas_zsymm(CblasColMajor,Side,CblasUpper,N,N,alpha,A,N,B,N,beta,C,N);
    #endif
#else
    if (Side == CblasLeft) {
        phi_gemm(N,alpha,A,B,beta,C);
    }
    else if (Side == CblasRight) {
        phi_gemm(N,alpha,B,A,beta,C);
    }
#endif
}


void phi_gemv(const int N, const Complex *alpha, 
                 const Complex *A, const Complex *X,
                 const Complex *beta, Complex *Y){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_cgemv(CblasColMajor,CblasNoTrans,N,N,alpha,A,N,X,1,beta,Y,1);
    #else
    cblas_zgemv(CblasColMajor,CblasNoTrans,N,N,alpha,A,N,X,1,beta,Y,1);
    #endif
#else
    int i,n;
       
    //multiply beta:
    for(i = 0; i < N; ++i)
        Y[i] *= (*beta);

    for(i = 0; i < N; ++i) {
        for (n = 0; n < N; ++n)
            Y[i] += (*alpha)*A[n*N+i]*X[n];
    } 
#endif
}

void phi_hemv(const int N, const Complex *alpha, const Complex *A,
              const Complex *X,
              const Complex *beta, Complex *Y){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_cgemv(CblasColMajor,CblasNoTrans,N,N,alpha,A,N,X,1,beta,Y,1);
    #else
    cblas_zgemv(CblasColMajor,CblasNoTrans,N,N,alpha,A,N,X,1,beta,Y,1);
    #endif
#else
    phi_gemv(N,alpha,A,X,beta,Y);
#endif
}

void phi_her(const int N, const Float alpha, const Complex *X,
                Complex *A){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_cher(CblasRowMajor,CblasUpper,N,alpha,X,1,A,N);
    #else
    cblas_zher(CblasRowMajor,CblasUpper,N,alpha,X,1,A,N);
    #endif
#else
    int i,j;
    for(i = 0; i < N; ++i){
        for(j = 0; j < N; ++j){
            A[i*N+j] = alpha*X[i]*conj(X[j]);
        }
    }
#endif
}

void phi_gerc(const int N,
                 const Complex *alpha, const Complex *X, 
                 const Complex *Y, Complex *A){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION
    cblas_cgerc(CblasRowMajor,N,N,alpha,X,1,Y,1,A,N);
    #else
    cblas_zgerc(CblasRowMajor,N,N,alpha,X,1,Y,1,A,N);
    #endif
#else
    int i,j;
    for(i = 0; i < N; ++i){
        for(j = 0; j < N; ++j){
            A[i*N+j] = (*alpha)*X[i]*Y[j];
        }
    }
#endif

}
  
void phi_copy(const int N, const Complex *X, const int incX,
                 Complex *Y, const int incY){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION 
    cblas_ccopy(N,X,1,Y,1);
    #else
    cblas_zcopy(N,X,1,Y,1);
    #endif
#else
    int i;
    for(i = 0; i < N; ++i){
        Y[i] = X[i];
    }
#endif
}

void phi_axpy(const int N, const Complex *alpha, const Complex *X,
                 const int incX, Complex *Y, const int incY){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION 
    cblas_caxpy(N,alpha,X,1,Y,1);
    #else
    cblas_zaxpy(N,alpha,X,1,Y,1);
    #endif
#else
    int i;
    for(i = 0; i < N; ++i){
        Y[i] = (*alpha)*X[i]+Y[i];
    }
#endif
}

void phi_dotc_sub(const int N, const Complex *X, const int incX,
                      const Complex *Y, const int incY, Complex *dotc){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION 
    cblas_cdotc_sub(N,X,1,Y,1,dotc);
    #else
    cblas_zdotc_sub(N,X,1,Y,1,dotc);
    #endif
#else
    int i;
    *dotc = 0;
    for(i = 0; i < N; ++i, X+=incX, Y+=incY){
        *dotc += (*X)*conj(*Y);
    }
#endif

}

void phi_dotu_sub(const int N, const Complex *X, const int incX,
                       const Complex *Y, const int incY, Complex *dotu){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION 
    cblas_cdotu_sub(N,X,incX,Y,incY,dotu);
    #else
    cblas_zdotu_sub(N,X,incX,Y,incY,dotu);
    #endif
#else
    int i;
    *dotu = 0;
    for(i = 0; i < N; ++i, X+=incX, Y+=incY){
        *dotu += (*X)*(*Y);
    }
#endif
}

void phi_scal(const int N, const Complex *alpha, Complex *X, const int incX){
#ifndef NOBLAS
    #ifdef SINGLEPRECISION 
    cblas_cscal(N,alpha,X,incX);
    #else
    cblas_zscal(N,alpha,X,incX);
    #endif
#else
    int i;
    for(i = 0; i < N; ++i, X+=incX){
        *X *= (*alpha);
    }
#endif
}


// Computes I_N (*) A - A^H (*) I_N where
// (*) is the outer-product,
// I_N is the NxN identity and
// ^H is the Hermitian conjugate
Complex * toLiouville(const Complex *A, const int N){
    Complex *Result;
    Result = new Complex[N*N*N*N];
    int iiN = 0;
    int jjN = 0;
    for (int ii = 0; ii < N*N*N*N; ++ii)
        Result[ii] = 0;

    for (int ii = 0; ii < N; ++ii, iiN+=N) {
        jjN = 0;
        for (int jj = 0; jj < N; ++jj, jjN+=N) {
            for (int i = 0; i < N; ++i) {
                    Result[(iiN+i)*N*N+jjN+i] = A[iiN+jj];
            }
        }
    } 
    iiN = 0;
    for (int ii = 0; ii < N; ++ii, iiN+=N) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                Result[(iiN+i)*N*N+iiN+j] -= conj(A[j*N+i]);
            }
        }
    } 
    return Result;       
}


//Print a matrix to the screen 
void printMatrix(const Complex *A,const int N){
    for (int i = 0; i < N*N; i+=N) {
        for (int j = 0 ; j < N; ++j){
            cout << "(" << showpos << setw(5) << scientific << setprecision(2) << real(A[i+j]) << ",";
            //if (imag(A[i+j]) < 0) sgn = '-';
            //else sgn = '+';
            cout << setprecision(2) << showpos << scientific << setw(5) << imag(A[i+j]) << ") ";
        }
        cout << endl;
    }
    cout.unsetf(ios_base::fixed);
    cout.unsetf(ios_base::showpos);
};

void printMatrix(const Complex *A,const Index2D *I, const int N){
    
    for (int i = 0; i < N; ++i) {
            cout << setw(5) << I[i].i << " " << I[i].j << " "; 
            cout << "(" << showpos << setw(4) << fixed << setprecision(1) << real(A[i]) << ",";
            cout << setprecision(1) << showpos << fixed << setw(4) << imag(A[i]) << ")\n";
            cout.unsetf(ios_base::fixed);
            cout.unsetf(ios_base::showpos);    
    }

};

void writeMatrixToFile(const Complex *M, const int N, ofstream &outf){
    outf << "#Dimension = " << N << "\n";
    outf << "#Format: i j value\n";
    for(int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            if( abs(M[i*N+j]) > 1e-10 )
                outf << i << " " << j << " " << M[i*N+j] << "\n";
        }
    }
    
};

