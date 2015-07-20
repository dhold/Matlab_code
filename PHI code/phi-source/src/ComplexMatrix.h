#ifndef ComplexMatrix_h
#define ComplexMatrix_h
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifndef NOBLAS
    extern "C" {
        #include <cblas.h>
    }

//    #ifdef SINGLEPRECISION
//        #define phi_hemm  cblas_chemm
//        #define phi_symm  cblas_csymm
//        #define phi_her   cblas_cher
//        #define phi_gerc  cblas_cgerc
//        #define phi_copy  cblas_ccopy
//        #define phi_axpy  cblas_caxpy
//        #define phi_dotc_sub cblas_cdotc_sub
//        #define phi_dotu_sub cblas_cdotu_sub
//        #define phi_scal  cblas_zscal
//        #define phi_hemv  cblas_chemv
//        #define phi_gemv  cblas_cgemv
//    #else
//        #define phi_hemm  cblas_zhemm
//        #define phi_symm  cblas_zsymm
//        #define phi_her   cblas_zher
//        #define phi_gerc  cblas_zgerc
//        #define phi_copy  cblas_zcopy
//        #define phi_axpy  cblas_zaxpy
//        #define phi_dotc_sub cblas_zdotc_sub
//        #define phi_dotu_sub cblas_zdotu_sub
//        #define phi_scal  cblas_zscal
//        #define phi_hemv  cblas_zhemv
//        #define phi_gemv  cblas_zgemv
//    #endif

#else

    #ifndef CBLAS_ENUM_DEFINED_H
        #define CBLAS_ENUM_DEFINED_H
        enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
        enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
        enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
        enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
        enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};
    #endif

#endif

using namespace std;

#ifdef SINGLEPRECISION
typedef float Float;
typedef complex<float> Complex;
#else
typedef double Float;
typedef complex<double> Complex;
#endif

struct Index2D {
    int i;
    int j;
};


//#ifdef NOBLAS
void phi_gemm(const int N, const Complex *alpha, const Complex *A,
              const Complex *B, const Complex *beta, Complex *C);

void phi_hemm(const enum CBLAS_SIDE Side, const int N, const Complex *alpha, const Complex *A,
                 const Complex *B, const Complex *beta, Complex *C);

void phi_symm(const enum CBLAS_SIDE Side, const int N, const Complex *alpha, const Complex *A,
                 const Complex *B, const Complex *beta, Complex *C);

void phi_her(const int N, const Float alpha, const Complex *X,
                Complex *A);

void phi_gerc(const int N, const Complex *alpha, const Complex *X,
                 const Complex *Y, Complex *A);
  
void phi_copy(const int N, const Complex *X, const int incX,
                 Complex *Y, const int incY); 

void phi_axpy(const int N, const Complex *alpha, const Complex *X,
                 const int incX, Complex *Y, const int incY); 

void phi_dotc_sub(const int N, const Complex *X, const int incX,
                       const Complex *Y, const int incY, Complex *dotc);

void phi_dotu_sub(const int N, const Complex *X, const int incX,
                       const Complex *Y, const int incY, Complex *dotu);

void phi_scal(const int N, const Complex *alpha, Complex *X, const int incX);

void phi_hemv(const int N, const Complex *alpha, const Complex *A,
                 const Complex *X, const Complex *beta, Complex *Y);

void phi_gemv(const int N, const Complex *alpha, const Complex *A,
                 const Complex *X, const Complex *beta, Complex *Y);

//#endif

//Equal size square matrices only
Complex * outerProduct(const Complex *,const Complex *, const int);

Complex * toLiouville(const Complex *,const int);

void printMatrix(const Complex *, const int);

void printMatrix(const Complex *, const Index2D *, const int);

void writeMatrixToFile(const Complex *, const int, ofstream &);


#endif
