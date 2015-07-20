#ifndef HierarchyIntegrator_h
#define HierarchyIntegrator_h
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <ctime>
#include <sys/time.h>
#include <map>
#include <pthread.h>
#ifdef THREADAFFINITY
#include <sched.h>
#endif
#include "Barrier.h"
#include "ComplexMatrix.h"
#include "HierarchyNode.h"
#include "PhiParameters.h"

#define NDEBUG
#define BROKENINPUT -1
#define NONE 0
#define RK4 1
#define RKF45 2
#define RK4SPECTRUM 3
#define RKF45SPECTRUM 4
#define BICGSTABL 5
#define BICGSTAB  6
#define BICGSTABU 7
#define PRINTHIERARCHY 8
#define PRINTMEMORY 9

using namespace std;

class Heirarchy{
    protected:
	Complex *zeros;
	HierarchyNode *rho, *rho_mid;
	map< int, MultiIndex> I;
	map< MultiIndex, int, CompareMultiIndex> Irev; 
	int MatrixCount;
        int *LMatrixCount;
	int Lmax;

        //COMPILE TIME CONSTANT Ns & Kt
	int Kmax;
	int Ns;
	int Kt;
        int M;
	int Ns2;
        int MKt;
	int Completed;

        int * MDiagInd;

	int *L;
	int find(int *);
	bool TL_trunc;
        bool multiBath;

        bool AssignDesityMatrices;
    public:
	Heirarchy();
	Heirarchy(int, int, int, int);
    	void ConstructI();
	~Heirarchy();
	int verbose;
};



class HierarchyIntegrator : public Heirarchy{
    private:
        int integrationMethod;

	Complex *H, *Hdag, *rho0, *spectrum_rho0;
        Complex **H_th;
        Complex **Hdag_th;
        Complex **rho_matrices_th;
        Complex **rho_mid_matrices_th;
        Complex *TDM;
    
        int numrho0inputlines;

	Float *kappa, *gamma, *lambda, *Vbath_re;
	Float dt, t, t_start;
        Float hbar, T, kB, kBT;

	Complex L_prefactor, negL_prefactor;
	Complex negone,one,two,zero,dtover2, dtover6, dt_c;
	Complex *C;
	Float *absC, maxrho0;
	Float filter_tolerance;
	
	string outFileName,restartFileName,restartIn;
	ofstream outf, restartf;
       
        bool parametersOk; 

        bool spectrum; 
        bool corrBath;
        bool rho_normalized;
	bool rho0_restart;
	int restartStep;

        //Variables for TL Truncation
	Complex **Q_tl, **Q_tldag;
        Complex *Evec, *Eval;
        Complex **Vsum_tl, **Qsum_tl;
        Complex **dE_tl, **freq_tl;
        bool EvecsPresent;
        bool EvalsPresent;

        
        Complex neg_imag_const;
        Complex imag_const;

        Float *Maxddrho;
        Float *RKF45_Diff;

        Float *MaxRhoVal;

        Float *min_dtval;
        Float *max_dtval;
       
        //Tolerances used by RKF45 integration
        //and steady-state solvers
        Float rkf45_mindt;
        Float rkf45_tolerance;
        int N_tolavg_ALL;
        Float bicgstab_tolerance;


        //BiCGSTAB Variables
        Complex *omegaBtm, *omegaTop, *alpha;
        Complex *rho1;
        Float *diff;
        int bicgstabEll;
    
        //Housekeeping variables for threads
        pthread_cond_t MatupdateCond;
        pthread_mutex_t MatupdateMutex;
        pthread_mutex_t io_lock;
        int MatupdateLimit, MatupdateCount;
        pthread_cond_t WriterCond;
        pthread_mutex_t WriterMutex;
        int WriterLimit, WriterCount;
        barrier_t MatupdateBarrier, WriterBarrier;
        int ncores;
        int *cpuaffinity;

        bool stupid_hierarchy_partition;
        int ** MatrixIndices_th;
        int * MatrixCount_th;

    public:
        

        HierarchyIntegrator(int,char* []);

        void readInputParameters(string);
        bool checkInputParameters();

        void run(int, int, int);


        void WriteInputParams(int);
	void initialize(int,bool);

        //Memory assignment
        void printMemoryRequirements(bool);
	void printMemory(int);

	const int GetMatrixCount();
	int allocate_memory(int, int, int);
        void partitionHierarchy(int);
        void partitionHierarchySimple(int);
        void countInterThreadConnections();
  
	void ConstructRho(HierarchyNode *, int, int, Complex *);
        void prepForSteadyState(int, int, int);

	void prepConstructTL();
	void constructTL_BLAS_threaded(Float, int);

	void printHeirarchyConnections(int);

	void rk4_integrate(int, int);
        void rkf45_integrate(int, int);
        //void rkf45_integrate_multistep(int, int);
	void integration_step(const HierarchyNode &, 
                              Complex *,const int &, 
                              Complex *, Complex *, 
                              Complex *, Complex *);
	void integration_step_adaptTrunc(const HierarchyNode &, 
                              Complex *,const int &, 
                              Complex *, Complex *, 
                              Complex *, Complex *);
	void integration_step_fullBath(const HierarchyNode &, 
                                        Complex *,const int &, 
                                        Complex *, Complex *, 
                                        Complex *, Complex *);


	void vec_integration_step(const HierarchyNode &, 
                                  Complex *,const int &,	
                                  Complex *, Complex *, 
                                  Complex *, Complex *);

	void vec_integration_step_fullBath(const HierarchyNode &, 
                                  Complex *,const int &,	
                                  Complex *, Complex *, 
                                  Complex *, Complex *);


        void bicgstab_steadystate(int, int);
        void bicgstab_l_steadystate(int, int);
        void bicgstab_steadystate_unstable(int, int);
	void minimize_step(const HierarchyNode &, 
                              Complex *,const int &, 
                              Complex *, Complex *, 
                              Complex *, Complex *);


        void constructSparseOperator();
        
        void output_timing(int, int, timeval *, int *, int);


	~HierarchyIntegrator();
};



#endif   
