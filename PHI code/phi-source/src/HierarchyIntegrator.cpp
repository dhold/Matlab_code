#include "HierarchyIntegrator.h"


Heirarchy::Heirarchy(){
    Ns = 0;
    Kt = 0;
    M = 0;
    Lmax = 0;
    MatrixCount = 0;
    rho = NULL;
    L = NULL;
};

Heirarchy::Heirarchy(int mxin, int k, int L_in, int vb){
    Ns = mxin;
    Kt = k+1;
    Lmax = L_in;
    MatrixCount = 0;
    rho = NULL;
    L = NULL;
    ConstructI();
 };

void Heirarchy::ConstructI(){
    if (verbose > 0)
      cout << "Constructing Index graph.\n";
    //if (vebose > 0)
    //  cout << "Lmax = " << Lmax << "\n";
  
    MultiIndex I0(M*Kt), I1(M*Kt); //Dummy indices to populate I[.]

    L = new int[Lmax];
    LMatrixCount = new int[Lmax];
    for(int i = 0; i < Lmax; i++)
	L[i] = 0;
    
    L[0] = 1;
    LMatrixCount[0] = 1;

    Irev[I0] = 0;
    I[0].Create(M*Kt);
    I[0] = I0;
    if (verbose > 1) cout << "\nHeirarchy Level: " << 0 << "\n";
    if (verbose > 2) cout << "" << MatrixCount << ": ";
    if (verbose > 2) I[0].Print();

    MatrixCount++;
    if( Lmax >= 2){
        if (verbose > 1) cout << "\nHeirarchy Level: " << 1 << "\n";
	L[1] = M*Kt;
        LMatrixCount[1] = 1+M*Kt;
	for(int i = 0; i < M*Kt; i++){
	    I[i+1].Create(M*Kt);
	    I[i+1] = I[0];
	    I[i+1].inc(i+1);
            I[i+1].updateWeight();
	    Irev[I[i+1]] = i+1;
	    if (verbose > 2) cout << "" << MatrixCount << ": ";

	    if (verbose > 2) I[i+1].Print();
	    MatrixCount++;
	}
       if (verbose > 1) cout << "\t" << L[1] << " Matrices\n";
    }

    // FOLLOWING ORDERS
    if (Lmax > 2) {
        int I_prev_n = 0;
        for(int Ln = 2; Ln < Lmax; Ln++){
            if (verbose > 1) cout << "\nHeirarchy Level: " << Ln << "\n";
            //FOR EACH POSITION
            for(int k = 0; k < MKt; k++){	
                //FOR EACH ENTRY IN PREVIOUS LEVEL
                for(int n = 0; n < L[Ln-1]; n++){
                    I_prev_n = MatrixCount - L[Ln-1] + n;
                    //GET PREVIOUS I
                    I1 = I[I_prev_n];
                    //UPDATE I[n+k] by 1
                    I1.inc((n+k)%MKt+1);
                    I1.updateWeight();

    //		//CHECK FOR SAME I in CURRENT LEVEL OF HEIRARCHY
                    //IF NOT FOUND INSERT
                    if(Irev.count(I1) == 0){
                        I[MatrixCount+L[Ln]].Create(MKt); 
                        I[MatrixCount+L[Ln]] = I[I_prev_n];
                        I[MatrixCount+L[Ln]].inc((n+k)%MKt+1);
                        I[MatrixCount+L[Ln]].updateWeight();
                        if (verbose > 2) cout << "" << MatrixCount+L[Ln] << ""
                                              << ": ";
                        
                        Irev[I[MatrixCount+L[Ln]]] = MatrixCount+L[Ln];
                        if (verbose > 2) I[MatrixCount+L[Ln]].Print();
                        L[Ln]++;
                    }	
                }
            }
            if (verbose > 1) cout << "\t" << L[Ln] << " Matrices\n";
            MatrixCount += L[Ln];
            LMatrixCount[Ln] = LMatrixCount[Ln-1]+L[Ln];
        }
    }
    //I NOW CONTAINTS MATRIX HEIRARCHY.
    // The above looping will create e.g. a 2-hierarchy with numbering:
    //                           1
    //                        2     3
    //                     4     6     5
    //                  7     9     10    8 
    //               11   13     15    14   12 
    //            16   18    20     21   19   17 
    //         22   24    26    28     27   25   23

    if (verbose > 0)
      cout << "\nTotal: " <<  MatrixCount << " matrices\n\n";
}



Heirarchy::~Heirarchy(){
    delete []rho;
    delete []L;
};

bool fexists(const char *filename)
{
    ifstream ifile(filename);
    return ifile;
}

struct thread_data {
    HierarchyIntegrator *E;
    int id;
    int num_threads;
    int method;
};

void *launchThread(void *p){
    struct thread_data *Tdata;
    Tdata = (struct thread_data *) p;
    (static_cast<HierarchyIntegrator*>(Tdata->E)->run)(Tdata->id,Tdata->num_threads,Tdata->method);
    return 0; 
}

/////////////////////////////////
////  Begin HierarchyIntegrator class ////
/////////////////////////////////

//    Constructor         
//    ------------
//    Read in parameters. 
//    Set up some matrices and constants.

HierarchyIntegrator::HierarchyIntegrator(int num_args,char* argv[]) : Heirarchy() {
    //Some constants
    hbar = 5.3088354;
    kB = 0.6950344;
    L_prefactor = Complex(0,-1./hbar);
    negL_prefactor = Complex(0,1./hbar);
    one = Complex(1,0);
    negone = Complex(-1,0);
    neg_imag_const = Complex(0,-1);
    imag_const = Complex(0,1);
    two = Complex(2,0);
    zero = Complex(0,0);
    verbose=1;

    spectrum = false;
    integrationMethod = NONE;



    int num_threads = 1;
    stringstream strStream;    

    // Threading
    pthread_t *tg;
    size_t stacksize;

    cout << "PHI Parallel Hierarchy Integrator by Johan Strumpfer, 2009-2012.\n"
         << "Integration of density matrix using hierarchy equations of motion.\n\n"
         << "Please cite Strumpfer and Schulten, J. Chem. Theor. Comp. (2012)\n"
         << "in all publications reporting results obtained with PHI\n\n"
#ifdef NOBLAS
         << "Using internal matrix functions\n"
#else
         << "Using BLAS for matrix functions\n"
#endif        
#ifdef SINGLEPRECISION
        << "Using single precision floating point operations\n\n";
#else
        << "Using double precision floating point operations\n\n";
#endif
//         << "Using RK4 or RKF45 based adaptive timestepping algorithm.\n\n";
    if (num_args < 2){
	cout << "Usage: \n\n  phi paramfile {integrator [threads],print}\n"
             << "    integrator  - One of rk4, rkf45, rk4spectrum, rkf45spectrum, steadystate.\n"
             << "    threads     - Number of threads\n"
             << "    memory      - Prints memory requirements\n"
             << "    print       - Prints hierarchy connections \n";
    }
    else{
	cout << "Reading parameters from " << argv[1] << "\n\n";
        if (fexists(argv[1])) {        
            readInputParameters(argv[1]);
            if (num_args > 2) {
                if (num_args==4) {
                    strStream.str(argv[3]);
                    strStream >> num_threads;
                }
                if (num_threads < 1) {
                    cerr << "Error: threads must be greater than 0.\n";
                }
                else {
                    /* Set up integration if needed */
                    if ((strcmp(argv[2],"rk4spectrum")) == 0){
                        integrationMethod = RK4SPECTRUM;
                        initialize(num_threads,true);
                    }
                    else if ((strcmp(argv[2],"rkf45spectrum")) == 0){
                        integrationMethod = RKF45SPECTRUM;
                        initialize(num_threads,true);
                    }
                    else if ((strcmp(argv[2],"rk4")) == 0){
                        integrationMethod = RK4;
                        initialize(num_threads,true);
                    }
                    else if ((strcmp(argv[2],"rkf45")) == 0){
                        integrationMethod = RKF45;
                        initialize(num_threads,true);
                    }
                    else if(strcmp(argv[2],"print")==0){
                        integrationMethod = PRINTHIERARCHY;
                        initialize(1,false);
                        printHeirarchyConnections(num_threads);
                    }
                    else if(strcmp(argv[2],"memory")==0){
                         integrationMethod = PRINTMEMORY;
                         initialize(1,false);
                         if (verbose <= 0)
                           printMemoryRequirements(num_threads);
                    }
                    else if (strcmp(argv[2],"bicgstabl")==0){
                        integrationMethod = BICGSTABL;
                        initialize(num_threads,true);
                    }
                    else if (strcmp(argv[2],"bicgstab")==0){
                        integrationMethod = BICGSTABL;
                        initialize(num_threads,true);
                    }
                    else if (strcmp(argv[2],"bicgstabu")==0){
                        integrationMethod = BICGSTABU;
                        initialize(num_threads,true);
                    }
                    /* Now launch threads */
                    switch (integrationMethod) {
                     case RK4:
                     case RKF45:
                     case RK4SPECTRUM:
                     case RKF45SPECTRUM:
                     case BICGSTABL:
                     case BICGSTABU:
                     case BICGSTAB: 
                        {
                            if (verbose > 0)
                              cout << "\nLaunching threads...\n";
                            pthread_attr_t attr;
                            pthread_attr_init(&attr);
                            /* Set Max Stacksize =  8MB on os x */
                            stacksize = 8388608;
                            int err = pthread_attr_setstacksize(&attr,stacksize);
                            if (err) {
                                    cout << "\tSet Stacksize Error=" << err << endl;
                            }
                            /* Ensure Joinable; */
                            pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
                            thread_data *tgData;
                            /* num_threads integrator threads + 1 writer */
                            tg = new pthread_t[num_threads+1];
                            tgData = new thread_data[num_threads+1];
                            /* and launch */
                            for(int i = 0; i <= num_threads; i++){
                                    tgData[i].E = this;
                                    tgData[i].id = i;
                                    tgData[i].num_threads = num_threads;
                                    tgData[i].method = integrationMethod;
                                    pthread_create(&tg[i],&attr,launchThread,&tgData[i]);
                            }	   
                            /* Wait for all to finish */
                            void *status;
                            for(int i = 0; i <= num_threads; i++){
                                    pthread_join(tg[i],&status);
                            }
                        }
                        break;
                      default:
                        break;
                    }
                    cout << "Done.\n";
                } // end if more than 0 threads
            } // end if more than 2 parameters
        } // end if input file exists
    } // end if more than 1 parameter   
}


void HierarchyIntegrator::readInputParameters(string ParmFile)  {

    //DEFAULTS
    T = 300;
    verbose = 1;
    rho0_restart = false;
    maxrho0 = 0;	
    rho_normalized = true;
    multiBath = false;
    Ns = 0; //System Size
    M = 0;  //# of System-Bath coupling matrices
    Kt = 1;  //MatsubaraTerms
    Lmax = 0; //HierarchyTruncation
    rkf45_mindt = 1e-6;
    rkf45_tolerance = 1e-6;
    filter_tolerance = 0; //  Shi et al 2009 truncationa
    bicgstab_tolerance = 1e-6;
    dt = 0.001;
    t = 0;
    kappa = NULL;
    gamma = NULL;
    lambda = NULL;
    H = NULL;
    Hdag = NULL;
    rho0 = NULL;
    spectrum_rho0 = NULL;
    corrBath = false;
    TL_trunc = false;
    EvecsPresent=false;
    EvalsPresent=false;
    restartStep = 0;
    bicgstabEll = 2;
    stupid_hierarchy_partition=false;
    cpuaffinity = NULL;
  
    MatupdateLimit = 0;
    MatupdateCount = 0;
    WriterLimit = 0;
    WriterCount = 0;
    ncores = 0;

    parametersOk = true;

    ifstream inPRM;
    inPRM.open(ParmFile.c_str());
    string line, curvarName, curvar;
    stringstream strStream;
    int sepPos,Hi=0,Hj=0,iprev = 0,n=0;

    //Counters
    numrho0inputlines = 0;
    int nH = 0;
    int nVbath = 0;
    int nEigvec = 0;
    int nEigval = 0;
    int nLambda = 0;
    int nGamma = 0;
    int nKappa = 0;

    Complex cbuf3;
    getline(inPRM,line);

    while(!inPRM.eof() && parametersOk){
        while (line.find("#") != string::npos){
            curvar = "";
            getline(inPRM,line);        
            
        }
	if (line.find("=") != string::npos){
	    sepPos = line.find("=");
	    curvarName = line.substr(0,sepPos); 
	    curvar = line.substr(sepPos+1,line.length());
	}
	else if(line.find(":") != string::npos){
	    sepPos = line.find(":");
	    curvarName = line.substr(0,sepPos); 
	    getline(inPRM,curvar);
	    Hi = 0;
	    Hj = 0;
	}
	else if (line.size() > 0) {
	    curvar = line;
	}
        else {
          curvarName = "";
        }
	getline(inPRM,line);
	strStream.clear();	
	if(curvarName ==  "fout" || curvarName == "Output" || curvarName == "OutputFile"){
	    outFileName = curvar;
	    if (verbose > 0)
		cout << "Output being written to " << curvar << endl;
            //Add code to make backup of file if it exists
	}
	else if(curvarName ==  "Ns" || curvarName == "N" || curvarName == "NumStates"){
	    strStream.str(curvar);
	    strStream >> Ns;
            if (Ns <= 0) {
                cerr << "ERROR: NumStates must be positive.\n";
                parametersOk = false;
            }
            else {
              if (verbose > 0)
                cout << "Number of states = " << Ns << endl;
              Ns2 = Ns*Ns;
              H = new Complex[Ns2];
              Hdag = new Complex[Ns2];
              Evec = new Complex[Ns2];
              Eval = new Complex[Ns];
              rho0 = new Complex[Ns2];
              kappa = new Float[Ns];
              TDM = new Complex[Ns*3];
              for(int i=0; i < Ns; i++){
                  kappa[i] = 0;
                  for(int j = 0; j < Ns; j++){
                      H[i*Ns+j] = Complex(0,0);
                      rho0[i*Ns+j] = Complex(0,0);
                  }
                  for(int j = 0; j < 3; j++)
                      TDM[j] = Complex(0);
              }
              if (M == 0) {
                  M = Ns;
              } 
              Vbath_re = new Float[M*Ns];
              lambda = new Float[M];
              gamma = new Float[M];
              MDiagInd = new int[M];
              for(int i=0; i < M; i++){
                  MDiagInd[i] = i;
                  lambda[i] = 0;
                  gamma[i] = 0;
                  for(int j = 0; j < Ns; j++)
                      Vbath_re[i*Ns+j] = 0;
              } 
           }
	}
        else if (curvarName == "M" || curvarName == "NumCouplingTerms"){
    	    strStream.str(curvar);
	    strStream >> M;
            if (M <= 0) {
                cerr << "ERROR: NumCouplingTerms must be positive.\n";
                parametersOk = false;
            }
            else{
              if (verbose > 0)
                cout << "Number of system-bath coupling terms = " << M << endl;
              if (Ns != 0) {
                  delete []Vbath_re;
                  delete []lambda;
                  delete []gamma;
                  delete []MDiagInd;
                  Vbath_re = new Float[M*Ns];
                  lambda = new Float[M];
                  gamma = new Float[M];
                  MDiagInd = new int[M];
                  for(int i=0; i < M; i++){
                      MDiagInd[i] = i;
                      lambda[i] = 0;
                      gamma[i] = 0;
                      for(int j = 0; j < Ns; j++)
                          Vbath_re[i*Ns+j] = 0;
                  }
              }
           }
        }
	else if(curvarName == "Kt" || curvarName == "MatsubaraTerms"){
            parametersOk = getVal(curvar,curvarName,Kt,PARAMNONNEGATIVE);
            if (parametersOk) {
              if (verbose > 0)
                  cout << "Number of Matsubara terms = " << Kt << endl;
              Kt += 1;
            }
	}
	else if(curvarName ==  "Lh" || curvarName == "HierarchyTruncation"){
            parametersOk = getVal(curvar,curvarName,Lmax,PARAMNONE);
            if (parametersOk) {
              if (verbose > 0)
                cout << "Hierarchy truncation = " << Lmax << endl;
              if (Lmax < 1) {
                  cerr << "ERROR: HierarchyTruncation must be greater than 1.\n";
                  parametersOk = false;
                  break;
              }
            }
	}
	else if(curvarName ==  "dt" || curvarName == "Timestep"){
            parametersOk = getVal(curvar,curvarName,dt,PARAMPOSITIVE);
            if (parametersOk) {
              if (verbose > 0)
                cout << "Timestep = " << dt << " ps.\n";
                dtover2 = Complex(dt/2,0);
                dtover6 = Complex(dt/6,0);
                dt_c = Complex(dt,0);
            }
            if (dt < 1e-10) {
              cout << "Warning: dt less than 1e-10 ps.\n";
            }
	}
	else if(curvarName ==  "t" || curvarName ==  "Time"){
            parametersOk = getVal(curvar,curvarName,t,PARAMPOSITIVE);
            if (parametersOk) {
              if (verbose > 0) 
                cout << "Runlength = " << t << " ps.\n";
            }
	}
	//  Shi et al 2009 truncation
        else if(curvarName == "filter" || curvarName == "filtered" || curvarName == "FilterTolerance"){
            parametersOk = getVal(curvar,curvarName,filter_tolerance,PARAMPOSITIVE);
            if (parametersOk) {
              if (verbose > 0) 
                cout << "Setting ADM Filter tolerance to " << filter_tolerance << ".\n";
            }
	}
	else if(curvarName == "TL" || curvarName == "TimeLocal"){
            parametersOk = getVal(curvar,curvarName,TL_trunc);
            if (parametersOk) {
              if (verbose > 0 && TL_trunc) cout << "Using time local truncation.\n";
              if (verbose > 0 && !TL_trunc) cout << "Using time non-local truncation.\n";
            }
	}
	else if(curvarName == "rho_normalized" || curvarName == "rho_normalised" ){
	    strStream.str(curvar);
	    int v;
	    strStream >> v;
            if (v != 0 && v != 1) {
                cerr << "ERROR: rho_normalized must be 0 or 1.\n";
                parametersOk = false;
                break;
            }
            else {
              if (v==1) rho_normalized = true;
              if (v==0) rho_normalized = false;
              if (verbose > 0 && rho_normalized) cout << "Using rho hierarchy scaling.\n";
            }
	}
	else if(curvarName ==  "H" || curvarName == "Hamiltonian"){
            if (Ns <= 0) {
                cerr << "ERROR: Ns must be specified before Hamiltonian.\n";
                parametersOk = false;
                break;
            }
              
            if (nH < Ns){
              if(!getLineOfVal(curvar,(H+nH),Ns,Ns)) {
                cerr << "ERROR: incorrect Hamiltonian input format.\n";
                parametersOk = false;
                break;
              }
              nH++;
            }
            else {
              cerr << "ERROR: Hamiltonian already specified\n";
                parametersOk = false;
                break;
            }
	}
        else if(curvarName == "TDM"){
            if (Ns <= 0) {
                cerr << "ERROR: Ns must be specified before TDM.\n";
                parametersOk = false;
                break;
            }
              
            if (!getLineOfVal(curvar,TDM+Hi,3,3)) {
              cerr << "ERROR: incorrect TDM input format.\n";
                parametersOk = false;
                break;
            }
	    Hi++;
        }
	else if(curvarName ==  "rho0" || curvarName == "InitialDensityMatrix"){
            if (Ns <= 0) {
                cerr << "ERROR: Ns must be specified before InitialDensityMatrix.\n";
                parametersOk = false;
                break;
            }
              
            if (numrho0inputlines < Ns) {
              if (!getLineOfVal(curvar,(rho0+numrho0inputlines),Ns,Ns)) {
                cerr << "ERROR: incorrect InitialDensityMatrix input format.\n";
                parametersOk = false;
                break;
              }
              for(int j = 0; j < Ns; ++j)
                maxrho0 = ( abs(rho0[numrho0inputlines+Ns*j]) > maxrho0 ) ? abs(rho0[numrho0inputlines+Ns*j]) : maxrho0;
              numrho0inputlines++;
            }
            else {
              cerr << "ERROR: rho0 already specified.\n";
                parametersOk = false;
                break;
            }
            
	}
        else if (curvarName == "spectrum_rho0" && spectrum_rho0 == NULL) {

            if (Ns <= 0) {
                cerr << "ERROR: NumStates must be specified before spectrum_rho0.\n";
                parametersOk = false;
                break;
            }
            if (verbose > 0){
                cout << "Specified spectrum calculation initial state.\n";
            }
            spectrum_rho0 = new Complex[Ns];      
        
            if (!getLineOfVal(curvar,spectrum_rho0,Ns,1)) {
              cerr << "ERROR: incorrect spectrum_rho0 input format.\n";
              parametersOk = false;
            }
        }
        else if(curvarName == "Multibath" || curvarName == "multiBath" || curvarName == "DiagonalCouplingIndices") {
            if (Ns <= 0) {
                cerr << "ERROR: NumCouplingTerms must be specified before " << curvarName << ".\n";
                parametersOk = false;
                break;
            }
            if (!getLineOfVal(curvar,MDiagInd,M,1)){
              cerr << "ERROR: incorrect  " << curvarName << " input format.\n";
                parametersOk = false;
                break;
            }
            //Check indices less than Ns
            
        }
        else if(curvarName == "Vbath" || curvarName == "CorrelatedCouplingTerms"){ //correlated System-Bath coupling
            if (Ns <= 0) {
                cerr << "ERROR: NumCouplingTerms must be specified before " << curvarName << ".\n";
                parametersOk = false;
                break;
            }
            if (nVbath < M) {
              if (!getLineOfVal(curvar,Vbath_re+M*nVbath,Ns,1)) {
                cerr << "ERROR: incorrect " << curvarName << " input format.\n";
                parametersOk = false;
                break;
              }
              nVbath++;
            } 
            else {
                cerr << "ERROR: too many " << curvarName << " input lines.\n";
                parametersOk = false;
                break;
            }
	}
        // TIME-LOCAL TRUNCATION
	else if(curvarName ==  "eigenvectors" || curvarName ==  "Eigenvectors" || curvarName == "HamiltonianEigenvectors"){
            if (Ns <= 0) {
                cerr << "ERROR: Ns must be specified before Eigenvectors.\n";
                parametersOk = false;
                break;
            }
              
            if (nEigvec < Ns) {
              if (!getLineOfVal(curvar,(Evec+nEigvec),Ns,Ns)) {
                cerr << "ERROR: incorrect Eigenvectors input format.\n";
                parametersOk = false;
                break;
              }
              nEigvec++;
              EvecsPresent = true;
              if (verbose > 0 && nEigvec == Ns)
                  cout << "Read eigenvectors.\n";
            }
            else {
              cerr << "ERROR: Eigenvectors already specified.\n";
                parametersOk = false;
                break;
            }
	}
	else if(curvarName ==  "eigenvalues" || curvarName ==  "Eigenvalues" || curvarName == "HamiltonianEigenvalues" ){
            if (Ns <= 0) {
                cerr << "ERROR: Ns must be specified before Eigenvalues.\n";
                parametersOk = false;
                break;
            }
              
            if (nEigval < 1) {
              if (!getLineOfVal(curvar,(Eval),Ns,1)) {
                cerr << "ERROR: incorrect Eigenvalues input format.\n";
                parametersOk = false;
                break;
              }
              nEigval++;
              //CHECK FOR COMPLEX EIGENVALUES!!!
              EvalsPresent = true;
              if (verbose > 0)
                  cout << "Read eigenvalues.\n";
            }
            else {
              cerr << "ERROR: Eigenvalues already specified.\n";
              parametersOk = false;
              break;
            }
	}
	else if(curvarName == "RestartInput"){
	    strStream.str(curvar);
	    strStream >> restartIn;
	    rho0_restart = true;
	    if (verbose > 0) cout << "Taking initial state from " << restartIn << ".\n";
	}
	else if(curvarName ==  "lambda" || curvarName ==  "Lambda" ){
            if (M <= 0) {
                cerr << "ERROR: Ns must be specified before " << curvarName << ".\n";
                parametersOk = false;
                break;
            }
              
            if (nLambda < 1) {
              if (!getLineOfVal(curvar,lambda,M,1)){
                cerr << "ERROR: incorrect " << curvarName << " input format.\n";
                parametersOk = false;
                break;
              }
              nLambda++;
              if (verbose > 0)
                  cout << "Read reorganization energies lambda.\n";
              if (verbose > 1) {
                for(int i = 0; i < M; ++i) {
                  cout << "lambda["<<i<<"]="<<lambda[i] << endl;
                }
              }
            }
            else {
              cerr << "ERROR: " << curvarName << " already specified.\n";
                parametersOk = false;
                break;
            }
	}
	else if(curvarName ==  "gamma" || curvarName == "Gamma"){
            if (M <= 0) {
                cerr << "ERROR: Ns must be specified before " << curvarName << ".\n";
                parametersOk = false;
                break;
            }
              
            if (nGamma < 1) {
              if (!getLineOfVal(curvar,gamma,M,1)) {
                cerr << "ERROR: incorrect " << curvarName << " input format.\n";
                parametersOk = false;
                break;
              }
              nGamma++;
              if (verbose > 0)
                  cout << "Read response frequencies gamma.\n";
              if (verbose > 1) {
                for(int i = 0; i < M; ++i) {
                  cout << "gamma["<<i<<"]="<<gamma[i] << endl;
                }
              }
            }
            else {
              cerr << "ERROR: " << curvarName << " already specified.\n";
              parametersOk = false;
              break;
            }
	}
	else if(curvarName ==  "kappa"){
            if (M <= 0) {
                cerr << "ERROR: Ns must be specified before " << curvarName << ".\n";
                parametersOk = false;
                break;
            }
              
            if (nKappa < 1) {
              if (!getLineOfVal(curvar,kappa,M,1)){
                cerr << "ERROR: incorrect " << curvarName << " input format.\n";
                parametersOk = false;
                break;
              }
              nKappa++;
              if (verbose > 0)
                  cout << "Read Kappa.\n";
              if (verbose > 1) {
                for(int i = 0; i < M; ++i) {
                  cout << "kappa["<<i<<"]="<<kappa[i] << endl;
                }
              }
            }
            else {
              cerr << "ERROR: " << curvarName << " already specified.\n";
              parametersOk = false;
              break;
            }
	}
	else if(curvarName == "T" || curvarName == "Temperature" ){
	    strStream.str(curvar);
	    strStream >> T;
	    if (verbose > 0) cout << "Temperature = " << T << " Kelvin.\n";
	}
	else if(curvarName == "verbose"){
	    strStream.str(curvar);
	    strStream >> verbose;
	}
	else if(curvarName == "restartFile"){
	    strStream.str(curvar);
	    strStream >> restartFileName;
	    if (verbose > 0) cout << "Restart File = " << restartFileName << endl;
	}
	else if(curvarName == "restartStep"){
	    strStream.str(curvar);
	    strStream >> restartStep;
	    if (verbose > 0 && restartStep > 0) cout << "Writing restart file every " << restartStep << " steps.\n";
	}	    
	else if(curvarName == "tol_low" || curvarName == "RKF45mindt"){
	    strStream.str(curvar);
	    strStream >> rkf45_mindt;
	    if (verbose > 0) cout << "Minimum RKF45 timestep = " << rkf45_mindt << " ps\n";
            //compatibility
            if (curvarName == "tol_low") {
              bicgstab_tolerance = rkf45_mindt;
              if (verbose > 0) cout << "bicgstab tolerance = " << bicgstab_tolerance << "\n";
            }

	}	    
       	else if(curvarName == "tol_high" || curvarName == "RKF45tolerance"){
	    strStream.str(curvar);
	    strStream >> rkf45_tolerance;
	    if (verbose > 0) cout << "RKF45 integration tolerance = " << rkf45_tolerance << "\n";
	}
	else if(curvarName == "bicgstab_tolerance"){
	    strStream.str(curvar);
	    strStream >> bicgstab_tolerance;
	    if (verbose > 0) cout << "bicgstab tolerance = " << bicgstab_tolerance << "\n";
	}
        else if(curvarName == "bicgstab"){            
	    strStream.str(curvar);
            strStream >> bicgstabEll;
            if (bicgstabEll < 0) {
                cout << "bicgstab cannot be less than 0 assuming bicgstab="<< -bicgstabEll << endl;
                bicgstabEll *= -1 ;
            }
            if (verbose > 0) {
                if (bicgstabEll > 0)
                    cout << "BiCGSTAB(" << bicgstabEll << ") will be used for steadystate.\n";             
                else if (bicgstabEll == 0)
                    cout << "BiCGSTAB will be used for steadystate.\n";             
            }
        }
        else if (curvarName == "partition") {
            if (curvar == "stupid" || curvar == "no" || curvar == "stoopid")
                stupid_hierarchy_partition = true;
            if (verbose > 0 && stupid_hierarchy_partition){
                cout << "Using simple hierarchy partition (round-robin).\n";
            }
        }
#ifdef THREADAFFINITY
        else if (curvarName == "cpuaffinity" && cpuaffinity == NULL) {
            ncores = 1;
	    for(int i = curvar.find(",",0); i != string::npos; i = curvar.find(",",i+1)){
                ncores++;
            }
            cpuaffinity = new int[ncores];
	    n = 0;
	    iprev = 0;
	    for(int i = curvar.find(",",0); i != string::npos; i = curvar.find(",",i+1)){
		strStream.str(curvar.substr(iprev,i));
		strStream >> cpuaffinity[n];
		if (verbose > 1) cout << "cpuaffinity["<<n<<"]="<<cpuaffinity[n] << endl;
		n++;
		iprev = i+1;
		strStream.clear();
	    }
	    strStream.str(curvar.substr(iprev,curvar.length()-iprev));
	    strStream >> cpuaffinity[n];
	    if (verbose > 1) cout << "cpuaffinity["<<n<<"]="<<cpuaffinity[n] << endl ;            
        }
#endif
            
    }

    inPRM.close();

    if (verbose > 1)
        cout << "Input file read.\n";
}

bool HierarchyIntegrator::checkInputParameters() {

  string INPERR = "INPUT ERROR ";
  bool retval = true;
  if (Ns <= 0 ) {
    cerr << INPERR << "NumStates must be > 0\n";
    retval = false;
  }
  if ( Lmax < 1 ) {
    cerr << INPERR << "HierarchyTruncation must be > 1\n";
    retval = false;
  }
  if ( M <= 0 ) {
    cerr << INPERR << "NumCouplingTerms must be > 0\n";
    retval = false;
  }
  if ( Kt < 0 ) {
    cerr << INPERR << "NumMatsubaraTerms must be >= 0\n";
    retval = false;
  }
  if ( H == NULL ) {
    cerr << INPERR << "Hamiltonian must be specified\n";
    retval = false;
  }
  if ( gamma == NULL ) {
    cerr << INPERR << "Gamma values must be specified\n";
    retval = false;
  }
  if ( lambda == NULL ) {
    cerr << INPERR << "Lambda values must be specified\n";
    retval = false;
  }
  if (rho0 == NULL && spectrum_rho0 == NULL && rho0_restart == false) {
    cerr << INPERR << "Initial state must be specified\n";
    retval = false;
  }
  return retval;
}

void HierarchyIntegrator::printMemory(int num_threads){
    printMemoryRequirements(false);
};


void HierarchyIntegrator::printMemoryRequirements(bool fileout){
        float memreq = 0;
        float rk4 = 0;
        float rkf45 = 0;
        float spectrum_rk4 = 0;
        float spectrum_rkf45 = 0;
        float steadystate = 0;
        float matrices = 0;

        matrices =  (float)(sizeof(Complex)*Ns*Ns*MatrixCount)/1024/1024;
        //HierarchyNode rho and rho_mid

        memreq += sizeof(int)*(M+4+M+M+M+M+5)+sizeof(Complex)*(1+M+M+M+M+Ns*Ns);
        memreq *= 2*MatrixCount;
        //Integration Variables
        memreq += sizeof(Complex)*(3+Ns2)+sizeof(int)*6+sizeof(time_t)*2;
        //Make sure Vbath matrices can be taken into account
        // Kt+2 * M matrices sized Ns*Ns 
        if (corrBath) {
            memreq += sizeof(Complex)*M*Ns*Ns*(Kt+2);
        }
        rk4 = (float)(sizeof(Complex)*Ns*Ns*MatrixCount*4 + memreq)/1024/1024;

        rkf45 += sizeof(Complex)*Ns*Ns*MatrixCount*6;
        rkf45 += sizeof(Complex)*31;
        rkf45 += sizeof(Float)*7+sizeof(int)*2;
        rkf45 = ((float)(memreq+rkf45))/1024/1024;

        //Spectrum
        //HierarchyNode rho and rho_mid
        memreq = 0;
        memreq += sizeof(int)*(M+4+M+M+M+M+5)+sizeof(Complex)*(1+M+M+M+M+Ns);
        memreq *= 2*MatrixCount;
        //Integration Variables
        memreq += sizeof(Complex)*(3+Ns)+sizeof(int)*6+sizeof(time_t)*2;

        spectrum_rk4 += sizeof(Complex)*(Ns+1)*MatrixCount*4; 
        spectrum_rk4 = ((float)(spectrum_rk4+memreq))/1024/1024;

        spectrum_rkf45 += sizeof(Complex)*3*Ns*MatrixCount*6;
        spectrum_rkf45 += sizeof(Complex)*31;
        spectrum_rkf45 += sizeof(Float)*7+sizeof(int)*2;
        spectrum_rkf45 = ((float)(spectrum_rkf45+memreq))/1024/1024;
        
        //SteadyState
        steadystate = 0;
        steadystate += sizeof(int)*(M+4+M+M+M+M+5)+sizeof(Complex)*(1+M+M+M+M); // Hierarchy
        steadystate += 6*(sizeof(Complex)*(Ns*Ns)); // BiCGSTAB vectors & SteadyStateRho
        steadystate *= MatrixCount;
        steadystate += 9*sizeof(Complex) + 3*sizeof(int) + 2*sizeof(Float);
        steadystate = steadystate/1024/1024;

        if(!fileout) {
            cout << "Memory to store matrices for calculations:\n";
            cout << "\thierarchy: " << matrices << " Mb.\n";
            cout << "\trk4: " << rk4 << " Mb.\n";
            cout << "\trkf45: " <<  rkf45 << " Mb.\n";
            cout << "\trk4spectrum: " << spectrum_rk4 << " Mb.\n";
            cout << "\trkf45spectrum: " << spectrum_rkf45 << " Mb.\n";
            cout << "\tsteadystate: " << steadystate << " Mb.\n";
        }
        else {
            outf << "#Memory to store matrices for calculations:\n";
            outf << "#hierarchy: " << matrices << " Mb.\n";
            outf << "#rk4: " << rk4 << " Mb.\n";
            outf << "#rkf45: " <<  rkf45 << " Mb.\n";
            outf << "#rk4spectrum: " << spectrum_rk4 << " Mb.\n";
            outf << "#rkf45spectrum: " << spectrum_rkf45 << " Mb.\n";
            outf << "#steadystate: " << steadystate << " Mb.\n";
        }
}



void HierarchyIntegrator::WriteInputParams(int threads){
    outf.open((outFileName).c_str());
    double version = 1.0;
//    hbar = 5.29;	    
    outf << "#PHI " << version <<"\n";
    outf << "#Parameters:\n";
    outf << "#-----------\n";
    outf << "#NumStates=" << Ns << "\n";
    outf << "#MatsubaraTerms=" << Kt-1 << "\n";
    outf << "#CouplingTerms=" << M << "\n";
    outf << "#Temperature=" << T << "\n";
    outf << "#HierarchyTruncation=" << Lmax << "\n";
    outf << "#TL="<< TL_trunc<<"\n";
    outf << "#dt=" << dt << "\n#t=" << t << endl;
    outf << "#filter=" << filter_tolerance << endl;
    outf << "#RKF45tolerance = " << rkf45_tolerance << endl;
    outf << "#RKF45mindt = " << rkf45_mindt << endl;
    outf << "#gamma:\n";
    outf << "#";
    for( int i = 0; i < M-1; i++)
	outf << gamma[i] << ",";
    outf << gamma[M-1] <<endl;
    outf << "#lambda:\n";
    outf << "#";
    for( int i = 0; i < M-1; i++)
	outf << lambda[i] << ",";
    outf << lambda[M-1] << endl;
    outf << "#kappa:\n";
    outf << "#";
    for( int i = 0; i < Ns-1; i++)
	outf << kappa[i] << ",";
    outf << kappa[Ns-1] << endl;
    outf << "#H:\n";
    for( int i = 0; i < Ns; i++){
	outf << "#";
	for(int j = 0; j < Ns-1; j++)
	    outf << H[i+j*Ns] << ",";
	outf << H[i+Ns*(Ns-1)] << endl;
    }
    if (TL_trunc) {
      outf << "#Eigenvalues:\n";
      outf << "#";
      for(int j = 0; j < Ns-1; j++)
          outf << Eval[j] << ",";
      outf << Eval[(Ns-1)] << endl;
      outf << "#Eigenvectors:\n";
      for( int i = 0; i < Ns; i++){
          outf << "#";
          for(int j = 0; j < Ns-1; j++)
              outf << Evec[i+j*Ns] << ",";
          outf << Evec[i+Ns*(Ns-1)] << endl;
      }
    }
    //outf << "#TDM:\n";
    //for (int i = 0; i < Ns; i++){
    //    outf << "#";
    //    for (int j = 0; j < 2; ++j)
    //        outf << TDM[i*Ns+j] << ",";
    //    outf << TDM[i*Ns+2] << "\n";
    //}

    if (restartStep > 0){
    	restartf.open(restartFileName.c_str());
	//outf << "#Restart File=" << restartFileName << "\n";
	//outf << "#Restart file being written every " << restartStep << " steps.\n";
    }
    /*if(!rho0_restart){
	outf << "#rho0:\n";
        for( int i = 0; i < Ns; i++){
	    outf << "#";
	    for( int j = 0; j < Ns-1; j++)
		outf << rho0[i+Ns*j] << ",";
	    outf << rho0[i+Ns*(Ns-1)] << endl;
	}
    }*/
    if (multiBath) {
        outf << "#multiBath:\n";
        outf << "#";
        for( int i = 0; i < M-1; i++){
                outf << MDiagInd[i] << ",";
        }
        outf << MDiagInd[M-1] << endl;
    }        
    if (corrBath){
        outf << "#Vbath:\n";
        for( int i = 0; i < Ns; i++){
            outf << "#";
            for(int j = 0; j < Ns-1; j++)
                outf << Vbath_re[i*Ns+j] << ",";
            outf << Vbath_re[i*Ns+(Ns-1)] << endl;
        }


    }
    outf << "#-------------\n";
    outf << "#Total Number of matrices= " << MatrixCount << endl;
    printMemoryRequirements(true);
    if (integrationMethod == RK4 || integrationMethod == RKF45 || integrationMethod == RK4SPECTRUM || integrationMethod == RKF45SPECTRUM) {
        outf << "#Threads=" << threads << endl;
        outf << "#Integrator=";
        if (integrationMethod==RK4) outf << "rk4\n";
        else if (integrationMethod==RKF45) outf << "rkf45\n";
        else if (integrationMethod==RK4SPECTRUM) outf << "rk4spectrum\n";
        else if (integrationMethod==RKF45SPECTRUM) outf << "rkf45spectrum\n";
    } 
    else if (integrationMethod == BICGSTAB || integrationMethod == BICGSTABU || integrationMethod == BICGSTABL) {
        outf << "#Threads=" << threads << endl;
        if (integrationMethod == BICGSTABL)
            outf << "#Calculating steady-state using BiCGSTAB("<< bicgstabEll <<")\n";
        else
            outf << "#Calculating steady-state using BiCGSTAB\n";
    }

};

void HierarchyIntegrator::printHeirarchyConnections(int num_threads){
    int MCminus = 0;
    int vb = verbose;
    verbose = -1;
    for (int i = 0; i < num_threads; ++i){
        ConstructRho(rho,i,1,NULL);
    }
    verbose = vb;
    cout << "Hierarchy connections:\n";
    for(int i = Lmax - 1; i >= 0; i--){
	MCminus += L[i];
    	cout << "Level " << i << ":\n";
	for(int m = L[i]-1; m >= 0; m--){
	   cout << "\t" << MatrixCount-MCminus+m << " "; 
	   cout << "-> ";
	    for(int k = 0; k < rho[MatrixCount-MCminus+m].Nprev; k++)
		cout << rho[MatrixCount-MCminus+m].Prev[k]->id << " ";
	    cout << endl;
	}
    }

};
void HierarchyIntegrator::prepConstructTL(){
    if (verbose > 0) cout << "Prep TL truncation matrices (with threads).\n";
    Q_tl = new Complex*[M];
    Q_tldag = new Complex*[M];
    Vsum_tl = new Complex*[M];
    freq_tl = new Complex*[Ns*Ns];
    Qsum_tl = new Complex*[Ns*Ns];
    Complex dE, f ; //dE_tl   = new Complex[Ns*Ns];
    Complex *Evec_alpha_ptr, *Evec_beta_ptr;
    Float *Vbath_ptr;
    Float nu = 2*PI*kBT/hbar;
    int ni;
    
    for (int i = 0; i < Ns*Ns; i++){
        freq_tl[i] = new Complex[MKt];
        Qsum_tl[i] = new Complex[MKt];
    }
        
    for(int i = 0; i < M; i++){ //M System-Bath matrices
        Q_tl[i] = new Complex[Ns*Ns]; // Matrix for site $i$  
        Q_tldag[i] = new Complex[Ns*Ns]; 
        Vsum_tl[i] = new Complex[Ns*Ns];
    }
    for(int i = 0; i < M; i++){ //M System-Bath matrices
        Vbath_ptr = (Vbath_re+i*M);
        ni = MDiagInd[i];
        for(int a = 0; a < Ns; a++){
            Evec_alpha_ptr = (Evec+a*Ns);
            for(int b = 0; b < Ns; b++){
                /* Prepare Vsum matrix */
                Q_tl[i][a*Ns+b] = 0.0;
                Evec_beta_ptr = Evec+b*Ns;
                if (corrBath){
                    Vsum_tl[i][a*Ns+b] = Complex(0);
                    for(int n = 0; n < M; ++n)
                        Vsum_tl[i][a*Ns+b] += Complex(Vbath_ptr[n])*conj(Evec_alpha_ptr[MDiagInd[n]])*Evec_beta_ptr[MDiagInd[n]];
                }
                else
                    Vsum_tl[i][a*Ns+b] = conj(Evec[ni+Ns*a])*(Evec[ni+Ns*b]);

                /* Prepare Freq matrix */
                dE = Complex(0,1/hbar)*(Eval[a]-Eval[b]);
                f = (gamma[i]+dE);
                freq_tl[a*Ns+b][i*Kt] = f;
                Qsum_tl[a*Ns+b][i*Kt] = C[i*Kt] / f ;
                for(int k = 1; k < Kt; k++){ //Matsuraba terms
                    f =  (Complex(k*nu)+dE);
                    freq_tl[a*Ns+b][i*Kt+k] =f;
                    Qsum_tl[a*Ns+b][i*Kt+k] = C[i*Kt+k] / f;
               }
            }
        }
    }
}


const int HierarchyIntegrator::GetMatrixCount(){
    return MatrixCount;
};

// Set up the integration with integration type Method:
// 1 = Runga-Kutta 4 fixed timesteps. Time Local Truncation and Arb. Diagonal
//     Bath Operators allowed. Shi-Adaptive Truncation NOT allowed.
// 2 = Runga-Kutta-Fehlberg 4/5 adaptive timestepping. Time-Location truncation
//     allowed. Arb. Diagonal Bath Operators NOT allowed. Shi-Adaptive Truncation allowed.
// 3 = Spectrum calculation using Runga-Kutta 4 fixed timestepping. 
// 4 = Spectrum calculation using Runga-Kutta-Fehlberg 4/5 adaptive
// timestepping
// 5 = BiCGSTAB Steady State Calculation Arb. Diagonal Bath Operators NOT
//     allowed. TL Truncation NOT allowed.
void HierarchyIntegrator::initialize(int num_threads, bool assignRhoMemory)
{


      parametersOk = checkInputParameters();

     

      kBT = kB*T;
      if (verbose > 0)
        cout << "Initializing.\n";

      if (parametersOk) {
        if (TL_trunc && (!EvecsPresent || !EvalsPresent)) {
            cout << "Eigenvectors or Eigenvalues for Hamiltonian not input. Turning off Time-Local Truncation.\n";
            TL_trunc = false;
        }
        for( int i = 0; i < Ns; i++){
            H[i+Ns*i] += Complex(0,-kappa[i]*hbar);
        }


        //check rho0
        if (numrho0inputlines == 1  && spectrum_rho0 == NULL) {
          spectrum_rho0 = new Complex[Ns];
          for (int i = 0; i < Ns; ++i)
            spectrum_rho0[i] = rho0[i*Ns];
        }
       
        if (gamma != NULL && lambda != NULL) { 
          MKt = M*Kt;
          C = new Complex[MKt];
          absC = new Float[MKt];
          if (verbose > 2)
              cout << "C:";
          Float MinAbsC;
          MinAbsC = 1e9;
          for(int i = 0; i < M; i++){
              if(gamma[i] > 0)
                  C[i*Kt] = Complex(gamma[i]*lambda[i]/hbar,0)*Complex(1/tan(gamma[i]*hbar/(2*kBT)),-1);
              else
                  C[i*Kt] = Complex(0);
              absC[i*Kt] = abs(C[i*Kt]);
              if(absC[i*Kt] < MinAbsC) MinAbsC = absC[i*Kt];
              if (verbose > 2) cout << C[i*Kt] << ",";
              for(int k =1; k < Kt; k++){
                  if(gamma[i] > 0)
                      C[i*Kt+k] = Complex(4*lambda[i]*gamma[i]*kBT/hbar/hbar*(2*PI*kBT*k/hbar)
                                          /(pow(2*PI*kBT*k/hbar,2)-pow(gamma[i],2)));
                  else
                      C[i*Kt+k] = Complex(0);
                  absC[i*Kt+k] = abs(C[i*Kt+k]);
                  if(absC[i*Kt+k] < MinAbsC) MinAbsC = absC[i*Kt+k];
                  if(verbose > 2) cout << C[i*Kt+k] << ",";
              }
          }
          if (verbose > 2) cout << endl;
        }
        if (H != NULL){
          for(int i = 0; i < Ns; i++)
              for (int j = 0; j < Ns; j++)
                  Hdag[j+i*Ns] = conj(H[i+j*Ns]);

        if(corrBath)
            cout << "Using correlated diagonal fluctuations in bath.\n";
      }


      if (Ns > 0 && Lmax > 0 && M > 0 && Kt > 0) {
        int vb = verbose;
        if (integrationMethod == PRINTHIERARCHY)
          verbose = 3;
        ConstructI();
        if (integrationMethod == PRINTHIERARCHY)
          verbose = vb;
        if ((verbose > 0 && integrationMethod != PRINTHIERARCHY)  || integrationMethod == PRINTMEMORY)
          printMemoryRequirements(false);
      }
    } //end if parameters OK
    else 
        cout << "Parameter error.\n";    

    if (!parametersOk && integrationMethod != NONE) {
       integrationMethod = BROKENINPUT;
       return;
    }

    if ((integrationMethod == RKF45 || integrationMethod == RKF45SPECTRUM) && !rho_normalized){
        if (verbose > 0) cout << "Using rho hierarchy scaling for RKF45 integration.\n";
        rho_normalized=true;
    }

    if (integrationMethod == RK4 || 
        integrationMethod == RKF45 || 
        integrationMethod == RK4SPECTRUM || 
        integrationMethod == RKF45SPECTRUM || 
        integrationMethod == BICGSTAB  ||
        integrationMethod == BICGSTABL ||
        integrationMethod == BICGSTABU) 
        WriteInputParams(num_threads);

    if (integrationMethod == BICGSTAB 
        || integrationMethod == BICGSTABL 
        || integrationMethod == BICGSTABU)
          TL_trunc = false; // Can't do steady state with TL truncation
    if ((integrationMethod == BICGSTAB 
        || integrationMethod == BICGSTABL 
        || integrationMethod == BICGSTABU) 
        && !rho_normalized) {
          if (verbose > 0) cout << "Using rho hierarchy scaling for steady state calculation.\n";
            rho_normalized=true;
    }

    if (integrationMethod == BICGSTABL && bicgstabEll < 1){
        integrationMethod = BICGSTAB;
    }

    if ( TL_trunc && (integrationMethod == RK4 || 
        integrationMethod == RKF45 || 
        integrationMethod == RK4SPECTRUM || 
        integrationMethod == RKF45SPECTRUM || 
        integrationMethod == BICGSTAB  ||
        integrationMethod == BICGSTABL ||
        integrationMethod == BICGSTABU) )
	prepConstructTL();

    if (verbose > 2)
        cout << "Creating rho[MatrixCount]\n";
    rho = new HierarchyNode[MatrixCount];

    MatrixCount_th = new int [num_threads+1];
    MatrixIndices_th = new int* [num_threads+1];
    H_th = new Complex * [num_threads];
    Hdag_th = new Complex * [num_threads];
    rho_matrices_th = new Complex * [num_threads];
    rho_mid_matrices_th = new Complex * [num_threads];
    for (int t = 0; t < num_threads; ++t){
        rho_matrices_th[t] = NULL;
        rho_mid_matrices_th[t] = NULL;
    }

    spectrum=false;
    if (integrationMethod == RK4SPECTRUM 
      || integrationMethod == RKF45SPECTRUM) 
        spectrum = true;   

    // check to see if we're doing emission spectrum:
    if (spectrum && rho0_restart) {
        L_prefactor *= Complex(-1);

    }

    if (verbose > 2)
        cout << "Creating rho_mid[MatrixCount]\n";
    if (integrationMethod == RK4 || integrationMethod == RKF45 || integrationMethod == RK4SPECTRUM || integrationMethod == RKF45SPECTRUM) 
        rho_mid = new HierarchyNode[MatrixCount];

    if (integrationMethod == BICGSTAB || integrationMethod == BICGSTABU || integrationMethod == BICGSTABL) {
        diff = new Float[num_threads+1];
        alpha = new Complex[num_threads];
        omegaTop = new Complex[num_threads];
        omegaBtm = new Complex[num_threads];
        rho1 = new Complex[num_threads];
    }

    MaxRhoVal = new Float[MatrixCount];
    for(int i = 0; i < MatrixCount; i++)
        MaxRhoVal[i] = 0;

    //Maxddrho = new Float[num_threads];

    RKF45_Diff = new Float[num_threads+1];
    min_dtval = new Float[num_threads+1];
    max_dtval = new Float[num_threads+1];

    /* Set up integrator and writer barriers */
    barrier_init(&MatupdateBarrier,num_threads);
    barrier_init(&WriterBarrier,num_threads+1);

    /* Set up mutexes */
    if (verbose > 1)
        cout << "Initializing io mutex\n";
    pthread_mutex_init(&io_lock,NULL);    


#ifdef THREADAFFINITY
    /* Check cpu affinities */
    if (verbose > 2)
        cout << "Checking cpu affinities\n";
    int ncpus = 0;
    for(int i = 0; i < ncores; ++i){
        if (cpuaffinity[i] > ncpus)
            ncpus = cpuaffinity[i];
    } 
    if (ncores < num_threads && cpuaffinity != NULL){
        cout << "WARNING: Trying to run " << num_threads << " thread on a computer with " << ncores << " cores.\n";
    }
#endif
    if (verbose > 2)
        cout << "Creating Maxddrho\n";
    Maxddrho = new Float[MatrixCount];
    for(int i = 0; i < MatrixCount; ++i)
        Maxddrho[i] = 0;        

    for(int i = 0; i < num_threads+1; ++i){
        //if (i < num_threads)
            //Maxddrho[i] = 0;
        RKF45_Diff[i] = 1E30;
        max_dtval[i] = 0.;
        min_dtval[i] = 1;
    }

    
    if (!stupid_hierarchy_partition){
        partitionHierarchy(num_threads);
    } 
    else {
        partitionHierarchySimple(num_threads);
    }
    
    countInterThreadConnections();

};

void HierarchyIntegrator::partitionHierarchySimple(int num_threads){
    cout << "Partitioning hierarchy into " << num_threads << " sets.\n";
    for(int id = 0; id < num_threads+1; ++id){
        if (id < num_threads){
            if (id < MatrixCount%num_threads)
                MatrixCount_th[id] = ceil((Float)MatrixCount/num_threads);
            else
                MatrixCount_th[id] = floor((Float)MatrixCount/num_threads);
        } 
        else
            MatrixCount_th[id] = MatrixCount;
        MatrixIndices_th[id] = new int[MatrixCount_th[id]];
    
    }
    
    if (verbose > 2){
        cout << "Matrices per thread:\n";
        for(int i = 0; i < num_threads; ++i)
            cout << i << ": " << MatrixCount_th[i] << endl;
    }
    // Cyclic partitioning
    int n = 0;
    for (int t = 0; t < num_threads; ++t){
        n = 0;
        for(int i = t; i < MatrixCount; i += num_threads) {
            MatrixIndices_th[t][n] = i;
            I[i].t = t;
            ++n;
        }
    }
    /*// Block partitioning
    int n = 0;
    for(int t = 0; t < num_threads; ++t) {
      for (int i = 0; i < MatrixCount_th[t]; ++i) {
        MatrixIndices_th[t][i] = n;
        I[n].t = t;
        ++n;
      }
    }*/
    for(int i = 0; i < MatrixCount; ++i) {
        MatrixIndices_th[num_threads][i] = i;
    }
}

void HierarchyIntegrator::partitionHierarchy(int num_threads){
    if (verbose > 1)
        cout << "Partitioning hierarchy into " << num_threads << " sets.\n";
    MultiIndex In_2(MKt);
    map< MultiIndex, int, CompareMultiIndex>::iterator fnd_result;
   
    
    for(int id = 0; id < num_threads+1; ++id){
        if (id < num_threads) {
            if (id < MatrixCount%num_threads)
                MatrixCount_th[id] = ceil((Float)MatrixCount/num_threads);
            else
                MatrixCount_th[id] = floor((Float)MatrixCount/num_threads);
        }
        else
            MatrixCount_th[id] = MatrixCount;
        MatrixIndices_th[id] = new int[MatrixCount_th[id]];
    }

    // Writer thread:
    for(int i = 0; i < MatrixCount; ++i) {
        MatrixIndices_th[num_threads][i] = i;
    }
   
    
 
    if (verbose > 2){
        cout << "Matrices per thread:\n";
        for(int i = 0; i < num_threads; ++i)
            cout << i << ": " << MatrixCount_th[i] << endl;
    }

    
    int n = 0;
    int *m_as;
    m_as = new int[num_threads];
    for(int i = 0; i < num_threads; ++i){
        m_as[i] = 0;
    }
    int nl = 0;
    
    // Do corners - then update rest by parent thread weight, if equal
    // weight then, use least assigned thread, if equal assigned and less
    // than matrix_count_th, then assigned lowest number, else assigned
    // lowest thread weight overall

    //debug
    //cout << "Assigning L=0\n"
    //debug

    MatrixIndices_th[0][0] = 0;
    m_as[0] = 1;
    I[0].t = 0;
    n = 1;


    if (Lmax > 1) { 
      // assign L=1 and corners;
      for (int i = 0; i < L[1]; ++i){
          //debug
          //cout << "rho[" << i << "]: assigned to thread " << i%num_threads << "\n";
          //debug
          MatrixIndices_th[i%num_threads][m_as[i%num_threads]] = n+i;
          m_as[i%num_threads]+=1;
          I[n+i].t = i%num_threads;

          //debug
          //cout << "Stepping down hierarchy:\n";
          //debug
          

          /*
          for(int l = 2; l < Lmax; ++l){
              if(m_as[i%num_threads] < MatrixCount_th[i%num_threads]){
                  nl = LMatrixCount[l-1]+i;
                  MatrixIndices_th[i%num_threads][m_as[i%num_threads]] = nl;
                  m_as[i%num_threads]+=1;
                  I[nl].t = i%num_threads;
                  //debug
                  //cout << "\tl=" << l << ":\n";
                  //cout << "\trho[" << nl << "]: assigned to thread " << i%num_threads << "\n";
                  //debug
              }
          }
          */
              
      }


      //debug
      //cout << "Processing rest of hierarchy.\n";
      //ndebug
      int next_t = 0;
      if (L[1] < num_threads)
          next_t = L[1];

      int max_t = -1;

      n = 1+L[1];
      int asstype = 0;
      int prev_weights[num_threads];
      for (int l = 2; l < Lmax; ++l){

          //debug
          //cout << "\tl=" << l << " next_t=" << next_t << endl;
          //debug

          for(int i = 0; i < L[l]; ++i){
              if (I[n+i].t == -1 && next_t == 0){
                  for(int t = 0; t < num_threads; ++t)
                      prev_weights[t] = 0;

                  //debug
                  //cout << "\t\tgetting weights of rho[" << n+i << "]:\n";
                  //debug

                  //determine which thread will call I[n+i]
                  for(int j = 0; j < MKt; ++j){ // run through prev connections
                      In_2 = I[n+i];
                      if (In_2[j] > 0){
                          In_2.dec(j+1);
                          fnd_result = Irev.find(In_2);
                          nl = fnd_result->second;
                          prev_weights[I[nl].t] += 1; 
                      }
                  }
                  //debug
                  //cout << "\t\t";
                  //for(int t = 0; t < num_threads; ++t)
                  //    cout << " " << prev_weights[t];
                  //cout << endl;
                  //debug

                  //find the thread that will call it the most
                  //if two or more call it equally then assign least 
                  //occupied thread.
                  max_t = 0;
                  for(int t = 0; t < num_threads; ++t){
                      if (prev_weights[t] > prev_weights[max_t]) {
                          max_t = t;                    
                          asstype = 0;
                      }
                      else if (prev_weights[t] == prev_weights[max_t] 
                                  && m_as[t] < m_as[max_t] 
                                  && prev_weights[t] > 0) {
                          max_t = t;
                          asstype = 1;
                      }
                  }
                  //debug
                  //cout << "\tmax_t=" << max_t << " at=" << asstype << endl;
                  //debug

                  //if assigned thread is already maximally occupied
                  //assign least occupied thread.
  //                if ( m_as[max_t] >= MatrixCount_th[max_t] ){
                  if ( m_as[max_t] >= ceil((Float)MatrixCount/num_threads) ) {
                      max_t = 0;
                      for(int t = 0; t < num_threads; ++t){
                          if (m_as[max_t] > m_as[t]) {
                              max_t = t;
                              asstype = 3;
                          }
                      }
                  }
                  //debug
                  //cout << "\trho[" << n+i << "] assigned to thread " << max_t << " (" << asstype << ")\n";
                  //debug

                  //now that we have identified an appropriate thread,
                  //do assignment
                  MatrixIndices_th[max_t][m_as[max_t]] = n+i;
                  m_as[max_t]+=1;
                  I[n+i].t = max_t;

                  //debug
                  //cout << "\tm_as:";
                  //for(int t = 0; t < num_threads; ++t)
                  //    cout << " " << m_as[t];
                  //cout << endl;
                  //debug
              }
              else if (I[n+i].t == -1 && next_t != 0){
                  //debug
                  //cout << "\trho[" << n+i << "] assigned to thread " << next_t << " (4)\n";
                  //debug

                  MatrixIndices_th[next_t][m_as[next_t]] = n+i;
                  m_as[next_t]+=1;
                  I[n+i].t = next_t;
                  next_t++;
                  next_t = next_t%num_threads;
                  asstype = 4;
                  //debug
                  //cout << "\tm_as:";
                  //for(int t = 0; t < num_threads; ++t)
                  //    cout << " " << m_as[t];
                  //cout << endl;
                  //debug

                  //debug
                  //cout << "\trho[" << n+i << "] assigned to thread " << next_t << " (" << asstype << ")\n";
                  //debug
              }
          }
          n += L[l];
      }
    }  
    //Now do output+testing
    //debug
    //for (int i = 0; i < MatrixCount; ++i){
    //    cout << i << "(t=" << I[i].t << "): ";
    //    I[i].Print();
    //}

    for (int t = 0; t < num_threads; ++t){
        if (m_as[t] != MatrixCount_th[t]){
            //cout << "ERROR: m_as[" << t << "] != MatrixCount_th[" << t << "] !\n";
            MatrixCount_th[t] = m_as[t];
        }
    }

    //debug
    //for (int t = 0; t < num_threads; ++t){
    //        cout << "m_as[" << t << "]=" << m_as[t] << " MatrixCount_th[" << t << "]=" << MatrixCount_th[t] << "\n";
    //}
    //debug

    if (verbose > 2) {
        for (int t = 0; t < num_threads; ++t){
            cout << "[" << t << "]:\t";
            for (int m = 0; m < m_as[t]; m++)
                cout << MatrixIndices_th[t][m] << " ";
            cout << endl;   
        }
    }
   
}


void HierarchyIntegrator::countInterThreadConnections() {
  int vertexID, prevVertexID;
  int j;
  int numInterThreadEdges=0;
  map< MultiIndex, int, CompareMultiIndex>::iterator fnd_result;
  MultiIndex In_2(MKt);
  for( vertexID = 1; vertexID < MatrixCount; ++vertexID ) {
     for ( j = 0; j < MKt; ++j) {
        In_2 = I[vertexID];
        if (In_2[j] > 0){
           In_2.dec(j+1);
           fnd_result = Irev.find(In_2);
           prevVertexID = fnd_result->second;
           if ( I[prevVertexID].t != I[vertexID].t )
             numInterThreadEdges++;
        }
     }
  }
  if ( (integrationMethod == RK4 || 
        integrationMethod == RKF45 || 
        integrationMethod == RK4SPECTRUM || 
        integrationMethod == RKF45SPECTRUM || 
        integrationMethod == BICGSTAB  ||
        integrationMethod == BICGSTABL ||
        integrationMethod == BICGSTABU  
     ) && verbose > 1 )
    cout << "Number of inter-thread connections: " << numInterThreadEdges << endl;


}
/////////////////////////////////////////////////////////////////////////////
///////                                                             /////////
/////// THREADED CODE BELOW HERE --------- THREADED CODE BELOW HERE /////////
///////                                                             /////////
/////////////////////////////////////////////////////////////////////////////

void HierarchyIntegrator::run(int id, int num_threads, int method){
#ifdef THREADAFFINITY
    cpu_set_t cpuset;
    pthread_t thread;
    int s;
    thread = pthread_self();
    if (id < num_threads && cpuaffinity != NULL){
        CPU_ZERO(&cpuset);
        CPU_SET(cpuaffinity[id],&cpuset);
        s = pthread_setaffinity_np(thread,sizeof(cpu_set_t),&cpuset);
        if ( s != 0 ){
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "]:\tFailed to set cpu affinity (errno=" << s << ").\n";
            pthread_mutex_unlock(&io_lock);
        }
        else {
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "]:\tRunning on cpu " << cpuaffinity[id] << ".\n";
            pthread_mutex_unlock(&io_lock);
        }
    }
#endif
    allocate_memory(id, num_threads, method);
    barrier_wait(&WriterBarrier);

    switch (method){
        case RK4:
            rk4_integrate(id, num_threads);
            break;
        case RKF45:
//            rkf45_integrate_multistep(id, num_threads);
            rkf45_integrate(id, num_threads);
            break;
        case RK4SPECTRUM:
            rk4_integrate(id, num_threads);
            break; 
        case RKF45SPECTRUM:
            rkf45_integrate(id, num_threads);
            break;
        case BICGSTABL:
            bicgstab_l_steadystate(id, num_threads);
            break;
        case BICGSTAB:
            bicgstab_steadystate(id, num_threads);
            break;
        case BICGSTABU:
            bicgstab_steadystate_unstable(id, num_threads);
            break;
        default:
            break;
    }
    pthread_exit(0);
};


int HierarchyIntegrator::allocate_memory(int id, int num_threads, int Method){
    //bool spectrum=false;
    //if (Method == RK4SPECTRUM || Method == RKF45SPECTRUM) spectrum = true;   
    //NOTE THIS IS FOR THE PREFACTORS OF THE DYNAMICS
    //NOTE 2: FOR PREV LINKING TERMS IN DYNAMICS USE +[]rho.row +[]rho.col
    
    //debug
    //if (id == 0 ) { cout << "allocating memory\n"; }
    //debug
    int nm = 0;
    Float nu = 2*PI*kBT/hbar;
    barrier_wait(&WriterBarrier);
    if (id < num_threads){
        H_th[id] = new Complex[Ns2];
        Hdag_th[id] = new Complex[Ns2];

        for(int i = 0; i < Ns2; ++i){
            H_th[id][i] = H[i];
            Hdag_th[id][i] = Hdag[i];
        }
        if (!spectrum) {
            rho_matrices_th[id] = new Complex[MatrixCount_th[id]*Ns2];
            for (int i = 0; i < MatrixCount_th[id]; ++i){
                nm = MatrixIndices_th[id][i];
                rho[nm].rho = rho_matrices_th[id]+i*Ns2;
                for(int n = 0; n < Ns2; ++n)
                    rho[nm].rho[n] = 0;
            }
           
        }
        else {
            rho_matrices_th[id] = new Complex[MatrixCount_th[id]*Ns];
            for (int i = 0; i < MatrixCount_th[id]; ++i){
                nm =  MatrixIndices_th[id][i];
                rho[nm].rho = rho_matrices_th[id]+i*Ns;
                for(int n = 0; n < Ns; ++n)
                    rho[nm].rho[n] = 0;
            }
        }

        //Primary hierarchy
        ConstructRho(rho,id,num_threads,rho_matrices_th[id]);
    }
    barrier_wait(&WriterBarrier);
    int vb = verbose;
    if (id == 0)    
        verbose = -1;
    barrier_wait(&WriterBarrier);
    if (id < num_threads){    
        if (Method == RK4 || Method == RKF45 || Method == RK4SPECTRUM || Method == RKF45SPECTRUM){
            if (!spectrum) {
                rho_mid_matrices_th[id] = new Complex[MatrixCount_th[id]*Ns2];
                for (int i = 0; i < MatrixCount_th[id]; ++i){
                    nm =  MatrixIndices_th[id][i];
                    rho_mid[nm].rho = rho_mid_matrices_th[id]+i*Ns2;
                    for(int n = 0; n < Ns2; ++n)
                        rho_mid[nm].rho[n] = 0;
                }   
            }
            else {
                rho_mid_matrices_th[id] = new Complex[MatrixCount_th[id]*Ns];
                for (int i = 0; i < MatrixCount_th[id]; ++i){
                    nm =  MatrixIndices_th[id][i];
                    rho_mid[nm].rho = rho_mid_matrices_th[id]+i*Ns;
                    for(int n = 0; n < Ns; ++n)
                        rho_mid[nm].rho[n] = 0;
                }   
            }
            //Secondary hierarchy for storing integration steps
            ConstructRho(rho_mid,id,num_threads,rho_mid_matrices_th[id]);

        }
    }
    barrier_wait(&WriterBarrier);
    if (id == 0)
      verbose = vb;
    barrier_wait(&WriterBarrier);
    if (id < num_threads){    
        int i=0;
        for(int ii = 0; ii < MatrixCount_th[id]; ++ii){
            i = MatrixIndices_th[id][ii];
            rho[i].NZI_prefactor = Complex(0,0); // 
            rho[i].Same_prefactor= new Complex[M];
            rho[i].Next_prefactor= new Complex[rho[i].Nnext];
            for(int j = 0; j < M; j++){
                rho[i].Same_prefactor[j] = Complex(0,0);
            }
            rho[i].Prev_prefactor_row = new Complex[rho[i].Nprev];
            rho[i].Prev_prefactor_col = new Complex[rho[i].Nprev];
            for(int j = 0; j < rho[i].Nprev; j++){
                rho[i].Prev_prefactor_row[j] = Complex(0,0);
                rho[i].Prev_prefactor_col[j] = Complex(0,0);
            }
            
           //NEXT n+ 
            for(int j = 0; j < rho[i].Nnext; j++){
                if(rho_normalized){
                    nm = rho[i].NextInd[j];
                    rho[i].Next_prefactor[j] = Complex(0,sqrt((rho[i].I[nm]+1)*absC[nm]));
                }
                else {
                    rho[i].Next_prefactor[j] = Complex(0,1);
                 }
                //check for emission spectrum:
                if (spectrum && rho0_restart) rho[i].Next_prefactor[j] *= Complex(-1);
             }   
            //SAME - AI Temperature Fi
            //SAME n
            for(int j = 0; j < M; j++){
                if (gamma[j] > 0)
                    rho[i].Same_prefactor[j] = Complex(2*lambda[j]*kBT/hbar/hbar/gamma[j]-lambda[j]/hbar/tan(hbar*gamma[j]/(2*kBT)),0);
                else
                    rho[i].Same_prefactor[j] = Complex(0);
                for(int k = 1; k < Kt; k++){
                    if (gamma[j] > 0)
                        rho[i].Same_prefactor[j] -= Complex(4*lambda[j]*kBT/hbar/hbar*gamma[j]/(pow(nu*k,2)-pow(gamma[j],2)),0);
                    else
                        rho[i].Same_prefactor[j] += Complex(0);
                }
                //same_prefactor stays the same for emission and absorption
                //if (spectrum && rho0_restart) rho[i].Same_prefactor[j] *= Complex(-1);
             }	 

            //PREV AND NZI
             for(int j = 0; j < rho[i].Nprev; j++){   
                nm = rho[i].PrevInd[j];
                if(rho_normalized){
                    rho[i].Prev_prefactor_row[j] =      C[nm] *Complex(0,sqrt(rho[i].I[nm]/absC[nm]));
                    rho[i].Prev_prefactor_col[j] = conj(C[nm])*Complex(0,sqrt(rho[i].I[nm]/absC[nm]));
                }
                if (nm%Kt == 0){
                    rho[i].NZI_prefactor += Complex(rho[i].I[nm]*gamma[nm/Kt],0);
                    if(!rho_normalized){
                        if (gamma[j] > 0 ) {           // i(I[nm]*lambda*gamma/hbar(cot[hbar*gamma/2kT] -i) )
                                                       // i*I[nm]*C[mn]
                            rho[i].Prev_prefactor_row[j] = Complex( rho[i].I[nm] * lambda[nm/Kt] * gamma[nm/Kt]/hbar,
                                                            rho[i].I[nm] * lambda[nm/Kt]*gamma[nm/Kt]/hbar
                                                            /tan(hbar*gamma[nm/Kt]/(2*kBT)));
                                                       // i(I[nm]*lambda*gamma/hbar(cot[hbar*gamma/2kT]+i) )
                                                       // i*I[mn]*conj(C[mn])
                            rho[i].Prev_prefactor_col[j] = Complex(-rho[i].I[nm] * lambda[nm/Kt] * gamma[nm/Kt]/hbar,
                                                            rho[i].I[nm] * lambda[nm/Kt] * gamma[nm/Kt]/hbar
                                                            /tan(hbar*gamma[nm/Kt]/(2*kBT)));
                        }
                        else {
                            rho[i].Prev_prefactor_row[j] = Complex(0);
                            rho[i].Prev_prefactor_col[j] = Complex(0);
                        }
                    }
                }
                else{
                    rho[i].NZI_prefactor += Complex(rho[i].I[nm]*nu*(nm%Kt),0);
                    if(!rho_normalized){
                        if (gamma[j] > 0) {
                            rho[i].Prev_prefactor_row[j] = Complex(0,rho[i].I[nm]*(nm%Kt)*nu*4*lambda[(nm-nm%Kt)/Kt]
                                                                *gamma[(nm-nm%Kt)/Kt]/pow(hbar,2)*kBT
                                                                /(pow((nm%Kt)*nu,2)-pow(gamma[(nm-nm%Kt)/Kt],2)));
                            rho[i].Prev_prefactor_col[j] = Complex(0,rho[i].I[nm]*(nm%Kt)*nu*4*lambda[(nm-nm%Kt)/Kt]
                                                                *gamma[(nm-nm%Kt)/Kt]/pow(hbar,2)*kBT
                                                                /(pow((nm%Kt)*nu,2)-pow(gamma[(nm-nm%Kt)/Kt],2)));
                        }
                        else {
                            rho[i].Prev_prefactor_row[j] = Complex(0);
                            rho[i].Prev_prefactor_col[j] = Complex(0);
                        }
                    }
                }
                //check for emission spectrum:
                //if (spectrum && rho0_restart) 
                //    rho[i].Prev_prefactor_row[j] = Complex(-1)*conj(rho[i].Prev_prefactor_row[j]);
                if (spectrum && rho0_restart) 
                    rho[i].Prev_prefactor_row[j] = Complex(-1)*rho[i].Prev_prefactor_col[j];
            }
            rho[i].H = H_th[id];
            rho[i].Hdag = Hdag_th[id];
            if (Method == RK4 || Method == RKF45 || Method == RK4SPECTRUM || Method == RKF45SPECTRUM) {
                rho_mid[i].NZI_prefactor = rho[i].NZI_prefactor;
                rho_mid[i].Next_prefactor = rho[i].Next_prefactor;//new Complex[rho[i].Nnext];
                rho_mid[i].Same_prefactor = rho[i].Same_prefactor; //new Complex[rho[i].Ns];
                rho_mid[i].Prev_prefactor_row = rho[i].Prev_prefactor_row; //new Complex[rho[i].Nprev];
                rho_mid[i].Prev_prefactor_col = rho[i].Prev_prefactor_col;// new Complex[rho[i].Nprev];
                rho_mid[i].H = rho[i].H;
                rho_mid[i].Hdag = rho[i].Hdag;
            }
        }
    }
    barrier_wait(&WriterBarrier);
    Float rhoScale;
    if((Method == RK4 || Method == RKF45 || Method == RK4SPECTRUM || Method == RKF45SPECTRUM)) {
        if (id < num_threads) {
            // set up rho(t=0)
            if(!spectrum){
                if(rho0_restart){
                    if (id == 0) { 
                        pthread_mutex_lock(&io_lock);
                        ifstream rhoin(restartIn.c_str());
                        rhoin >> t_start;
                        for(int m = 0; m < MatrixCount; m++)
                            for(int j = 0; j < Ns; j++)
                                for( int k = 0; k < Ns; k++)
                                    rhoin >> rho[m].rho[j+Ns*k];
                        rhoin.close();
                        pthread_mutex_unlock(&io_lock);
                    }
                }
                else{
                    t_start = 0.0;
                    for( int m = id; m < MatrixCount; m+=num_threads){  
                        rhoScale = 1;
                        if(rho_normalized)
                            for (int i = 0; i < M*Kt; i++)
                                rhoScale *= sqrt(pow(absC[i],rho[m].I[i])*iter_factorial(rho[m].I[i]));
                        for( int i = 0; i < Ns; i++){
                            for( int j = 0; j < Ns; j++){
                                rho[m].rho[j+Ns*i] = rho0[j+Ns*i]/rhoScale;
                            }
                        }
                    }
                }
            }
            else{
                if(rho0_restart){
                    //This is for calculation of emission spectra - 
                    //basically calculate the dynamics starting from the
                    //the equilibrium density matrices.
                    if (id == 0) { 
                        pthread_mutex_lock(&io_lock);
                        ifstream rhoin(restartIn.c_str());
                        rhoin >> t_start;
                        Complex val;
                        for(int m = 0; m < MatrixCount; m++){
                            for(int j = 0; j < Ns; j++)
                                rho[m].rho[j] = 0;
                            for(int j = 0; j < Ns; j++){
                                for( int k = 0; k < Ns; k++){
                                    rhoin >> val; 
                                    // all dipoles aligned
                                    if (spectrum_rho0 == NULL) 
                                        rho[m].rho[k] += val;
                                    // calculate for specific component of
                                    // dipole moments 
                                    else
                                        rho[m].rho[k] += val*spectrum_rho0[k];
                                }
                            }
                        }
                        rhoin.close();
                        pthread_mutex_unlock(&io_lock);
                    }
                }
                else {
                    //if (id == 0) cout << "Initializing rho0 for spectrum calc" << endl;
                    for( int m = id; m < MatrixCount; m+=num_threads){  
                        rhoScale = 1;
                        if(rho_normalized){
                            for (int i = 0; i < M*Kt; i++)
                                rhoScale *= sqrt(iter_factorial(rho[m].I[i])*pow(absC[i],rho[m].I[i]));
                        }
                        if (spectrum_rho0 == NULL) {
                            // all dipoles aligned
                            for( int i = 0; i < Ns; i++){
                                rho[m].rho[i] = 1/rhoScale;
                            }
                        }
                        else {
                            // calculate for specific component of
                            // dipole moments 
                            for( int i = 0; i < Ns; i++){
                                rho[m].rho[i] = spectrum_rho0[i]/rhoScale;
                                //if (m == 0) cout << i << " " << spectrum_rho0[i] << endl;
                            }
                        }
                    }
                }
            }
        }
        barrier_wait(&WriterBarrier);
        if (id == num_threads) {
            pthread_mutex_lock(&io_lock);
            outf << 0.000 << " ";
            int j = 0;	
            if(!spectrum)
                for(j = 0; j < Ns; j++)
                    for( int k = 0; k < Ns; k++)
                        outf << rho[0].rho[k+Ns*j] << " ";
            else
                for( int k = 0; k < Ns; k++)
                    outf << rho[0].rho[k] << " ";
            outf << endl;
            pthread_mutex_unlock(&io_lock);
        }
    }
    if (Method == BICGSTAB || Method == BICGSTABU || Method == BICGSTABL) {
        if (id < num_threads) {
            prepForSteadyState(id, num_threads, Method);
        }
    }
    barrier_wait(&WriterBarrier);
    return 0;

};

///////////////////////////////////////////////////////////
//////////  Construct Pointer-Linked Hierarchy  ///////////
///////////////////////////////////////////////////////////
void HierarchyIntegrator::ConstructRho(HierarchyNode* sigma, int id, int num_threads, Complex * rho_matrices){
    if (verbose > 1){
        pthread_mutex_lock(&io_lock);
        cout << "[" << id << "]:\tConstructing density matrices.\n";
        pthread_mutex_unlock(&io_lock);
    }
    MultiIndex In_2(MKt);
    int PrevCount = 0;

    int n_fnd;
    int prev_assigned = 0;
    map< MultiIndex, int, CompareMultiIndex>::iterator fnd_result;
    HierarchyNode **MatrPtr;
    int *Ind;
    
    Ind = new int[MKt];
    MatrPtr = new HierarchyNode*[MKt];
    int i = 0;
    for(int ii = 0; ii < MatrixCount_th[id]; ++ii){
        i = MatrixIndices_th[id][ii];
        //debug
        //cout << "i: " << i << " ii: " << ii << " id: " << id << " MatrixCount_th[id] " <<  MatrixCount_th[id] << endl;
        //cout << "\t";  I[i].Print(); 
        //debug
        if (rho_matrices != NULL) {
            if(spectrum) sigma[i].createVec(Ns,M,Kt,I[i].I,&(rho_matrices[i*Ns]));
            else sigma[i].create(Ns,M,Kt,I[i].I,&(rho_matrices[i*Ns2]));
        }
        else {
            if(spectrum) sigma[i].createVec(Ns,M,Kt,I[i].I,NULL);
            else sigma[i].create(Ns,M,Kt,I[i].I,NULL);
        }
	sigma[i].id = i;
	sigma[i].Ns2 = Ns*Ns;
	sigma[i].TL_truncation = false;
        if (multiBath) {
            for (int m = 0; m < M; ++m) 
                sigma[i].SameInd[m] = MDiagInd[m];
        }
        else {
            for (int m = 0; m < M; ++m)
                sigma[i].SameInd[m] = m;
        }
    }
    barrier_wait(&MatupdateBarrier);
    for(int ii = 0; ii < MatrixCount_th[id]; ++ii){
        i = MatrixIndices_th[id][ii];        
	PrevCount = 0;
	if (verbose > 2) {
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "] ";
            sigma[i].printI();
            cout << ":\n";
        }
	for(int j=0; j < MKt; j++){
	    if (sigma[i].I[j] > 0)
		PrevCount++;
	}


	//PREV:
	if (PrevCount > 0){
	    if (verbose > 2) {
                cout << "[" << id << "]\tLinking Prev Matrices.\n";
            }
	    prev_assigned = 0;
    	    for(int j=0; j < MKt; j++){
	    	In_2 = I[i];    		
		if (In_2[j] > 0){
		    In_2.dec(j+1);
		    fnd_result = Irev.find(In_2);
		    if(fnd_result != Irev.end()){
			n_fnd = fnd_result->second;
			MatrPtr[prev_assigned] = &sigma[n_fnd];
			Ind[prev_assigned] = j;
			if (verbose > 2){
			    cout << "\t" << j << ": ";
			    sigma[n_fnd].printI();
			    cout << " - ";
			    MatrPtr[prev_assigned]->printI();
			    cout << "\n";
			}
			prev_assigned++;
		    }
		}
    		    
    	    }	
	    sigma[i].Nprev = prev_assigned;
	    sigma[i].Prev = new HierarchyNode*[prev_assigned];
	    sigma[i].PrevInd = new int[prev_assigned];
	    sigma[i].PrevIndMatrix = new int[prev_assigned];
	    for(int j = 0; j < prev_assigned; j++){
		sigma[i].Prev[j] = MatrPtr[j];
		sigma[i].PrevInd[j] = Ind[j];
		sigma[i].PrevIndMatrix[j] = MDiagInd[(Ind[j]-Ind[j]%Kt)/Kt];
	    }
	}	

	// NEXT
	if (sigma[i].L < Lmax-1)
	{
	    prev_assigned = 0;
	    if (verbose > 2) {
                cout << "[" << id << "]\tLinking Next Matrices.\n";
            }
	    for(int j=0; j < MKt; j++){
		In_2 = I[i];
		In_2.inc(j+1);
		Ind[prev_assigned] = j;
		fnd_result = Irev.find(In_2);
		if(fnd_result != Irev.end()){
		    n_fnd = fnd_result->second;
		    MatrPtr[j] = &sigma[n_fnd];
		    if (verbose > 2){
                        cout << "\t" << j << ": ";
                        sigma[n_fnd].printI();
		        cout << " - ";
		        MatrPtr[prev_assigned]->printI();
		        cout << "\n";
		    }
    		    prev_assigned++;
		}
                else{
                    sigma[i].TL_truncation = TL_trunc;
                }
	    }
            if (prev_assigned > 0) {
                sigma[i].Nnext = prev_assigned;
                sigma[i].Next = new HierarchyNode*[prev_assigned];
                sigma[i].NextInd = new int[prev_assigned];
                sigma[i].NextIndMatrix = new int[prev_assigned];
                for(int j = 0; j < prev_assigned; j++){
                    sigma[i].Next[j] = MatrPtr[j];
                    sigma[i].NextInd[j] = Ind[j];
                    sigma[i].NextIndMatrix[j] = MDiagInd[(Ind[j]-Ind[j]%Kt)/Kt];
                }
                sigma[i].TL_truncation = false;
            }
            else{
                sigma[i].TL_truncation = TL_trunc;
            }

	}
        if (verbose > 2)
               pthread_mutex_unlock(&io_lock);

	else if(sigma[i].L == Lmax-1){
	    sigma[i].TL_truncation = TL_trunc;
	}
    }
    barrier_wait(&MatupdateBarrier);
 
};

/////////////////////////////////////////////////////////
///////  Prepare Steady Matrices and Parameters  ////////
/////////////////////////////////////////////////////////
void HierarchyIntegrator::prepForSteadyState(int id, int num_threads, int Method){
    if (verbose > -1){
        pthread_mutex_lock(&io_lock);
        cout << "[" << id << "] Initialized Steady State Calculation.\n";
        pthread_mutex_unlock(&io_lock);
    }
    for (int i = id; i < MatrixCount; i+=num_threads){
        rho[i].freeRho();
        if (Method == BICGSTAB || Method == BICGSTABU)
            rho[i].bicgstabInit();
        else if (Method == BICGSTABL)
            rho[i].bicgstablInit(bicgstabEll);
    }
    barrier_wait(&MatupdateBarrier);
    // Create Liouville Space operator
    // for system density matrix Hamiltonian
    if (id == 0) {
        Complex * H_l;
        H_l = new Complex[Ns*Ns*Ns*Ns];
        H_l = toLiouville(H,Ns);
        int m, n;
        m = n = 0;
        for (int i = 0; i < Ns*Ns*Ns*Ns; ++i){
                H_l[i] *= L_prefactor;
        }
        for (int i = 0; i < Ns; ++i){
            for (int ii = 0; ii < Ns; ++ii){
                for (int j = 0; j < M; ++j){
                    n = MDiagInd[j];
                    if ((i == n && ii != n) || (i != n && ii == n)){ 
                        H_l[(i*Ns+ii)*Ns*Ns+i*Ns+ii] -= rho[0].Same_prefactor[j];
                    }
                }
            }
        }
        rho[0].SameLiouville = new Complex[Ns*Ns*Ns*Ns];
        for(int i = 0; i < Ns*Ns*Ns*Ns; ++i)
            rho[0].SameLiouville[i] = H_l[i];
        //for(int i = 1; i < Ns2; ++i)
        //    cout << i << " " << -rho[0].SameLiouville[i] << endl;

        /* Reduce Storage + Access Times - Need to modify minimize_step(...);
        // then
        NumHentries = Ns*Ns*(2*Ns-1);
        rho[0].SameLiouville = new Complex[NumHentries];
        m = 0; 
        n = 0;
        for (int i = 0; i < Ns; ++i){
            for(int ii=0; ii < Ns; ++ii) {
                for (int j = 0; j < i; ++j) {
                    rho[0].SameLiouville[m+j] = H_l[n+j*Ns+ii];                    
                }
                for (int j = 0; j < Ns; ++j){
                    rho[0].SameLiouville[m+i+j] = H_l[n+i*Ns+j]; 
                }
                for (int j = 0; j < Ns-1-i; ++j){
                    rho[0].SameLiouville[m+i+Ns+j] = H_l[n+(i+1)*Ns+ii+j*Ns]; 
                }
                m += 2*Ns-1;
                n += Ns*Ns;
            }
            
        }
        */
    }
    // Create Liouville space operator
    // for operations pointing to system
    // density matrix
    int NumPrevEntries = 2*Ns-1;
    int **PrevIndices;
    PrevIndices = new int*[M];
    int n, m;
    for (int nj = 0; nj < M; ++nj) {
        n = 0;
        m = MDiagInd[nj]; 
        PrevIndices[nj] = new int[NumPrevEntries];
        for (int i = 0; i < m; ++i) {
            PrevIndices[nj][n] = i*Ns+m;
            ++n;
        }
        for (int i = 0 ; i < m; ++i){
            PrevIndices[nj][n] = m*Ns+i;
            ++n;
        }
        PrevIndices[nj][n] = m*Ns+m;
        ++n;
        for (int i = m+1 ; i < Ns; ++i){
            PrevIndices[nj][n] = m*Ns+i;
            ++n;
        }
        for (int i = m+1; i < Ns; ++i) {
            PrevIndices[nj][n] = i*Ns+m;
            ++n;
        }
    }
    for (int i = id; i <= L[1]; i+=num_threads) {
        rho[i].createPrevLiouvilleOperator(NumPrevEntries);
        rho[i].PrevIndices = PrevIndices;
        
    }

    //
    // Guess a steady-state matrix
    //
    if (id == 0) {
        if (verbose > 1){
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "] Guessing steady-state rho.\n";
            pthread_mutex_unlock(&io_lock);
        }
        if (EvalsPresent && EvecsPresent) {
            if (verbose > 1){
                pthread_mutex_lock(&io_lock);
                cout << "[" << id << "] Eigenvalues and vectors present.\n";
                pthread_mutex_unlock(&io_lock);
            }
            Complex EvalPops[Ns];
            Complex Z = 0;
            for (int i = 0; i < Ns; ++i) {
                EvalPops[i] = exp(-Eval[i]/kBT);
                Z += EvalPops[i];
            }
            for (int i = 0; i < Ns; ++i) {
                EvalPops[i] /= Z;
            }
            Complex * A;
            A = new Complex[Ns*Ns];
            for(int i = 0; i < Ns*Ns; ++i)
                A[i] = 0;
            for (int i = 0; i < Ns; ++i) {
                //phi_her(CblasRowMajor,CblasUpper,Ns,real(EvalPops[i]),&(Evec[i*Ns]),1,A,Ns);
                phi_her(Ns,real(EvalPops[i]),&(Evec[i*Ns]),A);
            }

    //        for (int m = 0; m < MatrixCount; m += 1)
                int m = 0;
                for (int i = 0; i < Ns; ++i){
                    rho[m].SteadyStateRho[i*Ns+i] = A[i*Ns+i];
                    for (int j = i+1; j < Ns; ++j) {
                        rho[m].SteadyStateRho[j*Ns+i] = conj(A[i*Ns+j]);
                        rho[m].SteadyStateRho[i*Ns+j] = A[i*Ns+j];
                    }
                }
        }   
        else {
            //for(int m = 0; m < MatrixCount; m += 1)
                int m = 0;
                for (int i = 0; i < Ns; ++i)
                    rho[m].SteadyStateRho[i*Ns+i] = Complex(1./Ns);
        } 
            
        // Write out initial guess
        if (verbose > 2) {
            pthread_mutex_lock(&io_lock);            
            cout << "Initial Guess of Steady State:\n";
            printMatrix(rho[0].SteadyStateRho,Ns);
            pthread_mutex_unlock(&io_lock);

        }
        pthread_mutex_lock(&io_lock);        
        outf << "#Initial rho:\n#";
        for (int i = 0; i < Ns*Ns; ++i)            
            outf << rho[0].SteadyStateRho[i] << " ";
        outf << endl;
        pthread_mutex_unlock(&io_lock);
        
    }
    barrier_wait(&MatupdateBarrier);
};


////////////////////////////////////////////////////////////
//    Prepare system for Emission Spectrum Calculation    //
////////////////////////////////////////////////////////////
//void HierarchyIntegrator::prepE

////////////////////////////////////////////////////////////
///////////////  Time Local Truncation    //////////////////
////////////////////////////////////////////////////////////
/* TL model. At each timestep, recalculate (with threads) new
   Q_tl and Q_tldag
*/

void HierarchyIntegrator::constructTL_BLAS_threaded(Float t, int i){
    Complex * Q_tl_ptr, * Q_tldag_ptr, * Evec_alpha_ptr, * Evec_beta_ptr;
    Complex Q_sum;
    Complex tnow;   
    Complex f; 
    Float *Vbath_ptr;
//    time_t t1,t2;
    tnow = Complex(t);
    // i = Current system-bath coupling term:

    // NOT CURRENTLY USED!!! NEED TO REDO DiagBath
    Vbath_ptr = (Vbath_re+i*Ns);


    Q_tl_ptr = Q_tl[i];
    Q_tldag_ptr = Q_tldag[i];
    for (int a =0; a < Ns*Ns; a++)
        Q_tl_ptr[a] = zero;

    for(int a = 0; a < Ns; a++){
        Evec_alpha_ptr = (Evec+a*Ns);
        for(int b = 0; b < Ns; b++){
            Evec_beta_ptr = Evec+b*Ns;

            f = freq_tl[a*Ns+b][i*Kt];
            Q_sum = Qsum_tl[a*Ns+b][i*Kt] *(one-exp(-f*tnow));
            for(int k = 1; k < Kt; k++){ //Matsuraba terms
                f = freq_tl[a*Ns+b][i*Kt+k];
                Q_sum += Qsum_tl[a*Ns+b][i*Kt+k] *(one-exp(-f*tnow));
            }
            
            Q_sum = Q_sum*Vsum_tl[i][a*Ns+b]; 
            // CHECK IF Evec_beta_ptr does not have to be conjugated!
            //phi_gerc(CblasColMajor,Ns,Ns,(&Q_sum),Evec_alpha_ptr,1,Evec_beta_ptr,1,Q_tl_ptr,Ns);
            phi_gerc(Ns,&Q_sum,Evec_alpha_ptr,Evec_beta_ptr,Q_tl_ptr);
        }
    }

    for(int a = 0; a < Ns; a++){
        Q_tl_ptr = (Q_tl[i]+a);
        for(int b = 0; b < Ns; b++){
            *Q_tldag_ptr = conj(*Q_tl_ptr);
             Q_tldag_ptr++;
             Q_tl_ptr += Ns;
        }
    }
}
//////////////////////////////////////////
        // INTEGRATION STEP //
//////////////////////////////////////////

void HierarchyIntegrator::integration_step(const HierarchyNode &sigma, Complex * drho,const int &time_step, Complex * Qrho_rhoQ, Complex * Vnext, Complex * Vprev, Complex * Vsame){
    int nj;
    int m;
    int j;

    Complex * drho_col_elem, *drho_row_elem, * rho_row_elem, *rho_col_elem;

    phi_symm(CblasLeft, Ns,&L_prefactor,   sigma.H,   sigma.rho,&zero,drho);
    phi_symm(CblasRight,Ns,&negL_prefactor,sigma.Hdag,sigma.rho,&one, drho);

    //phi_symm(cblascolmajor,cblasleft,cblasupper,ns,ns,&l_prefactor,sigma.h,ns,sigma.rho,ns,&zero,drho,ns);
    //phi_symm(cblascolmajor,cblasright,cblasupper,ns,ns,&negl_prefactor,sigma.hdag,ns,sigma.rho,ns,&one,drho,ns);

    rho_col_elem = sigma.rho;
    drho_col_elem = drho;
    for(nj = 0; nj < Ns2;nj++, drho_col_elem++, rho_col_elem++){
        *drho_col_elem -= sigma.NZI_prefactor*(*rho_col_elem); 
    }

    if (!multiBath) {
        for(m = 0; m < Ns; m++){
            drho_row_elem = &(drho[m*Ns]);
            rho_row_elem = &(sigma.rho[m*Ns]);
            for(nj = 0; nj < Ns; nj++){
                drho_row_elem[nj] -= (sigma.Same_prefactor[nj]+sigma.Same_prefactor[m])*rho_row_elem[nj];
            }
            nj = (Ns+1)*m;
            drho[nj] += two*sigma.Same_prefactor[m]*sigma.rho[nj];
        }
    } 
    else {
      for(j = 0; j < M; j++){
        nj = MDiagInd[j];
        drho_row_elem = &(drho[nj]);
        drho_col_elem = &(drho[Ns*nj]);
        rho_row_elem = &(sigma.rho[nj]);
        rho_col_elem = &(sigma.rho[Ns*nj]);
        for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
            if (m != nj ) {
              *drho_row_elem -= sigma.Same_prefactor[j]*(*rho_row_elem);
              *drho_col_elem -= sigma.Same_prefactor[j]*(*rho_col_elem);	
            }
        }
      }
    }
    


    if(sigma.TL_truncation){
        //- sum_j=1^M [ |j><j| (Q_j rho - rho Q_j^H) - (Q_j rho - rho Q_j^H) |j><j| ]
        if (multiBath) {
            for(j = 0; j < M; j++){
                nj = MDiagInd[j];
                phi_hemm(CblasRight,Ns,&one,   sigma.rho,Q_tl[j],   &zero,Qrho_rhoQ);
                phi_hemm(CblasLeft, Ns,&negone,sigma.rho,Q_tldag[j],&one, Qrho_rhoQ);
                //phi_hemm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&one,sigma.rho,Ns,Q_tl[j],Ns,&zero,Qrho_rhoQ,Ns);
                //phi_hemm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&negone,sigma.rho,Ns,Q_tldag[j],Ns,&one,Qrho_rhoQ,Ns);
                drho_row_elem = &(drho[nj]);
                drho_col_elem = &(drho[Ns*nj]);
                rho_row_elem = &(Qrho_rhoQ[nj]);
                rho_col_elem = &(Qrho_rhoQ[Ns*nj]);
                for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                    *drho_row_elem -= (*rho_row_elem);
                    *drho_col_elem += (*rho_col_elem);	
                }
            }
        } 
        else {
            for(nj = 0; nj < Ns; nj++){
                phi_hemm(CblasRight,Ns,&one,   sigma.rho,Q_tl[nj],   &zero,Qrho_rhoQ);
                phi_hemm(CblasLeft, Ns,&negone,sigma.rho,Q_tldag[nj],&one, Qrho_rhoQ);
                //phi_hemm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&one,sigma.rho,Ns,Q_tl[nj],Ns,&zero,Qrho_rhoQ,Ns);
                //phi_hemm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&negone,sigma.rho,Ns,Q_tldag[nj],Ns,&one,Qrho_rhoQ,Ns);
                drho_row_elem = &(drho[nj]);
                drho_col_elem = &(drho[Ns*nj]);
                rho_row_elem = &(Qrho_rhoQ[nj]);
                rho_col_elem = &(Qrho_rhoQ[Ns*nj]);
                for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                    *drho_row_elem -= (*rho_row_elem);
                    *drho_col_elem += (*rho_col_elem);	
                }
            }
        }
    }
    else{
	for(j = 0; j < sigma.Nnext; j++){ // this is sum over Ns & Kt
            nj = sigma.NextIndMatrix[j];
            drho_col_elem = &(drho[Ns*nj]);
            rho_col_elem = &(sigma.Next[j]->rho[Ns*nj]);
            for(m=0; m < Ns; m++){
                drho_col_elem[m] += (sigma.Next_prefactor[j])*rho_col_elem[m];	
            }
	}
	for(m = 0; m < Ns; m++){ //, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
	    drho_row_elem = &(drho[Ns*m]);
	    for(j = 0; j < sigma.Nnext; j++){
                nj = sigma.NextIndMatrix[j];
                rho_row_elem = &(sigma.Next[j]->rho[Ns*m]);
                drho_row_elem[nj] -= sigma.Next_prefactor[j]*rho_row_elem[nj];
            }
	}
    }

    for(j = 0; j < sigma.Nprev; j++){
        m = sigma.PrevIndMatrix[j];
        drho_col_elem = &drho[Ns*m];
        rho_col_elem = &(sigma.Prev[j]->rho[Ns*m]);
        for(m = 0; m < Ns; m++) //, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
            drho_col_elem[m] += sigma.Prev_prefactor_col[j]*rho_col_elem[m];	
    }
    
    for(m = 0; m < Ns; m++){// rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
	drho_row_elem = &drho[m*Ns];
	for(j = 0; j < sigma.Nprev; j++){
            nj = sigma.PrevIndMatrix[j];  //nj runs from 0 to Ns (e.g. 0,0,0,1,1,1,2,2,2,...,Ns,Ns,Ns  for Kt = 2)
            rho_row_elem = &(sigma.Prev[j]->rho[m*Ns]);
            drho_row_elem[nj] -= sigma.Prev_prefactor_row[j]*rho_row_elem[nj];
	}
    }
};

//////////////////////////////////////////
// ADAPTIVE TRUNCATION INTEGRATION STEP //
//////////////////////////////////////////

void HierarchyIntegrator::integration_step_adaptTrunc(const HierarchyNode &sigma, Complex * drho,const int &time_step, Complex * Qrho_rhoQ, Complex * Vnext, Complex * Vprev, Complex * Vsame){
    int nj;
    int m;
    int j;
    Complex * drho_col_elem, *drho_row_elem, * rho_row_elem, *rho_col_elem;
    if ( sigma.active ){ 
        phi_symm(CblasLeft,Ns,&L_prefactor,sigma.H,sigma.rho,&zero,drho);
        phi_symm(CblasRight,Ns,&negL_prefactor,sigma.Hdag,sigma.rho,&one,drho);

        //phi_symm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&L_prefactor,sigma.H,Ns,sigma.rho,Ns,&zero,drho,Ns);
        //phi_symm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&negL_prefactor,sigma.Hdag,Ns,sigma.rho,Ns,&one,drho,Ns);


        rho_col_elem = sigma.rho;
        drho_col_elem = drho;
        for(nj = 0; nj < Ns2;nj++, drho_col_elem++, rho_col_elem++){
            *drho_col_elem -= sigma.NZI_prefactor*(*rho_col_elem); 
        }

        if (!multiBath) {
            for(m = 0; m < Ns; m++){
                drho_row_elem = &(drho[m*Ns]);
                rho_row_elem = &(sigma.rho[m*Ns]);
                for(nj = 0; nj < Ns; nj++){
                    drho_row_elem[nj] -= (sigma.Same_prefactor[nj]+sigma.Same_prefactor[m])*rho_row_elem[nj];
                }
                nj = (Ns+1)*m;
                drho[nj] += two*sigma.Same_prefactor[m]*sigma.rho[nj];
            }
        } 
        else {
          for(j = 0; j < M; j++){
            nj = MDiagInd[j];
            drho_row_elem = &(drho[nj]);
            drho_col_elem = &(drho[Ns*nj]);
            rho_row_elem = &(sigma.rho[nj]);
            rho_col_elem = &(sigma.rho[Ns*nj]);
            for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                if (m != nj ) {
                  *drho_row_elem -= sigma.Same_prefactor[j]*(*rho_row_elem);
                  *drho_col_elem -= sigma.Same_prefactor[j]*(*rho_col_elem);	
                }
            }
          }
        }

        
    }
    


    if(sigma.TL_truncation){
        //- sum_j=1^M [ |j><j| (Q_j rho - rho Q_j^H) - (Q_j rho - rho Q_j^H) |j><j| ]
        if (sigma.active) {
            if (multiBath) {
                for(j = 0; j < M; j++){
                    nj = MDiagInd[j];
                    phi_hemm(CblasRight,Ns,&one,   sigma.rho,Q_tl[j],   &zero,Qrho_rhoQ);
                    phi_hemm(CblasLeft, Ns,&negone,sigma.rho,Q_tldag[j],&one, Qrho_rhoQ);
                    //phi_hemm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&one,sigma.rho,Ns,Q_tl[j],Ns,&zero,Qrho_rhoQ,Ns);
                    //phi_hemm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&negone,sigma.rho,Ns,Q_tldag[j],Ns,&one,Qrho_rhoQ,Ns);
                    drho_row_elem = &(drho[nj]);
                    drho_col_elem = &(drho[Ns*nj]);
                    rho_row_elem = &(Qrho_rhoQ[nj]);
                    rho_col_elem = &(Qrho_rhoQ[Ns*nj]);
                    for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                        *drho_row_elem -= (*rho_row_elem);
                        *drho_col_elem += (*rho_col_elem);	
                    }
                }
            } 
            else {
                for(nj = 0; nj < Ns; nj++){
                    phi_hemm(CblasRight,Ns,&one,   sigma.rho,Q_tl[nj],   &zero,Qrho_rhoQ);
                    phi_hemm(CblasLeft, Ns,&negone,sigma.rho,Q_tldag[nj],&one, Qrho_rhoQ);
                    //phi_hemm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&one,sigma.rho,Ns,Q_tl[nj],Ns,&zero,Qrho_rhoQ,Ns);
                    //phi_hemm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&negone,sigma.rho,Ns,Q_tldag[nj],Ns,&one,Qrho_rhoQ,Ns);
                    drho_row_elem = &(drho[nj]);
                    drho_col_elem = &(drho[Ns*nj]);
                    rho_row_elem = &(Qrho_rhoQ[nj]);
                    rho_col_elem = &(Qrho_rhoQ[Ns*nj]);
                    for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                        *drho_row_elem -= (*rho_row_elem);
                        *drho_col_elem += (*rho_col_elem);	
                    }
                }
            }
        }
    }
    else{
	for(j = 0; j < sigma.Nnext; j++){ // this is sum over Ns & Kt
            if (sigma.Next[j]->active) {
                nj = sigma.NextIndMatrix[j];
                drho_col_elem = &(drho[Ns*nj]);
                rho_col_elem = &(sigma.Next[j]->rho[Ns*nj]);
                for(m=0; m < Ns; m++){
                    drho_col_elem[m] += (sigma.Next_prefactor[j])*rho_col_elem[m];	
                }
            }
	}
	for(m = 0; m < Ns; m++){ //, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
	    drho_row_elem = &(drho[Ns*m]);
	    for(j = 0; j < sigma.Nnext; j++){
                if (sigma.Next[j]->active) {
                    nj = sigma.NextIndMatrix[j];
                    rho_row_elem = &(sigma.Next[j]->rho[Ns*m]);
                    drho_row_elem[nj] -= sigma.Next_prefactor[j]*rho_row_elem[nj];
                }
            }
	}
    }

    for(j = 0; j < sigma.Nprev; j++){
        if (sigma.Prev[j]->active) {
            m = sigma.PrevIndMatrix[j];
            drho_col_elem = &drho[Ns*m];
            rho_col_elem = &(sigma.Prev[j]->rho[Ns*m]);
            for(m = 0; m < Ns; m++) //, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                drho_col_elem[m] += sigma.Prev_prefactor_col[j]*rho_col_elem[m];	
        }
    }
    
    for(m = 0; m < Ns; m++){// rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
	drho_row_elem = &drho[m*Ns];
	for(j = 0; j < sigma.Nprev; j++){
            if (sigma.Prev[j]->active) {            
                nj = sigma.PrevIndMatrix[j];  //nj runs from 0 to Ns (e.g. 0,0,0,1,1,1,2,2,2,...,Ns,Ns,Ns  for Kt = 2)
                rho_row_elem = &(sigma.Prev[j]->rho[m*Ns]);
                drho_row_elem[nj] -= sigma.Prev_prefactor_row[j]*rho_row_elem[nj];
            }
	}
    }
};


////////////////////////////////////////// For full-diagonal system-bath
/////////// INTEGRATION STEP ///////////// coupling.
//////////////////////////////////////////
void HierarchyIntegrator::integration_step_fullBath(const HierarchyNode &sigma, Complex * drho,const int &time_step, Complex * Qrho_rhoQ, Complex * Vnext, Complex * Vprev, Complex *Vsame){
    int nj;
    int m;
    int j;
    Complex * drho_col_elem, * rho_row_elem, *rho_col_elem;

    phi_symm(CblasLeft, Ns,&L_prefactor   ,sigma.H,   sigma.rho,&zero,drho);
    phi_symm(CblasRight,Ns,&negL_prefactor,sigma.Hdag,sigma.rho,&one, drho);
    //phi_symm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&L_prefactor,sigma.H,Ns,sigma.rho,Ns,&zero,drho,Ns);
    //phi_symm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&negL_prefactor,sigma.Hdag,Ns,sigma.rho,Ns,&one,drho,Ns);


    rho_col_elem = sigma.rho;
    drho_col_elem = drho;
    for(nj = 0; nj < Ns2;nj++, drho_col_elem++, rho_col_elem++){
    	*drho_col_elem -= sigma.NZI_prefactor*(*rho_col_elem); 
    }

    rho_col_elem = sigma.rho;
    for(m = 0; m < Ns; m++){
        rho_row_elem = &(Vsame[m*Ns2]);
        for(nj = 0; nj < Ns2; nj++){
    	    drho[nj] -= rho_row_elem[nj]*rho_col_elem[nj];
	}
    }

    if(sigma.TL_truncation){
        //- sum_m [ V_m (Q_m rho - rho Q_m) - (Q_m rho - rho Q_m) V_m ]
        // Note Qrho_rhoQ includes -i, and V_m includes -i 
	for(nj = 0; nj < Ns; nj++){
	    phi_hemm(CblasRight,Ns,&neg_imag_const,
                 sigma.rho,Q_tl[nj],&zero,Qrho_rhoQ);

	    phi_hemm(CblasLeft,Ns,&imag_const,
                 sigma.rho,Q_tldag[nj],&one,Qrho_rhoQ);

	    //phi_hemm(CblasColMajor,CblasRight,CblasUpper,Ns,Ns,&neg_imag_const,
            //     sigma.rho,Ns,Q_tl[nj],Ns,&zero,Qrho_rhoQ,Ns);

	    //phi_hemm(CblasColMajor,CblasLeft,CblasUpper,Ns,Ns,&imag_const,
            //            sigma.rho,Ns,Q_tldag[nj],Ns,&one,Qrho_rhoQ,Ns);

            rho_row_elem = &(Vnext[nj*Ns2]);
            for(m = 0; m < Ns2; m++)
                drho[m] += rho_row_elem[m]*Qrho_rhoQ[m];
	    }
    }
    else{
       	for(j = 0; j < sigma.Nnext; j++){ //Nnext runs from 0 to NsKt
	    nj = sigma.NextIndMatrix[j];//runs from 0 to Ns
	    rho_col_elem = sigma.Next[j]->rho;
            rho_row_elem = &(Vnext[nj*Ns2]);
	    for(m=0; m < Ns2; m++){
                drho[m] += rho_row_elem[m]*rho_col_elem[m];
	    }
	}
    }
    for(j = 0; j < sigma.Nprev; j++){
        //nj = sigma.PrevIndMatrix[j]; //runs from 0 to Ns
        nj = sigma.PrevInd[j]; //runs from 0 to Ns*Kt-1
        rho_col_elem = sigma.Prev[j]->rho;
        rho_row_elem = &(Vprev[nj*Ns2]);
        for(m=0; m < Ns2; m++){
            drho[m] += Complex(sigma.I[nj])*rho_row_elem[m]*rho_col_elem[m];
        }
    }
    
};


/////////////////////////////////////////////////////////////////////////
////////////           ADS INTEGRATOR FUNCTION               ////////////
/////////////////////////////////////////////////////////////////////////
/*
void HierarchyIntegrator::rkf45_integrate_multistep(int id, int num_threads){
    void (HierarchyIntegrator::*take_step)(const HierarchyNode &, Complex *,const int &, 
                                        Complex *, Complex *, Complex *, Complex *);

    if (!spectrum) {
        if (corrBath)
            take_step = &HierarchyIntegrator::integration_step_fullBath;
        else {
            if (!filter_tolerance) {
                take_step = &HierarchyIntegrator::integration_step;
            } 
            else
                take_step = &HierarchyIntegrator::integration_step_adaptTrunc;
            }
    }
    else{
        if (corrBath)
            take_step = &HierarchyIntegrator::vec_integration_step_fullBath;
        else 
            take_step = &HierarchyIntegrator::vec_integration_step; 
    }
    bool restart=false;
    if(restartStep > 0)
	restart=true;

    // Time Checkers
    time_t t1,t2;
  
  
    int *m_start; 
    int *m_end;  
    int *m;
    // Memory assigned by thread. 
    if (id != num_threads){
        m_start = MatrixIndices_th[id]; //pointer to first element of MatrixIndices
        m_end = MatrixIndices_th[id]+MatrixCount_th[id]; //pointer to last element of MatrixIndices
    }
    else{
        m_start = new int;
        *m_start = 0;
        m_end = new int;
        *m_end = MatrixCount;
    }
    int m_cnt = 0;
    Complex *k0,*k1,*k2,*k3,*k4,*k5, *rho_mts;
    Complex *Qrho_rhoQ = NULL;
    
    Complex *Vbath_next = NULL;
    Complex *Vbath_prev = NULL;
    Complex *Vbath_same = NULL;

    //ADAPTIVE TS VARIABLES
    Float mydt = dt;
    Float Diff = 0;
    Complex yRKF4;
    Float Tnow = 0;
    Float Tprev = 0;
    Float nextWriteT = 0;
    Float dT = 0.0001;
    Complex dtup = 1.5 ;
    Complex dtdn = 0.5;

    RKF45_Diff[id] = 0.0;
    Float tolerance = rkf45_tolerance;
    Float mindt = rkf45_mindt;

    ofstream maxRho;

    bool AdapTrunc = false;
    if (filter_tolerance > 0) 
        AdapTrunc = true;
    barrier_wait(&WriterBarrier);    

    int Ns2_th = Ns2;
    if (spectrum) Ns2_th = Ns;
    int Ns2M_th = Ns2_th*MatrixCount_th[id];

    int **cycle_m;
    int **cycle_start_m;
    int **cycle_end_m;

    if (id < num_threads){
        //cycle_start_m[0] = cycle_m[0];
        //cycle_end_m[0] = cycle_m[0]+MatrixCount_th[id];
	k0 = new Complex [Ns2M_th];
	k1 = new Complex [Ns2M_th];
	k2 = new Complex [Ns2M_th];
	k3 = new Complex [Ns2M_th];
	k4 = new Complex [Ns2M_th];
	k5 = new Complex [Ns2M_th];
        rho_mts = new Complex [Ns2M_th];
        //debug
        //cout << "[" << id << "] initialising dt's\n";
        //debug
        for (m = m_start; m != m_end; ++m) {
            rho[*m].dt = Complex(mydt);
            rho[*m].step = 1;
        }
        //debug
        //cout << "[" << id << "] initialised dt's\n";
        //debug     
	for(int i = 0; i < MatrixCount_th[id]; i++){
            //cycle_m[0][i] = MatrixIndices_th[id][i]; 
            rho[MatrixIndices_th[id][i]].k0 = k0+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k1 = k1+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k2 = k2+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k3 = k3+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k4 = k4+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k5 = k5+i*Ns2_th;
	    for(int j = 0; j < Ns2_th; j++){
                k0[i*Ns2_th+j] = Complex(0);
                k1[i*Ns2_th+j] = Complex(0);
                k2[i*Ns2_th+j] = Complex(0);
                k3[i*Ns2_th+j] = Complex(0);
                k4[i*Ns2_th+j] = Complex(0);
                k5[i*Ns2_th+j] = Complex(0);
                rho_mts[i*Ns2_th+j] = rho[MatrixIndices_th[id][i]].rho[j];
            }
	}
        if (TL_trunc)
            Qrho_rhoQ = new Complex[Ns2_th];
        if(corrBath){
            Vbath_same = new Complex[Ns*Ns2_th];
            Vbath_next = new Complex[Ns*Ns2_th];
            Vbath_prev = new Complex[Ns*Kt*Ns2_th];
            for(int n = 0; n < Ns; n++){ //n is BChl index (aka coupling term index)
                if (!spectrum) {
                    for(int i=0; i < Ns; i++)
                        for(int j=0; j < Ns; j++){
                                Vbath_next[n*Ns*Ns+i*Ns+j] = Complex(0,-1*(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i]));
                            Vbath_same[n*Ns*Ns+i*Ns+j] = Complex(pow(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i],2))
                                                            *rho[id].Same_prefactor[n];
                        }
                    for(int k=0; k < Kt; k++)
                        for(int i=0; i < Ns; i++)
                            for(int j=0; j < Ns; j++)
                                Vbath_prev[n*Kt*Ns*Ns+k*Ns*Ns+i*Ns+j] = 
                                        Complex(Vbath_re[n*Ns+j])*C[n*Kt+k]*Complex(0,-1)
                                    - Complex(Vbath_re[n*Ns+i])*conj(C[n*Kt+k])*Complex(0,-1);
                }
                else {
                    for(int i=0; i < Ns; i++)
                        for(int j=0; j < Ns; j++){
                                Vbath_next[n*Ns*Ns+i*Ns+j] = Complex(0,-1*(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i]));
                            Vbath_same[n*Ns*Ns+i*Ns+j] = Complex(pow(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i],2))
                                                                *rho[id].Same_prefactor[n];
                        }
                    for(int k=0; k < Kt; k++)
                        for(int i=0; i < Ns; i++)
                            for(int j=0; j < Ns; j++)
                                Vbath_prev[n*Kt*Ns*Ns+k*Ns*Ns+i*Ns+j] = 
                                        Complex(Vbath_re[n*Ns+j])*C[n*Kt+k]*Complex(0,-1)
                                      - Complex(Vbath_re[n*Ns+i])*conj(C[n*Kt+k])*Complex(0,-1);
                }
            }
        }	
        if (verbose > 0) {
          pthread_mutex_lock(&io_lock);
          cout << "[" << id << "]:\tMemory for " << MatrixCount_th[id] << " density matrices assigned with (" << Ns2M_th << " elements)\n";
          pthread_mutex_unlock(&io_lock);
        }
	barrier_wait(&MatupdateBarrier);
        if (verbose > -1) {
          pthread_mutex_lock(&io_lock);
          cout << "[" << id << "]:\tIntegrating using adaptive timestepping RKF45.\n"; 
          pthread_mutex_unlock(&io_lock);
        }
    }
    else{ //thread #: numthreads
        pthread_mutex_lock(&io_lock);
        if (verbose > 1)
	cout << "[" << id << "]:\tWriter Waiting.\n"; 
        //if (AdapTrunc) {
        //    maxRho.open("MaxRhoVals.txt");
        //    maxRho << "# time ";
        //    for(int m = 0; m < MatrixCount; ++m){
        //        for(int i = 0; i < MKt-1; i++)
        //                maxRho << rho[m].I[i] << ",";
        //        maxRho << rho[m].I[MKt-1] << " | ";
        //            
        //    }
        //    maxRho << endl;
        //}

        //debug
        maxRho.open("maxDiffVals.txt");
        maxRho << "# time ";
        for(int n = 0; n < MatrixCount; ++n){
            for(int i = 0; i < MKt-1; i++)
                    maxRho << rho[n].I[i] << ",";
            maxRho << rho[n].I[MKt-1] << " | ";
                
        }
        maxRho << endl;
        //debug
        pthread_mutex_unlock(&io_lock);

    }
    barrier_wait(&WriterBarrier);

    time(&t1);
    int i = -1;
    int iAccepted = 0;
    int inext = 0;
    Complex b10(2./9.);
    Complex b20(1./12.), b21(1./4.);
    Complex b30(69./128.), b31(-243./128.), b32(135./64.);
    Complex b40(-17./12.), b41(27./4.), b42(-27./5), b43(16./15);
    Complex b50(65./432.), b51(-5./16), b52(13./16.), b53(4./27.), b54(5./144.);
    //Complex a1(2./9), a2(1./3), a3(3./4), a4(1.0), a5(5./6);
    Complex c50(1./9), c52(9./20), c53(16./45), c54(1./12);
    Complex c60(47./450), c62(12./25), c63(32./225), c64(1./30), c65(6./25);
#ifdef TIMERS
    
    timeval *state_time;
    int *state;
    state_time = new timeval[10000];
    state = new int[10000];
    int t_i = 0;
#endif  

    int stepspercycle = 1;

    while (Tnow < t){
        //debug
        //cout << "[" << id << "] next timestep\n";
        //debug
        i += 1; 

	if(id < num_threads){
            //debug
            //cout << "[" << id << "] creating next cycle\n";
            //debug
           //each thread re-create new cycle
            cycle_m = new int *[stepspercycle];
            cycle_start_m = new int *[stepspercycle];
            cycle_end_m = new int *[stepspercycle];
            for(int step = 0; step < stepspercycle; ++step){ 
                cycle_m[step] = new int[MatrixCount_th[id]];
                m_cnt = 0;
                for(m = m_start; m != m_end; ++m){
 //                   if ( !((step+stepspercycle/2)%rho[*m].step) ) {
                    if ( !((step+stepspercycle-1)%rho[*m].step) ) {
                        cycle_m[step][m_cnt] = *m;
                        ++m_cnt;
                    }
                }
                cycle_start_m[step] = cycle_m[step];
                cycle_end_m[step] = cycle_m[step]+m_cnt;
                //debug
                //cout << "[" << id << "] T=" << Tnow << " | step=" << step << " :";              
                //for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) 
                //    cout << " " << *m;
                //cout << endl;             
                //debug
            }
 //         for (int step = (int)(stepspercycle/2); step < (int)(stepspercycle*1.5); ++step) {
            if (Tnow > Tprev && TL_trunc) {
                //Iterate over system-bath coupling terms
                for (int n = id; n < M; n+=num_threads){
                       
                    constructTL_BLAS_threaded(Tnow,n);
                }
                
            }
          for(int step = 0; step < stepspercycle; ++step) {
#ifdef TIMERS
            state[t_i] = 0;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 1;
            gettimeofday(state_time+(t_i++),NULL);
#endif
            //debug
            //cout << "[" << id << "] k0\n";
            //debug
	    /////////////////////////	
	    //          k0         //
	    /////////////////////////
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) 
//	    for(m = m_start; m != m_end; ++m)
//                if ( !(step%rho[*m].step) )
                    (this->*take_step)(rho[*m], rho[*m].k0,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
 
	    //WAIT FOR k0 UPDATED
#ifdef TIMERS
            state[t_i] = -1;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 11;
            gettimeofday(state_time+(t_i++),NULL);
#endif


           //m = m_start-1;       
           //for(int j = 0; j < Ns2M_th; ++j){
           //     if (!(j%Ns2_th)) ++m;
           //     if ( !(step%rho[*m].step) )
           //     rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*b10*rho[*m].dt;
           // }

            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
                for(int j = 0; j < Ns2_th; ++j)
                    rho_mid[*m].rho[j] = rho[*m].rho[j] + rho[*m].k0[j]*b10*rho[*m].dt;
            }

#ifdef TIMERS
            state[t_i] = -11;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 2;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////


            //debug
            //cout << "[" << id << "] k1\n";
            //debug
	    /////////////////////////	
	    //          k1         //
	    /////////////////////////
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
//	    for(m = m_start; m != m_end; ++m)
//                if ( !(step%rho[*m].step) )
		(this->*take_step)(rho_mid[*m], rho[*m].k1,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            }
	    //WAIT FOR k1 UPDATED
#ifdef TIMERS
            state[t_i] = -2;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 12;
            gettimeofday(state_time+(t_i++),NULL);
#endif

           //m = m_start-1;       
           //for(int j = 0; j < Ns2M_th; ++j){
           //     if (!(j%Ns2_th)) ++m;                
           //     if ( !(step%rho[*m].step) )
           //     rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + (k0[j]*b20
           //                                  + k1[j]*b21)*rho[*m].dt;
           // }

            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
                for(int j = 0; j < Ns2_th; ++j)
                    rho_mid[*m].rho[j] = rho[*m].rho[j] + (rho[*m].k0[j]*b20
                                            + rho[*m].k1[j]*b21)*rho[*m].dt;
            }

            
#ifdef TIMERS
            state[t_i] = -12;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 3;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////

            //debug
            //cout << "[" << id << "] k2\n";
            //debug
	    /////////////////////////	
	    //         k2          //
	    /////////////////////////
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
//	    for(m = m_start; m != m_end; ++m)
//                if ( !(step%rho[*m].step) )
		(this->*take_step)(rho_mid[*m], rho[*m].k2,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            }
	    //WAIT FOR k2 UPDATED
#ifdef TIMERS
            state[t_i] = -3;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 13;
            gettimeofday(state_time+(t_i++),NULL);
#endif

           //m = m_start-1;
           //for(int j = 0; j < Ns2M_th; ++j){
           //     if (!(j%Ns2_th)) ++m;                
           //     if ( !(step%rho[*m].step) )
           //    rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + (k0[j]*b30+
           //                                  + k1[j]*b31 + k2[j]*b32)*rho[*m].dt ;
           // }

            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
                for(int j = 0; j < Ns2_th; ++j)
                    rho_mid[*m].rho[j] = rho[*m].rho[j] + (rho[*m].k0[j]*b30
                                     + rho[*m].k1[j]*b31+rho[*m].k2[j]*b32)*rho[*m].dt;
            }

#ifdef TIMERS
            state[t_i] = -13;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 4;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    ///////////////////////////

            //debug
            //cout << "[" << id << "] k3\n";
            //debug
	    //////////////////////////	
	    //          k3          //
	    //////////////////////////
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
	    //for(m = m_start; m != m_end; ++m)
            //    if ( !(step%rho[*m].step) )
		(this->*take_step)(rho_mid[*m], rho[*m].k3,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            }
	    //WAIT FOR k3 UPDATED
#ifdef TIMERS
            state[t_i] = -4;
            gettimeofday(state_time+(t_i++),NULL);
#endif
            barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 14;
            gettimeofday(state_time+(t_i++),NULL);
#endif

            //m = m_start-1;
            //for(int j = 0; j < Ns2M_th; ++j){
            //    if (!(j%Ns2_th)) ++m;                
            //    if ( !(step%rho[*m].step) )
            //    rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + (k0[j]*b40
            //                                 + k1[j]*b41 + k2[j]*b42 + k3[j]*b43)*rho[*m].dt;
            //}
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
                for(int j = 0; j < Ns2_th; ++j)
                    rho_mid[*m].rho[j] = rho[*m].rho[j] + (rho[*m].k0[j]*b40
                                     + rho[*m].k1[j]*b41+rho[*m].k2[j]*b42
                                     + rho[*m].k3[j]*b43)*rho[*m].dt;
            }

#ifdef TIMERS
            state[t_i] = -14;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 5;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////	

            //debug
            //cout << "[" << id << "] k4\n";
            //debug
            //////////////////////////	
	    //          k4          //
	    //////////////////////////
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
	    //for(m = m_start; m != m_end; ++m)
            //    if ( !(step%rho[*m].step) )
		(this->*take_step)(rho_mid[*m], rho[*m].k4,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            }

	    //WAIT FOR k4 UPDATED
#ifdef TIMERS
            state[t_i] = -5;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 15;
            gettimeofday(state_time+(t_i++),NULL);
#endif

            //m = m_start-1;
            //for(int j = 0; j < Ns2M_th; ++j){
            //    if (!(j%Ns2_th)) ++m;                
            //    if ( !(step%rho[*m].step) )
            //    rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + (k0[j]*b50
            //                                 + k1[j]*b51 + k2[j]*b52 + k3[j]*b53
            //                                 + k4[j]*b54)*rho[*m].dt;
            //}
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
                for(int j = 0; j < Ns2_th; ++j)
                    rho_mid[*m].rho[j] = rho[*m].rho[j] + (rho[*m].k0[j]*b50
                                     + rho[*m].k1[j]*b51+rho[*m].k2[j]*b52
                                     + rho[*m].k3[j]*b53+rho[*m].k4[j]*b54)*rho[*m].dt;
            }

#ifdef TIMERS
            state[t_i] = -15;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 6;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////	

            //////////////////////////	
	    //          k5          //
	    //////////////////////////
            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
	    //for(m = m_start; m != m_end; ++m)
            //    if ( !(step%rho[*m].step) )
		(this->*take_step)(rho_mid[*m], rho[*m].k5,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            }
            
	    //WAIT FOR k5 UPDATED
#ifdef TIMERS
            state[t_i] = -6;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 16;
            gettimeofday(state_time+(t_i++),NULL);
#endif
            // Now we have k0,...,k5. Use yRKF4 to 
            // calculate RKF_4 update and use rho_mid for RKF_5 update. Then
            // compare MAX|RKF_4-RKF_5| and make sure this is less than TOL

           // MaxDiff = 0;
            //Note: Ns2M_th = Ns*Ns*MatrixCount_th
           // m_cnt = -1;
           // m = m_start-1;
           // for(int j = 0; j < Ns2M_th; ++j){
            //    if (!(j%Ns2_th)) ++m;                
            //    if ( !(step%rho[*m].step) ) {
            //    yRKF4           = rho_matrices_th[id][j] 
            //                      + (k0[j]*c50+k2[j]*c52
            //                      + k3[j]*c53+k4[j]*c54)*rho[*m].dt;

            //    rho_matrices_th[id][j] = rho_matrices_th[id][j] + (k0[j]*c60
            //                                 + k2[j]*c62 + k3[j]*c63
            //                                 + k4[j]*c64 + k5[j]*c65)*rho[*m].dt;

            //    Diff = abs((yRKF4-rho_matrices_th[id][j]));
            //    if (!(j%Ns2_th)) {
            //    } 
            //    else {
            //        if ( Diff > Maxddrho[*m] )
            //             Maxddrho[*m] = Diff;
            //    }
            //    }
            //}

            for(m = cycle_start_m[step]; m != cycle_end_m[step]; ++m) {
                Maxddrho[*m] = 0;
                for(int j = 0; j < Ns2_th; ++j) {
                       yRKF4  = rho[*m].rho[j] 
                                  + (rho[*m].k0[j]*c50+rho[*m].k2[j]*c52
                                  + rho[*m].k3[j]*c53+rho[*m].k4[j]*c54)*rho[*m].dt; 
                    
                        rho[*m].rho[j] +=  (rho[*m].k0[j]*c60
                                             + rho[*m].k2[j]*c62 + rho[*m].k3[j]*c63
                                             + rho[*m].k4[j]*c64 + rho[*m].k5[j]*c65)*rho[*m].dt;
                        Diff = abs((yRKF4-rho[*m].rho[j]));
                        Maxddrho[*m] = ( Diff > Maxddrho[*m] ) ? Diff : Maxddrho[*m];
                }
            }

#ifdef TIMERS
            state[t_i] = -16;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 7;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////

            //debug
            //cout << "[" << id << "] deleting cycle vars\n";
            //debug
            delete[] cycle_m[step];
          }  // end cycle      

          //debug
          //cout << "[" << id << "] end cycles." << endl;
          //debug
          delete [] cycle_m;
          delete [] cycle_start_m;
          delete [] cycle_end_m;
        } // end if id < numthreads
        
        //debug
        //cout << "wait: " << id << endl;
        //debug
        //Get maximum error (ddrho)
        RKF45_Diff[id] = 0;
        for(int n = 0; n < MatrixCount; ++n) {
            RKF45_Diff[id] = (RKF45_Diff[id] < Maxddrho[n]) ? Maxddrho[n] : RKF45_Diff[id];
        }

        if ( id < num_threads ) { 
           // update step size based on diff in last step
          max_dtval[id] = 0.00;
          min_dtval[id] = 1e30;
          for (m = m_start; m != m_end; ++m){
                // If we're globally within tolerance and locally much below
                // tolerance, increase dt
                if ( Maxddrho[*m]  < tolerance/5  && RKF45_Diff[id] < tolerance ) {
                    //debug
                    //cout << "[" << id << "]@" << Tnow << " increment " << *m << " from " << rho[*m].dt << " to " << rho[*m].dt*dtup << " (min=" << mindt << ")\n";
                    //debug
                    //increment dt
                    rho[*m].dt *= dtup;
                    //rho[*m].dt *= 1.1;
                }
                else if ( Maxddrho[*m]  > tolerance ) {
                    //debug
                    //cout << "[" << id << "]@" << Tnow << " decrement " << rho[*m].dt << " to " << rho[*m].dt*dtdn << " (min=" << mindt << ")\n";
                    //debug
                    //decrement dt
                    if (real(rho[*m].dt) > mindt ) {
                        //rho[*m].dt *= 0.5;
                        rho[*m].dt *= dtdn;
                        if (real(rho[*m].dt) < mindt ) 
                            rho[*m].dt = Complex(mindt);
                    }
                }
                min_dtval[id] = (real(rho[*m].dt) < min_dtval[id]) ? real(rho[*m].dt) : min_dtval[id];
                max_dtval[id] = (real(rho[*m].dt) > max_dtval[id]) ? real(rho[*m].dt) : max_dtval[id];
          }

          // If we're globally within tolerance, make backup of current state
          if ( RKF45_Diff[id] < tolerance  ) {
            for(int j = 0; j < Ns2M_th; ++j){ 
                rho_mts[j] =  rho_matrices_th[id][j];      
            }
          }
          // else we restore old state
          else {
            for(int j = 0; j < Ns2M_th; ++j){ 
                rho_matrices_th[id][j] = rho_mts[j];
            }
          }
        }

	//////////////////////////
	// WAIT FOR rho UPDATED //
	//////////////////////////
#ifdef TIMERS
        state[t_i] = -7;
        gettimeofday(state_time+(t_i++),NULL);
#endif
      	barrier_wait(&WriterBarrier);  //Barrier for all num_threads+1 threads
#ifdef TIMERS
        state[t_i] = 20;
        gettimeofday(state_time+(t_i++),NULL);
#endif
        // Use Writer Thread to distribute minimum dt value;
        if(id == num_threads) {
            //debug
            //maxRho << Tnow;
            //for(int n = 0; n < MatrixCount; ++n)
            //     maxRho << " " << rho[n].step;
            //maxRho << endl;
            //maxRho << Tnow;
            //for(int n = 0; n < MatrixCount; ++n)
            //     maxRho << " " << Maxddrho[n];
            //maxRho << endl;
            //debug
            max_dtval[id] = -1.00;
            min_dtval[id] = 1e30;
            for(int n = 0; n  < num_threads; n++) {
                min_dtval[id] = ( min_dtval[id] > min_dtval[n] ) ? min_dtval[n] : min_dtval[id];
                max_dtval[id] = ( max_dtval[id] < max_dtval[n] ) ? max_dtval[n] : max_dtval[id];
            }
            for(int n = 0; n  < num_threads; n++) {
                min_dtval[n] = min_dtval[id];
                max_dtval[n] = max_dtval[id];
            }
            RKF45_Diff[id] = RKF45_Diff[0];

        }
#ifdef TIMERS
        state[t_i] = -20;
        gettimeofday(state_time+(t_i++),NULL);
#endif
	barrier_wait(&WriterBarrier);  //Barrier for all num_threads+1 threads
#ifdef TIMERS
        state[t_i] = 21;
        gettimeofday(state_time+(t_i++),NULL);
#endif
    
        if ( id < num_threads) {
            //update dts to b integer values of min_dtval
            for (m=m_start;m < m_end; ++m) {
                rho[*m].step = (int)floor(real(rho[*m].dt)/min_dtval[id]);
                if  ( rho[*m].step > 5 )
                       rho[*m].step = 5; 
                rho[*m].dt   = Complex(min_dtval[id]*rho[*m].step);                
            } 
        }

        //update stepspercycle
        stepspercycle = (int)floor(max_dtval[id]/min_dtval[id]);
        if (stepspercycle > 5) stepspercycle = 5;
        //debug
        //if ( id == num_threads ) {
            //cout << "[" << id << "] " << Tnow << " stepspercycle=" << stepspercycle <<  " rkf45diff=" << RKF45_Diff[id] << endl;
            //for (int n = 0; n < MatrixCount; ++n) {
            //    cout << "\t" << n << ": dt=" << real(rho[n].dt) << endl;
            //}
        //}
        //debug
        //debug
        //pthread_mutex_lock(&io_lock);
        //for (int n = 0; n <= num_threads; ++n)
        //    if (n == id)
        //        cout << "[" << id << "] " << Tnow << " stepspercycle=" << stepspercycle
        //             <<  " rkf45diff=" << RKF45_Diff[id] << " mydt=" << mydt 
        //             <<   " mindt[id]=" << min_dtval[id] << " maxdt[id]=" << max_dtval[id] << endl;
        //pthread_mutex_unlock(&io_lock);
	//barrier_wait(&WriterBarrier);
        //if (id == num_threads)
        //    cout << endl;
        //debug
        if ( RKF45_Diff[id] < tolerance ) {
          Tprev = Tnow;
          Tnow += mydt;
          mydt = stepspercycle*min_dtval[id];
          iAccepted += 1;
        }
#ifdef TIMERS
        state[t_i] = -21;
        gettimeofday(state_time+(t_i++),NULL);
#endif
	barrier_wait(&WriterBarrier);
#ifdef TIMERS
        state[t_i] = 22;
        gettimeofday(state_time+(t_i++),NULL);
#endif
        if(id == num_threads && RKF45_Diff[id] < tolerance) { // && ( MaxDiff < tolerance || dt < mindt))
            dt = mydt;
            if (i >= inext){ //Write out some timing info
                time(&t2);
                pthread_mutex_lock(&io_lock);
                outf << "#Info: Time for single timestep: " << difftime(t2,t1)/(i+1) << " s.\n";
                outf << "#Info: Estimated time for 1 ps: " << difftime(t2,t1)/Tnow/60/60 << " hrs.\n";
                outf << "#Info: Total running time: " << difftime(t2,t1)*t/Tnow/60/60/24 << " days.\n";
                pthread_mutex_unlock(&io_lock);
                inext += 100;
            }       
            if ( Tnow > nextWriteT ) {         
                pthread_mutex_lock(&io_lock);
                //if (AdapTrunc) {
                //    maxRho << Tnow;
                //    for (int m = 0; m < MatrixCount; m++){
                //        maxRho << " " << MaxRhoVal[m];
                //    }
                //    maxRho << endl;
                //}

                // NOTE: WRITING OUTPUT IN COLUMN MAJOR ORDERING
                outf << Tnow << " ";
                for( int j = 0; j < Ns2_th; j++)
                   outf << rho_mid[0].rho[j] << " ";  
                outf << endl;

                // Write out restart file if neccessary
                if(restart && i%restartStep == 0 && i > 0){
                    restartf.seekp(0, ios_base::beg);
                    restartf << i*dt << endl;
                    for(int n = 0; n < MatrixCount; n++)
                        for(int mi = 0; mi < Ns2_th; mi++)
                            restartf << rho[n].rho[mi] << "\n";
                }
                pthread_mutex_unlock(&io_lock);
                //  End of restart loop
                nextWriteT = Tnow + dT;
            }
        }
#ifdef TIMERS
        state[t_i] = -22;
        gettimeofday(state_time+(t_i++),NULL);
        if(t_i >= 10000)
            t_i = 0;
#endif

    } // End of integration loop
   

#ifdef TIMERS
    output_timing(id, num_threads, state_time, state, t_i);
#endif 
    int nActive = 0;
    for(int j = 0; j < MatrixCount; j++)
        if (rho[j].active)
            nActive++;
    pthread_mutex_lock(&io_lock);
    if (AdapTrunc && id == num_threads && verbose > 0)
        cout << "Total number of active matrices = " << nActive << "/" << MatrixCount << endl;
    if (id == num_threads && verbose > -1) {
        cout << "Total number of integration steps = "<< i+1 << endl;
        cout << "Total number of accepted integration steps = " << iAccepted << endl;
        time(&t2);
        cout << "Total integration time = " << difftime(t2,t1)/60 << " mins.\n";
    }
    
    if(restart && id == num_threads){
	restartf.seekp(0, ios_base::beg);
	restartf << t << endl << setw(40) << setprecision(20);
	for(int mc = 0; mc < MatrixCount; mc++){
	    for(int mi = 0; mi < Ns2_th; mi++)
                restartf << rho[mc].rho[mi] << "\n";
	}
    }

//    if (id == num_threads)
//        maxRho.close();
    pthread_mutex_unlock(&io_lock);
};
*/

/////////////////////////////////////////////////////////////////////////
////////////           ADS INTEGRATOR FUNCTION               ////////////
/////////////////////////////////////////////////////////////////////////
void HierarchyIntegrator::rkf45_integrate(int id, int num_threads){
    void (HierarchyIntegrator::*take_step)(const HierarchyNode &, Complex *,const int &, 
                                        Complex *, Complex *, Complex *, Complex *);

    if (!spectrum) {
        if (corrBath)
            take_step = &HierarchyIntegrator::integration_step_fullBath;
        else {
            if (!filter_tolerance) {
                take_step = &HierarchyIntegrator::integration_step;
            } 
            else
                take_step = &HierarchyIntegrator::integration_step_adaptTrunc;
            }
    }
    else{
        if (corrBath)
            take_step = &HierarchyIntegrator::vec_integration_step_fullBath;
        else 
            take_step = &HierarchyIntegrator::vec_integration_step; 
    }

    bool restart=false;
    if(restartStep > 0)
	restart=true;

    // Time Checkers
    time_t t1,t2;
  
  
    int *m_start; 
    int *m_end;  
    int *m;
    // Memory assigned by thread. 
    if (id != num_threads){
        m_start = MatrixIndices_th[id]; //pointer to first element of MatrixIndices
        m_end = MatrixIndices_th[id]+MatrixCount_th[id]; //pointer to last element of MatrixIndices
    }
    else{
        m_start = new int;
        *m_start = 0;
        m_end = new int;
        *m_end = MatrixCount;
    }
    int m_cnt = 0;
    Complex *k0,*k1,*k2,*k3,*k4,*k5;
    Complex *Qrho_rhoQ = NULL;
    
    Complex *Vbath_next = NULL;
    Complex *Vbath_prev = NULL;
    Complex *Vbath_same = NULL;

    //ADAPTIVE TS VARIABLES
    Float mydt = dt;
    Float Diff = 0, MaxDiff = 0;
    Complex yRKF4;
    
    //Keep track of time
    Float Tnow = 0;
    Float Tprev = 0;
    Float nextWriteT = 0;
    //0.1 fs time resolution in output file
    Float dT = 0.0001;
    
    Complex dtup = 1.1;
    Complex dtdn = 0.5;

    RKF45_Diff[id] = 0.0;
    Float tolerance = rkf45_tolerance;
    Float mindt = rkf45_mindt;
    ofstream maxRho;

    bool AdapTrunc = false;
    if (filter_tolerance > 0) 
        AdapTrunc = true;
    barrier_wait(&WriterBarrier);    

    int Ns2_th = Ns2;
    if (spectrum) Ns2_th = Ns;
    int Ns2M_th = Ns2_th*MatrixCount_th[id];
#ifdef TESTING
//    Float NodeMaxDiff[MatrixCount_th[id]];
#endif
    if (id < num_threads){
	k0 = new Complex [Ns2M_th];
	k1 = new Complex [Ns2M_th];
	k2 = new Complex [Ns2M_th];
	k3 = new Complex [Ns2M_th];
	k4 = new Complex [Ns2M_th];
	k5 = new Complex [Ns2M_th];
     
	for(int i = 0; i < MatrixCount_th[id]; i++){
            rho[MatrixIndices_th[id][i]].k0 = k0+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k1 = k1+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k2 = k2+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k3 = k3+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k4 = k4+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k5 = k5+i*Ns2_th;
	    for(int j = 0; j < Ns2_th; j++){
                k0[i*Ns2_th+j] = Complex(0);
                k1[i*Ns2_th+j] = Complex(0);
                k2[i*Ns2_th+j] = Complex(0);
                k3[i*Ns2_th+j] = Complex(0);
                k4[i*Ns2_th+j] = Complex(0);
                k5[i*Ns2_th+j] = Complex(0);
            }
	}
        if (TL_trunc)
            Qrho_rhoQ = new Complex[Ns2_th];
        if(corrBath){
            Vbath_same = new Complex[Ns*Ns2_th];
            Vbath_next = new Complex[Ns*Ns2_th];
            Vbath_prev = new Complex[Ns*Kt*Ns2_th];
            for(int n = 0; n < Ns; n++){ //n is BChl index (aka coupling term index)
                if (!spectrum) {
                    for(int i=0; i < Ns; i++)
                        for(int j=0; j < Ns; j++){
                                Vbath_next[n*Ns*Ns+i*Ns+j] = Complex(0,-1*(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i]));
                            Vbath_same[n*Ns*Ns+i*Ns+j] = Complex(pow(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i],2))
                                                            *rho[id].Same_prefactor[n];
                        }
                    for(int k=0; k < Kt; k++)
                        for(int i=0; i < Ns; i++)
                            for(int j=0; j < Ns; j++)
                                Vbath_prev[n*Kt*Ns*Ns+k*Ns*Ns+i*Ns+j] = 
                                        Complex(Vbath_re[n*Ns+j])*C[n*Kt+k]*Complex(0,-1)
                                    - Complex(Vbath_re[n*Ns+i])*conj(C[n*Kt+k])*Complex(0,-1);
                }
                else {
                    for(int i=0; i < Ns; i++)
                        for(int j=0; j < Ns; j++){
                                Vbath_next[n*Ns*Ns+i*Ns+j] = Complex(0,-1*(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i]));
                            Vbath_same[n*Ns*Ns+i*Ns+j] = Complex(pow(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i],2))
                                                                *rho[id].Same_prefactor[n];
                        }
                    for(int k=0; k < Kt; k++)
                        for(int i=0; i < Ns; i++)
                            for(int j=0; j < Ns; j++)
                                Vbath_prev[n*Kt*Ns*Ns+k*Ns*Ns+i*Ns+j] = 
                                        Complex(Vbath_re[n*Ns+j])*C[n*Kt+k]*Complex(0,-1)
                                      - Complex(Vbath_re[n*Ns+i])*conj(C[n*Kt+k])*Complex(0,-1);
                }
            }
        }	
        if (verbose > 0) {
          pthread_mutex_lock(&io_lock);
          cout << "[" << id << "]:\tMemory for " << MatrixCount_th[id] 
               << " density matrices assigned (" << Ns2M_th << " elements)\n" << flush;
          pthread_mutex_unlock(&io_lock);
        }
	barrier_wait(&MatupdateBarrier);
        //debug
        //cout << " verbose=" << verbose << "\n";
        // debug
        if (verbose > -1) {
          pthread_mutex_lock(&io_lock);
          cout << "[" << id << "]:\tIntegrating using adaptive timestepping RKF45.\n" << flush; 
          pthread_mutex_unlock(&io_lock);
        }
    }
    else{ //thread #: numthreads
          if (verbose > 1) {
              pthread_mutex_lock(&io_lock);
              cout << "[" << id << "]:\tWriter Waiting.\n" << flush; 
              pthread_mutex_unlock(&io_lock);
          }
        //if (AdapTrunc) {
        //    maxRho.open("MaxRhoVals.txt");
        //    maxRho << "# time ";
        //    for(int m = 0; m < MatrixCount; ++m){
        //        for(int i = 0; i < MKt-1; i++)
        //                maxRho << rho[m].I[i] << ",";
        //        maxRho << rho[m].I[MKt-1] << " | ";
        //            
        //    }
        //    maxRho << endl;
        //}

    }
    barrier_wait(&WriterBarrier);
    //if (id == num_threads && verbose > -1) {
    //  pthread_mutex_lock(&io_lock);
    //  cout << "Integrating using adaptive timestepping RKF45.\n" << flush;
    //  pthread_mutex_unlock(&io_lock);
    //}

    time(&t1);
    int i = -1;
    int iAccepted = 0;
    int inext = 0;
    Complex b10(2./9.);
    Complex b20(1./12.), b21(1./4.);
    Complex b30(69./128.), b31(-243./128.), b32(135./64.);
    Complex b40(-17./12.), b41(27./4.), b42(-27./5), b43(16./15);
    Complex b50(65./432.), b51(-5./16), b52(13./16.), b53(4./27.), b54(5./144.);
    //Complex a1(2./9), a2(1./3), a3(3./4), a4(1.0), a5(5./6);
    Complex c50(1./9), c52(9./20), c53(16./45), c54(1./12);
    Complex c60(47./450), c62(12./25), c63(32./225), c64(1./30), c65(6./25);
    b10 *= dt_c;
    b20 *= dt_c; b21 *= dt_c;
    b30 *= dt_c; b31 *= dt_c; b32 *= dt_c;
    b40 *= dt_c; b41 *= dt_c; b42 *= dt_c; b43 *= dt_c;
    b50 *= dt_c; b51 *= dt_c; b52 *= dt_c; b53 *= dt_c; b54 *= dt_c;
    c50 *= dt_c;              c52 *= dt_c; c53 *= dt_c; c54 *= dt_c;
    c60 *= dt_c;              c62 *= dt_c; c63 *= dt_c; c64 *= dt_c; c65 *= dt_c;
#ifdef TIMERS
    
    timeval *state_time;
    int *state;
    state_time = new timeval[10000];
    state = new int[10000];
    int t_i = 0;
#endif  


    while (Tnow < t){
        i += 1;
	if(id < num_threads){
            if (Tnow > Tprev && TL_trunc) {
                //Iterate over system-bath coupling terms
                for (int idx = id; idx < M; idx+=num_threads){
                       
                    constructTL_BLAS_threaded(Tnow,idx);
                }
                
            }
#ifdef TIMERS
            state[t_i] = 0;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 1;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////	
	    //          k0         //
	    /////////////////////////
	    for(m = m_start; m != m_end; ++m)
		(this->*take_step)(rho[*m], rho[*m].k0,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
 
	    //WAIT FOR k0 UPDATED
#ifdef TIMERS
            state[t_i] = -1;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 11;
            gettimeofday(state_time+(t_i++),NULL);
#endif

#ifdef BLASUPDATE
            phi_copy(Ns2M_th,rho_matrices_th[id],1,rho_mid_matrices_th[id],1);           // p <- r
            phi_axpy(Ns2M_th,&b10,k0,1,rho_mid_matrices_th[id],1);
#else
           for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*b10;
            }
#endif

#ifdef TIMERS
            state[t_i] = -11;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 2;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////


	    /////////////////////////	
	    //          k1         //
	    /////////////////////////
            //cout << "k1\n";
	    for(m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k1,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
	    //WAIT FOR k1 UPDATED
#ifdef TIMERS
            state[t_i] = -2;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 12;
            gettimeofday(state_time+(t_i++),NULL);
#endif

#ifdef BLASUPDATE
            phi_copy(Ns2M_th,rho_matrices_th[id],1,rho_mid_matrices_th[id],1);           // p <- r
            phi_axpy(Ns2M_th,&b20,k0,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b21,k1,1,rho_mid_matrices_th[id],1);
#else                        
           for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*b20
                                             + k1[j]*b21;
            }
           
#endif
            
#ifdef TIMERS
            state[t_i] = -12;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 3;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////

	    /////////////////////////	
	    //         k2          //
	    /////////////////////////
            //cout << "k2\n";
	    for(m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k2,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            
	    //WAIT FOR k2 UPDATED
#ifdef TIMERS
            state[t_i] = -3;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 13;
            gettimeofday(state_time+(t_i++),NULL);
#endif

#ifdef BLASUPDATE
            phi_copy(Ns2M_th,rho_matrices_th[id],1,rho_mid_matrices_th[id],1);           // p <- r
            phi_axpy(Ns2M_th,&b30,k0,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b31,k1,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b32,k2,1,rho_mid_matrices_th[id],1);
#else           
           for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*b30+
                                             + k1[j]*b31 + k2[j]*b32 ;
            }
#endif

#ifdef TIMERS
            state[t_i] = -13;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 4;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    ///////////////////////////

	    //////////////////////////	
	    //          k3          //
	    //////////////////////////
	    for(m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k3,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
	    //WAIT FOR k3 UPDATED
#ifdef TIMERS
            state[t_i] = -4;
            gettimeofday(state_time+(t_i++),NULL);
#endif
            barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 14;
            gettimeofday(state_time+(t_i++),NULL);
#endif

#ifdef BLASUPDATE
            phi_copy(Ns2M_th,rho_matrices_th[id],1,rho_mid_matrices_th[id],1);           // p <- r
            phi_axpy(Ns2M_th,&b40,k0,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b41,k1,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b42,k2,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b43,k3,1,rho_mid_matrices_th[id],1);
#else    
            for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*b40
                                             + k1[j]*b41 + k2[j]*b42 + k3[j]*b43;
            }
#endif

#ifdef TIMERS
            state[t_i] = -14;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 5;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////	

            //////////////////////////	
	    //          k4          //
	    //////////////////////////
	    for(m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k4,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
           

	    //WAIT FOR k4 UPDATED
#ifdef TIMERS
            state[t_i] = -5;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 15;
            gettimeofday(state_time+(t_i++),NULL);
#endif

#ifdef BLASUPDATE
            phi_copy(Ns2M_th,rho_matrices_th[id],1,rho_mid_matrices_th[id],1);           // p <- r
            phi_axpy(Ns2M_th,&b50,k0,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b51,k1,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b52,k2,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b53,k3,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&b54,k4,1,rho_mid_matrices_th[id],1);
#else    
            for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*b50
                                             + k1[j]*b51 + k2[j]*b52 + k3[j]*b53
                                             + k4[j]*b54;
            }
#endif

#ifdef TIMERS
            state[t_i] = -15;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 6;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////	

            //////////////////////////	
	    //          k5          //
	    //////////////////////////
	    for(m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k5,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
            
	    //WAIT FOR k5 UPDATED
#ifdef TIMERS
            state[t_i] = -6;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 16;
            gettimeofday(state_time+(t_i++),NULL);
#endif
            // Now we have k0,...,k5. Use yRKF4 to 
            // calculate RKF_4 update and use rho_mid for RKF_5 update. Then
            // compare MAX|RKF_4-RKF_5| and make sure this is less than TOL

#ifdef BLASUPDATE
            phi_copy(Ns2M_th,rho_matrices_th[id],1,rho_mid_matrices_th[id],1);           // p <- r
            phi_axpy(Ns2M_th,&c60,k0,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&c62,k2,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&c63,k3,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&c64,k4,1,rho_mid_matrices_th[id],1);
            phi_axpy(Ns2M_th,&c65,k5,1,rho_mid_matrices_th[id],1);
            MaxDiff = 0;
	    for(m = m_start; m != m_end; ++m){
		rho_mid_elem = rho_mid[*m].rho;
                rho_elem = rho[*m].rho;
                MaxRhoVal[*m] = 0;
		for(int j = 0; j < Ns2_th; ++j) {
                    yRKF4           = rho_elem[j] 
                                      + rho[*m].k0[j]*c50+rho[*m].k2[j]*c52
                                      + rho[*m].k3[j]*c53+rho[*m].k4[j]*c54;
                    Diff = abs((yRKF4-rho_mid_elem[j]));
                    //if(Diff > MaxDiff ) { MaxDiff = Diff; }
                    MaxDiff = ( Diff > MaxDiff ) ? Diff : MaxDiff;
                    if(AdapTrunc && abs(rho_mid_elem[j]) > MaxRhoVal[*m]) MaxRhoVal[*m] = abs(rho_mid_elem[j]);
                }
	    }
            RKF45_Diff[id] = MaxDiff;
#else  
            MaxDiff = 0;
            //Note: Ns2M_th = Ns*Ns*MatrixCount_th
            m_cnt = -1;
            for(int j = 0; j < Ns2M_th; ++j){
                yRKF4           = rho_matrices_th[id][j] 
                                  + k0[j]*c50+k2[j]*c52
                                  + k3[j]*c53+k4[j]*c54;
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k0[j]*c60
                                             + k2[j]*c62 + k3[j]*c63
                                             + k4[j]*c64 + k5[j]*c65;
                Diff = abs((yRKF4-rho_mid_matrices_th[id][j]));
//                if(Diff > MaxDiff)  MaxDiff = Diff; 
                MaxDiff = ( Diff > MaxDiff ) ? Diff : MaxDiff;

            }
            RKF45_Diff[id] = MaxDiff;
#endif

#ifdef TIMERS
            state[t_i] = -16;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    barrier_wait(&MatupdateBarrier);
#ifdef TIMERS
            state[t_i] = 7;
            gettimeofday(state_time+(t_i++),NULL);
#endif
	    /////////////////////////	
	}
        //cout << "wait: " << id << endl;
	//////////////////////////
	// WAIT FOR rho UPDATED //
	//////////////////////////
#ifdef TIMERS
        state[t_i] = -7;
        gettimeofday(state_time+(t_i++),NULL);
#endif
      	barrier_wait(&WriterBarrier);  //Barrier for all num_threads+1 threads
#ifdef TIMERS
        state[t_i] = 20;
        gettimeofday(state_time+(t_i++),NULL);
#endif
        // Use Writer Thread to distribute MaxDiff
        if(id == num_threads) {
            //cout << "Distrib\n";
            MaxDiff=0;
            for(int n = 0; n  < num_threads; n++)
                if (MaxDiff < RKF45_Diff[n])
                    MaxDiff = RKF45_Diff[n];
            for(int n = 0; n  < num_threads+1; n++)
                RKF45_Diff[n] = MaxDiff;
        }
#ifdef TIMERS
        state[t_i] = -20;
        gettimeofday(state_time+(t_i++),NULL);
#endif
        //////////////////////////////////////////////////////////
        // Each thread compare MaxDiff with tolerance and       //
        // update mydt & dt_c. Calc Threads update rho and      //
        // writer thread writes out rho_mid (with stored RKF_5  //
        // result) if MaxDiff < tol                             //
        //////////////////////////////////////////////////////////
	barrier_wait(&WriterBarrier);  //Barrier for all num_threads+1 threads
#ifdef TIMERS
        state[t_i] = 21;
        gettimeofday(state_time+(t_i++),NULL);
#endif
        if ( RKF45_Diff[id] < tolerance || dt < mindt ) {
            // Increment time
            Tprev = Tnow;
            Tnow += mydt;
            iAccepted += 1;
            // Increase timestep if possible
            if (RKF45_Diff[id] < tolerance/10) {
                mydt *= real(dtup);
                b10 *= dtup;
                b20 *= dtup; b21 *= dtup;
                b30 *= dtup; b31 *= dtup; b32 *= dtup;
                b40 *= dtup; b41 *= dtup; b42 *= dtup; b43 *= dtup;
                b50 *= dtup; b51 *= dtup; b52 *= dtup; b53 *= dtup; b54 *= dtup;
                c50 *= dtup;              c52 *= dtup; c53 *= dtup; c54 *= dtup;
                c60 *= dtup;              c62 *= dtup; c63 *= dtup; c64 *= dtup; c65 *= dtup;                
            
               // if (id == num_threads)
               //     cout << "INCRSE @T="<<Tnow<<" dt: " << mydt/1.1 << "->" << mydt 
               //         << " Diff=" << RKF45_Diff[id] << endl;
                
            }
            // copy rho_mid -> rho.
            if (id < num_threads) {
                for(m = m_start; m != m_end; ++m){
                    if (AdapTrunc) MaxRhoVal[*m] = 0;
                    for(int j = 0; j < Ns2_th; ++j) {
                        rho[*m].rho[j] = rho_mid[*m].rho[j];
                        if ( AdapTrunc && abs(rho[*m].rho[j]) > MaxRhoVal[*m] )
                            MaxRhoVal[*m] = abs(rho[*m].rho[j]);
                    }
                    if (AdapTrunc) {
                        if ( MaxRhoVal[*m] < filter_tolerance && rho[*m].active ) {
                            rho[*m].active = false;
                            rho_mid[*m].active = false;
                        }
                        else if ( MaxRhoVal[*m] > filter_tolerance && rho[*m].active == false ) {
                            rho[*m].active = true;
                            rho_mid[*m].active = true;
                        }
                    }
                }
            } 
            else {
                dt = mydt;
                dt_c = Complex(mydt);
            }
        }
        else {
            // Decrease timestep
            mydt *= real(dtdn);
            b10 *= dtdn;
            b20 *= dtdn; b21 *= dtdn;
            b30 *= dtdn; b31 *= dtdn; b32 *= dtdn;
            b40 *= dtdn; b41 *= dtdn; b42 *= dtdn; b43 *= dtdn;
            b50 *= dtdn; b51 *= dtdn; b52 *= dtdn; b53 *= dtdn; b54 *= dtdn;
            c50 *= dtdn;              c52 *= dtdn; c53 *= dtdn; c54 *= dtdn;
            c60 *= dtdn;              c62 *= dtdn; c63 *= dtdn; c64 *= dtdn; c65 *= dtdn;                

            if (id == num_threads) {
               dt = mydt;            //Updates dt for all threads
               dt_c = Complex(mydt); //Updates dt_c for all threads
               if (dt < 1e-10) cout << "WARNING!!! dt < 1e-10\n";
            }
        }
#ifdef TIMERS
        state[t_i] = -21;
        gettimeofday(state_time+(t_i++),NULL);
#endif
	barrier_wait(&WriterBarrier);
#ifdef TIMERS
        state[t_i] = 22;
        gettimeofday(state_time+(t_i++),NULL);
#endif
        if(id == num_threads && ( MaxDiff < tolerance || dt < mindt)) {


            if (i >= inext){ //Write out some timing info
                time(&t2);
                pthread_mutex_lock(&io_lock);
                outf << "#Info: Time for single timestep: " << difftime(t2,t1)/(i+1) << " s.\n";
                outf << "#Info: Estimated time for 1 ps: " << difftime(t2,t1)/Tnow/60/60 << " hrs.\n";
                outf << "#Info: Total running time: " << difftime(t2,t1)*t/Tnow/60/60/24 << " days.\n";
                pthread_mutex_unlock(&io_lock);
                inext += 100;
            }       
            if ( Tnow > nextWriteT ) {         
                pthread_mutex_lock(&io_lock);
                //if (AdapTrunc) {
                //    maxRho << Tnow;
                //    for (int m = 0; m < MatrixCount; m++){
                //        maxRho << " " << MaxRhoVal[m];
                //    }
                //    maxRho << endl;
                //}

                // NOTE: WRITING OUTPUT IN COLUMN MAJOR ORDERING
                outf << Tnow << " ";
                for( int j = 0; j < Ns2_th; j++)
                   outf << rho_mid[0].rho[j] << " ";  
                outf << endl;

                // Write out restart file if neccessary
                if(restart && i%restartStep == 0 && i > 0){
                    restartf.seekp(0, ios_base::beg);
                    restartf << i*dt << endl;
                    for(int mat = 0; mat < MatrixCount; mat++)
                        for(int mi = 0; mi < Ns2_th; mi++)
                            restartf << rho[mat].rho[mi] << "\n";
                }
                pthread_mutex_unlock(&io_lock);
                //  End of restart loop
                nextWriteT = Tnow + dT;
            }
        }
#ifdef TIMERS
        state[t_i] = -22;
        gettimeofday(state_time+(t_i++),NULL);
        if(t_i >= 10000)
            t_i = 0;
#endif

    } // End of integration loop
   

#ifdef TIMERS
    output_timing(id, num_threads, state_time, state, t_i);
#endif 

    time(&t2);

    int nActive = 0;
    for(int j = 0; j < MatrixCount; j++)
        if (rho[j].active)
            nActive++;
    if (id == num_threads  && verbose > -1) {
        pthread_mutex_lock(&io_lock);
        if (AdapTrunc && verbose > 0)
          cout << "Total number of active matrices = " << nActive << "/" << MatrixCount << endl;
        
        cout << "Total number of integration steps = "<< i+1 << endl;
        cout << "Total number of accepted integration steps = " << iAccepted << endl;
        cout << "Total integration time = " << difftime(t2,t1)/60 << " mins.\n";
        pthread_mutex_unlock(&io_lock);
    }
    
    if(restart && id == num_threads){
	restartf.seekp(0, ios_base::beg);
	restartf << t << endl << setw(40) << setprecision(20);
	for(int mj = 0; mj < MatrixCount; mj++){
	    for(int mi = 0; mi < Ns2_th; mi++)
                restartf << rho[mj].rho[mi] << "\n";
	}
    }

    barrier_wait(&WriterBarrier);  //Barrier for all num_threads+1 threads
    //pthread_exit(0);

};


void HierarchyIntegrator::output_timing(int id, int num_threads, timeval * state_time, int * state, int Nt){
    float integration_time = 0;
    float vecadd_time = 0;
    float waiting_time = 0;
    float writing_time = 0;
    float maxdiff_time = 0;
    float timestepvar_time = 0;
    float bstart_t, bstep_t, bdiff_t, bats_t;
    float b00_t, b01_t, b10_t, b11_t;
    float b20_t, b21_t, b30_t, b31_t;
    float b40_t, b41_t, b50_t, b51_t;
    bstart_t = bstep_t =  bdiff_t = bats_t = 0;
    b00_t = b01_t = b10_t = b11_t = b20_t = b21_t =  0;
    b30_t = b31_t = b40_t = b41_t = b50_t = b51_t =  0;



    ofstream timing_out;
    stringstream ss;
    string fname;
    ss << id;
    fname = "timing_t-"+ss.str()+".txt";
    timing_out.open(fname.c_str());

    float td=0;
    float total=(writing_time+integration_time+vecadd_time+waiting_time+timestepvar_time);
    pthread_mutex_lock(&io_lock);    
    timing_out << "#state_time state\n";
    timing_out << id << " " << 0 << " " << state[0] << "\n";
    for (int j = 1; j < Nt; ++j){
        td = state_time[j].tv_sec-state_time[j-1].tv_sec;
        td += (state_time[j].tv_usec-state_time[j-1].tv_usec)*1e-6;
        switch(state[j-1]){
            case 0:
                bstart_t += td;
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
                integration_time += td;
                break;
            case 11:
            case 12:
            case 13:
            case 14:
            case 15:
            case 16:
                vecadd_time += td;
                break;
            case -1:
                b00_t += td; break;
            case -2:
                b10_t += td; break;
            case -3:
                b20_t += td; break;
            case -4:
                b30_t += td; break;
            case -5:
                b40_t += td; break;
            case -6:
                b50_t += dt; break;
            case -11:
                b01_t += td; break;
            case -12:
                b11_t += td; break;
            case -13:
                b21_t += td; break;
            case -14:
                b31_t += td; break;
            case -15:
                b41_t += td; break;
            case -16:
                b51_t += td; break;
            case -7:  //waiting time for full step
                bstep_t += td; break;
            case 20:  //calculate max difference between rkf4 and rkf5
                maxdiff_time += td; 
            case -20: //wait for max difference between rkf4 and rkf5
                bdiff_t += td; break;
            case 21:  //update adaptive timestep
                timestepvar_time += td; break;
            case -21: //adaptive_ts_update_wait
                bats_t += td; break;
            case 22:
                writing_time += td; break;
            case -22:
                break;
                
        
        } 
        td = state_time[j].tv_sec-state_time[0].tv_sec;
        td += (state_time[j].tv_usec-state_time[0].tv_usec)*1e-6;
        timing_out << id << " " << td*1000 << " " << state[j-1] << "\n";
        timing_out << id << " " << td*1000 << " " << state[j] << "\n";
    }
    timing_out << "\n";


    timing_out << (float)id+0.5 << " " << 0 << " " << state[0] << "\n";
    for (int j = 1; j < Nt; ++j){
        td = state_time[j].tv_sec-state_time[0].tv_sec;
        td += (state_time[j].tv_usec-state_time[0].tv_usec)*1e-6; 
        timing_out << (float)id+0.5 << " " << td*1000 << " " << state[j-1] << "\n";
        timing_out << (float)id+0.5 << " " << td*1000 << " " << state[j] << "\n";
    }
    timing_out << "\n";
    timing_out.close();

    waiting_time = bstart_t + bstep_t + bdiff_t + bats_t;
    waiting_time += b00_t+b01_t+b10_t+b11_t+b20_t+b21_t;
    waiting_time += b30_t+b31_t+b40_t+b41_t+b50_t+b51_t;
    total = integration_time + vecadd_time + waiting_time + writing_time + timestepvar_time + maxdiff_time;
    if (id == 0) {
        cout << "Total time counts (total = " << total << " s):\n"
             << "     Integrating |   Updating  |   Waiting   |   Writing  |  Adaptive Update \n"; 
    }
    pthread_mutex_unlock(&io_lock);
    barrier_wait(&WriterBarrier);
    pthread_mutex_lock(&io_lock);    
    if (id < num_threads) {
        cout << "[" << id << "] " << setprecision(4) << fixed
                                << setw(9) << integration_time << "    |" 
                                << setw(9) << vecadd_time << "    |"
                                << setw(9) << waiting_time << "    | "
                                << setw(9) << writing_time << "  | "
                                << setw(9) << timestepvar_time << "\n";
    }
    pthread_mutex_unlock(&io_lock);
    barrier_wait(&WriterBarrier);
    pthread_mutex_lock(&io_lock);    
    if (id == num_threads) {
        cout << "[W] "  
                                << setw(9) << integration_time << "    |" 
                                << setw(9) << vecadd_time << "    |"
                                << setw(9) << waiting_time << "    | "
                                << setw(9) << writing_time << "  | "
                                << setw(9) << timestepvar_time << "\n";
    }
    pthread_mutex_unlock(&io_lock);
    barrier_wait(&WriterBarrier);
    pthread_mutex_lock(&io_lock);    
    if (id == 0) {
        cout << "Barrier Timing counts:\n"
             << "   | bstart |   b00   |   b01   |   b10   |   b11   |   b20   |   b21   |"
               << "   b30   |   b31   |   b40   |   b41   |   b50   |   b51   |  bdiff  |  bats  | \n"; 
    }
    pthread_mutex_unlock(&io_lock);
    barrier_wait(&WriterBarrier);
    pthread_mutex_lock(&io_lock);    
    cout << setprecision(2);
    if (id < num_threads) {
        cout << "[" << id << "] " << setw(7) << bstart_t << " | "
                                << setw(7) << b00_t << " | " << setw(7) << b01_t << " | " 
                                << setw(7) << b10_t << " | " << setw(7) << b11_t << " | " 
                                << setw(7) << b20_t << " | " << setw(7) << b21_t << " | " 
                                << setw(7) << b30_t << " | " << setw(7) << b31_t << " | " 
                                << setw(7) << b40_t << " | " << setw(7) << b41_t << " | " 
                                << setw(7) << b50_t << " | " << setw(7) << b51_t << " | " 
                                << setw(7) << bdiff_t << " | " << setw(7) << bats_t << "\n"; 
    }
    pthread_mutex_unlock(&io_lock);
    barrier_wait(&WriterBarrier);
    pthread_mutex_lock(&io_lock);    
    if (id == num_threads) {
        cout << "[W] " << setw(7) << bstart_t << " | "
                                << setw(7) << b00_t << " | " << setw(7) << b01_t << " | " 
                                << setw(7) << b10_t << " | " << setw(7) << b11_t << " | " 
                                << setw(7) << b20_t << " | " << setw(7) << b21_t << " | " 
                                << setw(7) << b30_t << " | " << setw(7) << b31_t << " | " 
                                << setw(7) << b40_t << " | " << setw(7) << b41_t << " | " 
                                << setw(7) << b50_t << " | " << setw(7) << b51_t << " | " 
                                << setw(7) << bdiff_t << " | " << setw(7) << bats_t << "\n"; 
    }
    pthread_mutex_unlock(&io_lock);
    barrier_wait(&WriterBarrier);
}

/////////////////////////////////////////////////////////////////////////
/////////////             INTEGRATOR FUNCTION               /////////////
/////////////////////////////////////////////////////////////////////////
void HierarchyIntegrator::rk4_integrate(int id, int num_threads){
    void (HierarchyIntegrator::*take_step)(const HierarchyNode &, Complex *,const int &, 
                                        Complex *, Complex *, Complex *, Complex *);


    if (!spectrum) {
        if (corrBath)
            take_step = &HierarchyIntegrator::integration_step_fullBath;
        else {
            if (!filter_tolerance) {
                take_step = &HierarchyIntegrator::integration_step;
            } 
            else
                take_step = &HierarchyIntegrator::integration_step_adaptTrunc;
            }
    }
    else{
        if (corrBath)
            take_step = &HierarchyIntegrator::vec_integration_step_fullBath;
        else 
            take_step = &HierarchyIntegrator::vec_integration_step; 
    }

    int N = t/dt;
    bool restart=false;
    if(restartStep > 0)
	restart=true;
    time_t t1,t2;

    int *m_start; 
    int *m_end;  
    if (id != num_threads){
        m_start = MatrixIndices_th[id];
        m_end = MatrixIndices_th[id]+MatrixCount_th[id];
    }
    else{
        m_start = new int;
        *m_start = 0;
        m_end = new int;
        *m_end = MatrixCount;
    }
    int m_cnt = 0;
    int Ns2_th = Ns2;
    if (spectrum) Ns2_th = Ns;
    int Ns2M_th = Ns2_th*MatrixCount_th[id];
    Complex *k1,*k2,*k3,*k4;
    Complex *Qrho_rhoQ;
    Qrho_rhoQ = new Complex[Ns2_th];
    Complex *Vbath_next = NULL;
    Complex *Vbath_prev = NULL;
    Complex *Vbath_same = NULL;
    
    //ADAPTIVE TRUNCATION VARIABLES
    bool AdapTrunc = false;
    if (filter_tolerance > 0) 
        AdapTrunc = true;

    Float Tnow = 0;

    ofstream maxRho;

    ofstream maxf;
    barrier_wait(&WriterBarrier);    
    if (id < num_threads){
	k1 = new Complex [Ns2M_th];
	k2 = new Complex [Ns2M_th];
	k3 = new Complex [Ns2M_th];
	k4 = new Complex [Ns2M_th];
     
	for(int i = 0; i < MatrixCount_th[id]; i++){
            rho[MatrixIndices_th[id][i]].k1 = k1+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k2 = k2+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k3 = k3+i*Ns2_th;
            rho[MatrixIndices_th[id][i]].k4 = k4+i*Ns2_th;
	    for(int j = 0; j < Ns2_th; j++){
                k1[i*Ns2_th+j] = Complex(0);
                k2[i*Ns2_th+j] = Complex(0);
                k3[i*Ns2_th+j] = Complex(0);
                k4[i*Ns2_th+j] = Complex(0);
            }
	}
        if(corrBath){
            Vbath_same = new Complex[Ns*Ns2_th];
            Vbath_next = new Complex[Ns*Ns2_th];
            Vbath_prev = new Complex[Ns*Kt*Ns2_th];
            for(int n = 0; n < Ns; n++){ //n is BChl index (aka coupling term index)
                if(!spectrum) {
                    for(int i=0; i < Ns; i++)
                        for(int j=0; j < Ns; j++){
                            Vbath_next[n*Ns*Ns+i*Ns+j] = Complex(0,-1*(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i]));
    
                            Vbath_same[n*Ns*Ns+i*Ns+j] = Complex(pow(Vbath_re[n*Ns+j]-Vbath_re[n*Ns+i],2))
                                                                *rho[id].Same_prefactor[n];
                        }
                    for(int k=0; k < Kt; k++)
                        for(int i=0; i < Ns; i++)
                            for(int j=0; j < Ns; j++)
                                Vbath_prev[n*Kt*Ns*Ns+k*Ns*Ns+i*Ns+j] = 
                                        Complex(Vbath_re[n*Ns+j])*C[n*Kt+k]*Complex(0,-1)
                                    - Complex(Vbath_re[n*Ns+i])*conj(C[n*Kt+k])*Complex(0,-1);
                }
                else {
                    for(int i=0; i < Ns; i++){
                            Vbath_next[n*Ns+i] = Complex(0,-Vbath_re[n*Ns+i]);
                            Vbath_same[n*Ns+i] = Complex(pow(Vbath_re[n*Ns+i],2))*rho[id].Same_prefactor[n];
                    }
                    for(int k=0; k < Kt; k++)
                        for(int i=0; i < Ns; i++)
                            Vbath_prev[n*Kt*Ns+k*Ns+i] = Complex(0,-Vbath_re[n*Ns+i])*C[n*Kt+k];

                }

            }
        }	
        if (verbose > 0) {
            pthread_mutex_lock(&io_lock);
          
            cout << "[" << id << "]:\tMemory for " << MatrixCount_th[id] << " density matrices" 
                << " (each with " <<  Ns2_th << " elements) assigned\n";
            pthread_mutex_unlock(&io_lock);
        }
	barrier_wait(&MatupdateBarrier);
        if (verbose > -1) {
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "]:\tIntegrating using RK4.\n"; 
            pthread_mutex_unlock(&io_lock);
        }
    }
    else{
        pthread_mutex_lock(&io_lock);
        if (verbose > 1)
            cout << "[" << id << "]:\tWriter Waiting.\n"; 

        //if (AdapTrunc) {
        //    maxRho.open("MaxRhoVals.txt");
        //    maxRho << "# time ";
        //    for(int m = 0; m < MatrixCount; ++m){
        //        for(int i = 0; i < MKt-1; i++)
        //                maxRho << rho[m].I[i] << ",";
        //        maxRho << rho[m].I[MKt-1] << " | ";
        //           
        //    }
        //    maxRho << endl;
        //}
        pthread_mutex_unlock(&io_lock);
    }
    barrier_wait(&WriterBarrier);
    time(&t1);
    int i = -1;
    while (Tnow < t){
        i += 1;
//    for(int i =0; i < N; i++){
	if(id < num_threads){
            if (TL_trunc) {
                //Iterate over system-bath coupling terms
                for (int m = id; m < M; m+=num_threads)
                    constructTL_BLAS_threaded(Tnow,m);
                barrier_wait(&MatupdateBarrier);
            }
	    /////////////////////////	
	    //          k1         //
	    /////////////////////////
	    for(int *m = m_start; m != m_end; ++m)
		(this->*take_step)(rho[*m], rho[*m].k1,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
	    //WAIT FOR k1 UPDATED
	    barrier_wait(&MatupdateBarrier);
            for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k1[j]*dtover2;
            }
	    barrier_wait(&MatupdateBarrier);
	    ///////////////////////

	    /////////////////////////	
	    //         k2          //
	    /////////////////////////
	    for(int *m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k2,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
	    //WAIT FOR k2 UPDATED
	    barrier_wait(&MatupdateBarrier);
            for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k2[j]*dtover2;
            }
	    barrier_wait(&MatupdateBarrier);
	    ///////////////////////////

	    //////////////////////////	
	    //          k3          //
	    //////////////////////////
	    for(int *m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k3,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
	    //WAIT FOR k3 UPDATED
	    barrier_wait(&MatupdateBarrier);
            for(int j = 0; j < Ns2M_th; ++j){
                rho_mid_matrices_th[id][j] = rho_matrices_th[id][j] + k3[j]*dt_c;
            }
	    barrier_wait(&MatupdateBarrier);
	    /////////////////////////	

	    /////////////////////////	
	    //          k4         //
	    /////////////////////////
	    for(int *m = m_start; m != m_end; ++m)
		(this->*take_step)(rho_mid[*m], rho[*m].k4,i,Qrho_rhoQ,Vbath_next,Vbath_prev,Vbath_same);
	    //WAIT FOR drho UPDATED
	    barrier_wait(&MatupdateBarrier);
            Maxddrho[*m_start] = 0; 
            m_cnt = 0;
            for(int j = 0; j < Ns2M_th; ++j){
                rho_matrices_th[id][j] += dtover6*(k1[j]+two*k2[j]+two*k3[j]+k4[j]);
            }

            if (AdapTrunc) {
                for(int * m = m_start; m != m_end; ++m){
                    MaxRhoVal[*m] = 0;
                    for (int j = 0; j < Ns2_th; ++j){
                        if ( abs(rho[*m].rho[j]) > MaxRhoVal[*m] )
                            MaxRhoVal[*m] = abs(rho[*m].rho[j]);
                    }
                    if ( MaxRhoVal[*m] < filter_tolerance && rho[*m].active ) {
                        rho[*m].active = false;
                        rho_mid[*m].active = false;

                    }
                    else if ( MaxRhoVal[*m] > filter_tolerance && rho[*m].active == false ) {
                        rho[*m].active = true;
                        rho_mid[*m].active = true;
                    }
                }
	    }
    	    /////////////////////////
	}
	//////////////////////////
	// WAIT FOR rho UPDATED //
	//////////////////////////
	barrier_wait(&WriterBarrier);
        Tnow += dt;
       	if(id == num_threads){	    
	    if(((int)(i*dt*1000))%(1) == 0){
                if (i==0 or i==99){
                    time(&t2);
                    outf << "#Time for 1st timestep: " << difftime(t2,t1)/(i+1) << " s.\n";
                    outf << "#Estimated time for 1 ps: " << difftime(t2,t1)/dt/60/60/(i+1) << " hrs.\n";
                    outf << "#Total running time: " << difftime(t2,t1)*N/60/60/24/(i+1) << " days.\n";
                }

                //if (AdapTrunc) {
                //    maxRho << Tnow;
                //    for (int m = 0; m < MatrixCount; m++){
                //        maxRho << " " << MaxRhoVal[m];
                //    }
                //    maxRho << endl;
                //}

                outf << Tnow << " ";
                for( int j = 0; j < Ns2_th; j++)
                        outf << rho[0].rho[j] << " ";  // NOTE WRITING OUTPUT IN COLUMN MAJOR ORDERING!!!!!!!
                outf << endl;
            }
	    if(restart && i%restartStep == 0 && i > 0){
		restartf.seekp(0, ios_base::beg);
		restartf << i*dt << endl;
		for(int m = 0; m < MatrixCount; m++){
		    for(int mi = 0; mi < Ns2_th; mi++)
			    restartf << rho[m].rho[mi] << "\n";
		}
	    }
	} 
    }
    int nActive = 0;
    for(int j = 0; j < MatrixCount; j++)
        if (rho[j].active)
            nActive++;
    if (AdapTrunc && id == num_threads && verbose > 0)
        cout << "Total number of active matrices = " << nActive << "/" << MatrixCount << endl;
    if (id == num_threads && verbose > -1 ) {
        cout << "Total number of integration steps = "<< i+1 << endl;
        time(&t2);
        cout << "Total integration time = " << difftime(t2,t1)/60 << " mins\n";
    }
    
    if(restart && id == 0){
	restartf.seekp(0, ios_base::beg);
	restartf << t << endl;
	for(int m = 0; m < MatrixCount; m++){
	    for(int mi = 0; mi < Ns; mi++)
		for (int mj=0; mj < Ns; mj++)
		    restartf << rho[m].rho[mi+Ns*mj] << "\n";
	}
    }

    //if (id == num_threads)
        //maxRho.close();

    //pthread_exit(0);

};

//////////////////////////////////////////
//////// VECTOR INTEGRATION STEP ///////// CORRELATED BATH
//////////////////////////////////////////
void HierarchyIntegrator::vec_integration_step_fullBath(const HierarchyNode &sigma, Complex *drho,
                                               const int &time_step, Complex *Qrho, 
                                               Complex * Vnext, Complex * Vprev, Complex * Vsame){
    int j = 0;
    int nj = 0;
    int m = 0;
    Complex * drho_col_elem, *rho_col_elem, *rho_row_elem;

    phi_hemv(Ns,&L_prefactor,sigma.H,sigma.rho,&zero,drho);
    //phi_hemv(CblasColMajor,CblasUpper,Ns,&L_prefactor,sigma.H,Ns,sigma.rho,1,&zero,drho,1);


//    cout << "NZI\n";
    rho_col_elem = sigma.rho;
    drho_col_elem = drho;
    for(nj = 0; nj < Ns;nj++, drho_col_elem++, rho_col_elem++){
    	*drho_col_elem -= sigma.NZI_prefactor*(*rho_col_elem); 
    }


    //SAME
//    cout << "Temp Corr\n";
    rho_col_elem = sigma.rho;
    for(m = 0; m < Ns; m++){
        rho_row_elem = &(Vsame[m*Ns]);
        for(nj = 0; nj < Ns; nj++){
    	    drho[nj] -= rho_row_elem[nj]*rho_col_elem[nj];
	}
    }


    //NEXT
    if(sigma.TL_truncation){
        //- sum_m [ V_m Q_m rho]
        // Note Qrho_rhoQ includes no -i, and V_m includes -i*-i (different to
        // no TL_truncation)
	for(nj = 0; nj < Ns; nj++){
	    phi_gemv(Ns,&neg_imag_const,Q_tl[nj],sigma.rho,&zero,Qrho);
	    //phi_gemv(CblasColMajor,CblasNoTrans,Ns,Ns,&neg_imag_const,Q_tl[nj],Ns,sigma.rho,1,&zero,Qrho,1);
            rho_row_elem = &(Vnext[nj*Ns]);
            for(m = 0; m < Ns; m++)
                drho[m] += rho_row_elem[m]*Qrho[m];
	}

    }
    else{
       	for(j = 0; j < sigma.Nnext; j++){ //Nnext runs from 0 to NsKt
	    nj = sigma.NextIndMatrix[j];//runs from 0 to Ns
	    rho_col_elem = sigma.Next[j]->rho;
            rho_row_elem = &(Vnext[nj*Ns]);
	    for(m=0; m < Ns; m++){
                drho[m] += rho_row_elem[m]*rho_col_elem[m];
	    }
	}

    }

    //PREV
    for(j = 0; j < sigma.Nprev; j++){
        //nj = sigma.PrevIndMatrix[j]; //runs from 0 to Ns
        nj = sigma.PrevInd[j]; //runs from 0 to Ns*Kt-1
        rho_col_elem = sigma.Prev[j]->rho;
        rho_row_elem = &(Vprev[nj*Ns]);
        for(m=0; m < Ns; m++){
            drho[m] += Complex(sigma.I[nj])*rho_row_elem[m]*rho_col_elem[m];
        }
    }
};


////////////////////////////////////////// 
//////// VECTOR INTEGRATION STEP /////////
//////////////////////////////////////////
void HierarchyIntegrator::vec_integration_step(const HierarchyNode &sigma, Complex *drho,const int &time_step, Complex *Qrho, Complex * Vnext, Complex * Vprev, Complex * Vsame){
    int nj = 0;
    int m = 0;
    Complex * drho_col_elem, *rho_col_elem;

    phi_hemv(Ns,&L_prefactor,sigma.H,sigma.rho,&zero,drho);
    /*for(nj = 0; nj < Ns; ++nj){
        drho[nj] = 0;
        for (m = 0; m < Ns; ++m)
            drho[nj] += sigma.H[nj*Ns+m];
        drho[nj]*=L_prefactor;
    }*/
//    cout << "NZI\n";
    rho_col_elem = sigma.rho;
    drho_col_elem = drho;
    for(nj = 0; nj < Ns;nj++, drho_col_elem++, rho_col_elem++){
    	*drho_col_elem -= sigma.NZI_prefactor*(*rho_col_elem); 
    }


//    cout << "Temp Corr\n";
    drho_col_elem = drho;
    rho_col_elem = sigma.rho;
    if (!multiBath)
        for(nj = 0; nj < Ns; nj++,rho_col_elem++,drho_col_elem++){
            *drho_col_elem -= sigma.Same_prefactor[nj]*(*rho_col_elem);	
        }
    else
        for(nj = 0; nj < M; nj++){
            m = MDiagInd[nj];
            drho[m] -= sigma.Same_prefactor[nj]*sigma.rho[m];	
        }
        
//    cout << "Next\n";
    if(sigma.TL_truncation){
        if (!multiBath) {
            if (!rho0_restart) {
                for(nj = 0; nj < Ns; nj++){
                    phi_gemv(Ns,&one,Q_tl[nj],sigma.rho,&zero,Qrho);
                    //phi_gemv(CblasColMajor,CblasNoTrans,Ns,Ns,&one,Q_tl[nj],Ns,sigma.rho,1,&zero,Qrho,1);
                    drho[nj] -= Qrho[nj];
                }
            }
            else {
                for(nj = 0; nj < Ns; nj++){
                    phi_gemv(Ns,&one,Q_tldag[nj],sigma.rho,&zero,Qrho);
                    //phi_gemv(CblasColMajor,CblasTrans,Ns,Ns,&one,Q_tldag[nj],Ns,sigma.rho,1,&zero,Qrho,1);
                    drho[nj] -= Qrho[nj];
                }
            }
        }
        else {
            if (!rho0_restart) {
                for(nj = 0; nj < M; nj++){
                    m = MDiagInd[nj];
                    phi_gemv(Ns,&one,Q_tl[nj],sigma.rho,&zero,Qrho);
                    //phi_gemv(CblasColMajor,CblasNoTrans,Ns,Ns,&one,Q_tl[nj],Ns,sigma.rho,1,&zero,Qrho,1);
                    drho[m] -= Qrho[m];
                }
            }
            else {
               for(nj = 0; nj < M; nj++){
                    m = MDiagInd[nj];
                    phi_gemv(Ns,&one,Q_tldag[nj],sigma.rho,&zero,Qrho);
                    //phi_gemv(CblasColMajor,CblasTrans,Ns,Ns,&one,Q_tldag[nj],Ns,sigma.rho,1,&zero,Qrho,1);
                    drho[m] -= Qrho[m];
                }
             }
        }
    }
    else{
	for(nj = 0; nj < sigma.Nnext; nj++){ // this is sum over M & Kt
            if (sigma.Next[nj]->active) {
                m = sigma.NextIndMatrix[nj];
                drho[m] -= (sigma.Next_prefactor[nj])*sigma.Next[nj]->rho[m];
            }
	}
    }

    for(nj = 0; nj < sigma.Nprev; nj++){
        if (sigma.Prev[nj]->active) {
            m = sigma.PrevIndMatrix[nj];
            drho[m] -= sigma.Prev_prefactor_row[nj]*sigma.Prev[nj]->rho[m];
        }        
    }
};



//
// Treats rho[0].rho and appropriate drho as having only Ns*Ns-1 entries,
// starting from index 1.
// 
//
//
void HierarchyIntegrator::minimize_step(const HierarchyNode &sigma, Complex * drho,const int &time_step, Complex * Qrho_rhoQ, Complex * Vnext, Complex * Vprev, Complex * Vsame){
    int nj;
    int m;
    int j;
    Complex * drho_col_elem, *drho_row_elem, * rho_row_elem, *rho_col_elem;
//    if ( sigma.active ){ 
        if (sigma.id != 0) {
            phi_symm(CblasLeft, Ns,&L_prefactor,   H,   sigma.rho,&zero,drho);
            phi_symm(CblasRight,Ns,&negL_prefactor,Hdag,sigma.rho,&one, drho);
            //phi_symm(CblasRowMajor,CblasLeft,CblasUpper,Ns,Ns,&L_prefactor,H,Ns,sigma.rho,Ns,&zero,drho,Ns);
            //phi_symm(CblasRowMajor,CblasRight,CblasUpper,Ns,Ns,&negL_prefactor,Hdag,Ns,sigma.rho,Ns,&one,drho,Ns);


            rho_col_elem = sigma.rho;
            drho_col_elem = drho;
            for(nj = 0; nj < Ns2;nj++, drho_col_elem++, rho_col_elem++){
                *drho_col_elem -= sigma.NZI_prefactor*(*rho_col_elem); 
            }

            if (!multiBath) {
                for(m = 0; m < Ns; m++){
                    drho_row_elem = &(drho[m*Ns]);
                    rho_row_elem = &(sigma.rho[m*Ns]);
                    for(nj = 0; nj < Ns; nj++){
                        drho_row_elem[nj] -= (sigma.Same_prefactor[nj]+sigma.Same_prefactor[m])*rho_row_elem[nj];
                    }
                    nj = (Ns+1)*m;
                    drho[nj] += two*sigma.Same_prefactor[m]*sigma.rho[nj];
                }
            } 
            else {
              for(j = 0; j < M; j++){
                nj = MDiagInd[j];
                drho_row_elem = &(drho[nj]);
                drho_col_elem = &(drho[Ns*nj]);
                rho_row_elem = &(sigma.rho[nj]);
                rho_col_elem = &(sigma.rho[Ns*nj]);
                for(m = 0; m < Ns; m++, rho_row_elem+=Ns, rho_col_elem++, drho_row_elem+=Ns, drho_col_elem++){
                    if (m != nj ) {
                      *drho_row_elem -= sigma.Same_prefactor[j]*(*rho_row_elem);
                      *drho_col_elem -= sigma.Same_prefactor[j]*(*rho_col_elem);	
                    }
                }
              }
            }
        }
        else {
            Complex dot1, dot2, dot3;
            int Nmj = Ns-1;
            m = 0;
            // j = 0 case: //First Ns Rows of SameLiouville
            rho_row_elem = sigma.SameLiouville+Ns*Ns;
            drho_row_elem = drho+1;
            for (nj = 1; nj < Ns; ++nj){
                rho_col_elem = sigma.rho;
                phi_dotu_sub(Ns,rho_row_elem,1,rho_col_elem,1,&dot2);
                rho_row_elem += Ns;
                rho_col_elem += Ns;
                phi_dotu_sub(Nmj,rho_row_elem+nj,Ns,rho_col_elem+nj,Ns,&dot3);
                rho_row_elem += Nmj*Ns;
//                rho_row_elem++;
                *drho_row_elem = dot2+dot3;
                drho_row_elem++;
            }
            --Nmj;
 
            // j = 1 -> Ns-1 : //Next Ns*(Ns-2) Rows of SameLiouville
            for (j = 1; j < Ns-1; ++j){
                for (nj = 0; nj < Ns; ++nj){
                    rho_col_elem = sigma.rho;
                    phi_dotu_sub(j,rho_row_elem+nj,Ns,rho_col_elem+nj,Ns,&dot1);
                    rho_row_elem += j*Ns;
                    rho_col_elem += j*Ns;
                    phi_dotu_sub(Ns,rho_row_elem,1,rho_col_elem,1,&dot2);
                    rho_row_elem += Ns;
                    rho_col_elem += Ns;  //Nmj = Ns-1-j
                    phi_dotu_sub(Nmj,rho_row_elem+nj,Ns,rho_col_elem+nj,Ns,&dot3);
                    rho_row_elem += Nmj*Ns;
                    *drho_row_elem = dot1+dot2+dot3;
                    drho_row_elem++;
                    
                }
                --Nmj;    
            }

            // Final Ns Rows of SameLiouville
            j = Ns-1; 
            for (nj = 0; nj < Ns; ++nj){
                rho_col_elem = sigma.rho;
                phi_dotu_sub(j,rho_row_elem+nj,Ns,rho_col_elem+nj,Ns,&dot1);
                rho_row_elem += j*Ns;
                rho_col_elem += j*Ns;
                phi_dotu_sub(Ns,rho_row_elem,1,rho_col_elem,1,&dot2);
                rho_row_elem += Ns;
                rho_col_elem += Ns;  //Nmj = Ns-1-j
                *drho_row_elem = dot1+dot2;
                drho_row_elem++;
                
            }
        }
    
//    }
    /*
    if (sigma.id == 3) {
        cout << "!2: " << endl;
        printMatrix(drho,Ns);
    }
    */
    ////// NEXT HIERARCHY MATRICES
    for(j = 0; j < sigma.Nnext; j++){ // this is sum over Ns & Kt
//        if (sigma.Next[j]->active) {
            nj = sigma.NextIndMatrix[j];
            drho_col_elem = &(drho[Ns*nj]);
            rho_col_elem = &(sigma.Next[j]->rho[Ns*nj]);
            for(m=0; m < Ns; m++){
                drho_col_elem[m] -= (sigma.Next_prefactor[j])*rho_col_elem[m];	
            }
//        }
    }
    for(m = 0; m < Ns; m++){ 
        drho_row_elem = &(drho[Ns*m]);
        for(j = 0; j < sigma.Nnext; j++){
//            if (sigma.Next[j]->active) {
                nj = sigma.NextIndMatrix[j];
                rho_row_elem = &(sigma.Next[j]->rho[Ns*m]);
                drho_row_elem[nj] += sigma.Next_prefactor[j]*rho_row_elem[nj];
 //           }
        }
    }
    /*
    if (sigma.id == 3) {
        cout << "!3: " << endl;
        printMatrix(drho,Ns);
    }
    */
    ////// PREV HIERARCHY MATRICES
    for(j = 0; j < sigma.Nprev; j++){
        if (sigma.Prev[j]->id != 0 || sigma.PrevInd[j] != 0) {
            rho_col_elem = sigma.Prev[j]->rho;
            m = sigma.PrevIndMatrix[j];
            nj = m;
            for (int i = 0; i < m; ++i) {
               drho[nj] += sigma.Prev_prefactor_col[j]*rho_col_elem[nj];
               nj += Ns;
            }

            nj -= m;
            for (int i = 0; i < m; ++i) {
               drho[nj] -= sigma.Prev_prefactor_row[j]*rho_col_elem[nj];
                ++nj;
            }
            drho[nj] += (sigma.Prev_prefactor_col[j]-sigma.Prev_prefactor_row[j])*rho_col_elem[nj];
            ++nj;
            for (int i = m+1; i < Ns; ++i) {
               drho[nj] -= sigma.Prev_prefactor_row[j]*rho_col_elem[nj];
               ++nj;
            }

            nj += m;
            for (int i = m+1; i < Ns; ++i) {
               drho[nj] += sigma.Prev_prefactor_col[j]*rho_col_elem[nj];
               nj += Ns;
            }
        }
        else {
                //rho_col_elem = sigma.Prev[j]->rho+Ns+1; 
                //for (nj = 1; nj < Ns; ++nj, rho_col_elem+=Ns+1){
                //    drho[0] -= sigma.PrevLiouville[j][0]*(*rho_col_elem);
                //}
                for (nj = 1; nj < 2*Ns-1; ++nj) {
                    m = sigma.PrevIndices[j][nj];
                    drho[m] += sigma.PrevLiouville[j][nj]*sigma.Prev[j]->rho[m];
                    
                } 
        }
    }
    /*
    if (sigma.id == 3) {
        cout << "!4: " << endl;
        printMatrix(drho,Ns);
    }
    */
    //if (sigma.id == 0)
    //    drho[0] = 0;
};


/////////////////////////////////////////////
//             STEADY STATE                //
// Bi-Conjugate Gradient Stabilized method //
//          including l-GMRES steps        //
/////////////////////////////////////////////
void HierarchyIntegrator::bicgstab_l_steadystate(int id, int num_threads){
    int m_start = id;
    time_t t1,t2;
    time(&t1);
    int l = bicgstabEll;
    int lp1 = l+1;
    Complex *b;
    b = new Complex[Ns*Ns];


    if (id == 0) {
        // Calculate target vector:
        b = new Complex[Ns*Ns];
        for (int j = 0; j < Ns*Ns; ++j)
            b[j] = -rho[0].SameLiouville[j*Ns*Ns];

        // Set top row to zero:
        for (int j = 0; j < Ns*Ns; ++j)
            rho[0].SameLiouville[j] = 0;
    }

    Complex result;
    Complex rho0, rho1_th, omegaBtm_th, omega, beta, alpha_th;
    Complex *tau;
    Complex sigma;
    Complex *gamma;
    Complex *gammap;
    Complex *gammapp;
    
    tau = new Complex[lp1*lp1];
    gamma = new Complex[lp1];
    gammap = new Complex[lp1];
    gammapp = new Complex[lp1];
    for (int i = 0; i < lp1*lp1; ++i)
        tau[i] = 0;
    for (int i = 0; i < lp1; ++i){
        gamma[i] = 0;
        gammap[i] = 0;
        gammapp[i] = 0;

    }

    rho0 = 1;
    alpha_th = 0;
    omega = 1;
    result = 0;

    int MaxSteps = 999999;
    Float diff_th = 1e30;
    Float eps = bicgstab_tolerance;
    int step = 0;

    if (verbose > 2 && id == num_threads){
        pthread_mutex_lock(&io_lock);   
        cout << "Getting initial residuals.\n";
        pthread_mutex_unlock(&io_lock);
    }

    barrier_wait(&WriterBarrier);
    //for (int l = 0; l < Lmax; ++l) {

    // Calculate initial residuals
    if (id < num_threads){ 
        for(int m = m_start; m < MatrixCount; m+=num_threads){
            rho[m].rho = rho[m].SteadyStateRho;
            for (int i = 0; i < Ns2; ++i)
                rho[m].rho[i] = 0;
        }
        if (id == 0){
            rho[0].rho[0] = 0;
        }
        barrier_wait(&MatupdateBarrier);
        for(int m = m_start; m < MatrixCount; m+=num_threads)
            minimize_step(rho[m], rho[m].r,step, NULL, NULL, NULL, NULL);    
        diff[id] = 0;
        for(int m = m_start; m < MatrixCount; m+=num_threads) {
            if (m >= 1+L[1]) {
                for (int i = 0; i < Ns*Ns; ++i){
                    rho[m].r[i] = -rho[m].r[i];
                    rho[m].r0[i] = rho[m].r[i];
                }
            } 
            else if (m > 0){
                rho[m].r[0] = -rho[m].r[0];
                if (rho[m].PrevInd[0] == 0)
                    rho[m].r[0] -= rho[m].PrevLiouville[0][0];
                rho[m].r0[0] = rho[m].r[0];
                for (int i = 1; i < Ns*Ns; ++i){
                    rho[m].r[i] = -rho[m].r[i];
                    rho[m].r0[i] = rho[m].r[i];
                } 
            }
            else {
                for (int i = 0; i < Ns*Ns; ++i){
                    rho[m].r[i] = b[i]-rho[m].r[i];
                    rho[m].r0[i] = rho[m].r[i];
                }            
                rho[m].r[0] = 0;
                rho[m].r0[0] = 0;       
            }
            phi_dotc_sub(Ns2,rho[m].r,1,rho[m].r,1,&result);
            diff[id] += real(result);
        
        }
        barrier_wait(&MatupdateBarrier);   
        diff_th = 0; 
        for(int m = 0; m < num_threads; ++m){
                diff_th += diff[m];
        } 
        diff_th = sqrt(diff_th);
        result = Complex(1./diff_th);
        for (int m = m_start; m < MatrixCount; m+=num_threads) {
            phi_scal(Ns2,&result,rho[m].r0,1);
        }
    }
    if (verbose > 1) {
        if (id < num_threads) {
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "]  Calculating steady state density matrix using BiCGSTAB(" << l << ").\n";
            pthread_mutex_unlock(&io_lock);
        }
    }
    else if (verbose > 0 && id == num_threads){
        pthread_mutex_lock(&io_lock);
        cout << "Calculating steady state density matrix using  BiCGSTAB(" << l << ").\n";
        pthread_mutex_unlock(&io_lock);
    }

    barrier_wait(&WriterBarrier);    
    int n = 0;
    //MaxSteps = 6;
    while (diff_th > eps && step < MaxSteps) {
        if ( id < num_threads ) {
            // Take l BiCG steps
            rho0 = -omega*rho0;
            //cout << step << " Doing BiCG Calc:\n";
            for(int j = 0; j < l; ++j){
                rho1[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r0,1,rho[m].r+j*Ns2,1,&result);
                    rho1[id] += result;
                }

                ///////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                rho1_th = 0;
                for(int m = 0; m < num_threads; ++m){
                        rho1_th += rho1[m];
                }
                beta = -alpha_th*rho1_th/rho0; 
                //////////////////////////////////////////
                rho0 = rho1_th;

                for(int m = m_start; m < MatrixCount; m+=num_threads){
                    phi_scal(Ns2*(j+1),&beta,rho[m].v,1);
                    phi_axpy(Ns2*(j+1),&one,rho[m].r,1,rho[m].v,1);
                    rho[m].rho = rho[m].v+j*Ns2;
                }

                ///////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                for(int m = m_start; m < MatrixCount; m += num_threads)                
                    minimize_step(rho[m], rho[m].v+(j+1)*Ns2,step,NULL,NULL,NULL,NULL);
                if (id == 0)
                    rho[0].v[(j+1)*Ns2] = 0;
                ///////////////////////////////////////////

                alpha[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r0,1,rho[m].v+(j+1)*Ns2,1,&result);
                    alpha[id] += result;            
                }
                /////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                alpha_th = 0;
                for(int m = 0; m < num_threads; ++m){
                        alpha_th += alpha[m];
                }
                alpha_th = -rho1_th/alpha_th; 
                ///////////////////////////////////////// 
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2*(j+1),&alpha_th,rho[m].v+Ns2,1,rho[m].r,1);
                    rho[m].rho = rho[m].r+j*Ns2;
                }                

                alpha_th = -alpha_th;
                barrier_wait(&MatupdateBarrier);
                for(int m = m_start; m < MatrixCount; m+=num_threads) 
                    minimize_step(rho[m],rho[m].r+(j+1)*Ns2,step,NULL,NULL,NULL,NULL);
                if (id == 0)
                   rho[0].r[(j+1)*Ns2] = 0;

                for(int m = m_start; m < MatrixCount; m+=num_threads) 
                    phi_axpy(Ns2,&alpha_th,rho[m].v,1,rho[m].SteadyStateRho,1);                    
                //cout << step << " alpha= " << alpha_th << endl;

            }
            //// GMRES SECTION
            for(int j = 1; j <= l; ++j){
                for(int i = 1; i < j; ++i){
                    alpha[id] = 0;
                    for(int m = m_start; m < MatrixCount; m+=num_threads) {
                        phi_dotc_sub(Ns2,rho[m].r+j*Ns2,1,rho[m].r+i*Ns2,1,&result);
                        alpha[id] += result;            
                    }
                    /////////////////////////////////////////
                    barrier_wait(&MatupdateBarrier);
                    n = i*lp1+j;
                    tau[n] = 0;
                    for(int m = 0; m < num_threads; ++m){
                        tau[n] += alpha[m];
                    }
                    tau[n] /= -gamma[i];
                    /////////////////////////////////////////
                    for(int m = m_start; m < MatrixCount; m += num_threads)
                        phi_axpy(Ns2,tau+n,rho[m].r+i*Ns2,1,rho[m].r+j*Ns2,1);                    
                    tau[n] *= -1;
                }

                omegaTop[id] = 0;
                omegaBtm[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r+j*Ns2,1,rho[m].r+j*Ns2,1,&result);
                    omegaBtm[id] += result;            
                    phi_dotc_sub(Ns2,rho[m].r+j*Ns2,1,rho[m].r,1,&result);
                    omegaTop[id] += result;      
                }
                /////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                gamma[j] = 0;
                gammap[j] = 0;
                for(int m = 0; m < num_threads; ++m){
                    gamma[j] += omegaBtm[m];
                    gammap[j] += omegaTop[m];
                }
                /////////////////////////////////////////
                gammap[j] /= gamma[j];
                //cout << step << " " << j << " " << gammap[j] << endl;

            }


            gamma[l] = gammap[l];
            omega = gammap[l];
            for(int j = l-1; j > 0; --j){
                sigma = 0;
                for(int i = j+1; i <= l; ++i){
                    sigma += tau[j*lp1+i]*gamma[i];
                    //cout << step << " " << j << " " << i << " " << tau[j*lp1+i] << endl;
                }
                gamma[j] = gammap[j]- sigma;
                //cout << step << " " << j << " " << gamma[j] << endl;
            }

            for(int j = 1; j < l; ++j){
                sigma = 0;
                for(int i  = j+1; i < l; ++i)
                    sigma += tau[j*lp1+i]*gamma[i+1];
                gammapp[j] = gamma[j+1]+sigma;
            }

            for (int m = m_start; m < MatrixCount; m+=num_threads) {
               phi_axpy(Ns2,gamma+1,rho[m].r,1,rho[m].SteadyStateRho,1);
               for(int j = 1; j < l; ++j)
                   phi_axpy(Ns2,gammapp+j,rho[m].r+j*Ns2,1,rho[m].SteadyStateRho,1);     
            }

            for (int j = 0; j <= l; ++j){
                gamma[j] = -gamma[j];
                gammap[j] = -gammap[j];
            }
            
            /// Update residuals and search directions
            for (int m = m_start; m < MatrixCount; m+=num_threads){
                for (int j = l; j > 0; --j)
                   phi_axpy(Ns2,gammap+j,rho[m].r+j*Ns2,1,rho[m].r,1);
                for (int j = l; j > 0; --j)
                   phi_axpy(Ns2,gamma+j,rho[m].v+j*Ns2,1,rho[m].v,1);
            }

            // Recalculate ||r||
            ////////////////////////////////////////
            barrier_wait(&MatupdateBarrier);
            diff[id] = 0;
            for (int m = m_start; m < MatrixCount; m+=num_threads) {
                phi_dotc_sub(Ns2,rho[m].r,1,rho[m].r,1,&result);
                 diff[id] += real(result);
            }      
            //////////////////////////////////////
            

        }
        else{
            diff[id] = 0;
        }
        
        /////////////////////////////////////
        barrier_wait(&WriterBarrier);   
        diff_th = 0; 
        for(int m = 0; m < num_threads; ++m){
           //if (diff_th < diff[m])
                diff_th += diff[m];
        } 
        diff_th = sqrt(diff_th);
        /////////////////////////////////////
        if (id == num_threads) {
            //cout << step << " diff= " << diff_th << endl;
            if (step%1 == 0) {
                //rho[0].SteadyStateRho[0] = 1;
                result = 1;
                for(int i = 1; i < Ns; ++i)
                    result += rho[0].SteadyStateRho[(1+Ns)*i];
                outf << diff_th << " " << 1./real(result);
                for (int i = 0; i < Ns*Ns; ++i)
                    outf << rho[0].SteadyStateRho[i]/real(result) << " ";
                outf << endl;
            }
        }
        step += l;

    }// end while loop

    //}
    barrier_wait(&WriterBarrier);   
    time(&t2);
    // Finished, now write out the results
    if (id == num_threads) {
        pthread_mutex_lock(&io_lock);
        if (step == MaxSteps) {
            cout << "Convergence failed! (steps = " << step << ", diff = " << diff_th << ")\n";
            outf << "#Info: Total steps = " << step << ". No convergence.\n";
        }
        else{
            cout << "Converged to diff = " << diff_th << " in " << step << " steps.\n";
            outf << "#Info: Converged to diff = " << diff_th << " in " << step << " steps.\n";
        }

        if (verbose > 0)
            cout << "Writing final steady-state density matrix.\n";
    
        outf << "#Info: Total running time: " << difftime(t2,t1)/60/60/24 << " days.\n";
        rho[0].SteadyStateRho[0] = 1;
        result = 1;
        for(int i = 1; i < Ns; ++i)
            result += rho[0].SteadyStateRho[(1+Ns)*i];
        for (int i = 0; i < Ns*Ns; ++i)
            rho[0].SteadyStateRho[i] /= result;
        outf << diff_th << " ";
        for (int i = 0; i < Ns*Ns; ++i){
                outf << rho[0].SteadyStateRho[i] << " ";
        }
        outf << endl;

        if (verbose > 1)
            printMatrix(rho[0].SteadyStateRho,Ns);
        if (verbose > 2)
            for (int m = 1; m < MatrixCount; ++m) {
                cout << m << ":\n";
                printMatrix(rho[m].SteadyStateRho,Ns);   
            }
        cout << "Total running time: " <<  difftime(t2,t1)/60/60 << " hours.\n";
        if(restartStep > 0){
            restartf << step << endl;
            for(int m = 0; m < MatrixCount; m++){
                for(int mi = 0; mi < Ns2; mi++)
                        restartf << rho[m].SteadyStateRho[mi] << "\n";
            }
        }   
        pthread_mutex_unlock(&io_lock);
    }
    barrier_wait(&WriterBarrier);   
    //pthread_exit(0);

};

/////////////////////////////////////////////
//             STEADY STATE                //
// Bi-Conjugate Gradient Stabilized method //
/////////////////////////////////////////////
void HierarchyIntegrator::bicgstab_steadystate(int id, int num_threads){
    int m_start = id;
    time_t t1,t2;
    time(&t1);

    Complex *b;
    b = new Complex[Ns*Ns];

    if (id == 0) {
        // Correction for improved stability
        b = new Complex[Ns*Ns];
        for (int j = 0; j < Ns*Ns; ++j){
            b[j] = -rho[0].SameLiouville[j*Ns*Ns];
            rho[0].SameLiouville[j*Ns*Ns] = 0;
        }
        // Set top row to zero:
        for (int j = 0; j < Ns*Ns; ++j)
            rho[0].SameLiouville[j] = 0;

    }

    Complex result;
    Complex rho0, rho1_th, omegaBtm_th, omega, beta, alpha_th;
    Complex omegaInv;
    rho0 = rho1_th = alpha_th = beta = 1;
    omega = -1;
    result = 0;
    int MaxSteps = 9999;
    Float diff_th = 1e30;
    Complex diff_th_c = 0;
    Float eps = bicgstab_tolerance;
    int step = 0;

    if (verbose > 2 && id == num_threads){
        pthread_mutex_lock(&io_lock);   
        cout << "Getting initial residuals.\n";
        pthread_mutex_unlock(&io_lock);
    }

    barrier_wait(&WriterBarrier);
    //for (int l = 0; l < Lmax; ++l) {
        // Calculate initial residuals
        if (id < num_threads){ 
            for(int m = m_start; m < MatrixCount; m+=num_threads)
                rho[m].rho = rho[m].SteadyStateRho;
            if (id == 0){
                rho[0].rho[0] = 0;
            }
            barrier_wait(&MatupdateBarrier);
            //for(int m = m_start; m < MatrixCount; m+=num_threads)
            //    minimize_step(rho[m], rho[m].r,step, NULL, NULL, NULL, NULL);    
            diff[id] = 0;
            for(int m = m_start; m < MatrixCount; m+=num_threads) {
                if (m >= 1+L[1]) {
                    for (int i = 0; i < Ns*Ns; ++i){
                        rho[m].r[i] = -rho[m].r[i];
                        rho[m].r0[i] = rho[m].r[i];
                    }
                } 
                else if (m > 0){
                    rho[m].r[0] = -rho[m].r[0];
                    if (rho[m].PrevInd[0] == 0)
                        rho[m].r[0] -= rho[m].PrevLiouville[0][0];
                    rho[m].r0[0] = rho[m].r[0];
                    for (int i = 1; i < Ns*Ns; ++i){
                        rho[m].r[i] = -rho[m].r[i];
                        rho[m].r0[i] = rho[m].r[i];
                    } 
                }
                else {
                    for (int i = 0; i < Ns*Ns; ++i){
                        rho[m].r[i] = b[i]-rho[m].r[i];
                        rho[m].r0[i] = rho[m].r[i];
                    }            
                    rho[m].r[0] = 0;
                    rho[m].r0[0] = 0;        
                }
                phi_dotc_sub(Ns2,rho[m].r,1,rho[m].r,1,&result);
                diff[id] += real(result);
            }

            barrier_wait(&MatupdateBarrier);   
            diff_th = 0; 
            for(int m = 0; m < num_threads; ++m){
                if (diff_th < diff[m])
                    diff_th = diff[m];
            } 
            diff_th_c = Complex(1./diff_th);
            for (int m = m_start; m < MatrixCount; m+=num_threads) {
                phi_scal(Ns2,&diff_th_c,rho[m].r0,1);
            }
        
        }
        if (verbose > 1) {
            if(id < num_threads) {
                pthread_mutex_lock(&io_lock);
                cout << "[" << id << "]  Calculating steady state density matrix using BiCGSTAB.\n";
                pthread_mutex_unlock(&io_lock);
            }
        }
        else if (verbose > 0 && id == num_threads) {
            pthread_mutex_lock(&io_lock);
            cout << "Calculating steady state density matrix using BiCGSTAB.\n";
            pthread_mutex_unlock(&io_lock);
        }

        barrier_wait(&WriterBarrier);    

        while (diff_th > eps && step < MaxSteps) {
            ++step;
            if ( id < num_threads ) {
                rho1[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r0,1,rho[m].r,1,&result);
                    rho1[id] += result;
                }

                ///////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                rho1_th = 0;
                for(int m = 0; m < num_threads; ++m){
                        rho1_th += rho1[m];
                }
                beta = -rho1_th/rho0*alpha_th; 
                rho0 = rho1_th;
                //////////////////////////////////////////


                omegaInv = 1./real(omega);
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2,&omegaInv,rho[m].p,1,rho[m].v,1); // v <- v-p/omega
                    phi_copy(Ns2,rho[m].r,1,rho[m].p,1);           // p <- r
                    phi_axpy(Ns2,&beta,rho[m].v,1,rho[m].p,1);     // p <- p-beta*v
                    rho[m].rho = rho[m].p;
                }

                ///////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                for(int m = m_start; m < MatrixCount; m+=num_threads)
                    minimize_step(rho[m], rho[m].v,step, NULL, NULL, NULL, NULL); // v <- A*p
                if (id == 0)
                    rho[0].v[0] = 0;
                ///////////////////////////////////////////

                alpha[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r0,1,rho[m].v,1,&result);
                    alpha[id] += result;            
                }

                /////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                alpha_th = 0;
                for(int m = 0; m < num_threads; ++m){
                        alpha_th += alpha[m];
                }
                alpha_th = -rho1_th/alpha_th; 
                /////////////////////////////////////////

                for (int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2,&alpha_th,rho[m].v,1,rho[m].r,1); // r <- r-alpha*v
                    rho[m].rho = rho[m].r;                              // ready A*r
                }
                alpha_th = -alpha_th; 

                ////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                for(int m = m_start; m < MatrixCount; m+=num_threads)
                    minimize_step(rho[m], rho[m].t,step, NULL, NULL, NULL, NULL); // t <- A*r
                if (id == 0)
                    rho[0].t[0] = 0;
                //////////////////////////////////////

                omegaBtm[id] = 0;
                omegaTop[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads){
                   phi_dotc_sub(Ns2,rho[m].t,1,rho[m].r,1,&result); // dot(t,r)
                   omegaTop[id] += result;
                   phi_dotc_sub(Ns2,rho[m].t,1,rho[m].t,1,&result); // dot(t,t)
                   omegaBtm[id] += result;
                }

                /////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                omegaBtm_th = 0;
                omega = 0;
                for(int m = 0; m < num_threads; ++m){
                        omega += omegaTop[m];
                        omegaBtm_th += omegaBtm[m];
                }
                omega /= omegaBtm_th;  
                /////////////////////////////////////////

                // Update rho(\infty)
                for (int m = m_start; m < MatrixCount; m+=num_threads) {
                   phi_axpy(Ns2,&alpha_th,rho[m].p,1,rho[m].SteadyStateRho,1);  //rho <- rho + alpha*p
                   phi_axpy(Ns2,&omega,rho[m].r,1,rho[m].SteadyStateRho,1);      //rho <- rho + omega*s
                }

                omega = -omega;

                // Recalculate residuals            
                ////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                diff[id] = 0;
                for (int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2,&omega,rho[m].t,1,rho[m].r,1);   // r <- r-omega*t
                    phi_dotc_sub(Ns2,rho[m].r,1,rho[m].r,1,&result);
                    diff[id] += real(result);
                }      
                //////////////////////////////////////

            }
            else{
                diff[id] = 0;
            }
            
            /////////////////////////////////////
            barrier_wait(&WriterBarrier);   
            diff_th = 0; 
            for(int m = 0; m < num_threads; ++m){
                if (diff_th < diff[m])
                    diff_th = diff[m];
            } 
            /////////////////////////////////////
            if (id == num_threads) {
                if (step%1 == 0) {
                    //rho[0].SteadyStateRho[0] = 1;
                    result = 1;
                    for(int i = 1; i < Ns; ++i)
                        result += rho[0].SteadyStateRho[(1+Ns)*i];
                    outf << diff_th << " " << 1./real(result) << " ";
                    for (int i = 1; i < Ns*Ns; ++i)
                        outf << rho[0].SteadyStateRho[i]/result << " ";
                    outf << endl;
                }
            }
        }// end while loop

    //}
    barrier_wait(&WriterBarrier);   
    time(&t2);
    // Finished, now write out the results
    if (id == num_threads) {
        pthread_mutex_lock(&io_lock);
        if (step == MaxSteps) {
            cout << "Convergence failed! (steps = " << step << ", diff = " << diff_th << ")\n";
            outf << "#Info: Total steps = " << step << ". No convergence.\n";
        }
        else{
            cout << "Converged to diff = " << diff_th << " in " << step << " steps.\n";
            outf << "#Info: Converged to diff = " << diff_th << " in " << step << " steps.\n";
        }

        if (verbose > 0)
            cout << "Writing final steady-state density matrix.\n";
    
        outf << "#Info: Total running time: " << difftime(t2,t1)/60/60/24 << " days.\n";

        rho[0].SteadyStateRho[0] = 1;
        result = 1;
        for(int i = 1; i < Ns; ++i)
            result += rho[0].SteadyStateRho[(1+Ns)*i];
        for (int i = 0; i < Ns*Ns; ++i)
            rho[0].SteadyStateRho[i] /= result;
        outf << step << " ";
        for (int i = 0; i < Ns*Ns; ++i){
                outf << rho[0].SteadyStateRho[i] << " ";
        }
        outf << endl;

        if (verbose > 1)
            printMatrix(rho[0].SteadyStateRho,Ns);
        if (verbose > 2)
            for (int m = 1; m < MatrixCount; ++m) {
                cout << m << ":\n";
                printMatrix(rho[m].SteadyStateRho,Ns);   
            }
        cout << "Total running time: " <<  difftime(t2,t1)/60/60 << " hours.\n";
        if(restartStep > 0){
            restartf << step << endl;
            for(int m = 0; m < MatrixCount; m++){
                for(int mi = 0; mi < Ns2; mi++)
                        restartf << rho[m].SteadyStateRho[mi] << "\n";
            }
        }   
        pthread_mutex_unlock(&io_lock);
    }   
    barrier_wait(&WriterBarrier);   
    //pthread_exit(0);

};


void HierarchyIntegrator::bicgstab_steadystate_unstable(int id, int num_threads){
    int m_start = id;
    time_t t1,t2;
    time(&t1);

    Complex result;
    Complex rho0, rho1_th, omegaBtm_th, omega, beta, alpha_th;
    Complex omegaInv;
    rho0 = rho1_th = alpha_th = beta = 1;
    omega = -1;
    result = 0;
    int MaxSteps = 9999;
    Float diff_th = 1;
    Float eps = bicgstab_tolerance;
    int step = 0;
    int Ns2 = Ns*Ns;

    if (verbose > 1){ 
        if (id < num_threads) {
            pthread_mutex_lock(&io_lock);
            cout << "[" << id << "]  Calculating steady state density matrix using BiCGSTAB.\n";
            pthread_mutex_unlock(&io_lock);
        }
    }
    else if (verbose > 0 && id == num_threads) {
        pthread_mutex_lock(&io_lock);
        cout << "Calculating steady state density matrix using BiCGSTAB.\n";
        pthread_mutex_unlock(&io_lock);
    }


    barrier_wait(&WriterBarrier);
    //for (int l = 0; l < Lmax; ++l) {
        // Calculate initial residuals
        if (id < num_threads){ 
            for(int m = m_start; m < MatrixCount; m+=num_threads)
                rho[m].rho = rho[m].SteadyStateRho;

            barrier_wait(&MatupdateBarrier);
            for(int m = m_start; m < MatrixCount; m+=num_threads)
                integration_step(rho[m], rho[m].r,step, NULL, NULL, NULL, NULL);        

            for(int m = m_start; m < MatrixCount; m+=num_threads) {
                for (int i = 0; i < Ns2; ++i){
                    rho[m].r[i] = -rho[m].r[i];
                    rho[m].r0[i] = rho[m].r[i];
                }
            }
        }

        barrier_wait(&WriterBarrier);    
        while (diff_th > eps && step < MaxSteps) {
            ++step;
            if ( id < num_threads ) {
                rho1[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r0,1,rho[m].r,1,&result);
                    rho1[id] += result;
                }

                ///////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                rho1_th = 0;
                for(int m = 0; m < num_threads; ++m){
                        rho1_th += rho1[m];
                }
                beta = -rho1_th/rho0*alpha_th; 
                //////////////////////////////////////////

                omegaInv = 1./real(omega);
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2,&omegaInv,rho[m].p,1,rho[m].v,1); // v <- v-p/omega
                    phi_copy(Ns2,rho[m].r,1,rho[m].p,1);           // p <- r
                    phi_axpy(Ns2,&beta,rho[m].v,1,rho[m].p,1);     // p <- p-beta*v
                    rho[m].rho = rho[m].p;
                }

                ///////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                for(int m = m_start; m < MatrixCount; m+=num_threads)
                    integration_step(rho[m], rho[m].v,step, NULL, NULL, NULL, NULL); // v <- A*p
                ///////////////////////////////////////////

                alpha[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_dotc_sub(Ns2,rho[m].r0,1,rho[m].v,1,&result);
                    alpha[id] += result;            
                }

                /////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                alpha_th = 0;
                for(int m = 0; m < num_threads; ++m){
                        alpha_th += alpha[m];
                }
                alpha_th = -rho1_th/alpha_th; 
                /////////////////////////////////////////

                for (int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2,&alpha_th,rho[m].v,1,rho[m].r,1); // r <- r-alpha*v
                    rho[m].rho = rho[m].r;                              // ready A*r
                }
                alpha_th = -alpha_th; 

                ////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                for(int m = m_start; m < MatrixCount; m+=num_threads)
                    integration_step(rho[m], rho[m].t,step, NULL, NULL, NULL, NULL); // t <- A*r
                //////////////////////////////////////

                omegaBtm[id] = 0;
                omegaTop[id] = 0;
                for(int m = m_start; m < MatrixCount; m+=num_threads){
                   phi_dotc_sub(Ns2,rho[m].t,1,rho[m].r,1,&result); // dot(t,r)
                   omegaTop[id] += result;
                   phi_dotc_sub(Ns2,rho[m].t,1,rho[m].t,1,&result); // dot(t,t)
                   omegaBtm[id] += result;
                }

                /////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                omegaBtm_th = 0;
                omega = 0;
                for(int m = 0; m < num_threads; ++m){
                        omega += omegaTop[m];
                        omegaBtm_th += omegaBtm[m];
                }
                omega /= omegaBtm_th;  
                /////////////////////////////////////////

                // Update rho(\infty)
                for (int m = m_start; m < MatrixCount; m+=num_threads) {
                   phi_axpy(Ns2,&alpha_th,rho[m].p,1,rho[m].SteadyStateRho,1);  //rho <- rho + alpha*p
                   phi_axpy(Ns2,&omega,rho[m].r,1,rho[m].SteadyStateRho,1);      //rho <- rho + omega*s
                }

                omega = -omega;

                // Recalculate residuals            
                ////////////////////////////////////////
                barrier_wait(&MatupdateBarrier);
                diff[id] = 0;
                for (int m = m_start; m < MatrixCount; m+=num_threads) {
                    phi_axpy(Ns2,&omega,rho[m].t,1,rho[m].r,1);   // r <- r-omega*t
                }
                if (id == 0) {
                    phi_dotc_sub(Ns2,rho[0].r,1,rho[0].r,1,&result);
                    diff[id] += real(result);
                }      
                //////////////////////////////////////

                rho0 = rho1_th;
            }
            else{

            }
            
            /////////////////////////////////////
            barrier_wait(&WriterBarrier);   
            diff_th = 0; 
            for(int m = 0; m < num_threads; ++m){
                diff_th += diff[m];
            } 
            /////////////////////////////////////
            if (id == num_threads) {
                if (step%10 == 0) {
                    outf << diff_th << " ";
                    for (int i = 0; i < Ns*Ns; ++i)
                        outf << rho[0].SteadyStateRho[i] << " ";
                    outf << endl;
                }
            }
        }// end while loop

    //}
    barrier_wait(&WriterBarrier);   
    time(&t2);
    if (id == num_threads) {
        pthread_mutex_lock(&io_lock);
        if (step == MaxSteps) {
            cout << "Convergence failed! (steps = " << step << ", diff = " << diff_th << ")\n";
            outf << "#Info: Total steps = " << step << ". No convergence.\n";
        }
        else{
            cout << "Converged to diff = " << diff_th << " in " << step << " steps.\n";
            outf << "#Info: Converged to diff = " << diff_th << " in " << step << " steps.\n";
        }

        if (verbose > 0)
            cout << "Writing final steady-state density matrix.\n";

        outf << "#Info: Total running time: " << difftime(t2,t1)/60/60/24 << " days.\n";

        outf << diff_th << " ";
        for (int i = 0; i < Ns*Ns; ++i){
                outf << rho[0].SteadyStateRho[i] << " ";
        }
        outf << endl;

        if (verbose > 1)
            printMatrix(rho[0].SteadyStateRho,Ns);
        if (verbose > 2)
            for (int m = 1; m < MatrixCount; ++m) {
                cout << m << ":\n";
                printMatrix(rho[m].SteadyStateRho,Ns);   
            }
        cout << "Total running time: " <<  difftime(t2,t1)/60/60 << " hours.\n";
        if(restartStep > 0){
            restartf << step << endl;
            for(int m = 0; m < MatrixCount; m++){
                for(int mi = 0; mi < Ns2; mi++)
                        restartf << rho[m].SteadyStateRho[mi] << "\n";
            }
        }   
        pthread_mutex_unlock(&io_lock);
    }   
    barrier_wait(&WriterBarrier);   
    //pthread_exit(0);

};


////////////////////////////////////////
//                                    //
// Set Up Sparse Hamiltonian Operator //
//                                    //
////////////////////////////////////////
void HierarchyIntegrator::constructSparseOperator(){
    cout << "Preparing steady state Liouville operator indices.\n";
    Complex MinusOneImagOverHbar(0,-1./5.29);
    long totalEntries = 0;
    Float sameTime, nextTime, prevTime;
    sameTime = nextTime = prevTime = 0;
    time_t t1,t2;


    ////////// Create Sparse Hamiltonian Matrix ///////////
    cout << "\tCreating sparse Hamiltonian operator.\n";
    Complex * H_l;
    H_l = new Complex[Ns*Ns*Ns*Ns];
    H_l = toLiouville(H,Ns);
    int NumHentries = 0;
    int m, n;
    m = n = 0;
    for (int i = 0; i < Ns*Ns; ++i){
        for (int j = 0; j < Ns*Ns; ++j){
            H_l[n] *= MinusOneImagOverHbar;
            if (abs(H_l[n]) > 1e-10 || i == j)
                NumHentries += 1;
            ++n;
         }
    }
    int n_diag = 0;
    Index2D * Hindices;
    Hindices = new Index2D[NumHentries];
    Complex * LiouvilleH;
    LiouvilleH = new Complex[NumHentries];
    int *DiagonalIndices;
    DiagonalIndices = new int[Ns*Ns];
    Index2D * HBlockIndices;
    HBlockIndices = new Index2D[Ns*Ns];
    m = n = 0;
    for (int i = 0; i < Ns*Ns; ++i){
      for (int j = 0; j < Ns*Ns; ++j){
        if (abs(H_l[n]) > 1e-10 || i == j){
            LiouvilleH[m] = H_l[n];
            Hindices[m].i = i;
            Hindices[m].j = j;
            if ( i == j ) {
                DiagonalIndices[n_diag] = m;
                HBlockIndices[n_diag].j = i%Ns;
                HBlockIndices[n_diag].i = (i-i%Ns)/Ns;
                ++n_diag;
            }
            ++m;
        }
        ++n;
      }
    }
    cout << "\tAssigning memory for SameLiouville operator.\n";
    ////////////////////////////////////////////////////
    cout << "\tNext and Prev matrix indices.\n";
    ////// Indices for Next and Prev Matrices //////
    int **NextIndices;
    int **PrevIndices;
    int NumNextEntries = 2*Ns-2;
    int NumPrevEntries = 2*Ns-1;
    NextIndices = new int*[M];
    PrevIndices = new int*[M];
    Index2D **NextI2D;
    NextI2D = new Index2D*[M];
    Index2D **PrevI2D;
    PrevI2D = new Index2D*[M];
    for (int nj = 0; nj < M; ++nj) {
        n = 0;
        m = MDiagInd[nj]; 
        NextIndices[nj] = new int[NumNextEntries];
        NextI2D[nj] = new Index2D[NumNextEntries];
        PrevIndices[nj] = new int[NumPrevEntries];
        PrevI2D[nj] = new Index2D[NumPrevEntries];
        for (int i = 0; i < m; ++i) {
            NextIndices[nj][n] = i*Ns+m;
            NextI2D[nj][n].i = NextIndices[nj][n];
            NextI2D[nj][n].j = NextIndices[nj][n];
            PrevIndices[nj][n] = i*Ns+m;
            PrevI2D[nj][n].i = PrevIndices[nj][n];
            PrevI2D[nj][n].j = PrevIndices[nj][n];
            ++n;
        }
        for (int i = 0 ; i < m; ++i){
            NextIndices[nj][n] = m*Ns+i;
            NextI2D[nj][n].i = NextIndices[nj][n];
            NextI2D[nj][n].j = NextIndices[nj][n];
            PrevIndices[nj][n] = m*Ns+i;
            PrevI2D[nj][n].i = PrevIndices[nj][n];
            PrevI2D[nj][n].j = PrevIndices[nj][n];
            ++n;
        }
        PrevIndices[nj][n] = m*Ns+m;
        PrevI2D[nj][n].i = PrevIndices[nj][n];
        PrevI2D[nj][n].j = PrevIndices[nj][n];
        for (int i = m+1 ; i < Ns; ++i){
            NextIndices[nj][n] = m*Ns+i;
            NextI2D[nj][n].i = NextIndices[nj][n];
            NextI2D[nj][n].j = NextIndices[nj][n];
            PrevIndices[nj][n+1] = m*Ns+i;
            PrevI2D[nj][n+1].i = PrevIndices[nj][n+1];
            PrevI2D[nj][n+1].j = PrevIndices[nj][n+1];
            ++n;
        }
        for (int i = m+1; i < Ns; ++i) {
            NextIndices[nj][n] = i*Ns+m;
            NextI2D[nj][n].i = NextIndices[nj][n];
            NextI2D[nj][n].j = NextIndices[nj][n];
            PrevIndices[nj][n+1] = i*Ns+m;
            PrevI2D[nj][n+1].i = PrevIndices[nj][n+1];
            PrevI2D[nj][n+1].j = PrevIndices[nj][n+1];
            ++n;
        }
    }
    ////////////////////////////////////////////////////////

    n = 0;
    //long averageSameEntries = 0;
    cout << "Constructing Individual Liouville operators.\n";
    for(int i = 0; i < MatrixCount; ++i){
        time(&t1);
        rho[i].createSameLiouvilleOperator(LiouvilleH,DiagonalIndices,HBlockIndices,NumHentries);
        time(&t2);
        sameTime += difftime(t2,t1)/60;
        time(&t1);
        rho[i].createNextLiouvilleOperator(NumNextEntries);
        time(&t2);
        nextTime += difftime(t2,t1)/60;
        time(&t1);
        rho[i].createPrevLiouvilleOperator(NumPrevEntries);
        time(&t2);
        prevTime += difftime(t2,t1)/60;
        for (int j = 0; j < rho[i].Nnext; ++j)
            totalEntries += NumNextEntries;
        for (int j = 0; j < rho[i].Nprev; ++j)
            totalEntries += NumPrevEntries;            
        n += NumHentries;
    }
    totalEntries += MatrixCount*NumHentries;
    //totalEntries += averageSameEntries;
    //cout << "Average Same Entries = " << ((Float)averageSameEntries)/MatrixCount << "\n"
    //     << "Sparse H Entries = " << NumHentries << "\n";
    cout << "\n" 
         << "Time for same: " << sameTime << " mins.\n" 
         << "Time for next: " << nextTime << " mins.\n" 
         << "Time for prev: " << prevTime << " mins.\n";

    // For Debugging 
    if (verbose > 2) {
        cout << "Ls[0]=\n";
        printMatrix(rho[0].SameLiouville,Hindices,NumHentries);
        /*
        cout << "Ls[1]=\n";
        printMatrix(rho[1].SameLiouville,rho[1].SameLiouvilleIndices,rho[1].NumSameLiouvilleEntries);
        cout << "\n";
        */
    
        cout << "Vmx[0]=\n";
        printMatrix(rho[0].NextLiouville[0],NextI2D[rho[0].NextIndMatrix[0]],NumNextEntries);
        cout << "Vmx[2]=\n";
        printMatrix(rho[0].NextLiouville[2],NextI2D[rho[0].NextIndMatrix[2]],NumNextEntries);
        
        cout << "cVmkx[0][0]=\n";
        cout << rho[1].Nprev << " " << rho[1].PrevIndMatrix[0] << endl;
        printMatrix(rho[1].PrevLiouville[0],PrevI2D[rho[1].PrevIndMatrix[0]],NumPrevEntries);
        cout << "cVmkx[0][1]=\n";
        cout << rho[2].Nprev << " " << rho[2].PrevIndMatrix[0] << endl;
        printMatrix(rho[2].PrevLiouville[0],PrevI2D[rho[2].PrevIndMatrix[0]],NumPrevEntries);
    }

    long lda = MatrixCount*Ns*Ns;
    Float MatrixMem = lda*lda*sizeof(Complex); //bytes
    MatrixMem = MatrixMem/1024/1024/1024;      //Gbytes
    Float ReducedMem = totalEntries*sizeof(Complex);
    ReducedMem = ReducedMem/1024/1024/1024;
    cout << "Writing total operator to file.\n";
    if (verbose > 0) {
        cout << "Total operator matrix size: " << lda << "x" << lda
             << " (" << MatrixMem << " Gb)\n";
        cout << "Total operator sparse matrix size: " << totalEntries 
             << " (" << ReducedMem << " Gb)\n";
        if (verbose > 1)
            cout << "\tsizeof(Complex)=" <<sizeof(Complex) << " bytes\n";
    }
};

HierarchyIntegrator::~HierarchyIntegrator(){

/*    if (Ns > 0) {
        if (H != NULL)
          delete []H;
        if (Hdag != NULL)
          delete []Hdag;
        if (Evec != NULL)
          delete []Evec;
        if (Eval != NULL)
          delete []Eval;
        if (rho0 != NULL)
          delete []rho0;
        if (kappa != NULL)
          delete []kappa;
        if (TDM != NULL)
          delete []TDM;
        if (Vbath_re != NULL)
          delete []Vbath_re;
        if (lambda != NULL)
          delete []lambda;
        if (gamma != NULL)
          delete []gamma;
    }

    if (spectrum_rho0 != NULL)
      delete []spectrum_rho0;
  */
}



