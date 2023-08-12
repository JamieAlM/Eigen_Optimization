#define _USE_MATH_DEFINES // for C++
#define FFTW_ESTIMATE (1U << 6)
//#define EIGEN_DONT_PARALLELIZE
//#define NDEBUG // for assert, when including this line it will disable the error message 
#include <cmath>
#include<math.h>
#include<stdio.h>
#include "fftw3.h"
#include <cstring>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <iostream> 
#include <vector>
#include <fstream> //for file operations
#include <iomanip>
#include <numeric>
#include <assert.h>
#include<valarray>
#include <complex>

static const int nx = 256;//128;
static const int ny = 256;//128; 
static const int nz = 256;//128;
static const int ncomp = 2;
static const int nzk = nz/2 + 1;
using namespace std;
void cross_product(double vector_a[], double vector_b[], double temp[]);
double mag1DArray(double arr[]);
double makeSpatialMesh3D(double dx, double dy, double dz,double xarr[], double yarr[],double zarr[]);
//void FourierMesh3D(double Lx,double Ly,double Lz, double KX[], double KY[],double KZ[], double ksqu[], double ninvksqu[], int FourierMeshType);
void FourierMesh3D(double Lx,double Ly,double Lz, double kx[], double ky[], double kz[], double KX[], double KY[],double KZ[], double ksqu[], double ninvksqu[], int FourierMeshType);
void c2r3D(double cArr[][ncomp], double rArr[]);
void r2c3D(double rArr[], double cArr[][ncomp]);
void Derivk3D(double vark[][ncomp], double K[], double derivative[][ncomp]);
void calcV_ExBk3D(double dphidxk[][ncomp], double dphidyk[][ncomp],double dphidzk[][ncomp], double B[], double B2, double vexbkx[][ncomp], double vexbky[][ncomp],double vexbkz[][ncomp]);
void fourierDivision3D(double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]);
void calc_diamag3D(double dpdxk[][ncomp], double dpdyk[][ncomp], double dpdzk[][ncomp],double B[], double B2, double qa, double nak[][ncomp], double diamagxk[][ncomp], double diamagyk[][ncomp], double diamagzk[][ncomp]);
void CollFreqk_inertia3D(double nek[][ncomp], double Tik[][ncomp], double Tek[][ncomp], double kb, double eps0, double mi, double me, double ri, double rn, double nn, double Oci, double Oce, double e,double nuink[][ncomp], double nuiek[][ncomp] , double nuiik[][ncomp] , double nuenk[][ncomp] ,double nueek[][ncomp] ,double nueik[][ncomp]  , double isigPk[][ncomp] , double invnk[][ncomp] , double hallIk[][ncomp] ,double hallEk[][ncomp] );
double Arr1DMax(double arr[], int arrLength);
double max3Dk(double arr2D[]);
double calc_dt3D(double U[], double vexbx[], double vexby[],double vexbz[] ,double diamagxi[], double diamagyi[], double diamagzi[],double diamagxe[], double diamagye[], double diamagze[],double cfl, double kmax, double maxdt);
void Convolve3D( double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]);
void Arr2DArr3DMult(double arrIn0[], double arrIn1[], double arrOut[]);
void calc_residualn3D(double vexbxk[][ncomp], double vexbyk[][ncomp], double vexbzk[][ncomp], double nink[][ncomp], double residnoutk[][ncomp], double kx[], double ky[], double kz[]);
void calc_residualt3D(double voxk[][ncomp], double voyk[][ncomp],double vozk[][ncomp], double tempink[][ncomp], double tempoutk[][ncomp], double kx[], double ky[],double kz[]);
void calcPotSourcek_inertia3D(double ne0k[][ncomp], double nek[][ncomp],double dndx0k[][ncomp], double dndy0k[][ncomp],double dndz0k[][ncomp], double dndxk [][ncomp], double dndyk [][ncomp],double dndzk [][ncomp], double dphidx0k [][ncomp], double dphidy0k [][ncomp], double dphidz0k [][ncomp], double dphidx1k [][ncomp], double dphidy1k [][ncomp], double dphidz1k [][ncomp],double Pi1k[][ncomp], double Pe1k[][ncomp], double uxB[], double e, double Cm, double hallEk [][ncomp], double hallIk [][ncomp], double vexbkx0[][ncomp], double vexbky0[][ncomp],double vexbkz0[][ncomp], double vexbkx[][ncomp],  double vexbky[][ncomp], double vexbkz[][ncomp],double kx[], double ky[], double kz[], double ksqu[], double potSourcek_inertia[][ncomp]);
int Potentialk3D(double invnk[][ncomp], double dndxk[][ncomp], double dndyk[][ncomp], double dndzk[][ncomp],double phik[][ncomp], double potSourcek[][ncomp], double kx[], double ky[], double kz[],double ninvksqu[], double err_max, int max_iter);
double max_absComp3D(double arr3D[][ncomp]);
void calc_sourcen3D(double ksqu[], double nk[][ncomp], double d, double sourcenk[][ncomp]);
void Laplaciank3D( double vark[][ncomp], double ksqu[], double derivative[][ncomp]);
void RK4(double f[][ncomp], double dt, double residual[][ncomp], double source[][ncomp], int stage, double fout[][ncomp]);

//g++ MemoryLeakTest.cpp -lfftw3 -lm -ffast-math -fno-math-errno -march=native -Ofast -ggdb3 -o ../../../Datasim1/test.out && valgrind --track-origins=yes ./../../../Datasim1/test.out
#define sizee 3
#define IMAG 1
#define REAL 0
int main(){	
	double saveFrequency = 1.0; //20.;  //1.0;    // Save data every this many time steps
	double dt_max = 20.;//20.; //0.01;//0.3; //20;       // Set max allowable time step
	double tend = 150;// 4000; //500.;     //1000.     // Set the end time
	double err_max = 1.e-8;    //1.e-6?  // Set max allowable error for potential solver
	double CFL = 1.;//0.05;//1.; //3.;             // Set CFL number. We can just leave this at 3.
	double Dart = 0.1;//0.05; //0.1;// OR 0.11e3;     //m^2/s // Set artifical diffusion constants
	int phi_iter_max = 500;      // Max allowable iterations for potential solver
	int phi_iter = 0;
	int saveNum = 1;
	// Calculated parameters	
	int iter_max = 100000;
    // Set physical constants
	double e = 1.6021766208E-19; //e = 1.6021766208e-19;      % C
	double kb = 1.380649E-23; //kb = 1.380649e-23;
	double me = 9.10938356E-31; //me = 9.10938356e-31;       % kg
	double eps0 = 8.854187817E-12; //eps0 = 8.854187817e-12;    % F/m
	double mi = 16. * 1.66053906660E-27;   // Mass O+
	double mn = mi;               // Mass O
	double ri = 152.E-12;          // Effective collision radius
	double rn = ri;
	double B[3] = {0,0,5E-5}; //{0,0,5E-5}; //{0,0,1};  this norm units so unrealistic regime
    double u[3] = {500, 0, 0}; //{10, 0, 0}; 
    double Lx = 640;//128; //2*128; //32;//2*EIGEN_PI;
    double Ly = 640;//128; //32;//2*EIGEN_PI;
    double Lz = 640;//128; //32;//4;//32;//2*EIGEN_PI;
    double x0 = Lx/2;//64.; ///4.;//2D Cheby/Fourier: Lx/2;  
	double y0 = Ly/2;//64.; ///4.; //2D Cheby/Fourier: 0.0;
    double z0 = Lz/2;//64.; ///4.; 
    double rx0 = 100.;
    double ry0 =  100.;
    double rz0 =  100.;
    double uxB[sizee];
	double Bmag = mag1DArray(B); 
    double B2 = Bmag * Bmag;
   cross_product(u, B, uxB);
    double Oci = e*Bmag/mi;
	double Oce = e*Bmag/me;
	double Cm = 1./(1 / Oci + 1/Oce);
	double m = 1.;
    double dx = Lx / nx;
	double dy = Ly / ny;
    double dz = Lz / nz;
    double o = 2.5; //outer region
	double a = 4.; //for plasma enhancement case a =4 and for plasma depletion case a = -0.05
	double d = 1.; 
	double p = 3.; 
	double nn = 1e14; 
	
	

 

		
	double *XX;
	XX = (double*) fftw_malloc((nx*ny*nz) *sizeof(double));
	memset(XX, 42, (nx*ny*nz) * sizeof(double)); 

    double *YY;
	YY = (double*) fftw_malloc((nx*ny*nz) *sizeof(double));
	memset(YY, 42, (nx*ny*nz) * sizeof(double)); //test: XX for some reason isn't the correct dimensions, maybe this is why
	
	double *ZZ;
	ZZ = (double*) fftw_malloc((nx*ny*nz) *sizeof(double));
	memset(ZZ, 42, (nx*ny*nz) * sizeof(double));

	double kmax = makeSpatialMesh3D(dx,dy,dz,XX,YY,ZZ); 

	double *kx;
	kx = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(kx, 42, (nx*ny*nz)* sizeof(double)); 
	double *ky;
	ky = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(ky, 42, (nx*ny*nz)* sizeof(double)); 

    double *kz;
    kz = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(kz, 42, (nx*ny*nz)* sizeof(double)); 
    //##############################
    double *KXX;
	KXX = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(KXX, 42, (nx*ny*nz)* sizeof(double)); 
	double *KYY;
	KYY = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(KYY, 42, (nx*ny*nz)* sizeof(double)); 

    double *KZZ;
    KZZ = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(KZZ, 42, (nx*ny*nz)* sizeof(double)); 
    //###################
	double *ksqu;
	ksqu = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(ksqu, 42, (nx*ny*nz)* sizeof(double)); //test
	
	double *ninvksqu;
	ninvksqu = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(ninvksqu, 42, (nx*ny*nz)* sizeof(double)); 

	int FourierMeshType = 1;
    //FourierMesh3D(Lx, Ly, Lz, KXX, KYY, KZZ, ksqu, ninvksqu, FourierMeshType);
	FourierMesh3D(Lx, Ly, Lz, kx, ky, kz, KXX, KYY, KZZ, ksqu, ninvksqu, FourierMeshType);

	double *Theta;
	Theta = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Theta, 42, (nx*ny*nz)* sizeof(double));

    double *Phiz;
	Phiz = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Phiz, 42, (nx*ny*nz)* sizeof(double)); 

    double *P0;
	P0 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(P0, 42, (nx*ny*nz)* sizeof(double)); 

    double *nep; //ne1
	nep = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nep, 42, (nx*ny*nz)* sizeof(double)); 

    double *nebg; //ne0
	nebg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nebg, 42, (nx*ny*nz)* sizeof(double)); 

	double *neTot; //ne
	neTot = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(neTot, 42, (nx*ny*nz)* sizeof(double)); 

    double *Tip; //Ti_1
	Tip = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Tip, 42, (nx*ny*nz)* sizeof(double)); 

    double *Tep; //Te_1
	Tep = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Tep, 42, (nx*ny*nz)* sizeof(double));

    double *Tibg; //Ti_0
	Tibg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Tibg, 42, (nx*ny*nz)* sizeof(double)); 

    double *Tebg; //Te_0
	Tebg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Tebg, 42, (nx*ny*nz)* sizeof(double));

    double *TiTot; //Ti
	TiTot = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(TiTot, 42, (nx*ny*nz)* sizeof(double));

    double *TeTot; //Ti
	TeTot = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(TeTot, 42, (nx*ny*nz)* sizeof(double));

    double *phibg; //phi_0
	phibg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(phibg, 42, (nx*ny*nz)* sizeof(double));

    double *phip; //phi1
	phip = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(phip, 42, (nx*ny*nz)* sizeof(double));

    double *phiTot; //phi
	phiTot = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(phiTot, 42, (nx*ny*nz)* sizeof(double));

    double *piTot; //phi
	piTot = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(piTot, 42, (nx*ny*nz)* sizeof(double));

    double *peTot; //phi
	peTot = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(peTot, 42, (nx*ny*nz)* sizeof(double));

    double *pibg; //pi0
	pibg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(pibg, 42, (nx*ny*nz)* sizeof(double));

    double *pebg; //pe0
	pebg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(pebg, 42, (nx*ny*nz)* sizeof(double));

    double *pip; //pi1
	pip = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(pip, 42, (nx*ny*nz)* sizeof(double));

    double *pep; //pe1
	pep = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(pep, 42, (nx*ny*nz)* sizeof(double));

	double n0 = 1e11;
	double Ampl_1 = -0.07;//min
	double Ampl_2 = 0.07; //max
	double c = 12.; 
	for(int i = 0; i< nx; i++){
			for(int j = 0; j< ny; j++){
				for(int k = 0; k< nz; k++){ 
                    Theta[k + nz * (j + ny * i)] = atan2(YY[k + nz * (j + ny * i)]  -y0,  XX[k + nz * (j + ny * i)]  -x0)  - (M_PI/c); 
					Phiz[k + nz * (j + ny * i)] = atan2(ZZ[k + nz * (j + ny * i)]  -z0,  YY[k + nz * (j + ny * i)]  -y0)  - (M_PI/c);
                    P0[k + nz * (j + ny * i)] = (pow(XX[k + nz * (j + ny * i)]-x0,2)/pow(rx0,2) + pow(YY[k + nz * (j + ny * i)]-y0,2)/pow(ry0,2) + pow(ZZ[k + nz * (j + ny * i)]-z0,2)/pow(rz0,2) );               
                    nep[k + nz * (j + ny * i)] =( d + a * exp(pow(-1. * P0[k + nz * (j + ny * i)],p)) + ((rand() % 101)/100. * (Ampl_2 - Ampl_1) + Ampl_1) ) * n0  ;
                    nebg[k + nz * (j + ny * i)] =0.;
                    neTot[k + nz * (j + ny * i)] = nep[k + nz * (j + ny * i)]+ nebg[k + nz * (j + ny * i)];
					Tip[k + nz * (j + ny * i)] = 0.; 
			    	Tep[k + nz * (j + ny * i)] = 0.; 
			    	Tibg[k + nz * (j + ny * i)] = 1000.;
			    	Tebg[k + nz * (j + ny * i)] = 1000.;
                    TiTot[k + nz * (j + ny * i)] = Tip[k + nz * (j + ny * i)] + Tibg[k + nz * (j + ny * i)];
                    TeTot[k + nz * (j + ny * i)] = Tep[k + nz * (j + ny * i)] + Tebg[k + nz * (j + ny * i)];
                    phibg[k + nz * (j + ny * i)] = 0.; 
                    phip[k + nz * (j + ny * i)] = 0.; 
                    phiTot[k + nz * (j + ny * i)] = phibg[k + nz * (j + ny * i)] + phip[k + nz * (j + ny * i)] ; 
                    piTot[k + nz * (j + ny * i)] = neTot[k + nz * (j + ny * i)] * TiTot[k + nz * (j + ny * i)] * kb;
                    peTot[k + nz * (j + ny * i)] = neTot[k + nz * (j + ny * i)] * TeTot[k + nz * (j + ny * i)] * kb;
                    pibg[k + nz * (j + ny * i)] = nebg[k + nz * (j + ny * i)] * Tibg[k + nz * (j + ny * i)] *kb;
                    pebg[k + nz * (j + ny * i)] = nebg[k + nz * (j + ny * i)] * Tebg[k + nz * (j + ny * i)] *kb;
                    pip[k + nz * (j + ny * i)] = piTot[k + nz * (j + ny * i)] - pibg[k + nz * (j + ny * i)] *kb;
                    pep[k + nz * (j + ny * i)] = peTot[k + nz * (j + ny * i)] - pebg[k + nz * (j + ny * i)] *kb;
                    
				}
			}
	}


    fftw_complex *neK; //nek
	neK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(neK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *TiK; //nek
	TiK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(TiK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *TeK; //nek
	TeK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(TeK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *phiK; //nek
	phiK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(phiK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *piK; //nek
	piK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(piK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *peK; //nek
	peK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(peK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *nebgK; //nek
	nebgK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nebgK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *nepK; //nek
	nepK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nepK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *TebgK; //nek
	TebgK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(TebgK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *TibgK; //nek
	TibgK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(TibgK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *phibgK; //nek
	phibgK =  (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(phibgK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *phipK; //nek
	phipK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(phipK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *pibgK; //nek
	pibgK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(pibgK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *pipK; //nek
	pipK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(pipK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *pebgK; //nek
	pebgK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(pebgK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *pepK; //nek
	pepK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(pepK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	
	/*double *TestOut; //ne
	TestOut = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(TestOut, 42, (nx*ny*nz)* sizeof(double)); */

	/*r2c3D(neTot, neK);
	c2r3D(neK,TestOut);

	std::ofstream file1("KXX.txt");
    if (file1.is_open())
    {
        for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        file1 << KXX[k + nzk * (j + ny * i)] << '\n'; 
                    }
                }
        }
    }*/
	r2c3D(neTot, neK);
    r2c3D(TiTot, TiK);
    r2c3D(TeTot, TeK);
    r2c3D(phiTot, phiK);
    r2c3D(piTot, piK);
    r2c3D(peTot, peK);

    r2c3D(nebg, nebgK);
    r2c3D(nep, nepK);
    r2c3D(Tebg, TebgK);
    r2c3D(Tibg, TibgK);

    r2c3D(phibg, phibgK);
	r2c3D(phip, phipK);
	r2c3D(pibg, pibgK);
	r2c3D(pip, pipK);
	r2c3D(pebg, pebgK);
	r2c3D(pep, pepK);

	fftw_complex *hallEK; //nek
	hallEK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(hallEK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *hallIK; //nek
	hallIK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(hallIK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *nuinK; //nek
	nuinK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nuinK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *nuieK; //nek
	nuieK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nuieK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *nuiiK; //nek
	nuiiK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nuiiK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *nuenK; 
	nuenK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nuenK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *nueeK; 
	nueeK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nueeK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *nueiK; 
	nueiK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(nueiK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *isigPK; 
	isigPK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(isigPK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *invnK; 
	invnK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(invnK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    CollFreqk_inertia3D( neK, TiK, TeK, kb, eps0, mi, me, ri, rn, nn, Oci, Oce, e, nuinK, nuieK , nuiiK , nuenK ,nueeK ,nueiK , isigPK , invnK , hallIK , hallEK );
  
	fftw_complex *potSourceK_inertia;
	potSourceK_inertia = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(potSourceK_inertia, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	

	fftw_complex *dndxK;
	dndxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dndyK;
	dndyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dndzK;
	dndzK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndzK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dndx0K;
	dndx0K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndx0K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dndy0K;
	dndy0K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndy0K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dndz0K;
	dndz0K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndz0K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dndx1K;
	dndx1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndx1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *dndy1K;
	dndy1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndy1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dndz1K;
	dndz1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dndz1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dphidxK;
	dphidxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dphidyK;
	dphidyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dphidzK;
	dphidzK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidzK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dphidx0K;
	dphidx0K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidx0K, 42, (nx*ny*nzk)* sizeof(fftw_complex));


    fftw_complex *dphidx1K;
	dphidx1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidx1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
    
    
    fftw_complex *dphidy0K;
	dphidy0K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidy0K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *dphidy1K;
	dphidy1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidy1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dphidz0K;
	dphidz0K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidz0K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dphidz1K;
	dphidz1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dphidz1K, 42,(nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *d2phidxK;
	d2phidxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(d2phidxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *d2phidyK;
	d2phidyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(d2phidyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *d2phidzK;
	d2phidzK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(d2phidzK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dpedxK;
	dpedxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpedxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dpedx1K;
	dpedx1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpedx1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
   
    fftw_complex *dpedyK;
	dpedyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpedyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *dpedy1K;
	dpedy1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpedy1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dpedzK;
	dpedzK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpedzK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dpedz1K;
	dpedz1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpedz1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
  
    fftw_complex *dpidxK;
	dpidxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpidxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *dpidx1K;
	dpidx1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpidx1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dpidyK;
	dpidyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpidyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    fftw_complex *dpidy1K;
	dpidy1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpidy1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
    
    fftw_complex *dpidzK;
	dpidzK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpidzK, 42, (nx*ny*nzk)* sizeof(fftw_complex));
	
    fftw_complex *dpidz1K;
	dpidz1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(dpidz1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));
    
    fftw_complex *dphiKdt;
	dphiKdt = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        for(int l = 0; l<ncomp; ++l){
                            dphiKdt[k + nzk * (j + ny * i)][l] = 1.; 
                        }
                        
                    }
                }
        }

	Derivk3D(neK, KXX, dndxK);
	Derivk3D(neK, KYY, dndyK);
	Derivk3D(neK, KZZ, dndzK);

	Derivk3D(nebgK, KXX, dndx0K);
	Derivk3D(nebgK, KYY, dndy0K);
    Derivk3D(nebgK, KZZ, dndz0K);
	// add ne1k //dndx1k dndy1k
	Derivk3D(nepK, KXX, dndx1K);
	Derivk3D(nepK, KYY, dndy1K);
    Derivk3D(nepK, KZZ, dndz1K);

	Derivk3D(phiK, KXX, dphidxK);
	Derivk3D(phiK, KYY, dphidyK);
    Derivk3D(phiK, KZZ, dphidzK);

    // add phi0k
	Derivk3D(phibgK, KXX, dphidx0K);
	Derivk3D(phibgK, KYY, dphidy0K);
    Derivk3D(phibgK, KZZ, dphidz0K);
	// add phi1k
	Derivk3D(phipK, KXX, dphidx1K);
	Derivk3D(phipK, KYY, dphidy1K);
    Derivk3D(phipK, KZZ, dphidy1K);

    
	Derivk3D(peK, KXX, dpedxK);
	Derivk3D(peK, KYY, dpedyK);
    Derivk3D(peK, KZZ, dpedzK);

	Derivk3D(piK, KXX, dpidxK);
	Derivk3D(piK, KYY, dpidyK);
    Derivk3D(piK, KZZ, dpidzK);
	// add Pe1k and Pi1k ONLY for now
	Derivk3D(pepK, KXX, dpedx1K); 
	Derivk3D(pepK, KYY, dpedy1K);
    Derivk3D(pepK, KZZ, dpedz1K);

	Derivk3D(pipK, KXX, dpidx1K);
	Derivk3D(pipK, KYY, dpidy1K);
    Derivk3D(pipK, KZZ, dpidz1K);


	/*c2r3D(dndxK,TestOut);

	std::ofstream file1("dndx.txt");
    if (file1.is_open())
    {
        for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        file1 << TestOut[k + nz * (j + ny * i)] << '\n'; 
                    }
                }
        }
    }*/

	fftw_complex *vexbKx;
	vexbKx = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKx, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbKy;
	vexbKy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	memset(vexbKy, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbKz;
	vexbKz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	memset(vexbKz, 42, (nx*ny*nzk)* sizeof(fftw_complex));

    calcV_ExBk3D(dphidxK, dphidyK,dphidzK, B, B2,  vexbKx,  vexbKy,vexbKz );


	double *vexbX;
	vexbX = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbX, 42, (nx*ny*nz)* sizeof(double)); 

	double *vexbY;
	vexbY = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbY, 42, (nx*ny*nz)* sizeof(double)); 

	double *vexbZ;
	vexbZ = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbZ, 42, (nx*ny*nz)* sizeof(double)); 

	c2r3D(vexbKx, vexbX);
	c2r3D(vexbKy, vexbY);
    c2r3D(vexbKz, vexbZ);


	fftw_complex *vexbKx0;
	vexbKx0 = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKx0, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbKy0;
	vexbKy0 = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKy0, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbKz0;
	vexbKz0 = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKz0, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	calcV_ExBk3D(dphidx0K, dphidy0K, dphidz0K, B, B2, vexbKx0,vexbKy0,vexbKz0);
	
	double *vexbX0;
	vexbX0 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbX0, 42, (nx*ny*nz)* sizeof(double)); 

	double *vexbY0;
	vexbY0 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbY0, 42, (nx*ny*nz)* sizeof(double)); 

	double *vexbZ0;
	vexbZ0 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbZ0, 42, (nx*ny*nz)* sizeof(double)); 

	c2r3D(vexbKx0,vexbX0);
	c2r3D(vexbKy0, vexbY0);
    c2r3D(vexbKz0, vexbZ0);

	fftw_complex *vexbKx1;
	vexbKx1 = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKx1, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbKy1;
	vexbKy1 = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKy1, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbKz1;
	vexbKz1 = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vexbKz1, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	calcV_ExBk3D(dphidx1K, dphidy1K, dphidz0K, B, B2, vexbKx1,vexbKy1,vexbKz1);

	double *vexbX1;
	vexbX1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbX1, 42, (nx*ny*nz)* sizeof(double)); 

	double *vexbY1;
	vexbY1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbY1, 42, (nx*ny*nz)* sizeof(double)); 

	double *vexbZ1;
	vexbZ1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vexbZ1, 42, (nx*ny*nz)* sizeof(double));


	c2r3D(vexbKx1, vexbX1);
	c2r3D(vexbKy1, vexbY1);
    c2r3D(vexbKz1, vexbZ1);

   	fftw_complex *vdmexK;
	vdmexK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmexK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmeyK;
	vdmeyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmeyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmezK;
	vdmezK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmezK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	calc_diamag3D(dpedxK, dpedyK, dpedzK, B, B2, -1 * e,  neK, vdmexK, vdmeyK, vdmezK);
	fftw_complex *vdmixK;
	vdmixK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmixK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmiyK;
	vdmiyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmiyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmizK;
	vdmizK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmizK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	calc_diamag3D(dpidxK, dpidyK, dpidzK, B, B2, -1 * e,  neK, vdmixK, vdmiyK, vdmizK);
	
	double *vdmeX;
	vdmeX = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmeX, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmeY;
	vdmeY = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmeY, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmeZ;
	vdmeZ = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmeZ, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmiX;
	vdmiX = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmiX, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmiY;
	vdmiY = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmiY, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmiZ;
	vdmiZ = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmiZ, 42, (nx*ny*nz)* sizeof(double)); 

	c2r3D(vdmexK, vdmeX);
	c2r3D(vdmeyK, vdmeY);
    c2r3D(vdmezK, vdmeZ);

	c2r3D(vdmixK, vdmiX);
	c2r3D(vdmiyK, vdmiY);
    c2r3D(vdmizK, vdmiZ);

    
	fftw_complex *vdmex1K;
	vdmex1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmex1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmey1K;
	vdmey1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmey1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmez1K;
	vdmez1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmez1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	calc_diamag3D(dpedx1K, dpedy1K, dpedz1K, B, B2, -1 * e,  nepK, vdmex1K, vdmey1K, vdmez1K);

	double *vdmeX1;
	vdmeX1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmeX1, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmeY1;
	vdmeY1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmeY1, 42, (nx*ny*nz)* sizeof(double)); 

	double *vdmeZ1;
	vdmeZ1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmeZ1, 42, (nx*ny*nz)* sizeof(double)); 

	c2r3D(vdmex1K, vdmeX1); 
	c2r3D(vdmey1K, vdmeY1);
    c2r3D(vdmez1K, vdmeZ1);
	
	fftw_complex *vdmix1K;
	vdmix1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmix1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmiy1K;
	vdmiy1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmiy1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vdmiz1K;
	vdmiz1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vdmiz1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	calc_diamag3D(dpidx1K, dpidy1K, dpidz1K, B, B2, e,  nepK, vdmix1K, vdmiy1K, vdmiz1K);

	double *vdmiX1;
	vdmiX1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmiX1, 42, (nx*ny*nz)* sizeof(double));

	double *vdmiY1;
	vdmiY1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmiY1, 42, (nx*ny*nz)* sizeof(double));  

	double *vdmiZ1;
	vdmiZ1 = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(vdmiZ1, 42, (nx*ny*nz)* sizeof(double)); 

	c2r3D(vdmix1K, vdmiX1);
	c2r3D(vdmiy1K, vdmiY1);
	c2r3D(vdmiz1K, vdmiZ1); 

	std::vector<double> time(iter_max + 1,0.0); //This will create a vector of size iter_max + 1 all initialized to 0.0. You could use memset as well 

	c2r3D(neK, neTot);
	c2r3D(TeK, TeTot);	
	c2r3D(TiK, TiTot);
	c2r3D(phiK, phiTot);

	fftw_complex *ne1K_old;
	ne1K_old = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(ne1K_old, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *TeK_old;
	TeK_old = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(TeK_old, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *TiK_old;
	TiK_old = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(TiK_old, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *phi1K_old;
	phi1K_old = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(phi1K_old, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *residualnK;
	residualnK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(residualnK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vioxK;
	vioxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vioxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vioyK;
	vioyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(vioyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *viozK;
	viozK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(viozK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *residualtiK;
	residualtiK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(residualtiK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *veoxK;
	veoxK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(veoxK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *veoyK;
	veoyK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(veoyK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *veozK;
	veozK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(veozK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *residualteK;
	residualteK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(residualteK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *sourcen1K;
	sourcen1K = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(sourcen1K, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *residualK_phi;
	residualK_phi = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(residualK_phi, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *sourcetK;
	sourcetK = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	memset(sourcetK, 42, (nx*ny*nzk)* sizeof(fftw_complex));

	
    //begin time step loop here
	for (int iter = 0; iter < iter_max; iter++){ //for (int iter = 0; iter < iter_max; iter++){
		// total veloc
		//***************
		c2r3D(vexbKx, vexbX);
		c2r3D(vexbKy, vexbY);
        c2r3D(vexbKz, vexbZ);

		c2r3D(vdmexK, vdmeX);
		c2r3D(vdmeyK, vdmeY);
        c2r3D(vdmezK, vdmeZ);

		c2r3D(vdmixK, vdmiX);
		c2r3D(vdmiyK, vdmiY);
        c2r3D(vdmizK, vdmiZ);

	
        double dt =  calc_dt3D(u, vexbX,  vexbY, vexbZ, vdmiX, vdmiY,  vdmiZ,  vdmeX, vdmeY, vdmeZ, CFL,  kmax,  dt_max);

       time[iter + 1] = time[iter] + dt; 
		for(int i = 0; i< nx; i++){
		    for(int j = 0; j< ny; j++){
			    for(int k = 0; k< nzk; k++){ 
					for (int l = 0; l < ncomp; l++){
						ne1K_old[k + nzk * (j + ny * i)][l] = nepK[k + nzk * (j + ny * i)][l] ; // change to perturbed all of them: ne1k
						TeK_old[k + nzk * (j + ny * i)][l] = TeK[k + nzk * (j + ny * i)][l] ; //Tek and Ti are pertrubed
						TiK_old[k + nzk * (j + ny * i)][l] = TiK[k + nzk * (j + ny * i)][l] ;
						phi1K_old[k + nzk * (j + ny * i)][l] = phipK[k + nzk * (j + ny * i)][l] ; 
					}
				}
			}
		}

       
		 

        for (int stage = 0; stage < 4; stage++){ //for (int stage = 0; stage < 4; stage++){
            

			calc_residualn3D(vexbKx, vexbKy, vexbKz, neK, residualnK, KXX,  KYY, KZZ);
			calc_residualt3D(vioxK, vioyK,viozK,TiK, residualtiK, KXX,KYY, KZZ);
			calc_residualt3D(veoxK, veoyK,veozK,TeK, residualteK, KXX,KYY, KZZ);
            
			calcPotSourcek_inertia3D(nebgK, neK, dndx0K, dndy0K,dndz0K,dndxK , dndyK, dndzK ,dphidx0K , dphidy0K, dphidz0K, dphidx1K,dphidy1K,dphidz1K, pipK, pepK, uxB, e, Cm,  hallEK, hallIK , vexbKx0, vexbKy0, vexbKz0, vexbKx, vexbKy, vexbKz, KXX, KYY, KZZ, ksqu, potSourceK_inertia);
			phi_iter = Potentialk3D(invnK, dndxK, dndyK, dndzK, dphiKdt, potSourceK_inertia, KXX, KYY, KZZ, ninvksqu, err_max, phi_iter_max);
			
			if (phi_iter > phi_iter_max){
				printf("Blew up");
			break;
		    }
			//***************
			calc_sourcen3D( ksqu, nepK, Dart, sourcen1K);
			RK4(ne1K_old, dt, residualnK, sourcen1K, stage, nepK);
			RK4(phi1K_old, dt, residualK_phi, dphiKdt, stage, phipK);
			RK4(TiK_old, dt, residualtiK, sourcetK, stage, TiK); 
            RK4(TeK_old, dt, residualteK, sourcetK, stage, TeK); 
            
			// tot ne 
			//Arr3DArr3DAdd(ne0k,ne1k,nek); //nek = ne0k + ne1k;
			// phi
			//Arr3DArr3DAdd(phi0k,phi1k,phik); //phik = phi0k + phi1k;
			for(int i = 0; i < nx; i++){
				for( int j = 0; j < ny; j++){	
					for( int k = 0; k < nzk; k++){	
						for( int l = 0; l < ncomp; l++){
							neK[k + nzk * (j + ny * i)][l] = nebgK[k + nzk * (j + ny * i)][l] + nepK[k + nzk * (j + ny * i)][l];
							phiK[k + nzk * (j + ny * i)][l] = phibgK[k + nzk * (j + ny * i)][l] + phipK[k + nzk * (j + ny * i)][l];
						}
					}			
				}
			}
			Convolve3D(neK, TeK, peK);
			//Pek = kb * Pek; //rscalArr3DMult(Pek, kb, Pek); 
			for(int i = 0; i < nx; i++){
				for( int j = 0; j < ny; j++){
					for (int k = 0; k < nzk; k++){
						for (int l = 0; l < ncomp; l++){
							peK[k + nzk * (j + ny * i)][l] = peK[k + nzk * (j + ny * i)][l] * (kb);	
						}			
					}
				}
			}	
            Convolve3D(neK, TiK, piK);
            //Pik = kb * Pik;//rscalArr3DMult(Pik, kb, Pik);
			for(int i = 0; i < nx; i++){
				for( int j = 0; j < ny; j++){
					for (int k = 0; k < nzk; k++){
						for (int l = 0; l < ncomp; l++){
							piK[k + nzk * (j + ny * i)][l] = piK[k + nzk * (j + ny * i)][l] * (kb);	
						}			
					}
				}
			}
		    // pert pressure
			for(int i = 0; i < nx; i++){
				for( int j = 0; j < ny; j++){	
					for( int k = 0; k < nzk; k++){	
						for( int l = 0; l < ncomp; l++){
							pepK[k + nzk * (j + ny * i)][l] = peK[k + nzk * (j + ny * i)][l] + pebgK[k + nzk * (j + ny * i)][l];
							pipK[k + nzk * (j + ny * i)][l] = piK[k + nzk * (j + ny * i)][l] + pibgK[k + nzk * (j + ny * i)][l];
						}
					}			
				}
			}
			// all derv for tot
			Derivk3D(neK, KXX, dndxK);
			Derivk3D(neK, KYY, dndyK);
            Derivk3D(neK, KZZ, dndzK);

            Derivk3D(phiK, KXX, dphidxK);
			Derivk3D(phiK, KYY, dphidyK);
            Derivk3D(phiK, KZZ, dphidzK);

			Derivk3D(peK, KXX, dpedxK);
			Derivk3D(peK, KYY, dpedyK);
            Derivk3D(peK, KZZ, dpedzK);

			Derivk3D(piK, KXX, dpidxK);
			Derivk3D(piK, KYY, dpidyK); 
            Derivk3D(piK, KZZ, dpidzK);

            Derivk3D(phipK, KXX, dphidx1K);
			Derivk3D(phipK, KYY, dphidy1K);
            Derivk3D(phipK, KZZ, dphidz1K);

			CollFreqk_inertia3D(neK, TiK, TeK , kb, eps0, mi, me, ri, rn, nn, Oci, Oce, e, nuinK, nuieK, nuiiK, nuenK, nueeK, nueiK, isigPK, invnK, hallIK, hallEK);
			// tot velo
			calcV_ExBk3D(dphidxK, dphidyK,dphidzK, B, B2, vexbKx, vexbKy,vexbKz);
			
			calc_diamag3D(dpedxK, dpedyK, dpedzK,B, B2, -1 * e, neK, vdmexK, vdmeyK,vdmezK); 
			calc_diamag3D(dpidxK, dpidyK,dpidzK ,B, B2, e, neK, vdmixK, vdmiyK,vdmizK);


			// Get total velocity: use for loop or matrix operations
			for(int i = 0; i < nx; i++){
				for( int j = 0; j < ny; j++){	
					for( int k = 0; k < nzk; k++){	
						for( int l = 0; l < ncomp; l++){
							veoxK[k + nzk * (j + ny * i)][l] = vdmexK[k + nzk * (j + ny * i)][l] + vexbKx[k + nzk * (j + ny * i)][l];
							veoyK[k + nzk * (j + ny * i)][l] = vdmeyK[k + nzk * (j + ny * i)][l] + vexbKy[k + nzk * (j + ny * i)][l];
							veozK[k + nzk * (j + ny * i)][l] = vdmezK[k + nzk * (j + ny * i)][l] + vexbKz[k + nzk * (j + ny * i)][l];
							
							vioxK[k + nzk * (j + ny * i)][l] = vdmixK[k + nzk * (j + ny * i)][l] + vexbKx[k + nzk * (j + ny * i)][l];
							vioyK[k + nzk * (j + ny * i)][l] = vdmiyK[k + nzk * (j + ny * i)][l] + vexbKy[k + nzk * (j + ny * i)][l];
							viozK[k + nzk * (j + ny * i)][l] = vdmizK[k + nzk * (j + ny * i)][l] + vexbKz[k + nzk * (j + ny * i)][l];


						}
					}			
				}
			}
            
			
            
        } 
        
        if (phi_iter > phi_iter_max){
            printf("Blew up");
	        break;
	    }
            
         
        
        printf("Iteration = %d    t = %.10f   phi_iter = %d\n", iter, time[iter+1], phi_iter);
        //fprintf("Iteration = %d    t = %.10f   phi_iter = %d\n", iter, time[iter+1], phi_iter);
		if ((iter/saveFrequency) - saveNum == 0){
			c2r3D(neK, neTot);
			c2r3D(TeK, TeTot);	
			c2r3D(TiK, TiTot);
			c2r3D(phiK, phiTot);


                //print stuff
            char save [16] = {0};
			snprintf(save, 16,  "%d", saveNum);
            const char *type = ".txt";

            char nefilename[16] = {0};
            strcat(nefilename, "ne");
            

    
                
		    
  		   
		
            saveNum++;
			
        }
      
         
    }






	fftw_free(XX);
	fftw_free(YY);
	fftw_free(ZZ);
	fftw_free(kx);
	fftw_free(ky);
	fftw_free(kz);
	fftw_free(KXX);
	fftw_free(KYY);
	fftw_free(KZZ);
	fftw_free(ksqu);
	fftw_free(ninvksqu);
	fftw_free(Theta);
	fftw_free(Phiz);
	fftw_free(P0);
	fftw_free(nep);
	fftw_free(nebg);
	fftw_free(neTot);
	fftw_free(neK);
	//fftw_free(TestOut);
	fftw_free(Tip);
	fftw_free(Tep);
	fftw_free(Tibg);
	fftw_free(Tebg);
	fftw_free(TiTot);
	fftw_free(TeTot);
	fftw_free(phibg);
	fftw_free(phip);
	fftw_free(phiTot);
	fftw_free(piTot);
	fftw_free(peTot);
	fftw_free(pibg); 
	//fftw_free(pibg); 
    fftw_free(pebg); 
    fftw_free(pip); 
    fftw_free(pep); 
	fftw_free(TiK);
	fftw_free(TeK);
    fftw_free(phiK); 
    fftw_free(piK); 
    fftw_free(peK); 
    fftw_free(nebgK); 
    fftw_free(nepK); 
    fftw_free(TebgK); 
    fftw_free(TibgK); 
    fftw_free(phibgK); 
    fftw_free(phipK); 
    fftw_free(pibgK); 
    fftw_free(pipK); 
    fftw_free(pebgK); 
    fftw_free(pepK); 
	fftw_free(dndxK);
	fftw_free(dndyK);
    fftw_free(dndzK);
    fftw_free(dndx0K);
    fftw_free(dndy0K);
    fftw_free(dndz0K);
	fftw_free(dndx1K);
    fftw_free(dndy1K);
    fftw_free(dndz1K);
    fftw_free(dphidxK);
    fftw_free(dphidyK);
    fftw_free(dphidzK);
    fftw_free(dphidx0K);
    fftw_free(dphidx1K);
    fftw_free(dphidy0K);
    fftw_free(dphidy1K);
    fftw_free(dphidz0K);
    fftw_free(dphidz1K);
    fftw_free(d2phidxK);
    fftw_free(d2phidyK);
    fftw_free(d2phidzK);
    fftw_free(dpedxK);
    fftw_free(dpedx1K);
    fftw_free(dpedyK);
	fftw_free(dpedy1K);
	fftw_free(dpedzK);
	fftw_free(dpedz1K);
    fftw_free(dpidxK);
    fftw_free(dpidx1K);
    fftw_free(dpidyK);
    fftw_free(dpidy1K);
    fftw_free(dpidzK);
    fftw_free(dpidz1K);
    fftw_free(dphiKdt);
	fftw_free(vexbKx);
	fftw_free(vexbKy);
	fftw_free(vexbKz);
	fftw_free(vexbX);
	fftw_free(vexbY);
	fftw_free(vexbZ);
	fftw_free(vexbKx0);
	fftw_free(vexbKy0);
	fftw_free(vexbKz0);
	fftw_free(vexbX0);
	fftw_free(vexbY0);
	fftw_free(vexbZ0);
	fftw_free(vexbKx1);
	fftw_free(vexbKy1);
	fftw_free(vexbKz1);
	fftw_free(vexbX1);
	fftw_free(vexbY1);
	fftw_free(vexbZ1);
   	fftw_free(vdmexK);
	fftw_free(vdmeyK);
	fftw_free(vdmezK);
	fftw_free(hallEK); 
    fftw_free(hallIK); 
    fftw_free(nuinK); 
    fftw_free(nuieK); 
    fftw_free(nuiiK); 
    fftw_free(nuenK); 
    fftw_free(nueeK); 
    fftw_free(nueiK); 
    fftw_free(isigPK); 
    fftw_free(invnK); 
	fftw_free(potSourceK_inertia);
	fftw_free(vdmixK);
	fftw_free(vdmiyK);
	fftw_free(vdmizK);
	fftw_free(vdmeX);
	fftw_free(vdmeY);
	fftw_free(vdmeZ);
	fftw_free(vdmiX);
	fftw_free(vdmiY);
	fftw_free(vdmiZ);
	fftw_free(vdmex1K);
	fftw_free(vdmey1K);
	fftw_free(vdmez1K);
	fftw_free(vdmeX1);
	fftw_free(vdmeY1);
	fftw_free(vdmeZ1);
	fftw_free(vdmix1K);
	fftw_free(vdmiy1K);
	fftw_free(vdmiz1K);
	fftw_free(vdmiX1);
	fftw_free(vdmiY1);
	fftw_free(vdmiZ1);
	fftw_free(ne1K_old);
	fftw_free(TeK_old);

	fftw_free(TiK_old);

	fftw_free(phi1K_old);

	fftw_free(residualnK);

	fftw_free(vioxK);

	fftw_free(vioyK);

	fftw_free(viozK);

	fftw_free(residualtiK);

	fftw_free(veoxK);
	fftw_free(veoyK);
	fftw_free(veozK);

	fftw_free(residualteK);

	fftw_free(sourcen1K);

	fftw_free(residualK_phi);

	fftw_free(sourcetK);

	//Data allocated by fftw_malloc must be deallocated by fftw_free and not by the ordinary free.


}
double mag1DArray(double arr[]){
	return sqrt(arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2]);
}	 

void cross_product(double vector_a[], double vector_b[], double temp[]) {

   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];//   temp[1] = vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0];
   temp[2] = vector_a[1] * vector_b[0] - vector_a[0] * vector_b[1];// temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

}
double makeSpatialMesh3D(double dx, double dy, double dz,double xarr[], double yarr[],double zarr[]){
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){ 
            for(int l =0; l<nz; l++){
                xarr[l + nz * (j + ny * i)] = i*dx; //this index returns same indexing as Eigen
                yarr[l + nz * (j + ny * i)] = l*dy;
                zarr[l + nz * (j + ny * i)] = j*dz;
            }		
 			
		}		
	}
		
		if (dx < dy){
			return 1/(2*dx);
		}else {
			return 1/(2*dy);
		}
		

}
//void FourierMesh3D(double Lx,double Ly,double Lz, double KX[], double KY[],double KZ[], double ksqu[], double ninvksqu[], int FourierMeshType){
void FourierMesh3D(double Lx,double Ly,double Lz, double kx[], double ky[], double kz[], double KX[], double KY[],double KZ[], double ksqu[], double ninvksqu[], int FourierMeshType){
   
    
	double dx = Lx/nx; 
	double dy = Ly/ny; 
	double dz = Lz/nz;

    int k_counter = 0;		
	// Make kx. kx corresponds to modes [0:n/2-1 , -n/2:-1]. This is why there is an additional step, just due to FFTW3's structuring
	for(int i = 0; i < nx ; i++){ //change to ny
			if (i < nx/2){ 
				kx[i] = 2.*M_PI*i/Lx;
				ky[i] = 2.*M_PI*i/Ly;
                //eKZ(i) = 2.*EIGEN_PI*i/Lz;
				}
				if( i >= nx/2){ // i >= nx/2 --> from -128 to -1
				kx[i] = 2.* M_PI * (-i + 2.*k_counter) / Lx;//	2.* M_PI * (-i + 2.*k_counter) / Lx;		
				ky[i] = 2.*M_PI * (-i + 2.*k_counter) / Ly;
                //eKZ(i) = 2.* EIGEN_PI * (-i + 2.*k_counter) / Lz;
				//std::cout << eKX[i] << endl;
				}
				if( i >= nx/2){
					k_counter++;
				}
		
	    }
    
    int k_counterz = 0;
	
     for(int k = 0; k < nz ; k++){ 
			if (k < nz/2){ 
                kz[k] = 2.*M_PI*k/Lz;
				}
				if( k >= nz/2){ // i >= nx/2 --> from -128 to -1
                kz[k] = 2.*M_PI * (-k + 2.*k_counterz) / Lz;
				}
				if( k >= nz/2){
					k_counterz++;
				}
		
	}
	
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nzk; k++){ 
				//KX(i,j,k) = eKX(j); //
                KX[k + nzk * (j + ny * i)] = kx[j]; //kx[i];
                KY[k + nzk * (j + ny * i)] = ky[k];
                KZ[k + nzk * (j + ny * i)] = kz[i];
                
                //std::cout<<KX[k + nz * (j + ny * i)] <<endl;
                //KY(i,j,k) = eKY(i);
                //KZ(i,j,k) = eKZ(k);
			}

		}		
	}  
   
	
    
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nzk; k++){ 
				ksqu[k + nzk * (j + ny * i)] = ( sin(KX[k + nzk * (j + ny * i)] * dx/2)/(dx/2) * sin(KX[k + nzk * (j + ny * i)] * dx/2)/(dx/2) + sin(KY[k + nzk * (j + ny * i)] * dy/2)/(dy/2) * sin(KY[k + nzk * (j + ny * i)] * dy/2)/(dy/2) + sin(KZ[k + nzk * (j + ny * i)] * dz/2)/(dz/2) * sin(KZ[k + nzk * (j + ny * i)] * dz/2)/(dz/2) );
				ninvksqu[k + nzk * (j + ny * i)] = -1./ ksqu[k + nzk * (j + ny * i)];
                ninvksqu[0] = -1.;

			}
		}

	}		
	
	for(int i = 0; i< nx; i++){
            kx[i] = sin(kx[i]* dx)/dx ; 
			ky[i] = sin(ky[i]* dy)/dy ; 
    }
    for(int k = 0; k< nz; k++){
            kz[k] = sin(kz[k]* dz)/dz ; 
    }
  
    for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nzk; k++){ 
                //KX(i,j,k) = eKX(j); //
                KX[k + nzk * (j + ny * i)] = kx[j];
                KY[k + nzk * (j + ny * i)] = ky[k];
                KZ[k + nzk * (j + ny * i)] = kz[i];
                //std::cout<<KX[k + nz * (j + ny * i)]<<endl;
                //KY(i,j,k) = eKY(i);
                //KZ(i,j,k) = eKZ(k);
			}

		}		
	}
	
	/*int k_counter = 0;	
	// Make kx. kx corresponds to modes [0:n/2-1 , -n/2:-1]. This is why there is an additional step, just due to FFTW3's structuring
	for(int i = 0; i < nx ; i++){
		for(int j = 0; j < ny; j++){	
			for (int k = 0; k < nzk; ++k) {			
				if (i < nx/2){ // i<nx/2 --> 0 to 127
					KX[k + nzk * (j + ny * i)] =2.*M_PI*i/Lx; //2.*M_PI*i/Lx; 
					//K = 2pi/L* [0:nx/2-1	, -n/2:-1]' : creates a coloumn vector with nx/2 elements starting from 0 ( from 0 to 127 then -128 to -1) 256/2 is 128	
					KY[k + nzk * (j + ny * i)] =2*M_PI*i/Ly; //2*M_PI*j/Ly;	 
				}
				if( i >= nx/2){ // i >= nx/2 --> from -128 to -1
					KX[k + nzk * (j + ny * i)] =  2.* M_PI * (-i + 2.*k_counter) / Lx; //2.* M_PI * (-i + 2.*k_counter) / Lx;
					KY[k + nzk * (j + ny * i)] =  2.* M_PI * (-i + 2.*k_counter) / Ly;		
				}		
			}
		}
		if( i >= nx/2){
			k_counter++;
		}
	}
	

	// Make ky. Because we are using special FFTs for purely real valued functions, the last dimension (y in this case) only need n/2 + 1 points due to symmetry about the imaginary axis.
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){	
			for (int k = 0; k < nzk; ++k) {		
				KZ[k + nzk * (j + ny * i)] =2*M_PI*k/Ly;
			}						
		}	
	}
	double dx = Lx/nx; 
	double dy = Ly/ny; 
	double dz = Lz/nz;
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for (int k = 0; k < nzk; ++k) {	
				if ( FourierMeshType == 0){	
					ksqu[k + nzk * (j + ny * i)] = KX[k + nzk * (j + ny * i)] * KX[k + nzk * (j + ny * i)] + KY[k + nzk * (j + ny * i)] * KY[k + nzk * (j + ny * i)] + KZ[k + nzk * (j + ny * i)] * KZ[k + nzk * (j + ny * i)];						
					ninvksqu[k + nzk * (j + ny * i)] = -1 / (KX[k + nzk * (j + ny * i)] * KX[k + nzk * (j + ny * i)] + KY[k + nzk * (j + ny * i)] * KY[k + nzk * (j + ny * i)] + KZ[k + nzk * (j + ny * i)] * KZ[k + nzk * (j + ny * i)] );
				}
				if (FourierMeshType == 1){
					ksqu[k + nzk * (j + ny * i)] = (pow(sin(KX[k + nzk * (j + ny * i)] * dx/2)/(dx/2),2) + pow(sin(KY[k + nzk * (j + ny * i)] * dy/2)/(dy/2),2) + pow(sin(KZ[k + nzk * (j + ny * i)] * dz/2)/(dz/2),2)) ; // Use approximations for Kx, Ky, K^2
					ninvksqu[k + nzk * (j + ny * i)] = -1 /(pow(sin(KX[k + nzk * (j + ny * i)] * dx/2)/(dx/2),2) + pow(sin(KY[k + nzk * (j + ny * i)] * dy/2)/(dy/2),2) + pow(sin(KZ[k + nzk * (j + ny * i)] * dz/2)/(dz/2),2)) ;
					KX[k + nzk * (j + ny * i)] = sin(KX[k + nzk * (j + ny * i)]* dx)/dx ; //overwrite the exact
					KY[k + nzk * (j + ny * i)] = sin(KY[k + nzk * (j + ny * i)] *dy)/dy;
					KZ[k + nzk * (j + ny * i)] = sin(KZ[k + nzk * (j + ny * i)] *dz)/dz;  
				}
				
			}
			// For the "approximation" case. the discrt k is not periodic so causes issues 
		}	
	}
	
	// Account for the [0][0] point since there is a divide by 0 issue.
	ninvksqu[0] = -1.;*/
}
void r2c3D(double rArr[], double cArr[][ncomp]){
    //fftw_complex *input_array;
	//input_array = (fftw_complex*) fftw_malloc((nx*ny*nz)*sizeof(fftw_complex));
		
		
	//memcpy(input_array, rArr, (nx*ny*nz)*sizeof(fftw_complex));

	//fftw_plan forward = fftw_plan_dft_3d(nx, ny, nz, input_array, cArr, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan forward = fftw_plan_dft_r2c_3d(nx, ny, nz, &rArr[0], &cArr[0], FFTW_ESTIMATE);
	//fftw_plan fftw_plan_dft_r2c_3d(int n0, int n1, int n2,double *in, fftw_complex *out,unsigned flags);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
	fftw_cleanup();

    //fftw_free(input_array);

}
void c2r3D(double cArr[][ncomp], double rArr[]){
    fftw_complex *dummy;
	dummy = (fftw_complex*) fftw_malloc((nx*ny*nzk) * sizeof(fftw_complex));
	//to not overwrite data set up dummy
	for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        for(int l = 0; l<ncomp; ++l){
                            dummy[k + nzk * (j + ny * i)][l] = cArr[k + nzk * (j + ny * i)][l]; 
                        }
                        
                    }
                }
        }	
	
	//fftw_plan backward = fftw_plan_dft_3d(nx, ny, nz, cArr, output_array, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan backward = fftw_plan_dft_c2r_3d(nx, ny, nz, &dummy[0], &rArr[0], FFTW_ESTIMATE);

	fftw_execute(backward);
    fftw_destroy_plan(backward);
	fftw_cleanup();

	/*memcpy(rArr,output_array, (nx*ny*nz) * sizeof(double)); //size of double fixes mem error
		
		for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        rArr[k + nz * (j + ny * i)] = rArr[k + nz * (j + ny * i)]/(nx*ny*nz) ;
                    }
                }
        }*/
    fftw_free(dummy);

}
void Derivk3D(double vark[][ncomp], double K[], double derivative[][ncomp]){
	// Multiply by the corresponding k
    //Arr3DArr2DMult(vark, k, derivative);
    for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        for (int l = 0; l < ncomp; ++l){
                            derivative[k + nzk * (j + ny * i)][l] = vark[k + nzk * (j + ny * i)][l] * K[k + nzk * (j + ny * i)]; 


                        }
                        
        
            		}
        	    }
    }
   
	// Multiply by i
	//iArr3DMult(derivative, derivative);	
    double dummy;   // This is used so that data are not accidentally overwritten
	for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        dummy = derivative[k + nzk * (j + ny * i)][REAL];
                        derivative[k + nzk * (j + ny * i)][REAL] = -derivative[k + nzk * (j + ny * i)][IMAG];
                        derivative[k + nzk * (j + ny * i)][IMAG] = dummy;
            		}
        	    }
    }
    
}
void calcV_ExBk3D(double dphidxk[][ncomp], double dphidyk[][ncomp],double dphidzk[][ncomp], double B[], double B2, double vexbkx[][ncomp], double vexbky[][ncomp],double vexbkz[][ncomp]){

		for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        for (int l = 0; l < ncomp; ++l){
                            vexbkx[k + nzk * (j + ny * i)][l] = (dphidzk[k + nzk * (j + ny * i)][l] * B[1] - dphidyk[k + nzk * (j + ny * i)][l] * B[2]) /B2 ; 
							vexbky[k + nzk * (j + ny * i)][l]  =  (dphidxk[k + nzk * (j + ny * i)][l]  * B[2] - dphidzk[k + nzk * (j + ny * i)][l]  * B[0]) /B2;
							vexbkz[k + nzk * (j + ny * i)][l]  =  (dphidyk[k + nzk * (j + ny * i)][l]  * B[0] - dphidxk[k + nzk * (j + ny * i)][l]  * B[1]) /B2;


                        }
                        
        
            		}
        	    }
    }
	
}
void fourierDivision3D(double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]){
	double *f;
	f = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
		
	double *g;
	g = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
		
	double *fg;
	fg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
	// Take inverse ffts
	c2r3D(fk, f);
	c2r3D(gk, g);
	
	// Multiply in real space
	//Arr2DArr2DDiv(f, g, fg);
	for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
				fg[k + nz * (j + ny * i)] = f[k + nz * (j + ny * i)]/g[k + nz * (j + ny * i)];
			}
		}
	}
	
	r2c3D(fg, fgk); //here
	fftw_free(f); 
	fftw_free(g); 
	fftw_free(fg); 

}
void calc_diamag3D(double dpdxk[][ncomp], double dpdyk[][ncomp], double dpdzk[][ncomp],double B[], double B2, double qa, double nak[][ncomp], double diamagxk[][ncomp], double diamagyk[][ncomp], double diamagzk[][ncomp]){
		
		fftw_complex *predivx;
		predivx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

		fftw_complex *predivy;
		predivy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
		
		fftw_complex *predivz;
		predivz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
		
		for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nzk; ++k) {
                        for (int l = 0; l < ncomp; ++l){
							predivx[k + nzk * (j + ny * i)][l] = (dpdyk[k + nzk * (j + ny * i)][l] * B[2] - dpdzk[k + nzk * (j + ny * i)][l] * B[1])/(B2*qa*-1);
							predivy[k + nzk * (j + ny * i)][l] = (dpdzk[k + nzk * (j + ny * i)][l] * B[0] - dpdxk[k + nzk * (j + ny * i)][l] * B[2])/(B2*qa*-1);
							predivz[k + nzk * (j + ny * i)][l] = (dpdxk[k + nzk * (j + ny * i)][l] * B[1] - dpdyk[k + nzk * (j + ny * i)][l] * B[0])/(B2*qa*-1);
						}
		        }
		    }
	    }
		
	fourierDivision3D(predivx, nak, diamagxk);
	fourierDivision3D(predivy, nak, diamagyk);
	fourierDivision3D(predivz, nak, diamagzk);

	fftw_free(predivx); 
	fftw_free(predivy); 
	fftw_free(predivz);
}
void Arr2DArr3DMult(double arrIn0[], double arrIn1[], double arrOut[]){
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){
			for( int k = 0; k < nz; k++){
				arrOut[k + nz * (j + ny * i)] = arrIn0[k + nz * (j + ny * i)] * arrIn1[k + nz * (j + ny * i)];
			}
		}
	}
}
void Convolve3D( double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]){
	
	double *f;
	f = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));	
	double *g;
	g = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));	
	double *fg;
	fg = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
	// Take inverse ffts
	c2r3D(fk, f);
	c2r3D(gk, g);
	
	// Multiply in real space
	Arr2DArr3DMult(f, g, fg);
	
	// Take fft of fg
	r2c3D(fg, fgk);
	
	fftw_free(f);
	fftw_free(g);
	fftw_free(fg);

}
void CollFreqk_inertia3D(double nek[][ncomp], double Tik[][ncomp], double Tek[][ncomp], double kb, double eps0, double mi, double me, double ri, double rn, double nn, double Oci, double Oce, double e,double nuink[][ncomp], double nuiek[][ncomp] , double nuiik[][ncomp] , double nuenk[][ncomp] ,double nueek[][ncomp] ,double nueik[][ncomp]  , double isigPk[][ncomp] , double invnk[][ncomp] , double hallIk[][ncomp] ,double hallEk[][ncomp] ){
    
    double *ne; //pe1
	ne = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(ne, 42, (nx*ny*nz)*sizeof(double));

    double *Ti; //pe1
	Ti = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Ti, 42, (nx*ny*nz)*sizeof(double));

    double *Te; //pe1
	Te = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(Te, 42, (nx*ny*nz)*sizeof(double));

    c2r3D(nek,ne);
    c2r3D(Tik,Ti);
    c2r3D(Tek,Te);

	// Set scalar doubles for variables that are needed in the loops
	double Vthi, Vthe, lambdaD, Lambda;

	double *nuin; //pe1
	nuin = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nuin, 42, (nx*ny*nz)*sizeof(double));

    double *nuii; //pe1
	nuii = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nuii, 42, (nx*ny*nz)* sizeof(double));

    double *nuie; //pe1
	nuie = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nuie, 42, (nx*ny*nz)* sizeof(double));

    double *nuen; //pe1
	nuen = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nuen, 42, (nx*ny*nz)* sizeof(double));

    double *nuei; //pe1
	nuei = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nuei, 42, (nx*ny*nz)* sizeof(double));

    double *nuee; //pe1
	nuee = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(nuee, 42, (nx*ny*nz)* sizeof(double));

    double *isigP; //
	isigP = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(isigP, 42, (nx*ny*nz)* sizeof(double));

    double *invn; //
	invn = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(invn, 42, (nx*ny*nz)* sizeof(double));

    double *hallE; //
	hallE = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(hallE, 42, (nx*ny*nz)* sizeof(double));

    double *hallI; //
	hallI = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(hallI, 42, (nx*ny*nz)* sizeof(double));


	// Begin big loop to calculating everything.
	
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 

				Vthi = sqrt( 2. * kb * Ti[k + nz * (j + ny * i)] / mi);
				Vthe = sqrt( 2. * kb * Te[k + nz * (j + ny * i)] / me);
				nuin[k + nz * (j + ny * i)] =  nn * Vthi* M_PI * (ri+rn) * (ri+rn);
				nuen[k + nz * (j + ny * i)] = nn * Vthe * M_PI * rn * rn;
				// Calculate "inverse Pedersen conductivity"
				isigP[k + nz * (j + ny * i)] = 1.0 / (e * ( nuin[k + nz * (j + ny * i)]/ Oci + nuen[k + nz * (j + ny * i)]/ Oce) );
				// Calculate Debye length
				lambdaD = sqrt( eps0 * kb * Te[k + nz * (j + ny * i)]/ (ne[k + nz * (j + ny * i)] * e * e) );					
				// Calculate plasma parameter
				Lambda = 12.0 * M_PI * ne[k + nz * (j + ny * i)] * lambdaD * lambdaD * lambdaD;
				// Calculate electron-electron collision frequency
				nuee[k + nz * (j + ny * i)] = ne[k + nz * (j + ny * i)]* e * e *e *e * log(Lambda/3.0) / ( 2. * M_PI * eps0 * eps0 * me * me * Vthe * Vthe * Vthe);

				// Calculate ion-ion collision frequency
				nuii[k + nz * (j + ny * i)]= nuee[k + nz * (j + ny * i)] * sqrt(me/mi);
				// Calculate ion-electron collision frequency
				nuie[k + nz * (j + ny * i)] = nuee[k + nz * (j + ny * i)] * 0.5 * me / mi;
			
				// Calculate electron-ion collision frequency
				nuei[k + nz * (j + ny * i)] = nuee[k + nz * (j + ny * i)];
				// Calculate the inverse of the density
				// inverse of ne in Fourier space (which is needed for several terms in the temperature equation )	
				invn[k + nz * (j + ny * i)] = 1.0 / ne[k + nz * (j + ny * i)];

				//*************************	Hall parameters: TEST forb the inertial function ****************************
				hallE[k + nz * (j + ny * i)] = nuen[k + nz * (j + ny * i)]/Oce;
				hallI[k + nz * (j + ny * i)] = nuin[k + nz * (j + ny * i)]/Oci;
                


			}
		}
	}

    r2c3D(nuin, nuink);
    r2c3D(nuie, nuiek);
    r2c3D(nuii, nuiik);
    r2c3D(nuen, nuenk);
    r2c3D(nuei, nueik);
    r2c3D(nuee, nueek);
    r2c3D(invn, invnk);
    r2c3D(isigP, isigPk);
    r2c3D(hallE, hallEk);
    r2c3D(hallI, hallIk);

    fftw_free(ne);
    fftw_free(Ti);
    fftw_free(Te);
    fftw_free(nuin);
    fftw_free(nuii);
    fftw_free(nuie);
    fftw_free(nuen);
    fftw_free(nuei);
    fftw_free(nuee);
    fftw_free(isigP);
    fftw_free(invn);
    fftw_free(hallE);
    fftw_free(hallI);
}
double Arr1DMax(double arr[], int arrLength){
	double max = 0.0;
	for (int i = 0; i < arrLength; i++){
		if (arr[i] > max){
			max = arr[i];
		}
	}
	return max;
}
double max3Dk(double arr2D[]){ 
	double maxVal = 0.;
	for (int i = 0; i < nx; i++){
		for (int j = 0; j < ny; j++){//for (int j = 0; j < nyk; j++){
			for (int k = 0; k < nz; k++){
				if (arr2D[k + nz * (j + ny * i)] > maxVal){ //if (n)--> n isn't initialized here at
					//double *arr2D= new double [j+nyk*i]; //test this isn't initialized

					maxVal = arr2D[k + nz * (j + ny * i)];
				}
			}
		}
	}	
	return maxVal;
}
double calc_dt3D(double U[], double vexbx[], double vexby[],double vexbz[] ,double diamagxi[], double diamagyi[], double diamagzi[],double diamagxe[], double diamagye[], double diamagze[],double cfl, double kmax, double maxdt){
	//double bar = absolute(U);
		double *absArr;
	    absArr = (double*) fftw_malloc((nx*ny*nz) *sizeof(double));
		memset(absArr, 42, (nx*ny*nz)* sizeof(double));
		/*
		absArr = (double*) fftw_malloc(nx*nyk *sizeof(double));
		memset(absArr, 42, nx*nyk* sizeof(double));
		*/
	//	int N = 3;
		for (int i=0; i < sizee; i++){
			if (U[i] < 0){
				absArr[i] = abs(U[i]);
				
			}
		}

		//print2DPhysArray(absArr);

		double vMaxArr[10];
		
		vMaxArr[0] = max3Dk((vexbx));// max2D(vexbx); 
		vMaxArr[1] = max3Dk(vexby);
		vMaxArr[2] = max3Dk(vexbz);
		vMaxArr[3] = max3Dk(diamagxi);
		vMaxArr[4] = max3Dk(diamagyi);
		vMaxArr[5] = max3Dk(diamagzi);
		vMaxArr[6] = max3Dk(diamagxe);
		vMaxArr[7] = max3Dk(diamagye);
		vMaxArr[8] = max3Dk(diamagze);
		vMaxArr[9] = Arr1DMax(U,3); 
		double max = Arr1DMax(vMaxArr, 7);
		double dt = cfl / (max * kmax); // added 
	
		fftw_free(absArr);
		if (dt < maxdt){
		return dt;
		}else {
			return maxdt;
		}
		
}
void calc_residualn3D(double vexbxk[][ncomp], double vexbyk[][ncomp], double vexbzk[][ncomp], double nink[][ncomp], double residnoutk[][ncomp], double kx[], double ky[], double kz[]){
	fftw_complex *dninxk;
	dninxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dninyk;
	dninyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 

	fftw_complex *dninzk;
	dninzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
	
	fftw_complex *mult1;
	mult1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
	
	fftw_complex *mult2;
	mult2 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 

	fftw_complex *mult3;
	mult3 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
	
	Derivk3D(nink, kx, dninxk);
	Derivk3D(nink, ky, dninyk);
	Derivk3D(nink, kz, dninzk);
	
	Convolve3D(vexbxk, dninxk, mult1);
	Convolve3D(vexbyk, dninyk, mult2);
	Convolve3D(vexbzk, dninzk, mult3);
	
	//Arr3DArr3DAdd(mult1, mult2,mult3 ,residnoutk);
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){	
			for( int k = 0; k < nzk; k++){	
				for( int l = 0; l < ncomp; l++){
					residnoutk[k + nzk * (j + ny * i)][l] = mult1[k + nzk * (j + ny * i)][l] + mult2[k + nzk * (j + ny * i)][l] + mult3[k + nzk * (j + ny * i)][l];
				}
			}			
		}
	}	

	fftw_free(dninxk); //test no mult1/2 here 
	fftw_free(dninyk);
	fftw_free(dninzk);
	fftw_free(mult1); 
	fftw_free(mult2); 
	fftw_free(mult3); 

}
void calc_residualt3D(double voxk[][ncomp], double voyk[][ncomp],double vozk[][ncomp], double tempink[][ncomp], double tempoutk[][ncomp], double kx[], double ky[],double kz[]){
	
	fftw_complex *dtempinxk;
	dtempinxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *dtempinyk;
	dtempinyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dtempinzk;
	dtempinzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *dvoxk;
	dvoxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *dvoyk;
	dvoyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dvozk;
	dvozk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *divvo;
	divvo = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *mult1;
	mult1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *mult2;
	mult2 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *mult3;
	mult3 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *mult4;
	mult4 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	fftw_complex *dninxk;
	dninxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	Derivk3D(tempink, kx, dtempinxk);
	Derivk3D(tempink, ky, dtempinyk);
	Derivk3D(tempink, kz, dtempinzk);
	Derivk3D(voxk, kx, dvoxk);
	Derivk3D(voyk, ky, dvoyk);
	Derivk3D(vozk, kz, dvozk);

	//Arr3DArr3DAdd(dvoxk, dvoyk, divvo);
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){	
			for( int k = 0; k < nzk; k++){	
				for( int l = 0; l < ncomp; l++){
					divvo[k + nzk * (j + ny * i)][l] = dvoxk[k + nzk * (j + ny * i)][l] + dvoyk[k + nzk * (j + ny * i)][l] + dvozk[k + nzk * (j + ny * i)][l];
				}
			}			
		}
	}	
	
	Convolve3D(dtempinxk, voxk, mult1);
	Convolve3D(dtempinyk, voyk, mult2);
	Convolve3D(dtempinzk, vozk, mult4);
	Convolve3D(tempink, divvo, mult3);
	
	//rscalArr3DMult(mult3, 2/3, mult3);
	//void rscalArr3DMult(double arr[][ncomp], double rscal, double arrOut[][ncomp]){
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){
			for (int k = 0; k < nzk; k++){
				for (int l = 0; l < ncomp; l++){
					mult3[k + nzk * (j + ny * i)][l] = mult3[k + nzk * (j + ny * i)][l] * 2/3;	
				}			
			}
		}
	}	


	
	for (int i = 0; i < nx; i++){
		for (int j = 0; j < ny; j++){
			for (int k = 0; k < nzk; k++){
				for (int l = 0; l < ncomp; l++){
					tempoutk[k + nzk * (j + ny * i)][l] = mult1[k + nzk * (j + ny * i)][l] + mult2[k + nzk * (j + ny * i)][l] + mult3[k + nzk * (j + ny * i)][l] + mult4[k + nzk * (j + ny * i)][l] ;
				}
			}
		}
	}
	fftw_free(dtempinzk);
	fftw_free(dtempinxk); 
	fftw_free(dtempinyk);
	fftw_free(dvoxk);
	fftw_free(dvoyk);
	fftw_free(dvozk);
	fftw_free(divvo);
	fftw_free(mult1);
	fftw_free(mult2);
	fftw_free(mult3);
	fftw_free(mult4);
	fftw_free(dninxk);
	
}
double max_absComp3D(double arr3D[][ncomp]){
	// Take the absolute value
	double *absArr;
	absArr = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(absArr, 42, (nx*ny*nz)* sizeof(double)); //test here, if you pass 2D array to func decays to pointer and sizeof doesn't give size of array


	//absComp(arr3D, absArr);
	//void absComp(double arr3D[][ncomp], double arrOut[]){
	for (int i = 0; i < nx ; i++){
		for (int j = 0 ; j < ny; j++){
			for (int k = 0; k < nzk; k++){
				absArr[k + nzk * (j + ny * i)] = sqrt(arr3D[k + nzk * (j + ny * i)][0]*arr3D[k + nzk * (j + ny * i)][0] + arr3D[k + nzk * (j + ny * i)][1]*arr3D[k + nzk * (j + ny * i)][1]);
			}
		}
	}
	
	// Calculate the max value
	double maxVal = max3Dk(absArr); //by
	fftw_free(absArr);
	
	return maxVal;
	
}
void calcPotSourcek_inertia3D(double ne0k[][ncomp], double nek[][ncomp],double dndx0k[][ncomp], double dndy0k[][ncomp],double dndz0k[][ncomp], double dndxk [][ncomp], double dndyk [][ncomp],double dndzk [][ncomp], double dphidx0k [][ncomp], double dphidy0k [][ncomp], double dphidz0k [][ncomp], double dphidx1k [][ncomp], double dphidy1k [][ncomp], double dphidz1k [][ncomp],double Pi1k[][ncomp], double Pe1k[][ncomp], double uxB[], double e, double Cm, double hallEk [][ncomp], double hallIk [][ncomp], double vexbkx0[][ncomp], double vexbky0[][ncomp],double vexbkz0[][ncomp], double vexbkx[][ncomp],  double vexbky[][ncomp], double vexbkz[][ncomp],double kx[], double ky[], double kz[], double ksqu[], double potSourcek_inertia[][ncomp]){

	fftw_complex *d2Pi1k;
	d2Pi1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2Pe1k;
	d2Pe1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));


	Laplaciank3D(Pi1k, ksqu, d2Pi1k); 
	Laplaciank3D(Pe1k, ksqu, d2Pe1k); 
	
	fftw_complex *pikTerm;
	pikTerm = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *pekTerm;
	pekTerm = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nzk; k++){
				for(int l = 0; l < ncomp; l++){
					pikTerm[k + nzk * (j + ny * i)][l] = -1.*(Cm/e)* hallIk[k + nzk * (j + ny * i)][l]; 
					pekTerm[k + nzk * (j + ny * i)][l] = (Cm/e)*hallEk[k + nzk * (j + ny * i)][l];

				}
			}
		}		
	}
	

	
 	// First convolution for ions to multiply by Laplacian of pressure
	Convolve3D(d2Pi1k, pikTerm, pikTerm); 
	//convolve3D(d2Pi1k, pikTerm, pikTerm); 

 	// First convolution for electrons to multiply by Laplacian of pressure
	//convolve3D(d2Pe1k, pekTerm, pekTerm);
	Convolve3D(d2Pe1k, pekTerm, pekTerm);

	fftw_complex *ne1k;
	ne1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dndx1k;
	dndx1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dndy1k;
	dndy1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dndz1k;
	dndz1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dphidxk;
	dphidxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dphidyk;
	dphidyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dphidzk;
	dphidzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbkx1;
	vexbkx1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbky1;
	vexbky1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *vexbkz1;
	vexbkz1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nzk; k++){
				for(int l = 0; l < ncomp; l++){
					ne1k[k + nzk * (j + ny * i)][l] =  nek[k + nzk * (j + ny * i)][l] - ne0k[k + nzk * (j + ny * i)][l]; 
					dndx1k[k + nzk * (j + ny * i)][l] = dndxk[k + nzk * (j + ny * i)][l] - dndx0k[k + nzk * (j + ny * i)][l];
					dndy1k[k + nzk * (j + ny * i)][l] = dndyk[k + nzk * (j + ny * i)][l] - dndy0k[k + nzk * (j + ny * i)][l];
					dndz1k[k + nzk * (j + ny * i)][l] = dndzk[k + nzk * (j + ny * i)][l] - dndz0k[k + nzk * (j + ny * i)][l];
					dphidxk[k + nzk * (j + ny * i)][l] = dphidx0k[k + nzk * (j + ny * i)][l] + dphidx1k[k + nzk * (j + ny * i)][l];
					dphidyk[k + nzk * (j + ny * i)][l] = dphidy0k[k + nzk * (j + ny * i)][l] + dphidy1k[k + nzk * (j + ny * i)][l];
					dphidzk[k + nzk * (j + ny * i)][l] = dphidz0k[k + nzk * (j + ny * i)][l] + dphidz1k[k + nzk * (j + ny * i)][l];
					vexbkx1[k + nzk * (j + ny * i)][l] = vexbkx[k + nzk * (j + ny * i)][l] - vexbkx0[k + nzk * (j + ny * i)][l];
					vexbky1[k + nzk * (j + ny * i)][l] = vexbky[k + nzk * (j + ny * i)][l] - vexbky0[k + nzk * (j + ny * i)][l];
					vexbkz1[k + nzk * (j + ny * i)][l] = vexbkz[k + nzk * (j + ny * i)][l] - vexbkz0[k + nzk * (j + ny * i)][l];

				}
			}
		}		
	}

	
	fftw_complex *div_n_nabphi_x;
	div_n_nabphi_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_n_nabphi_y;
	div_n_nabphi_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_n_nabphi_z;
	div_n_nabphi_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	
  // divergence (n nabla phi) term:
	fftw_complex *div_dummy_x1;
	div_dummy_x1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_y1;
	div_dummy_y1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_z1;
	div_dummy_z1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_x0;
	div_dummy_x0 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_y0;
	div_dummy_y0 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_z0;
	div_dummy_z0 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummyTerm;
	div_dummyTerm = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

 	//**********************************************************
	fftw_complex *div_dummy_dx1;
	div_dummy_dx1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_dy1;
	div_dummy_dy1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_dz1;
	div_dummy_dz1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_dx0;
	div_dummy_dx0 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_dy0;
	div_dummy_dy0 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *div_dummy_dz0;
	div_dummy_dz0 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *coeff;
	coeff = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	
 //Do the convolution of nek and dphi1k , to find the sum of final div dummy term seperate them like this?
    Convolve3D(nek, dphidx1k, div_dummy_x1);
	Convolve3D(nek, dphidy1k, div_dummy_y1);
	Convolve3D(nek, dphidz1k, div_dummy_z1); 

	Convolve3D(ne1k, dphidx0k, div_dummy_x0);
	Convolve3D(ne1k, dphidy0k, div_dummy_y0);
	Convolve3D(ne1k, dphidz0k, div_dummy_z0); 
	


 //Take derivatives of dummy variable:
	Derivk3D(div_dummy_x1, kx, div_dummy_dx1);  
	Derivk3D(div_dummy_y1, ky, div_dummy_dy1);
	Derivk3D(div_dummy_z1, kz, div_dummy_dz1);

	Derivk3D(div_dummy_x0, kx, div_dummy_dx0); 
	Derivk3D(div_dummy_y0, ky, div_dummy_dy0);
	Derivk3D(div_dummy_z0, kz, div_dummy_dz0);


	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nzk; k++){
				for(int l = 0; l < ncomp; l++){
					div_dummyTerm[k + nzk * (j + ny * i)][l] = div_dummy_x1[k + nzk * (j + ny * i)][l] + div_dummy_y1[k + nzk * (j + ny * i)][l] + div_dummy_z1[k + nzk * (j + ny * i)][l] + div_dummy_x0[k + nzk * (j + ny * i)][l] + div_dummy_y0[k + nzk * (j + ny * i)][l] + div_dummy_z0[k + nzk * (j + ny * i)][l] ;
					coeff[k + nzk * (j + ny * i)][l] = -1.*(Cm)*(hallEk[k + nzk * (j + ny * i)][l]+ hallIk[k + nzk * (j + ny * i)][l]);
					div_n_nabphi_x[k + nzk * (j + ny * i)][l] = div_dummy_dx1[k + nzk * (j + ny * i)][l] + div_dummy_dx0[k + nzk * (j + ny * i)][l]; 
					div_n_nabphi_y[k + nzk * (j + ny * i)][l] = div_dummy_dy1[k + nzk * (j + ny * i)][l] + div_dummy_dy0[k + nzk * (j + ny * i)][l];
					div_n_nabphi_z[k + nzk * (j + ny * i)][l] = div_dummy_dz1[k + nzk * (j + ny * i)][l] + div_dummy_dz0[k + nzk * (j + ny * i)][l];
				


				}
			}
		}		
	}
	Convolve3D(div_n_nabphi_x,coeff, div_n_nabphi_x); 
	Convolve3D(div_n_nabphi_y, coeff, div_n_nabphi_y); 
	Convolve3D(div_n_nabphi_z, coeff, div_n_nabphi_z); 

	fftw_complex *coeff1;
	coeff1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *uxBTerm;
	uxBTerm = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nzk; k++){
				for(int l = 0; l < ncomp; l++){
					uxBTerm[k + nzk * (j + ny * i)][l] = uxB[0] * dndx1k[k + nzk * (j + ny * i)][l] + uxB[1] * dndy1k[k + nzk * (j + ny * i)][l] + uxB[2] * dndz1k[k + nzk * (j + ny * i)][l];
					coeff1[k + nzk * (j + ny * i)][l] = (Cm)*(hallEk[k + nzk * (j + ny * i)][l]+ hallIk[k + nzk * (j + ny * i)][l]);

				}
			}
		}		
	}
	Convolve3D(uxBTerm, coeff1, uxBTerm);



	fftw_complex *v_nabnab_phiTerm;
	v_nabnab_phiTerm = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_x;
	nabnab_dummyTerm_n0v0phi1_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_y;
	nabnab_dummyTerm_n0v0phi1_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_z;
	nabnab_dummyTerm_n0v0phi1_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_yz;
	nabnab_dummyTerm_n0v0phi1_yz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1;
	nabnab_dummyTerm_n0v0phi1 = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	 
	fftw_complex *nabnab_dummyTerm_n1vphi_x;
	nabnab_dummyTerm_n1vphi_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_y;
	nabnab_dummyTerm_n1vphi_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_z;
	nabnab_dummyTerm_n1vphi_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_yz;
	nabnab_dummyTerm_n1vphi_yz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_xy;
	nabnab_dummyTerm_n1vphi_xy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_xz;
	nabnab_dummyTerm_n1vphi_xz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_yx;
	nabnab_dummyTerm_n1vphi_yx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_zx;
	nabnab_dummyTerm_n1vphi_zx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_zy;
	nabnab_dummyTerm_n1vphi_zy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi;
	nabnab_dummyTerm_n1vphi = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_x;
	nabnab_dummyTerm_n0v1phi_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_y;
	nabnab_dummyTerm_n0v1phi_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_z;
	nabnab_dummyTerm_n0v1phi_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi;
	nabnab_dummyTerm_n0v1phi = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));


	fftw_complex *d2phidx1k;
	d2phidx1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidy1k;
	d2phidy1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidz1k;
	d2phidz1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidxk;
	d2phidxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidyk;
	d2phidyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidzk;
	d2phidzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidxyk;
	d2phidxyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidxzk;
	d2phidxzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidzyk;
	d2phidzyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidyxk;
	d2phidyxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidyzk;
	d2phidyzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidzxk; 
	d2phidzxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidxy1k;
	d2phidxy1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidzx1k;
	d2phidzx1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidzy1k;
	d2phidzy1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidxz1k;
	d2phidxz1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidyx1k;
	d2phidyx1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *d2phidyz1k;
	d2phidyz1k = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));
	
	fftw_complex *nabnab_dummyTerm_n0v0phi1_xy;
	nabnab_dummyTerm_n0v0phi1_xy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_yx;
	nabnab_dummyTerm_n0v0phi1_yx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_xz;
	nabnab_dummyTerm_n0v0phi1_xz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_zy;
	nabnab_dummyTerm_n0v0phi1_zy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_zx;
	nabnab_dummyTerm_n0v0phi1_zx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_xy;
	nabnab_dummyTerm_n0v1phi_xy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_xz;
	nabnab_dummyTerm_n0v1phi_xz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_zy;
	nabnab_dummyTerm_n0v1phi_zy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_yx;
	nabnab_dummyTerm_n0v1phi_yx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_yz;
	nabnab_dummyTerm_n0v1phi_yz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_zx;
	nabnab_dummyTerm_n0v1phi_zx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_yy;
	nabnab_dummyTerm_n0v0phi1_yy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v0phi1_zz;
	nabnab_dummyTerm_n0v0phi1_zz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_yy;
	nabnab_dummyTerm_n1vphi_yy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n1vphi_zz;
	nabnab_dummyTerm_n1vphi_zz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_zz;
	nabnab_dummyTerm_n0v1phi_zz = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *nabnab_dummyTerm_n0v1phi_yy;
	nabnab_dummyTerm_n0v1phi_yy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

 //Derivatives:
	Derivk3D(dphidx1k, kx, d2phidx1k); 
	Derivk3D(dphidy1k, ky, d2phidy1k);
	Derivk3D(dphidz1k, kz, d2phidz1k);
	
	Derivk3D(dphidz1k, ky, d2phidzy1k); 
	Derivk3D(dphidz1k, kx, d2phidzx1k);

	//cross deriv of dphidxk
	Derivk3D(dphidx1k, ky, d2phidxy1k); 
	Derivk3D(dphidy1k, kx, d2phidyx1k);
	Derivk3D(dphidx1k, kz, d2phidxz1k);
	Derivk3D(dphidy1k, kz, d2phidyz1k);


 // convolutions:
	Convolve3D(vexbkx0, d2phidx1k, nabnab_dummyTerm_n0v0phi1_x);
	Convolve3D(vexbky0, d2phidy1k, nabnab_dummyTerm_n0v0phi1_y);
	Convolve3D(vexbkz0, d2phidyz1k, nabnab_dummyTerm_n0v0phi1_yz);


	//cross terms convo
	Convolve3D(vexbky0, d2phidxy1k, nabnab_dummyTerm_n0v0phi1_xy);
	Convolve3D(vexbkx0, d2phidyx1k, nabnab_dummyTerm_n0v0phi1_yx);
	Convolve3D(vexbkz0, d2phidxz1k,nabnab_dummyTerm_n0v0phi1_xz); 


	//cross terms convo
	Convolve3D(vexbky0, d2phidzy1k, nabnab_dummyTerm_n0v0phi1_zy);
	Convolve3D(vexbkx0, d2phidzx1k, nabnab_dummyTerm_n0v0phi1_zx);
	Convolve3D(vexbkz0, d2phidz1k,nabnab_dummyTerm_n0v0phi1_z); 

	
 	// Calculate other n1vphi term:
 	//Derivatives:
	Derivk3D(dphidxk, kx, d2phidxk); 
	Derivk3D(dphidyk, ky, d2phidyk);
	Derivk3D(dphidzk, kz, d2phidzk);
	

	//cross deriv of dphidxk
	Derivk3D(dphidxk, ky, d2phidxyk); 
	Derivk3D(dphidxk, kz, d2phidxzk);
	Derivk3D(dphidyk, kx, d2phidyxk);
	Derivk3D(dphidyk, kz, d2phidyzk);//
	Derivk3D(dphidzk, ky, d2phidzyk); 
	Derivk3D(dphidzk, kx, d2phidzxk);

	
 // convolutions:
	Convolve3D(vexbkx, d2phidxk, nabnab_dummyTerm_n1vphi_x);
	Convolve3D(vexbky, d2phidyk, nabnab_dummyTerm_n1vphi_y);
	Convolve3D(vexbkz, d2phidzk, nabnab_dummyTerm_n1vphi_z);

	
	//cross terms convo
	Convolve3D(vexbky, d2phidxyk, nabnab_dummyTerm_n1vphi_xy);
	Convolve3D(vexbkx, d2phidyxk, nabnab_dummyTerm_n1vphi_yx);
	Convolve3D(vexbkz, d2phidyzk, nabnab_dummyTerm_n1vphi_yz);
	Convolve3D(vexbkz, d2phidxzk, nabnab_dummyTerm_n1vphi_xz);
	Convolve3D(vexbky, d2phidzyk, nabnab_dummyTerm_n1vphi_zy);
	Convolve3D(vexbkx, d2phidzxk, nabnab_dummyTerm_n1vphi_zx);
	

 //  Calculate other n0v1phi term:
 // convolutions:
	Convolve3D(vexbkx1, d2phidxk, nabnab_dummyTerm_n0v1phi_x);
	Convolve3D(vexbky1, d2phidyk, nabnab_dummyTerm_n0v1phi_y);


	// cross terms convo 
	Convolve3D(vexbky1, d2phidxyk, nabnab_dummyTerm_n0v1phi_xy);
	Convolve3D(vexbkx1, d2phidyxk, nabnab_dummyTerm_n0v1phi_yx);
	Convolve3D(vexbkz1, d2phidxzk, nabnab_dummyTerm_n0v1phi_xz);
	Convolve3D(vexbkz1, d2phidyzk, 	nabnab_dummyTerm_n0v1phi_yz);

	
	Convolve3D(vexbkz1, d2phidzk,nabnab_dummyTerm_n0v1phi_z);
	Convolve3D(vexbky1, d2phidzyk, nabnab_dummyTerm_n0v1phi_zy);
	Convolve3D(vexbkx1, d2phidzxk, nabnab_dummyTerm_n0v1phi_zx);
	
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nzk; k++){
				for(int l = 0; l < ncomp; l++){
					nabnab_dummyTerm_n0v0phi1[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v0phi1_x[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_xy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_xz[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n0v0phi1_yy[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v0phi1_yx[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_y[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_yz[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n0v0phi1_zz[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v0phi1_zx[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_zy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_z[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n1vphi[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n1vphi_x[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n1vphi_xy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n1vphi_xz[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n1vphi_yy[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n1vphi_yx[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n1vphi_y[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n1vphi_yz[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n1vphi_zz[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n1vphi_zx[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n1vphi_zy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n1vphi_z[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n0v1phi[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v1phi_x[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v1phi_xy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v1phi_xz[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n0v1phi_yy[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v1phi_yx[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v1phi_y[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v1phi_yz[k + nzk * (j + ny * i)][l] );
					nabnab_dummyTerm_n0v1phi_zz[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v1phi_zx[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v1phi_zy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v1phi_z[k + nzk * (j + ny * i)][l] );
				
				}
			}
		}		
	}
	
 // Multiply by n0:
	Convolve3D(nabnab_dummyTerm_n0v0phi1, ne0k, nabnab_dummyTerm_n0v0phi1);
	Convolve3D(nabnab_dummyTerm_n0v0phi1_yy, ne0k, nabnab_dummyTerm_n0v0phi1_yy);
	Convolve3D(nabnab_dummyTerm_n0v0phi1_zz, ne0k, nabnab_dummyTerm_n0v0phi1_zz);


 // Multiply by n1:
	Convolve3D(nabnab_dummyTerm_n1vphi, ne1k, nabnab_dummyTerm_n1vphi);
	Convolve3D(nabnab_dummyTerm_n1vphi_yy, ne1k, nabnab_dummyTerm_n1vphi_yy);
	Convolve3D(nabnab_dummyTerm_n1vphi_zz, ne1k, nabnab_dummyTerm_n1vphi_zz);


	

 // Multiply by n0:
	Convolve3D(ne0k,nabnab_dummyTerm_n0v1phi ,nabnab_dummyTerm_n0v1phi);
	Convolve3D(nabnab_dummyTerm_n0v1phi_yy, ne0k, nabnab_dummyTerm_n0v1phi_yy);
	Convolve3D(nabnab_dummyTerm_n0v1phi_zz, ne0k, nabnab_dummyTerm_n0v1phi_zz);


	fftw_complex *dnabnab_dummyTerm_n0v0phi1_x;
	dnabnab_dummyTerm_n0v0phi1_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v0phi1_y;
	dnabnab_dummyTerm_n0v0phi1_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v0phi1_z;
	dnabnab_dummyTerm_n0v0phi1_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n1vphi_z;
	dnabnab_dummyTerm_n1vphi_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n1vphi_x;
	dnabnab_dummyTerm_n1vphi_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n1vphi_y;
	dnabnab_dummyTerm_n1vphi_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v1phi_x;
	dnabnab_dummyTerm_n0v1phi_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v1phi_y;
	dnabnab_dummyTerm_n0v1phi_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v1phi_z;
	dnabnab_dummyTerm_n0v1phi_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v0phi1_xy;
	dnabnab_dummyTerm_n0v0phi1_xy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v0phi1_yx;
	dnabnab_dummyTerm_n0v0phi1_yx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n1vphi_xy;
	dnabnab_dummyTerm_n1vphi_xy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n1vphi_yx;
	dnabnab_dummyTerm_n1vphi_yx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v1phi_xy;
	dnabnab_dummyTerm_n0v1phi_xy = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dnabnab_dummyTerm_n0v1phi_yx;
	dnabnab_dummyTerm_n0v1phi_yx = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

   // try nabnab_dummyTerm_n0v0phi1
   	Derivk3D(nabnab_dummyTerm_n0v0phi1, kx, dnabnab_dummyTerm_n0v0phi1_x);
	Derivk3D(nabnab_dummyTerm_n0v0phi1_yy, ky, dnabnab_dummyTerm_n0v0phi1_y);
	Derivk3D(nabnab_dummyTerm_n0v0phi1_zz, kz, dnabnab_dummyTerm_n0v0phi1_z);

	
	//cross terms
	Derivk3D(nabnab_dummyTerm_n1vphi, kx, dnabnab_dummyTerm_n1vphi_x);
	Derivk3D(nabnab_dummyTerm_n1vphi_yy, ky, dnabnab_dummyTerm_n1vphi_y);
	Derivk3D(nabnab_dummyTerm_n1vphi_zz, kz, dnabnab_dummyTerm_n1vphi_z);

	
	// cross terms
	Derivk3D(nabnab_dummyTerm_n0v1phi, kx, dnabnab_dummyTerm_n0v1phi_x); // it was kx
	
	Derivk3D(nabnab_dummyTerm_n0v1phi_yy, ky, dnabnab_dummyTerm_n0v1phi_y);
	Derivk3D(nabnab_dummyTerm_n0v1phi_zz, ky, dnabnab_dummyTerm_n0v1phi_z);

	
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nzk; k++){
				for(int l = 0; l < ncomp; l++){
					nabnab_dummyTerm_n0v0phi1[k + nzk * (j + ny * i)][l] = (nabnab_dummyTerm_n0v0phi1_x[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_xy[k + nzk * (j + ny * i)][l] ) + ( nabnab_dummyTerm_n0v0phi1_xz[k + nzk * (j + ny * i)][l] );
					v_nabnab_phiTerm[k + nzk * (j + ny * i)][l] = - (dnabnab_dummyTerm_n0v0phi1_x[k + nzk * (j + ny * i)][l] +  dnabnab_dummyTerm_n0v0phi1_y[k + nzk * (j + ny * i)][l]  +  dnabnab_dummyTerm_n0v0phi1_z[k + nzk * (j + ny * i)][l]
					+ dnabnab_dummyTerm_n1vphi_x[k + nzk * (j + ny * i)][l] +  dnabnab_dummyTerm_n1vphi_y[k + nzk * (j + ny * i)][l]  +  dnabnab_dummyTerm_n1vphi_z[k + nzk * (j + ny * i)][l] 
					+ dnabnab_dummyTerm_n0v1phi_x[k + nzk * (j + ny * i)][l] +  dnabnab_dummyTerm_n0v1phi_y[k + nzk * (j + ny * i)][l]  +  dnabnab_dummyTerm_n0v1phi_z[k + nzk * (j + ny * i)][l] );
				
					potSourcek_inertia[k + nzk * (j + ny * i)][l] = v_nabnab_phiTerm[k + nzk * (j + ny * i)][l] + div_n_nabphi_x[k + nzk * (j + ny * i)][l] + div_n_nabphi_y[k + nzk * (j + ny * i)][l] +  pikTerm[k + nzk * (j + ny * i)][l] + pekTerm[k + nzk * (j + ny * i)][l] + uxBTerm[k + nzk * (j + ny * i)][l];

				}
			}
		}		
	}

	fftw_free(d2Pi1k);
	fftw_free(d2Pe1k);
	fftw_free(pikTerm);
	fftw_free(pekTerm);
	fftw_free(ne1k);
	fftw_free(dndx1k);
	fftw_free(dndy1k);
	fftw_free(dndz1k);
	fftw_free(dphidxk);
	fftw_free(dphidyk);
	fftw_free(dphidzk);
	fftw_free(vexbkx1);
	fftw_free(vexbky1);
	fftw_free(vexbkz1);
	fftw_free(div_n_nabphi_x);
	fftw_free(div_n_nabphi_y);
	fftw_free(div_n_nabphi_z);
	fftw_free(div_dummy_x1);
	//fftw_free(div_n_nabphi_y);
	//fftw_free(div_n_nabphi_z);
	//fftw_free(div_dummy_x1);
	fftw_free(div_dummy_y1);
	fftw_free(div_dummy_z1);
	fftw_free(div_dummy_x0);
	fftw_free(div_dummy_y0);
	fftw_free(div_dummy_z0);
	fftw_free(div_dummyTerm);
	fftw_free(div_dummy_dx1);
	fftw_free(div_dummy_dy1);
	fftw_free(div_dummy_dz1);
	fftw_free(div_dummy_dx0);
	fftw_free(div_dummy_dy0);
	fftw_free(div_dummy_dz0);
	fftw_free(coeff);
	fftw_free(coeff1);
	fftw_free(uxBTerm);
	fftw_free(v_nabnab_phiTerm);
	fftw_free(nabnab_dummyTerm_n0v0phi1_x);
	fftw_free(nabnab_dummyTerm_n0v0phi1_y);
	fftw_free(nabnab_dummyTerm_n0v0phi1_z);
	fftw_free(nabnab_dummyTerm_n0v0phi1_yz);
	fftw_free(nabnab_dummyTerm_n0v0phi1);
	fftw_free(nabnab_dummyTerm_n1vphi_x);
	fftw_free(nabnab_dummyTerm_n1vphi_y);
	fftw_free(nabnab_dummyTerm_n1vphi_z);
	fftw_free(nabnab_dummyTerm_n1vphi_yz);
	fftw_free(nabnab_dummyTerm_n1vphi_xy);
	fftw_free(nabnab_dummyTerm_n1vphi_xz);
	fftw_free(nabnab_dummyTerm_n1vphi_yx);
	fftw_free(nabnab_dummyTerm_n1vphi_zx);
	fftw_free(nabnab_dummyTerm_n1vphi_zy);
	fftw_free(nabnab_dummyTerm_n1vphi);
	fftw_free(nabnab_dummyTerm_n0v1phi_x);
	fftw_free(nabnab_dummyTerm_n0v1phi_y);
	fftw_free(nabnab_dummyTerm_n0v1phi_z);
	fftw_free(nabnab_dummyTerm_n0v1phi);
	fftw_free(d2phidx1k);

	fftw_free(d2phidy1k);
	fftw_free(d2phidz1k);
	fftw_free(d2phidxk);
	fftw_free(d2phidyk);
	fftw_free(d2phidzk);
	fftw_free(d2phidxyk);
	fftw_free(d2phidxzk);
	fftw_free(d2phidzyk);
	fftw_free(d2phidyxk);
	fftw_free(d2phidyzk);
	fftw_free(d2phidzxk);
	fftw_free(d2phidxy1k);
	fftw_free(d2phidzx1k);
	fftw_free(d2phidzy1k);
	fftw_free(d2phidxz1k);
	fftw_free(d2phidyx1k);
	fftw_free(d2phidyz1k);

	fftw_free(nabnab_dummyTerm_n0v0phi1_xy);
	fftw_free(nabnab_dummyTerm_n0v0phi1_yx);

	fftw_free(nabnab_dummyTerm_n0v0phi1_xz);
	fftw_free(nabnab_dummyTerm_n0v0phi1_zy);
	fftw_free(nabnab_dummyTerm_n0v0phi1_zx);
	fftw_free(nabnab_dummyTerm_n0v1phi_xy);

	fftw_free(nabnab_dummyTerm_n0v1phi_xz);
	fftw_free(nabnab_dummyTerm_n0v1phi_zy);
	fftw_free(nabnab_dummyTerm_n0v1phi_yx);

	fftw_free(nabnab_dummyTerm_n0v1phi_yz);
	fftw_free(nabnab_dummyTerm_n0v1phi_zx);
	fftw_free(nabnab_dummyTerm_n0v0phi1_yy);
	fftw_free(nabnab_dummyTerm_n0v0phi1_zz);
	fftw_free(nabnab_dummyTerm_n1vphi_yy);
	fftw_free(nabnab_dummyTerm_n1vphi_zz);
	fftw_free(nabnab_dummyTerm_n0v1phi_zz);
	fftw_free(nabnab_dummyTerm_n0v1phi_yy);

	fftw_free(dnabnab_dummyTerm_n0v0phi1_x);
	fftw_free(dnabnab_dummyTerm_n0v0phi1_y);
	fftw_free(dnabnab_dummyTerm_n0v0phi1_z);
	fftw_free(dnabnab_dummyTerm_n1vphi_z);
	fftw_free(dnabnab_dummyTerm_n1vphi_x);
	fftw_free(dnabnab_dummyTerm_n1vphi_y);
	fftw_free(dnabnab_dummyTerm_n0v1phi_x);

	fftw_free(dnabnab_dummyTerm_n0v1phi_y);
	fftw_free(dnabnab_dummyTerm_n0v1phi_z);
	fftw_free(dnabnab_dummyTerm_n0v0phi1_xy);
	fftw_free(dnabnab_dummyTerm_n0v0phi1_yx);
	fftw_free(dnabnab_dummyTerm_n1vphi_xy);
	fftw_free(dnabnab_dummyTerm_n1vphi_yx);
	fftw_free(dnabnab_dummyTerm_n0v1phi_xy);
	fftw_free(dnabnab_dummyTerm_n0v1phi_yx);




}
int Potentialk3D(double invnk[][ncomp], double dndxk[][ncomp], double dndyk[][ncomp], double dndzk[][ncomp],double phik[][ncomp], double potSourcek[][ncomp], double kx[], double ky[], double kz[],double ninvksqu[], double err_max, int max_iter){
	 // Initialize variables used in the function
	fftw_complex *dphidxk;
	dphidxk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
	
	fftw_complex *dphidyk;
	dphidyk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));

	fftw_complex *dphidzk;
	dphidzk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex));  
	
	fftw_complex *gradNgradPhi_x;
	gradNgradPhi_x = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
	
	fftw_complex *gradNgradPhi_y;
	gradNgradPhi_y = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 

	fftw_complex *gradNgradPhi_z;
	gradNgradPhi_z = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 
	
	fftw_complex *RHS;
	RHS = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 

	
	double phik_max, phik_max_old;
	double it_error;

	// Begin counter for the number of iterations it takes to calculate phi
	int count = 0;
	
	// Begin while loop
	do{			
		// Calculate phi derivatives
		Derivk3D(phik, kx, dphidxk);
		Derivk3D(phik, ky, dphidyk);
		Derivk3D(phik, kz, dphidzk);		
		// Do convolutions for grad n dot grad phi term
		// gradNgradphi = [ikx phik] * [ikx nk] - [iky phik] * [iky nk], where * is the convolution
		Convolve3D(dndxk, dphidxk, gradNgradPhi_x);
		Convolve3D(dndyk, dphidyk, gradNgradPhi_y);
		Convolve3D(dndzk, dphidzk, gradNgradPhi_z);
		// Subtract gradNgradphi from the source term. Calculate RHS of the equation:
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
				for (int k = 0; k < nzk; k++){
					for (int l = 0; l < ncomp; l++){
					RHS[k + nzk * (j + ny * i)][l]= potSourcek[k + nzk * (j + ny * i)][l] - gradNgradPhi_x[k + nzk * (j + ny * i)][l] - gradNgradPhi_y[k + nzk * (j + ny * i)][l] - gradNgradPhi_z[k + nzk * (j + ny * i)][l]  ;
					}
				}
			}
		}
		
		// Convolve RHS with invnk
		Convolve3D(RHS, invnk, RHS);			
		
		// Calculate maximum of absolute value of previous phi
		phik_max_old = max_absComp3D(phik);		
		
		// Multiply by ninvksqu to get the updated phi: this will output potential in Fourier space, phik
		//Arr3DArr2DMult(RHS, ninvksqu, phik);
    	for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				for (int k = 0; k < nzk; ++k) {
					for (int l = 0; l < ncomp; ++l){
						phik[k + nzk * (j + ny * i)][l] = RHS[k + nzk * (j + ny * i)][l] * ninvksqu[k + nzk * (j + ny * i)]; 


                    }
                        
        
            	}
        	}
    	}

		// Calculate maximum of absolute value of updated phi(new phi)
		phik_max = max_absComp3D(phik); //by
		//cout<<phik<<endl;
		// Increase iteration count by 1
		count = count + 1;
		
		// Calculate error
		it_error = fabs((phik_max-phik_max_old)/phik_max);	//err_max is the error we want to converge to	
		// If error is too high and we haven't reached the max iterations yet, repeat iterations
	 }while( it_error > err_max && count  <= max_iter && phik_max > err_max); // and instead &&
	 
	
	fftw_free(dphidxk);
	fftw_free(dphidyk);
	fftw_free(dphidzk);
	fftw_free(gradNgradPhi_x);
	fftw_free(gradNgradPhi_y);
	fftw_free(gradNgradPhi_z);
	fftw_free(RHS);
	 return count;

	 
}
void Laplaciank3D( double vark[][ncomp], double ksqu[], double derivative[][ncomp]){	
	//Arr3DArr2DMult(vark, k, derivative);
    for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (int k = 0; k < nzk; ++k) {
				for (int l = 0; l < ncomp; ++l){
					derivative[k + nzk * (j + ny * i)][l] = vark[k + nzk * (j + ny * i)][l] * ksqu[k + nzk * (j + ny * i)]; 


                }
                        
        
            }
        }
    }
	//rscalArr3DMult(derivative, -1., derivative);
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){
			for (int k = 0; k < nzk; k++){
				for (int l = 0; l < ncomp; l++){
					derivative[k + nzk * (j + ny * i)][l] = derivative[k + nzk * (j + ny * i)][l] * (-1.);	
				}			
			}
		}
	}
}	
void calc_sourcen3D(double ksqu[], double nk[][ncomp], double d, double sourcenk[][ncomp]){
	fftw_complex *lapnk;
	lapnk = (fftw_complex*) fftw_malloc((nx*ny*nzk)* sizeof(fftw_complex)); 

	Laplaciank3D(nk, ksqu, lapnk);
	//rscalArr3DMult(lapnk, d, sourcenk);
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){
			for (int k = 0; k < nzk; k++){
				for (int l = 0; l < ncomp; l++){
					sourcenk[k + nzk * (j + ny * i)][l] = lapnk[k + nzk * (j + ny * i)][l] * (d);	
				}			
			}
		}
	}

	fftw_free(lapnk);
    

}
void RK4(double f[][ncomp], double dt, double residual[][ncomp], double source[][ncomp], int stage, double fout[][ncomp]){
double alpha[4] = {1./4, 1./3, 1./2, 1.};

	for(int i = 0; i < nx; i++){
		for( int j = 0; j < ny; j++){
			for (int k = 0; k < nzk; k++){
				for (int l = 0; l < ncomp; l++){
					fout[k + nzk * (j + ny * i)][l] =  f[k + nzk * (j + ny * i)][l] - (alpha[stage] * dt * (residual[k + nzk * (j + ny * i)][l] - source[k + nzk * (j + ny * i)][l]));
				}			
			}
		}
	}
		
}
