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

static const int nx = 128;
static const int ny = 128; 
static const int nz = 128;
static const int ncomp = 2;

using namespace std;
void cross_product(double vector_a[], double vector_b[], double temp[]);
double mag1DArray(double arr[]);
double Arr1DMax(double arr[], int arrLength);
void Arr2DArr2DDiv(double arrIn0[], double arrIn1[], double arrOut[]);
void fourierDivision3D(double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]);
void Arr2DArr3DMult(double arrIn0[], double arrIn1[], double arrOut[]);

void FourierMesh3D(double Lx,double Ly,double Lz, double kx[], double ky[], double kz[], double KX[], double KY[],double KZ[], double ksqu[], double ninvksqu[], int FourierMeshType);
double makeSpatialMesh3D(double dx, double dy, double dz,double xarr[], double yarr[],double zarr[]);
void c2r3D(double cArr[][ncomp], double rArr[]);
void r2c3D(double rArr[], double cArr[][ncomp]);
void Derivk3D(double vark[][ncomp], double K[], double derivative[][ncomp]);
void Convolve3D( double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]);
void CollFreqk_inertia3D(double nek[][ncomp], double Tik[][ncomp], double Tek[][ncomp], double kb, double eps0, double mi, double me, double ri, double rn, double nn, double Oci, double Oce, double e,double nuink[][ncomp], double nuiek[][ncomp] , double nuiik[][ncomp] , double nuenk[][ncomp] ,double nueek[][ncomp] ,double nueik[][ncomp]  , double isigPk[][ncomp] , double invnk[][ncomp] , double hallIk[][ncomp] ,double hallEk[][ncomp] );
//void calcPotSourcek_inertia3D(double ne0k[][ncomp], double nek[][ncomp],double dndx0k[][ncomp], double dndy0k[][ncomp],double dndz0k[][ncomp], double dndxk [][ncomp], double dndyk [][ncomp],double dndzk [][ncomp], double dphidx0k [][ncomp], double dphidy0k [][ncomp], double dphidz0k [][ncomp], double dphidx1k [][ncomp], double dphidy1k [][ncomp], double dphidz1k [][ncomp],double Pi1k[][ncomp], double Pe1k[][ncomp], double uxB[], double e, double Cm, double hallEk [][ncomp], double hallIk [][ncomp], double vexbkx0[][ncomp], double vexbky0[][ncomp],double vexbkz0[][ncomp], double vexbkx[][ncomp],  double vexbky[][ncomp], double vexbkz[][ncomp],double kx[], double ky[], double kz[], double ksqu[], double potSourcek_inertia[][ncomp]);
//int Potentialk3D(double invnk[][ncomp], double dndxk[][ncomp], double dndyk[][ncomp], double dndzk[][ncomp],double phik[][ncomp], double potSourcek[][ncomp], double kx[], double ky[], double kz[],double ninvksqu[], double err_max, int max_iter);
void calcV_ExBk3D(double dphidxk[][ncomp], double dphidyk[][ncomp],double dphidzk[][ncomp], double B[], double B2, double vexbkx[][ncomp], double vexbky[][ncomp],double vexbkz[][ncomp]);
void calc_diamag3D(double dpdxk[][ncomp], double dpdyk[][ncomp], double dpdzk[][ncomp],double B[], double B2, double qa, double nak[][ncomp], double diamagxk[][ncomp], double diamagyk[][ncomp],double diamagzk[][ncomp]);
void fourierDivision3D(double fk[][ncomp], double gk[][ncomp], double fgk[][ncomp]);
//double calc_dt3D(double U[], double vexbx[], double vexby[],double vexbz[] ,double diamagxi[], double diamagyi[], double diamagzi[],double diamagxe[], double diamagye[], double diamagze[],double cfl, double kmax, double maxdt);
//double max3Dk(double arr2D[]);
//double max_absComp3D(double arr3D[][ncomp]);

//g++ MemoryLeakTest.cpp -lfftw3 -lm -ffast-math -fno-math-errno -march=native -Ofast -o ../../../Datasim1/test.out 
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
    double rx0 = 100.;//20.; ///4.;
    double ry0 =  100.;//20.; ///4.;
    double rz0 =  100.;//20.;///4.;
    double uxB[sizee];
	double Bmag = mag1DArray(B); 
    double B2 = Bmag * Bmag;
   cross_product(u, B, uxB);
    double Oci = e*Bmag/mi;
	double Oce = e*Bmag/me;
	double lgN_right = 3e4;;// convert to km
	double lgN_left = 3e4; //3e4;
	double Cm = 1./(1 / Oci + 1/Oce);
	double m = 1.;
    double dx = Lx / nx;
	double dy = Ly / ny;
    double dz = Lz / nz;
    double o = 2.5; //outer region
	double a = 4.; //for plasma enhancement case a =4 and for plasma depletion case a = -0.05
	double d = 1.; 
	double p = 3.; //was a2 //aL
	double b = 0.015; //was d2
    double c = 12.; 
    double nuin = 0.1;
	double nuen = 0.;
	double nn = 1e14; 
    

		
	double *XX;
	XX = (double*) fftw_malloc((nx*ny*nz) *sizeof(double));
	memset(XX, 42, (nx*ny*nz) * sizeof(double)); 

    double *YY;
	YY = (double*) fftw_malloc((nx*ny*nz) *sizeof(double));
	memset(YY, 42, (nx*ny*nz) * sizeof(double)); 
	
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
   
    double *KXX;
	KXX = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(KXX, 42, (nx*ny*nz)* sizeof(double)); 
	double *KYY;
	KYY = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(KYY, 42, (nx*ny*nz)* sizeof(double)); 

    double *KZZ;
    KZZ = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(KZZ, 42, (nx*ny*nz)* sizeof(double)); 

	double *ksqu;
	ksqu = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(ksqu, 42, (nx*ny*nz)* sizeof(double)); //test
	
	double *ninvksqu;
	ninvksqu = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	memset(ninvksqu, 42, (nx*ny*nz)* sizeof(double)); 

	int FourierMeshType = 1;
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
    //Eigen::Tensor<double, 3> ne(nx,ny,nz); 
    double n0 = 1e11;
	double Ampl_1 = -0.07;//-0.2; //-1E-3; // min
	double Ampl_2 = 0.07; //0.2;//1E-3; // max 
	//ttest for rand
	unsigned int seed = time(NULL);

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
	neK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(neK, 42, (nx*ny*nz)* sizeof(fftw_complex));	

    fftw_complex *TiK; //nek
	TiK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(TiK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *TeK; //nek
	TeK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(TeK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *phiK; //nek
	phiK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(phiK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *piK; //nek
	piK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(piK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *peK; //nek
	peK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(peK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *nebgK; //nek
	nebgK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nebgK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *nepK; //nek
	nepK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nepK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *TebgK; //nek
	TebgK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(TebgK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *TibgK; //nek
	TibgK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(TibgK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *phibgK; //nek
	phibgK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(phibgK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *phipK; //nek
	phipK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(phipK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *pibgK; //nek
	pibgK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(pibgK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *pipK; //nek
	pipK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(pipK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *pebgK; //nek
	pebgK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(pebgK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *pepK; //nek
	pepK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(pepK, 42, (nx*ny*nz)* sizeof(fftw_complex));
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
	hallEK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(hallEK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *hallIK; //nek
	hallIK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(hallIK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *nuinK; //nek
	nuinK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nuinK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *nuieK; //nek
	nuieK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nuieK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *nuiiK; //nek
	nuiiK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nuiiK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *nuenK; //nek
	nuenK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nuenK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *nueeK; //nek
	nueeK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nueeK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *nueiK; //nek
	nueiK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(nueiK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *isigPK; //nek
	isigPK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(isigPK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *invnK; //nek
	invnK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(invnK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    CollFreqk_inertia3D( neK, TiK, TeK, kb, eps0, mi, me, ri, rn, nn, Oci, Oce, e, nuinK, nuieK , nuiiK , nuenK ,nueeK ,nueiK , isigPK , invnK , hallIK , hallEK );
  
	fftw_complex *potSourceK_inertia;
	potSourceK_inertia = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(potSourceK_inertia, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *dndxK;
	dndxK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndxK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dndyK;
	dndyK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndyK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dndzK;
	dndzK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndzK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dndx0K;
	dndx0K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndx0K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dndy0K;
	dndy0K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndy0K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dndz0K;
	dndz0K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndz0K, 42, (nx*ny*nz)* sizeof(fftw_complex));
	fftw_complex *dndx1K;
	dndx1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndx1K, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *dndy1K;
	dndy1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndy1K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dndz1K;
	dndz1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dndz1K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dphidxK;
	dphidxK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidxK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dphidyK;
	dphidyK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidyK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dphidzK;
	dphidzK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidzK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dphidx0K;
	dphidx0K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidx0K, 42, (nx*ny*nz)* sizeof(fftw_complex));


    fftw_complex *dphidx1K;
	dphidx1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidx1K, 42, (nx*ny*nz)* sizeof(fftw_complex));
	fftw_complex *dphidy0K;
	dphidy0K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidy0K, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *dphidy1K;
	dphidy1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidy1K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dphidz0K;
	dphidz0K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidz0K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dphidz1K;
	dphidz1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dphidz1K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *d2phidxK;
	d2phidxK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(d2phidxK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *d2phidyK;
	d2phidyK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(d2phidyK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *d2phidzK;
	d2phidzK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(d2phidzK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dpedxK;
	dpedxK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpedxK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dpedx1K;
	dpedx1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpedx1K, 42, (nx*ny*nz)* sizeof(fftw_complex));
 
    fftw_complex *dpedyK;
	dpedyK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpedyK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *dpedy1K;
	dpedy1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpedy1K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dpedzK;
	dpedzK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpedzK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dpedz1K;
	dpedz1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpedz1K, 42, (nx*ny*nz)* sizeof(fftw_complex));
 
    fftw_complex *dpidxK;
	dpidxK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpidxK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *dpidx1K;
	dpidx1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpidx1K, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dpidyK;
	dpidyK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpidyK, 42, (nx*ny*nz)* sizeof(fftw_complex));

    fftw_complex *dpidy1K;
	dpidy1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpidy1K, 42, (nx*ny*nz)* sizeof(fftw_complex));
    
    fftw_complex *dpidzK;
	dpidzK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpidzK, 42, (nx*ny*nz)* sizeof(fftw_complex));
	
    fftw_complex *dpidz1K;
	dpidz1K = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(dpidz1K, 42, (nx*ny*nz)* sizeof(fftw_complex));
    
    fftw_complex *dphiKdt;
	dphiKdt = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        for(int l = 0; l<ncomp; ++l){
                            dphiKdt[k + nz * (j + ny * i)][l] = 1.; 
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
	fftw_complex *vexbKx;
	vexbKx = (fftw_complex*) fftw_malloc((nx*ny*nz)*sizeof(fftw_complex));
	memset(vexbKx, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vexbKy;
	vexbKy = (fftw_complex*) fftw_malloc((nx*ny*nz)* sizeof(fftw_complex));
	memset(vexbKy, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vexbKz;
	vexbKz = (fftw_complex*) fftw_malloc((nx*ny*nz)* sizeof(fftw_complex));
	memset(vexbKz, 42, (nx*ny*nz)* sizeof(fftw_complex));
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
	vexbKx0 = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vexbKx0, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vexbKy0;
	vexbKy0 = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vexbKy0, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vexbKz0;
	vexbKz0 = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vexbKz0, 42, (nx*ny*nz)* sizeof(fftw_complex));


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
	vexbKx1 = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vexbKx1, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vexbKy1;
	vexbKy1 = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vexbKy1, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vexbKz1;
	vexbKz1 = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vexbKz1, 42, (nx*ny*nz)* sizeof(fftw_complex));


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
	vdmexK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vdmexK, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vdmeyK;
	vdmeyK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vdmeyK, 42, (nx*ny*nz)* sizeof(fftw_complex));

	fftw_complex *vdmezK;
	vdmezK = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
	memset(vdmezK, 42, (nx*ny*nz)* sizeof(fftw_complex));

	calc_diamag3D(dpedxK, dpedyK, dpedzK, B, B2, -1 * e,  neK, vdmexK, vdmeyK, vdmezK);
	


    

	//free(XX);
	//free(YY);
	//free(ZZ);
	//free(fk);
	//free(neK);
	//free(neTot);
	//free(gk);
	//free(fgk);


}
void cross_product(double vector_a[], double vector_b[], double temp[]) {

   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];//   temp[1] = vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0];
   temp[2] = vector_a[1] * vector_b[0] - vector_a[0] * vector_b[1];// temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

}
double mag1DArray(double arr[]){
	return sqrt(arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2]);
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
			for(int k = 0; k< nz; k++){ 
				//KX(i,j,k) = eKX(j); //
                KX[k + nz * (j + ny * i)] = kx[j]; //kx[i];
                KY[k + nz * (j + ny * i)] = ky[k];
                KZ[k + nz * (j + ny * i)] = kz[i];
                
                //std::cout<<KX[k + nz * (j + ny * i)] <<endl;
                //KY(i,j,k) = eKY(i);
                //KZ(i,j,k) = eKZ(k);
			}

		}		
	}  
   
	
    
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				ksqu[k + nz * (j + ny * i)] = ( sin(KX[k + nz * (j + ny * i)] * dx/2)/(dx/2) * sin(KX[k + nz * (j + ny * i)] * dx/2)/(dx/2) + sin(KY[k + nz * (j + ny * i)] * dy/2)/(dy/2) * sin(KY[k + nz * (j + ny * i)] * dy/2)/(dy/2) + sin(KZ[k + nz * (j + ny * i)] * dz/2)/(dz/2) * sin(KZ[k + nz * (j + ny * i)] * dz/2)/(dz/2) );
                //ksqu(i,j,k) = ( pow(sin(KX(i,j,k) * dx/2)/(dx/2),2) + pow(sin(KY(i,j,k) * dy/2)/(dy/2),2) + pow(sin(KZ(i,j,k) * dz/2)/(dz/2),2) );

				ninvksqu[k + nz * (j + ny * i)] = -1./ ksqu[k + nz * (j + ny * i)];
                ninvksqu[0] = -1.;
                //std::cout<<ninvksqu[k + nz * (j + ny * i)] <<endl;

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
			for(int k = 0; k< nz; k++){ 
                //KX(i,j,k) = eKX(j); //
                KX[k + nz * (j + ny * i)] = kx[j];
                KY[k + nz * (j + ny * i)] = ky[k];
                KZ[k + nz * (j + ny * i)] = kz[i];
              
			}

		}		
	}
}
double makeSpatialMesh3D(double dx, double dy, double dz,double xarr[], double yarr[],double zarr[]){
	// Make x coordinates. We want them to change as we change the zeroth index and be constant as we change the first index.
	
	
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){ 
            for(int l =0; l<nz; l++){
                xarr[l + nz * (j + ny * i)] = i*dx; //this index returns same indexing as Eigen
                //[l + nz * (j + ny * i)] = (k,i,j)
                yarr[l + nz * (j + ny * i)] = l*dy;//j*dy;
                //YY[j + nz * (i + ny * l)] = l*dy;
                zarr[l + nz * (j + ny * i)] = j*dz; //l*dz;
                //std::cout<<ZZ[l + nz * (j + ny * i)]<<endl;
            }		
 			
		}		
	}
		
		if (dx < dy){
			return 1/(2*dx);
		}else {
			return 1/(2*dy);
		}
		

}
void Derivk3D(double vark[][ncomp], double K[], double derivative[][ncomp]){
	// Multiply by the corresponding k
    //Arr3DArr2DMult(vark, k, derivative);
    for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        for (int l = 0; l < ncomp; ++l){
                            derivative[k + nz * (j + ny * i)][l] = vark[k + nz * (j + ny * i)][l] * K[k + nz * (j + ny * i)]; 


                        }
                        
        
            		}
        	    }
    }
   

	// Multiply by i
	//iArr3DMult(derivative, derivative);	
    double dummy;   // This is used so that data are not accidentally overwritten
	for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        dummy = derivative[k + nz * (j + ny * i)][REAL];
                        derivative[k + nz * (j + ny * i)][REAL] = -derivative[k + nz * (j + ny * i)][IMAG];
                        derivative[k + nz * (j + ny * i)][IMAG] = dummy;
            		}
        	    }
    }
    
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

    double *test; //
	test = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));

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

    free(ne);
    free(Ti);
    free(Te);
    free(nuin);
    free(nuii);
    free(nuie);
    free(nuen);
    free(nuei);
    free(nuee);
    free(isigP);
    free(invn);
    free(hallE);
    free(hallI);
}
void calcV_ExBk3D(double dphidxk[][ncomp], double dphidyk[][ncomp],double dphidzk[][ncomp], double B[], double B2, double vexbkx[][ncomp], double vexbky[][ncomp],double vexbkz[][ncomp]){

		for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        for (int l = 0; l < ncomp; ++l){
                            vexbkx[k + nz * (j + ny * i)][l] = (dphidzk[k + nz * (j + ny * i)][l] * B[1] - dphidyk[k + nz * (j + ny * i)][l] * B[2]) /B2 ; 
							vexbky[k + nz * (j + ny * i)][l]  =  (dphidxk[k + nz * (j + ny * i)][l]  * B[2] - dphidzk[k + nz * (j + ny * i)][l]  * B[0]) /B2;
							vexbkz[k + nz * (j + ny * i)][l]  =  (dphidyk[k + nz * (j + ny * i)][l]  * B[0] - dphidxk[k + nz * (j + ny * i)][l]  * B[1]) /B2;


                        }
                        
        
            		}
        	    }
    }
	
}
void calc_diamag3D(double dpdxk[][ncomp], double dpdyk[][ncomp], double dpdzk[][ncomp],double B[], double B2, double qa, double nak[][ncomp], double diamagxk[][ncomp], double diamagyk[][ncomp], double diamagzk[][ncomp]){
		
		fftw_complex *predivx;
		predivx = (fftw_complex*) fftw_malloc((nx*ny*nz)* sizeof(fftw_complex));

		fftw_complex *predivy;
		predivy = (fftw_complex*) fftw_malloc((nx*ny*nz)* sizeof(fftw_complex));
		
		fftw_complex *predivz;
		predivz = (fftw_complex*) fftw_malloc((nx*ny*nz)* sizeof(fftw_complex)); 
		
		for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        for (int l = 0; l < ncomp; ++l){
							predivx[k + nz * (j + ny * i)][l] = (dpdyk[k + nz * (j + ny * i)][l] * B[2] - dpdzk[k + nz * (j + ny * i)][l] * B[1])/(B2*qa*-1);
							predivy[k + nz * (j + ny * i)][l] = (dpdzk[k + nz * (j + ny * i)][l] * B[0] - dpdxk[k + nz * (j + ny * i)][l] * B[2])/(B2*qa*-1);
							predivz[k + nz * (j + ny * i)][l] = (dpdxk[k + nz * (j + ny * i)][l] * B[1] - dpdyk[k + nz * (j + ny * i)][l] * B[0])/(B2*qa*-1);
						}
		        }
		    }
	    }
		
	fourierDivision3D(predivx, nak, diamagxk);
	fourierDivision3D(predivy, nak, diamagyk);
	fourierDivision3D(predivz, nak, diamagzk);
	/*double *fpredivx;
	fpredivx = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
		
	double *na;
	na = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
		
	double *diamagx;
	diamagx = (double*) fftw_malloc((nx*ny*nz)*sizeof(double));
	
	// Take inverse ffts
	c2r3D(predivx, fpredivx);
	c2r3D(nak, na);
	
	// Multiply in real space
	//Arr2DArr2DDiv(f, g, fg);
	for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
						diamagx[k + nz * (j + ny * i)] = fpredivx[k + nz * (j + ny * i)]/na[k + nz * (j + ny * i)];
					}
				}
	}
	
	// Take fft of fg
	r2c3D(diamagx, diamagxk); //here
	fftw_free(fpredivx); 
	fftw_free(na); 
	fftw_free(diamagx); */
	

	fftw_free(predivx); 
	fftw_free(predivy); 
	fftw_free(predivz);
}
void r2c3D(double rArr[], double cArr[][ncomp]){
    fftw_complex *input_array;
	input_array = (fftw_complex*) fftw_malloc((nx*ny*nz)*sizeof(fftw_complex));
		
		
	memcpy(input_array, rArr, (nx*ny*nz)*sizeof(fftw_complex));

	fftw_plan forward = fftw_plan_dft_3d(nx, ny, nz, input_array, cArr, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(forward);
    fftw_destroy_plan(forward);
	fftw_cleanup();

    fftw_free(input_array);

}
void c2r3D(double cArr[][ncomp], double rArr[]){
    fftw_complex *output_array;
	output_array = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
		
	
	fftw_plan backward = fftw_plan_dft_3d(nx, ny, nz, cArr, output_array, FFTW_BACKWARD, FFTW_ESTIMATE);
		
	fftw_execute(backward);
    fftw_destroy_plan(backward);
	fftw_cleanup();

	memcpy(rArr,output_array, (nx*ny*nz) * sizeof(double)); //size of double fixes mem error
		
		for (int i = 0; i < nx; ++i) {
        		for (int j = 0; j < ny; ++j) {
            		for (int k = 0; k < nz; ++k) {
                        rArr[k + nz * (j + ny * i)] = rArr[k + nz * (j + ny * i)]/(nx*ny*nz) ;
                    }
                }
        }
    fftw_free(output_array);

}
void Arr2DArr2DDiv(double arrIn0[], double arrIn1[], double arrOut[]){
	//Arr2DArr2DDiv(f, g, fg);
	for (int i = 0; i < nx; ++i) {
        	for (int j = 0; j < ny; ++j) {
            	for (int k = 0; k < nz; ++k) {
					arrOut[k + nz * (j + ny * i)] = arrIn0[k + nz * (j + ny * i)]/arrIn1[k + nz * (j + ny * i)];
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
	Arr2DArr2DDiv(f, g, fg);
	/*for (int i = 0; i < nx; ++i) {
        	for (int j = 0; j < ny; ++j) {
            	for (int k = 0; k < nz; ++k) {
					fg[k + nz * (j + ny * i)] = f[k + nz * (j + ny * i)]/g[k + nz * (j + ny * i)];
				}
			}
		}*/
	
	// Take fft of fg
	r2c3D(fg, fgk); //here

	fftw_free(f); 
	fftw_free(g); 
	fftw_free(fg); 

}
