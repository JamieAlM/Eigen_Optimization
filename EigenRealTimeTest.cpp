/*#define EIGEN_RUNTIME_NO_MALLOC
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <Eigen/SparseCore> //SparseMatrix and SparseVector classes, matrix assembly, basic sparse linear algebra (including sparse triangular solvers)
#include <Eigen/Sparse>

   //This is a simple test to show dynamic memory allocations should be avoided in real-time applications 
   //g++ -I /mnt/c/Users/lujain/Documents/eigen-3.4.0/eigen-3.4.0/ EigenRealTimeTest.cpp -lm -O3 -Ofast -DEIGEN_NO_DEBUG -o ../../../Datasim1/test.out
    //-DEIGEN_NO_DEBUG: disables Eigen's assertions if defined.
    //instead use the command:
    //g++ -I /mnt/c/Users/lujain/Documents/eigen-3.4.0/eigen-3.4.0/ EigenRealTimeTest.cpp -lm -O3 -Ofast -o ../../../Datasim1/test.out
    using namespace std;
    using namespace Eigen;
    
    void init(MatrixXd& a, MatrixXd& b, MatrixXd& c, int size)
    {
      a = MatrixXd::Ones(size,size);
      b = MatrixXd::Ones(size,size);
      c = 0.00001*MatrixXd::Ones(size,size);
    }
    
    void update(MatrixXd& a, MatrixXd& b, MatrixXd& c)
    {
       //Eigen::internal::is_malloc_allowed();
      // Pretty random equations to illustrate dynamic memory allocation
      Eigen::internal::set_is_malloc_allowed(false);
      b += c;
      //Eigen::internal::set_is_malloc_allowed(true);
      a += b*c;
      //Eigen::internal::set_is_malloc_allowed(true);
       
      c += b*c;
       Eigen::internal::set_is_malloc_allowed(true);
    }
    
    int main(int n_args, char** args)
    {
      int size = 32;
      if (n_args>1)
        size = atoi(args[1]);
      
      // Initialization (not real-time)
      MatrixXd a,b,c;
      init(a,b,c,size);
      
      // Real-time loop
      for (float t=0.0; t<0.1; t+=0.01)
        //Eigen::internal::set_is_malloc_allowed(false);
        update(a,b,c);
        //Eigen::internal::set_is_malloc_allowed(true);
      
      return 0;
    }*/
    

/*#define EIGEN_RUNTIME_NO_MALLOC // Define this symbol to enable runtime tests for allocations
//#define EIGEN_NO_MALLOC //This causes your program to abort whenever a temporary is created.
#include <Eigen/Core>
#include <Eigen/Dense>
//check link: http://eigen.tuxfamily.org/index.php?title=FAQ#Where_in_my_program_are_temporary_objects_created.3F
 
int main(int argc, char** argv)
{
  // It's OK to allocate here
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(20, 20);
 
  //Eigen::internal::set_is_malloc_allowed(false);
  // It's NOT OK to allocate here
  // An assertion will be triggered if an Eigen-related heap allocation takes place
  
  //Eigen::internal::set_is_malloc_allowed(true);
  // It's OK to allocate again
}
*/
#define EIGEN_DONT_VECTORIZE
#define EIGEN_RUNTIME_NO_MALLOC
#define _USE_MATH_DEFINES // for C++
#define FFTW_ESTIMATE (1U << 6)
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
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <Eigen/SparseCore> //SparseMatrix and SparseVector classes, matrix assembly, basic sparse linear algebra (including sparse triangular solvers)
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor> 
static const int nx = 128; //128; //64; //32; //4; //10;//64;//4; //128; 
static const int ny = 128; //128; //64; //32; //4;//4; //10;//64;//4;//128; 
static const int nz = 128; //16;
static const int  mm = (nx* 3) / 2;
static const int nyk = ny/2 + 1;
static const int nxk = nx/2 + 1;
static const int ncomp = 2;

#define sizee 3
using namespace std;
using namespace Eigen;
#define IMAG 1
#define REAL 0
void print2DArrf(char filename[], Ref<MatrixXd> arr); //test
void cross_product(double vector_a[], double vector_b[], double temp[]);
double mag1DArray(double arr[]);
double Arr1DMax(double arr[], int arrLength);
void print2DB(char filename[], double Lx, Ref<MatrixXd> arr);
//*****************************************
//using my_fftw_complex = double[2];
//using my_std_complex = std::complex<double>;
//test
std::string file_name_for_iteration(int n);
std::string file_name_for_iteration1(int n);
std::string file_name_for_iteration2(int n);
std::string file_name_for_iteration3(int n);

void makeFourierMesh3D(double Lx,double Ly,double Lz,  Ref<VectorXd> eKX,Ref<VectorXd> eKY,Ref<VectorXd> eKZ, Eigen::Tensor<double, 3>& KX, Eigen::Tensor<double, 3>& KY, Eigen::Tensor<double, 3>& KZ, Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<double, 3>& ninvksqu, int FourierMeshType);
double SpatialMesh3D(double dx,double dy, double dz,Eigen::Tensor<double, 3>& xarr, Eigen::Tensor<double, 3>& yarr, Eigen::Tensor<double, 3>& zarr);

void c2rfft3d(Eigen::Tensor<std::complex<double>, 3>& cArr, Eigen::Tensor<double, 3>& rArr);
void r2cfft3d(Eigen::Tensor<double, 3>& rArr, Eigen::Tensor<std::complex<double>, 3>& cArr);
void derivk3D(Eigen::Tensor<std::complex<double>, 3>& vark, Eigen::Tensor<double, 3>& K, Eigen::Tensor<std::complex<double>, 3>& derivative);
void laplaciank3D(Eigen::Tensor<std::complex<double>, 3>& vark,Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<std::complex<double>, 3>& derivative);
void convolve3D(Eigen::Tensor<std::complex<double>, 3>& fk, Eigen::Tensor<std::complex<double>, 3>& gk, Eigen::Tensor<std::complex<double>, 3>& fgk);
//void print3DPhysArray(Eigen::Tensor<double, 3>& arr);
void calcCollFreqk_inertia3D( Eigen::Tensor<std::complex<double>, 3>& nek, Eigen::Tensor<std::complex<double>, 3>& Tik, Eigen::Tensor<std::complex<double>, 3>& Tek, double kb, double eps0, double mi, double me, double ri, double rn, double nn, double Oci, double Oce, double e,Eigen::Tensor<std::complex<double>, 3>& nuink, Eigen::Tensor<std::complex<double>, 3>&  nuiek, Eigen::Tensor<std::complex<double>, 3>&  nuiik, Eigen::Tensor<std::complex<double>, 3>& nuenk,Eigen::Tensor<std::complex<double>, 3>&  nueek,Eigen::Tensor<std::complex<double>, 3>&  nueik, Eigen::Tensor<std::complex<double>, 3>& isigPk, Eigen::Tensor<std::complex<double>, 3>& invnk,Eigen::Tensor<std::complex<double>, 3>&  hallIk,Eigen::Tensor<std::complex<double>, 3>& hallEk);
void EcalcPotSourcek_inertia3D(Eigen::Tensor<std::complex<double>, 3>& ne0k,Eigen::Tensor<std::complex<double>, 3>&  nek,Eigen::Tensor<std::complex<double>, 3>& dndx0k, Eigen::Tensor<std::complex<double>, 3>& dndy0k,Eigen::Tensor<std::complex<double>, 3>& dndz0k,Eigen::Tensor<std::complex<double>, 3>& dndxk , Eigen::Tensor<std::complex<double>, 3>& dndyk, Eigen::Tensor<std::complex<double>, 3>& dndzk ,Eigen::Tensor<std::complex<double>, 3>& dphidx0k , Eigen::Tensor<std::complex<double>, 3>& dphidy0k, Eigen::Tensor<std::complex<double>, 3>& dphidz0k,Eigen::Tensor<std::complex<double>, 3>& dphidx1k,Eigen::Tensor<std::complex<double>, 3>& dphidy1k, Eigen::Tensor<std::complex<double>, 3>& dphidz1k,Eigen::Tensor<std::complex<double>, 3>& Pi1k, Eigen::Tensor<std::complex<double>, 3>& Pe1k, double uxB[], double e, double Cm, Eigen::Tensor<std::complex<double>, 3>& hallEk, Eigen::Tensor<std::complex<double>, 3>& hallIk , Eigen::Tensor<std::complex<double>, 3>& vexbkx0,Eigen::Tensor<std::complex<double>, 3>& vexbky0,  Eigen::Tensor<std::complex<double>, 3>& vexbkz0, Eigen::Tensor<std::complex<double>, 3>& vexbkx,Eigen::Tensor<std::complex<double>, 3>& vexbky, Eigen::Tensor<std::complex<double>, 3>& vexbkz,Eigen::Tensor<double, 3>& kx,Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz,Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<std::complex<double>, 3>& potSourcek_inertia);
int potentialk3D(Eigen::Tensor<std::complex<double>, 3>& invnk, Eigen::Tensor<std::complex<double>, 3>& dndxk, Eigen::Tensor<std::complex<double>, 3>& dndyk, Eigen::Tensor<std::complex<double>, 3>& dndzk,Eigen::Tensor<std::complex<double>, 3>& phik, Eigen::Tensor<std::complex<double>, 3>& potSourcek,Eigen::Tensor<double, 3>& kx, Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz, Eigen::Tensor<double, 3>&  ninvksqu, double err_max, int max_iter);
void EcalcV_ExBk3D(Eigen::Tensor<std::complex<double>, 3>& dphidxk,Eigen::Tensor<std::complex<double>, 3>& dphidyk,Eigen::Tensor<std::complex<double>, 3>& dphidzk, double B[], double B2, Eigen::Tensor<std::complex<double>, 3>& vexbkx, Eigen::Tensor<std::complex<double>, 3>& vexbky,Eigen::Tensor<std::complex<double>, 3>& vexbkz);
void Ecalc_diamag3D(Eigen::Tensor<std::complex<double>, 3>& dpdxk,Eigen::Tensor<std::complex<double>, 3>& dpdyk,Eigen::Tensor<std::complex<double>, 3>& dpdzk, double B[], double B2, double qa, Eigen::Tensor<std::complex<double>, 3>& nak, Eigen::Tensor<std::complex<double>, 3>& diamagxk,Eigen::Tensor<std::complex<double>, 3>& diamagyk, Eigen::Tensor<std::complex<double>, 3>& diamagzk);
double Ecalc_dt3D(double U[], Eigen::Tensor<double, 3> vexbx, Eigen::Tensor<double, 3> vexby, Eigen::Tensor<double, 3> vexbz,Eigen::Tensor<double, 3> diamagxi,Eigen::Tensor<double, 3> diamagyi,  Eigen::Tensor<double, 3> diamagzi, Eigen::Tensor<double, 3> diamagxe, Eigen::Tensor<double, 3> diamagye, Eigen::Tensor<double, 3> diamagze,double cfl, double kmax, double maxdt);
void ERK4(Eigen::Tensor<std::complex<double>, 3>& f, double dt, Eigen::Tensor<std::complex<double>, 3>& residual, Eigen::Tensor<std::complex<double>, 3>& source, int stage, Eigen::Tensor<std::complex<double>, 3>& fout);
void Ecalc_residualt3D(Eigen::Tensor<std::complex<double>, 3>& voxk, Eigen::Tensor<std::complex<double>, 3>& voyk, Eigen::Tensor<std::complex<double>, 3>& vozk,Eigen::Tensor<std::complex<double>, 3>& tempink, Eigen::Tensor<std::complex<double>, 3>& tempoutk, Eigen::Tensor<double, 3>& kx, Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz);
void Ecalc_residualn3D(Eigen::Tensor<std::complex<double>, 3>& vexbxk, Eigen::Tensor<std::complex<double>, 3>& vexbyk, Eigen::Tensor<std::complex<double>, 3>& vexbzk, Eigen::Tensor<std::complex<double>, 3>& nink, Eigen::Tensor<std::complex<double>, 3>& residnoutk, Eigen::Tensor<double, 3>& kx, Eigen::Tensor<double, 3>& ky,Eigen::Tensor<double, 3>& kz);
void Ecalc_sourcen3D(Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<std::complex<double>, 3>& nk, double d,Eigen::Tensor<std::complex<double>, 3>& sourcenk);



// g++ -I /mnt/c/Users/lujain/Documents/eigen-3.4.0/eigen-3.4.0/ EigenRealTimeTest.cpp -lfftw3 -lm -O3 -Ofast -o ../../../Datasim1/test.out
//speed up: g++ -I /mnt/c/Users/lujain/Documents/eigen-3.4.0/eigen-3.4.0/ EigenRealTimeTest.cpp -lfftw3 -lm -ffast-math -fno-math-errno -DNDEBUG -march=native -Ofast -o ../../../Datasim1/test.out
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
	double B[3] = {0,0,1}; 
    double u[3] = {10, 0, 0}; 
    double Lx = 128; //32;//2*EIGEN_PI;
    double Ly = 128; //32;//2*EIGEN_PI;
    double Lz = 128; //32;//4;//32;//2*EIGEN_PI;
    double x0 = 64.; ///4.;//2D Cheby/Fourier: Lx/2;  
	double y0 = 64.; ///4.; //2D Cheby/Fourier: 0.0;
    double z0 = 64.; ///4.; 
    double rx0 = 20.; ///4.;
    double ry0 =  20.; ///4.;
    double rz0 =  20.;///4.;
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
	double d = 1.; //1.6; //o; 
	// a,d,c, a2, c2,d2 are the coefficients for tanh
	double p = 3.; //was a2 //aL
	double b = 0.015; //was d2
    double c = 12.; 
    double nuin = 0.1;
	double nuen = 0.;
	double nn = 1e14; 
    double A = (2 * EIGEN_PI)/Lx;
	double A1 = (2 * EIGEN_PI)/ Ly;

    Eigen::Tensor<double, 3> eXX(128,128,128); 
	eXX.setZero(); 
    Eigen::Tensor<double, 3> eYY(128,128,128); 
	eYY.setZero(); 
    Eigen::Tensor<double, 3> eZZ(128,128,128); 
	eZZ.setZero(); 

    
	double kmax = SpatialMesh3D(dx, dy, dz,eXX,eYY, eZZ);
    
    //double mydata[3 * 4]; // Manage your own memory as you see fit
    //double* data_ptr = mydata;

    //Eigen::Map<Eigen::MatrixXi, Eigen::Unaligned> mymatrix(data_ptr, 3, 4);

    //MatrixXd mat;
    //mat = Map<MatrixXd>(data, rows, cols);

    //Eigen::MatrixXd mat = Eigen::Map<const Eigen::MatrixXd>(x.data(), rows, cols);

    //Eigen::Tensor<double, 3> KXtest = Eigen::Map<const Eigen::Tensor<double, 3>>(mydata.data(), 3, 4);

    Eigen::MatrixXd eKX((1),nx);
    eKX.setZero();
    Eigen::MatrixXd eKY((1),ny);
    eKY.setZero();
    Eigen::MatrixXd eKZ((1),ny);
    //eKZ.setZero();
    //std::cout << eKX << endl;
    Eigen::Tensor<double, 3> KX(nx,ny,nz); 
	//KX.setZero();
    Eigen::Tensor<double, 3> KY(nx,ny,nz); 
	//KY.setZero();
    Eigen::Tensor<double, 3> KZ(nx,ny,nz); 
	//KZ.setZero();

    
    Eigen::MatrixXd eksqux(1,nx); //(128,128,128);
    //eksqux.setZero();
	Eigen::MatrixXd eksquy(1,ny);
    //eksquy.setZero();
    Eigen::MatrixXd eksquz(1,nz); 
    //eksquz.setZero();

    Eigen::Tensor<double, 3> eksqu(nx,ny,nz); 
	//eksqu.setZero();
    Eigen::Tensor<double, 3> eninvksqu(nx,ny,nz); 
	//eninvksqu.setZero();
    int FourierMeshType = 1;
    

    makeFourierMesh3D( Lx, Ly, Lz, eKX, eKY, eKZ,  KX,  KY,  KZ, eksqu,  eninvksqu, FourierMeshType);
    

    Eigen::Tensor<double, 3> theta(nx,ny,nz);  
	//theta.setZero();
    Eigen::Tensor<double, 3> p0(nx,ny,nz); 
	//p0.setZero(); 
    Eigen::Tensor<double, 3> ne1(nx,ny,nz); 
	//ne1.setZero(); 
    Eigen::Tensor<double, 3> ne0(nx,ny,nz);
	//ne0.setZero(); 
    Eigen::Tensor<double, 3> ne(nx,ny,nz); 
	//ne.setZero(); 
    Eigen::Tensor<double, 3> Ti_1(nx,ny,nz);
	//Ti_1.setZero(); 
    Eigen::Tensor<double, 3> Te_1(nx,ny,nz); 
	//Te_1.setZero(); 
    Eigen::Tensor<double, 3> Ti_0(nx,ny,nz); 
	//Ti_0.setZero();
    Eigen::Tensor<double, 3> Te_0(nx,ny,nz); 
	//Te_0.setZero();
    Eigen::Tensor<double, 3> Te(nx,ny,nz); 
	//Te.setZero();
    Eigen::Tensor<double, 3> Ti(nx,ny,nz); 
	//Ti.setZero();
    Eigen::Tensor<double, 3> phi_0(nx,ny,nz); 
	//phi_0.setZero();
    Eigen::Tensor<double, 3> phi1(nx,ny,nz);
	//phi1.setZero();
    Eigen::Tensor<double, 3> phi(nx,ny,nz);
	//phi.setZero();
    Eigen::Tensor<double, 3> Pi(nx,ny,nz);
	//Pi.setZero();
    Eigen::Tensor<double, 3> Pe(nx,ny,nz);
	//Pe.setZero();
    Eigen::Tensor<double, 3> Pi0(nx,ny,nz);
	//Pi0.setZero();
    Eigen::Tensor<double, 3> Pe0(nx,ny,nz);
	//Pe0.setZero();
    Eigen::Tensor<double, 3> Pi1(nx,ny,nz);
	//Pi1.setZero();
    Eigen::Tensor<double, 3> Pe1(nx,ny,nz);
	//Pe1.setZero();

    Eigen::Tensor<double, 3> phiZ(nx,ny,nz); 

    //Eigen::Tensor<double, 3> ne(nx,ny,nz); 
    double n0 = 1e11;
	double Ampl_1 = -0.07;//-0.2; //-1E-3; // min
	double Ampl_2 = 0.07; //0.2;//1E-3; // max 

   for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
                theta(k,i,j) = 	atan2(eYY(k,i,j)  -y0,  eXX(k,i,j)  -x0)  - (EIGEN_PI/c); 
              
                //theta(i + ny*j)  = atan2(eYY(i + ny*j)  -y0,  eXX(i + ny*j)  -x0)  - (EIGEN_PI/c); 
                phiZ(k,i,j)  = atan2(eZZ(k,i,j)  -z0,  eYY(k,i,j)  -y0)  - (EIGEN_PI/c); //in y-z plane for now
                //((rand() % 101)/100. * (Ampl_2 - Ampl_1) + Ampl_1)*n0;
                //
                //p0(k,i,j) =   (( (eXX(k,i,j)-x0) * (eXX(k,i,j)-x0)  )/( rx0 * rx0  ) + ( (eYY(k,i,j)-y0) * (eYY(k,i,j)-y0)  )/( ry0 * ry0  ) + ( (eZZ(k,i,j)-z0) * (eZZ(k,i,j)-z0)  )/( rz0 * rz0  ) );
                p0(k,i,j) = (pow(eXX(k,i,j)-x0,2)/pow(rx0,2) + pow(eYY(k,i,j)-y0,2)/pow(ry0,2) + pow(eZZ(k,i,j)-z0,2)/pow(rz0,2) );               
              
                //actual si case (not normalized)
                ne1(k,i,j) = ( d + a * exp(pow(-1. * p0(k,i,j),p)) + ((rand() % 101)/100. * (Ampl_2 - Ampl_1) + Ampl_1) ) * n0  ;
             
                ne0(k,i,j)= 0.;
                ne(k,i,j) =ne1(k,i,j) + ne0(k,i,j);
			    Ti_1(k,i,j) = 0.; 
			    Te_1(k,i,j) = 0.;
			    Ti_0(k,i,j) = 1000.;
			    Te_0(k,i,j) = 1000.; 
                Te(k,i,j) = Te_1(k,i,j) + Te_0(k,i,j);
			    Ti(k,i,j) = Ti_1(k,i,j) + Ti_0(k,i,j);

                phi_0(k,i,j) =0.;  
			    phi1(k,i,j) = 0.;
			    phi(k,i,j) = phi_0(k,i,j) + phi1(k,i,j);
                Pi(k,i,j)= ne(k,i,j)* Ti(k,i,j) * kb; 
			    Pe(k,i,j) = ne(k,i,j) * Te(k,i,j) * kb; 
			
			    Pi0(k,i,j) = ne0(k,i,j) * Ti_0(k,i,j) * kb;
			    Pe0(k,i,j) = ne0(k,i,j) * Te_0(k,i,j) * kb;
			    //Get perturbed pressures:
			    Pi1(k,i,j)= Pi(k,i,j) - Pi0(k,i,j);
			    Pe1(k,i,j) = Pe(k,i,j) - Pe0(k,i,j);
			}

		}		
	}
    
    

      
	

    Eigen::Tensor<std::complex<double>, 3> nek(nx,ny,nz);
    Eigen::Tensor<std::complex<double>, 3> Tik(nx,ny,nz); 
	Tik.setZero();
    Eigen::Tensor<std::complex<double>, 3> Tek(nx,ny,nz); 
	Tek.setZero();
    Eigen::Tensor<std::complex<double>, 3> phik(nx,ny,nz);  
	phik.setZero();
    Eigen::Tensor<std::complex<double>, 3> Pik(nx,ny,nz); 
	Pik.setZero();
    Eigen::Tensor<std::complex<double>, 3> Pek(nx,ny,nz);  
	Pek.setZero();
    Eigen::Tensor<std::complex<double>, 3> ne0k(nx,ny,nz);  
	ne0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> ne1k(nx,ny,nz); 
	ne1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> Te0k(nx,ny,nz); 
	Te0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> Ti0k(nx,ny,nz);  
	Ti0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> phi0k(nx,ny,nz);  
	phi0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> phi1k(nx,ny,nz); 
	phi1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> Pi0k(nx,ny,nz); 
	Pi0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> Pi1k(nx,ny,nz);  
	Pi1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> Pe0k(nx,ny,nz); 
	Pe0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> Pe1k(nx,ny,nz); 
	Pe1k.setZero();

    r2cfft3d(ne, nek);
    r2cfft3d(Ti, Tik);
	r2cfft3d(Te, Tek);
	r2cfft3d(phi, phik);
	r2cfft3d(Pi, Pik);
	r2cfft3d(Pe, Pek);
	//std::cout <<nek <<endl;
	//************** Inertial TEST ******************
	r2cfft3d(ne0, ne0k);
	r2cfft3d(ne1, ne1k);
	r2cfft3d(Te_0, Te0k);
    r2cfft3d(Ti_0, Ti0k);

	r2cfft3d(phi_0, phi0k);
	r2cfft3d(phi1, phi1k);
	r2cfft3d(Pi0, Pi0k);
	r2cfft3d(Pi1, Pi1k);
	r2cfft3d(Pe0, Pe0k);
	r2cfft3d(Pe1, Pe1k);

    Eigen::Tensor<std::complex<double>, 3> hallEk(nx,ny,nz); 
	hallEk.setZero();
    Eigen::Tensor<std::complex<double>, 3> hallIk(nx,ny,nz);  
	hallIk.setZero();
    Eigen::Tensor<std::complex<double>, 3> nuink(nx,ny,nz);  
	nuink.setZero();
    Eigen::Tensor<std::complex<double>, 3> nuiek(nx,ny,nz); 
	nuiek.setZero();
    Eigen::Tensor<std::complex<double>, 3> nuiik(nx,ny,nz); 
	nuiik.setZero();
    Eigen::Tensor<std::complex<double>, 3> nuenk(nx,ny,nz);  
	nuenk.setZero();
    Eigen::Tensor<std::complex<double>, 3> nueek(nx,ny,nz); 
	nueek.setZero();
    Eigen::Tensor<std::complex<double>, 3> nueik(nx,ny,nz);  
	nueik.setZero();
    Eigen::Tensor<std::complex<double>, 3> isigPk(nx,ny,nz); 
	isigPk.setZero();
    Eigen::Tensor<std::complex<double>, 3> invnk(nx,ny,nz); 
	invnk.setZero();
    
	
    calcCollFreqk_inertia3D(nek, Tik, Tek , kb, eps0, mi, me, ri, rn, nn, Oci, Oce, e, nuink, nuiek, nuiik, nuenk, nueek, nueik, isigPk, invnk, hallIk, hallEk);

    Eigen::Tensor<std::complex<double>, 3> potSourcek_inertia(nx,ny,nz);  
	potSourcek_inertia.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndxk(nx,ny,nz);  
	dndxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndyk(nx,ny,nz);  
	dndyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndzk(nx,ny,nz); 
	dndzk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndx0k(nx,ny,nz); 
	dndx0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndy0k(nx,ny,nz); 
	dndy0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndz0k(nx,ny,nz); 
	dndz0k.setZero();
    
    Eigen::Tensor<std::complex<double>, 3> dndx1k(nx,ny,nz); 
	dndx1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndy1k(nx,ny,nz);  
	dndy1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dndz1k(nx,ny,nz); 
	dndz1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidxk(nx,ny,nz);  
	dphidxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidyk(nx,ny,nz); 
	dphidyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidzk(nx,ny,nz);  
	dphidzk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidx0k(nx,ny,nz); 
	dphidx0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidx1k(nx,ny,nz);  
	dphidx1k.setZero();

    Eigen::Tensor<std::complex<double>, 3> dphidy0k(nx,ny,nz); 
	dphidy0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidy1k(nx,ny,nz); 
	dphidy1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidz0k(nx,ny,nz);  
	dphidz0k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dphidz1k(nx,ny,nz);  
	dphidz1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> d2phidxk(nx,ny,nz); 
	d2phidxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> d2phidyk(nx,ny,nz);
	d2phidyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> d2phidzk(nx,ny,nz); 
	d2phidzk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpedxk(nx,ny,nz); 
	dpedxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpedx1k(nx,ny,nz);  
	dpedx1k.setZero();

    Eigen::Tensor<std::complex<double>, 3> dpedyk(nx,ny,nz); 
	dpedyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpedy1k(nx,ny,nz); 
	dpedy1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpedzk(nx,ny,nz);  
	dpedzk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpedz1k(nx,ny,nz); 
	dpedz1k.setZero();

    Eigen::Tensor<std::complex<double>, 3> dpidxk(nx,ny,nz); 
	dpidxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpidx1k(nx,ny,nz);  
	dpidx1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpidyk(nx,ny,nz);  
	dpidyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpidy1k(nx,ny,nz); 
	dpidy1k.setZero();

    Eigen::Tensor<std::complex<double>, 3> dpidzk(nx,ny,nz); 
	dpidzk.setZero();
    Eigen::Tensor<std::complex<double>, 3> dpidz1k(nx,ny,nz);  
	dpidz1k.setZero();

    Eigen::Tensor<std::complex<double>, 3> dphikdt(nx,ny,nz);  
	dphikdt.setConstant(1.0f); 
    //std::cout<<dphikdt<<endl;
    
    derivk3D(nek, KX, dndxk);
    derivk3D(nek, KY, dndyk);
    derivk3D(nek, KZ, dndzk);

    
	// add ne0k
	derivk3D(ne0k, KX, dndx0k);
	derivk3D(ne0k, KY, dndy0k);
    derivk3D(ne0k, KZ, dndz0k);
	// add ne1k //dndx1k dndy1k
	derivk3D(ne1k, KX, dndx1k);
	derivk3D(ne1k, KY, dndy1k);
    derivk3D(ne1k, KZ, dndz1k);

    
	derivk3D(phik, KX, dphidxk);
	derivk3D(phik, KY, dphidyk);
    derivk3D(phik, KZ, dphidzk);

    // add phi0k
	derivk3D(phi0k, KX, dphidx0k);
	derivk3D(phi0k, KY, dphidy0k);
    derivk3D(phi0k, KZ, dphidz0k);
	// add phi1k
	derivk3D(phi1k, KX, dphidx1k);
	derivk3D(phi1k, KY, dphidy1k);
    derivk3D(phi1k, KZ, dphidy1k);

    
	derivk3D(Pek, KX, dpedxk);
	derivk3D(Pek, KY, dpedyk);
    derivk3D(Pek, KZ, dpedzk);

	derivk3D(Pik, KX, dpidxk);
	derivk3D(Pik, KY, dpidyk);
    derivk3D(Pik, KZ, dpidzk);
	// add Pe1k and Pi1k ONLY for now
	derivk3D(Pe1k, KX, dpedx1k); 
	derivk3D(Pe1k, KY, dpedy1k);
    derivk3D(Pe1k, KZ, dpedz1k);

	derivk3D(Pi1k, KX, dpidx1k);
	derivk3D(Pi1k, KY, dpidy1k);
    derivk3D(Pi1k, KZ, dpidz1k);
    
  
    Eigen::Tensor<std::complex<double>, 3> vexbkx(nx,ny,nz); 
	vexbkx.setZero();
    Eigen::Tensor<std::complex<double>, 3> vexbky(nx,ny,nz); 
	vexbky.setZero();
    Eigen::Tensor<std::complex<double>, 3> vexbkz(nx,ny,nz); 
	vexbkz.setZero();
    
    EcalcV_ExBk3D(dphidxk, dphidyk, dphidzk, B, B2, vexbkx, vexbky,vexbkz);
    

    Eigen::Tensor<double, 3> vexbx(nx,ny,nz); 
	vexbx.setZero(); 
    Eigen::Tensor<double, 3> vexby(nx,ny,nz); 
	vexby.setZero(); 
    Eigen::Tensor<double, 3> vexbz(nx,ny,nz); 
	vexbz.setZero(); 
    c2rfft3d(vexbkx, vexbx);
	c2rfft3d(vexbky, vexby);
    c2rfft3d(vexbkz, vexbz);
    Eigen::Tensor<std::complex<double>, 3> vexbkx0(nx,ny,nz); 
	vexbkx0.setZero();
    Eigen::Tensor<std::complex<double>, 3> vexbky0(nx,ny,nz); 
	vexbky0.setZero();
    Eigen::Tensor<std::complex<double>, 3> vexbkz0(nx,ny,nz);  
	vexbkz0.setZero();
    
    EcalcV_ExBk3D(dphidx0k, dphidy0k, dphidz0k, B, B2, vexbkx0,vexbky0,vexbkz0);
    
	Eigen::Tensor<double, 3> vexbx0(nx,ny,nz);
	//vexbx0.setZero(); 
    Eigen::Tensor<double, 3> vexby0(nx,ny,nz);
	//vexby0.setZero(); 
    Eigen::Tensor<double, 3> vexbz0(nx,ny,nz);
	//vexbz0.setZero(); 
    
    c2rfft3d(vexbkx0, vexbx0);
	c2rfft3d(vexbky0, vexby0);
    c2rfft3d(vexbkz0, vexbz0);
   

    Eigen::Tensor<std::complex<double>, 3> vexbkx1(nx,ny,nz);  
	//vexbkx1.setZero();
    Eigen::Tensor<std::complex<double>, 3> vexbky1(nx,ny,nz);  
	//vexbky1.setZero();
    Eigen::Tensor<std::complex<double>, 3> vexbkz1(nx,ny,nz);
	//vexbkz1.setZero();

    EcalcV_ExBk3D(dphidx1k, dphidy1k, dphidz0k, B, B2, vexbkx1,vexbky1,vexbkz1);
	
    Eigen::Tensor<double, 3> vexbx1(nx,ny,nz);
	//vexbx1.setZero(); 
    Eigen::Tensor<double, 3> vexby1(nx,ny,nz); 
	//vexby1.setZero(); 
    Eigen::Tensor<double, 3> vexbz1(nx,ny,nz);
	//vexbz1.setZero();

    c2rfft3d(vexbkx1, vexbx1);
	c2rfft3d(vexbky1, vexby1);
    c2rfft3d(vexbkz1, vexbz1);
	
    Eigen::Tensor<std::complex<double>, 3> vdmexk(nx,ny,nz);
	//vdmexk.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmeyk(nx,ny,nz);
	//vdmeyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmezk(nx,ny,nz); 
	//vdmezk.setZero();
     
    //Eigen::internal::set_is_malloc_allowed(false);
	Ecalc_diamag3D(dpedxk, dpedyk, dpedzk, B, B2, -1 * e,  nek, vdmexk, vdmeyk, vdmezk);
    //Eigen::internal::set_is_malloc_allowed(true);

    Eigen::Tensor<std::complex<double>, 3> vdmixk(nx,ny,nz);  
	//vdmixk.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmiyk(nx,ny,nz);  
	//vdmiyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmizk(nx,ny,nz); 
	//vdmizk.setZero();

    Ecalc_diamag3D(dpidxk, dpidyk, dpidzk, B, B2, e,  nek, vdmixk, vdmiyk, vdmizk);
    Eigen::Tensor<double, 3> vdmex(nx,ny,nz); 
	//vdmex.setZero(); 
    Eigen::Tensor<double, 3> vdmey(nx,ny,nz);
	//vdmey.setZero(); 
    Eigen::Tensor<double, 3> vdmez(nx,ny,nz);
	//vdmez.setZero();
    Eigen::Tensor<double, 3> vdmix(nx,ny,nz);
	//vdmix.setZero();
    Eigen::Tensor<double, 3> vdmiy(nx,ny,nz);
	//vdmiy.setZero();
    Eigen::Tensor<double, 3> vdmiz(nx,ny,nz);
	//vdmiz.setZero();

    c2rfft3d(vdmexk, vdmex);
	c2rfft3d(vdmeyk, vdmey);
    c2rfft3d(vdmezk, vdmez);

	c2rfft3d(vdmixk, vdmix);
	c2rfft3d(vdmiyk, vdmiy);
    c2rfft3d(vdmizk, vdmiz);

    Eigen::Tensor<std::complex<double>, 3> vdmex1k(nx,ny,nz); 
	//vdmex1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmey1k(nx,ny,nz); 
	//vdmey1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmez1k(nx,ny,nz);
	//vdmez1k.setZero();

    Ecalc_diamag3D(dpedx1k, dpedy1k, dpedz1k, B, B2, -1 * e,  ne1k, vdmex1k, vdmey1k, vdmez1k);
    
    Eigen::Tensor<double, 3> vdmex1(nx,ny,nz); 
	//vdmex1.setZero(); 
    Eigen::Tensor<double, 3> vdmey1(nx,ny,nz);
	//vdmey1.setZero(); 
    Eigen::Tensor<double, 3> vdmez1(nx,ny,nz); 
	//vdmez1.setZero();
    c2rfft3d(vdmex1k, vdmex1); 
	c2rfft3d(vdmey1k, vdmey1);
    c2rfft3d(vdmez1k, vdmez1);

    Eigen::Tensor<std::complex<double>, 3> vdmix1k(nx,ny,nz); 
	//vdmix1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmiy1k(nx,ny,nz);  
	//vdmiy1k.setZero();
    Eigen::Tensor<std::complex<double>, 3> vdmiz1k(nx,ny,nz);
	//vdmiz1k.setZero();

    Ecalc_diamag3D(dpidx1k, dpidy1k, dpidz1k, B, B2, e,  ne1k, vdmix1k, vdmiy1k, vdmiz1k);
    
    Eigen::Tensor<double, 3> vdmix1(nx,ny,nz);
	//vdmix1.setZero(); 
    Eigen::Tensor<double, 3> vdmiy1(nx,ny,nz); 
	//vdmiy1.setZero(); 
    Eigen::Tensor<double, 3> vdmiz1(nx,ny,nz);
	//vdmiz1.setZero();

	c2rfft3d(vdmix1k, vdmix1);
	c2rfft3d(vdmiy1k, vdmiy1); 
    c2rfft3d(vdmiz1k, vdmiz1); 

    std::vector<double> time(iter_max + 1,0.0); //This will create a vector of size iter_max + 1 all initialized to 0.0. You could use memset as well 
	
    c2rfft3d(nek, ne);
	c2rfft3d(Tek, Te);	
	c2rfft3d(Tik, Ti);
	c2rfft3d(phik, phi);

    //print ICs

    Eigen::Tensor<std::complex<double>, 3> ne1k_old(nx,ny,nz);  
	//ne1k_old.setZero();
    Eigen::Tensor<std::complex<double>, 3> Tek_old(nx,ny,nz);
	//Tek_old.setZero();
    Eigen::Tensor<std::complex<double>, 3> Tik_old(nx,ny,nz);  
	//Tik_old.setZero();
    Eigen::Tensor<std::complex<double>, 3> phi1k_old(nx,ny,nz); 
	//phi1k_old.setZero();

    
    Eigen::Tensor<std::complex<double>, 3> residualnk(nx,ny,nz);  
	//residualnk.setZero();
    Eigen::Tensor<std::complex<double>, 3> vioxk(nx,ny,nz); 
	//vioxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> vioyk(nx,ny,nz);
	//vioyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> viozk(nx,ny,nz);
	//viozk.setZero();
    Eigen::Tensor<std::complex<double>, 3> residualtik(nx,ny,nz);
	//residualtik.setZero();

     Eigen::Tensor<std::complex<double>, 3> veoxk(nx,ny,nz);
	//veoxk.setZero();
    Eigen::Tensor<std::complex<double>, 3> veoyk(nx,ny,nz); 
	//veoyk.setZero();
    Eigen::Tensor<std::complex<double>, 3> veozk(nx,ny,nz);
	//veozk.setZero();
    Eigen::Tensor<std::complex<double>, 3> residualtek(nx,ny,nz); 
	//residualtek.setZero();
    Eigen::Tensor<std::complex<double>, 3> sourcen1k(nx,ny,nz);  
	//sourcen1k.setZero();
    
    Eigen::Tensor<std::complex<double>, 3> residualk_phi(nx,ny,nz); 
	//residualk_phi.setZero();
    Eigen::Tensor<std::complex<double>, 3> sourcetk(nx,ny,nz);
	//sourcetk.setZero();

    //begin time step loop here
     

    for (int iter = 0; iter < iter_max; iter++){
		// total veloc
		c2rfft3d(vexbkx, vexbx);
		c2rfft3d(vexbky, vexby);
        c2rfft3d(vexbkz, vexbz);
        // diamg e, ion
        c2rfft3d(vdmexk, vdmex);
		c2rfft3d(vdmeyk, vdmey);
        c2rfft3d(vdmezk, vdmez);

		c2rfft3d(vdmixk, vdmix);
		c2rfft3d(vdmiyk, vdmiy);
        c2rfft3d(vdmizk, vdmiz);
        
        double dt =  Ecalc_dt3D(u, vexbx,  vexby, vexbz, vdmix, vdmiy,  vdmiz,  vdmex, vdmey, vdmez, CFL,  kmax,  dt_max);
        //std::cout<<dt<<endl;
        time[iter + 1] = time[iter] + dt; 
        for(int i = 0; i< nx; i++){
		    for(int j = 0; j< ny; j++){
			    for(int k = 0; k< nz; k++){ 
					ne1k_old(k,i,j) = ne1k(k,i,j); // change to perturbed all of them: ne1k
					Tek_old(k,i,j) = Tek(k,i,j); //Tek and Ti are pertrubed
					Tik_old(k,i,j) = Tik(k,i,j);
					phi1k_old(k,i,j) = phi1k(k,i,j); 
				}
			}
		}


        for (int stage = 0; stage < 4; stage++){
            
			Ecalc_residualn3D(vexbkx, vexbky, vexbkz, nek, residualnk, KX,  KY, KZ);
			
            Ecalc_residualt3D(vioxk, vioyk,viozk,Tik, residualtik, KX,KY, KZ);

			Ecalc_residualt3D(veoxk, veoyk,veozk,Tek, residualtek, KX,KY, KZ);
            
			// pot source	            
            EcalcPotSourcek_inertia3D(ne0k, nek, dndx0k, dndy0k,dndz0k,dndxk , dndyk, dndzk ,dphidx0k , dphidy0k, dphidz0k, dphidx1k,dphidy1k,dphidz1k, Pi1k, Pe1k, uxB, e, Cm,  hallEk, hallIk , vexbkx0, vexbky0, vexbkz0, vexbkx, vexbky, vexbkz, KX, KY, KZ, eksqu, potSourcek_inertia);
           
            phi_iter = potentialk3D(invnk, dndxk, dndyk, dndzk, dphikdt, potSourcek_inertia, KX, KY, KZ, eninvksqu, err_max,  phi_iter_max);
            //std::cout<<dphikdt<<endl;
			if (phi_iter > phi_iter_max){
				printf("Blew up");
			break;
		    }
            Ecalc_sourcen3D(eksqu, ne1k, Dart, sourcen1k); 
            ERK4(ne1k_old, dt, residualnk, sourcen1k, stage, ne1k); 

            ERK4(phi1k_old, dt, residualk_phi, dphikdt, stage, phi1k);

            ERK4(Tik_old, dt, residualtik, sourcetk, stage, Tik); 

            ERK4(Tek_old, dt, residualtek, sourcetk, stage, Tek); 
			
            //Total ne and phi: for loop or matrix operations
            nek = ne0k + ne1k;
		    // phi
            phik = phi0k + phi1k;
            convolve3D(nek, Tek, Pek);
            Pek = kb * Pek; 
            convolve3D(nek, Tik, Pik);
            Pik = kb * Pik;
		    // pert pressure
            Pe1k = Pek + Pe0k;
            Pi1k = Pik + Pi0k;

            // all derv for tot
			derivk3D(nek, KX, dndxk);
			derivk3D(nek, KY, dndyk);
            derivk3D(nek, KZ, dndzk);

            derivk3D(phik, KX, dphidxk);
			derivk3D(phik, KY, dphidyk);
            derivk3D(phik, KZ, dphidzk);

			derivk3D(Pek, KX, dpedxk);
			derivk3D(Pek, KY, dpedyk);
            derivk3D(Pek, KZ, dpedzk);

			derivk3D(Pik, KX, dpidxk);
			derivk3D(Pik, KY, dpidyk); 
            derivk3D(Pik, KZ, dpidzk);

            derivk3D(phi1k, KX, dphidx1k);
			derivk3D(phi1k, KY, dphidy1k);
            derivk3D(phi1k, KZ, dphidz1k);

			calcCollFreqk_inertia3D(nek, Tik, Tek , kb, eps0, mi, me, ri, rn, nn, Oci, Oce, e, nuink, nuiek, nuiik, nuenk, nueek, nueik, isigPk, invnk, hallIk, hallEk);
	
            // tot velo
			EcalcV_ExBk3D(dphidxk, dphidyk,dphidzk, B, B2, vexbkx, vexbky,vexbkz);
			
			Ecalc_diamag3D(dpedxk, dpedyk, dpedzk,B, B2, -1 * e, nek, vdmexk, vdmeyk,vdmezk); 
			Ecalc_diamag3D(dpidxk, dpidyk,dpidzk ,B, B2, e, nek, vdmixk, vdmiyk,vdmizk);

            // Get total velocity: use for loop or matrix operations
            veoxk = vdmexk + vexbkx;	
		    veoyk = vdmeyk + vexbky;
            veozk = vdmezk + vexbkz;

            vioxk = vdmixk + vexbkx;
            vioyk = vdmiyk + vexbky;
            viozk = vdmizk + vexbkz;
            
        }
        
        if (phi_iter > phi_iter_max){
            printf("Blew up");
	        break;
	    }
            
         
        
        printf("Iteration = %d    t = %.10f   phi_iter = %d\n", iter, time[iter+1], phi_iter);
        //fprintf("Iteration = %d    t = %.10f   phi_iter = %d\n", iter, time[iter+1], phi_iter);
		if ((iter/saveFrequency) - saveNum == 0){
            c2rfft3d(nek, ne);
			c2rfft3d(Tek, Te);	
			c2rfft3d(Tik, Ti);
			c2rfft3d(phik, phi);

                //print stuff
            char save [16] = {0};
			snprintf(save, 16,  "%d", saveNum);
            const char *type = ".txt";

            char nefilename[16] = {0};
            strcat(nefilename, "ne");
            

            //ofstream myfile;
  		    //myfile.open ("neTest.txt");
            
		    //for (int iter = 0; iter < 1; iter++){ //iter = -1
            //for (int iter = 0; iter < iter_max; iter++){ 
                std::ofstream my_file( file_name_for_iteration(iter) );
  		        my_file << ne << '\n';
                my_file.close();
		    //}
  		   
		
            saveNum++;
			
        }
        
	
			
			


			

			
	


        

       
         
    }
    

}
std::string file_name_for_iteration(int n) 
{ return "ne" + std::to_string(n) + ".txt"; }
std::string file_name_for_iteration1(int n) 
{ return "phi" + std::to_string(n) + ".txt"; }
std::string file_name_for_iteration2(int n) 
{ return "Ti" + std::to_string(n) + ".txt"; }
std::string file_name_for_iteration3(int n) 
{ return "Te" + std::to_string(n) + ".txt"; }

void Arr3DArr2DMult(double arr3D[][ncomp], double arr2D[], double arrOut[][ncomp]){
	for(int i = 0; i < nx; i++){
		for(int j = 0; j < nyk; j++){
			for (int k = 0; k < ncomp; k++){
				arrOut[j + nyk*i][k] = arr3D[j + nyk*i][k] * arr2D[j + nyk*i];
			}
		}
	}	
}

// This function multiplies a 3D Fourier array by i.
void iArr3DMult(double arr[][ncomp], double arrOut[][ncomp]){
	double dummy;   // This is used so that data are not accidentally overwritten
	for(int i = 0; i < nx; i++){
		for( int j = 0; j < nyk; j++){	
			dummy = arr[j + nyk*i][0];
			arrOut[j + nyk*i][0] = -arr[j + nyk*i][1];
			arrOut[j + nyk*i][1] =  dummy;
		}
	}

	
} 
void print2DArrf(char filename[], Ref<MatrixXd> arr){
 //void print2DArrf(char filename[], double arr[]){

	FILE *fptr;
	
	fptr = fopen(filename, "w");

	//Eigen::Matrix< double, (ny+1), (nx)> dummy; 
	Eigen::MatrixXd dummy((ny+1),nx);
    dummy.setZero();
	dummy = arr;
	
	for (int i = 0; i < nx; i++) {
		for (int j = 0 ; j < ny+1 ; j++) { //for (int j = 0 ; j < ny ; j++) {
			fprintf(fptr, "%+4.16le  ",dummy(j + (ny+1)*i)); //arr[j + ny*i]);
			//fprintf(fptr, "%+4.16le  ",arr(j + (ny+1)*i)); //arr[j + ny*i]);

		}
		fprintf(fptr, "\n");		
	}
	
	fclose(fptr);
}
double mag1DArray(double arr[]){
	return sqrt(arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2]);
}	 

void cross_product(double vector_a[], double vector_b[], double temp[]) {

   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];//   temp[1] = vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0];
   temp[2] = vector_a[1] * vector_b[0] - vector_a[0] * vector_b[1];// temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

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
void print2DB(char filename[], double Lx, Ref<MatrixXd> arr){ 

	Eigen::MatrixXd dummy((ny+1),nx);
	dummy.setZero();
	dummy = arr;
	
    // open file to write solution
    FILE *fp = fopen(filename, "wb"); if (!fp) return;

    // write grid
    uint64_t ndim = 2;
    uint64_t cells[] = { nx,ny };
    double lower[] = { 0,-1 }, upper[] = { Lx,1 }; //  test
	//double lower[] = { 0,0 }, upper[] = { 1,1 }; // from 0 to Lx and -1 to 1

	//test

	uint64_t real_type = 2;
    fwrite(&real_type, sizeof(uint64_t), 1, fp);

    fwrite(&ndim, sizeof(uint64_t), 1, fp);
    fwrite(cells, 2*sizeof(uint64_t), 1, fp);
    fwrite(lower, 2*sizeof(double), 1, fp);
    fwrite(upper, 2*sizeof(double), 1, fp);

    uint64_t esznc = sizeof(double), size =nx*ny;
    fwrite(&esznc, sizeof(uint64_t), 1, fp);
    fwrite(&size, sizeof(uint64_t), 1, fp);
 
    fwrite(&dummy(0), esznc*size, 1, fp);
	//fwrite(&arr[0], esznc*size, 1, fp);
    fclose(fp);
	
	
}
void makeFourierMesh3D(double Lx,double Ly,double Lz,  Ref<VectorXd> eKX,Ref<VectorXd> eKY,Ref<VectorXd> eKZ, Eigen::Tensor<double, 3>& KX, Eigen::Tensor<double, 3>& KY, Eigen::Tensor<double, 3>& KZ, Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<double, 3>& ninvksqu, int FourierMeshType){
 
	
	double dx = Lx/nx; 
	double dy = Ly/ny; 
	double dz = Lz/nz; 

	int k_counter = 0;		
	// Make kx. kx corresponds to modes [0:n/2-1 , -n/2:-1]. This is why there is an additional step, just due to FFTW3's structuring
	//#pragma omp parallel for schedule(dynamic) reduction(+ \
													 : k_counter)
	for(int i = 0; i < nx ; i++){ //change to ny
			if (i < nx/2){ 
				eKX(i) = 2.*EIGEN_PI*i/Lx;
				eKY(i) = 2.*EIGEN_PI*i/Ly;
                //eKZ(i) = 2.*EIGEN_PI*i/Lz;
				}
				if( i >= nx/2){ // i >= nx/2 --> from -128 to -1
				eKX(i) = 2.* EIGEN_PI * (-i + 2.*k_counter) / Lx;//	2.* M_PI * (-i + 2.*k_counter) / Lx;		
				eKY(i) = 2.* EIGEN_PI * (-i + 2.*k_counter) / Ly;
                //eKZ(i) = 2.* EIGEN_PI * (-i + 2.*k_counter) / Lz;
				//std::cout << eKX[i] << endl;
				}
				if( i >= nx/2){
					k_counter++;
				}
		
	    }
    
    int k_counterz = 0;
	//#pragma omp parallel for schedule(dynamic) reduction(+ \
													 : k_counterz)
     for(int k = 0; k < nz ; k++){ //change to ny
			if (k < nz/2){ 
                eKZ(k) = 2.*EIGEN_PI*k/Lz;
				}
				if( k >= nz/2){ // i >= nx/2 --> from -128 to -1
                eKZ(k) = 2.* EIGEN_PI * (-k + 2.*k_counterz) / Lz;
				}
				if( k >= nz/2){
					k_counterz++;
				}
		
	}
	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				KX(i,j,k) = eKX(j); //
                KY(i,j,k) = eKY(i);
                KZ(i,j,k) = eKZ(k);
			}

		}		
	}
	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				ksqu(i,j,k) = ( sin(KX(i,j,k) * dx/2)/(dx/2) * sin(KX(i,j,k) * dx/2)/(dx/2) + sin(KY(i,j,k) * dy/2)/(dy/2) * sin(KY(i,j,k) * dy/2)/(dy/2) + sin(KZ(i,j,k) * dz/2)/(dz/2) * sin(KZ(i,j,k) * dz/2)/(dz/2) );
                //ksqu(i,j,k) = ( pow(sin(KX(i,j,k) * dx/2)/(dx/2),2) + pow(sin(KY(i,j,k) * dy/2)/(dy/2),2) + pow(sin(KZ(i,j,k) * dz/2)/(dz/2),2) );

				ninvksqu(i,j,k) = -1./ ksqu(i,j,k);
                ninvksqu(0) = -1.;
			}
		}

	}		
	
	for(int i = 0; i< nx; i++){
            eKX(i) = sin(eKX(i)* dx)/dx ; 
			eKY(i) = sin(eKY(i)* dy)/dy ; 
            //eKZ(i) = sin(eKZ(i)* dz)/dz ; 
    }
    for(int k = 0; k< nz; k++){
            eKZ(k) = sin(eKZ(k)* dz)/dz ; 
    }
     
    //[KX,KY,KZ] = meshgrid(kx,ky,kz);
	#pragma omp parallel for schedule(static) collapse(3)
    for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
                KX(i,j,k) = eKX(j); //
                KY(i,j,k) = eKY(i);
                KZ(i,j,k) = eKZ(k);
			}

		}		
	}
		
		
	
}
double SpatialMesh3D(double dx,double dy, double dz,Eigen::Tensor<double, 3>& xarr, Eigen::Tensor<double, 3>& yarr, Eigen::Tensor<double, 3>& zarr){
	//test
    //Eigen::internal::set_is_malloc_allowed(false);
	//#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				xarr(k,i,j) = i*dx;
				yarr(j,i,k) = j*dy;
				zarr(j,i,k) = k*dz; 
			}

		}		
	}
		if (dx < dy){
			return 1/(2*dx);
		}else {
			return 1/(2*dy);
		}
		
}
void r2cfft3d(Eigen::Tensor<double, 3>& rArr, Eigen::Tensor<std::complex<double>, 3>& cArr){	
	   
	

		fftw_complex *input_array;
		input_array = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
		
		
		memcpy(input_array, rArr.data(), (nx*ny*nz) * sizeof(fftw_complex));

		fftw_complex *output_array;
		output_array = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
		

		fftw_plan forward = fftw_plan_dft_3d(nx, ny, nz, input_array, output_array, FFTW_FORWARD, FFTW_ESTIMATE);
			
		fftw_execute(forward);
    	fftw_destroy_plan(forward);
		fftw_cleanup();
		//fftw_cleanup_threads();
	
		
		memcpy(cArr.data(),output_array, (nx*ny*nz) * sizeof(fftw_complex));
		
		
		
		fftw_free(input_array);
		fftw_free(output_array);
		
	//}
	

}
void c2rfft3d(Eigen::Tensor<std::complex<double>, 3>& cArr, Eigen::Tensor<double, 3>& rArr){	
	
	
	//#pragma omp parallel
	//{
		fftw_complex *input_array;
		input_array = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
		memcpy(input_array, cArr.data(), (nx*ny*nz) * sizeof(fftw_complex));
		
		fftw_complex *output_array;
		output_array = (fftw_complex*) fftw_malloc((nx*ny*nz) * sizeof(fftw_complex));
		
		fftw_plan backward = fftw_plan_dft_3d(nx, ny, nz, input_array, output_array, FFTW_BACKWARD, FFTW_ESTIMATE);
		
		fftw_execute(backward);
    	fftw_destroy_plan(backward);
		fftw_cleanup();
		//fftw_cleanup_threads();

		memcpy(rArr.data(),output_array, (nx*ny*nz) * sizeof(double)); //size of double fixes mem error
		
		Eigen::Tensor<double, 3> dummy(nx,ny,nz); 
        //Eigen::internal::set_is_malloc_allowed(false);
		dummy = 1.0/(nx*ny*nz) * rArr;
		rArr = dummy;
        //Eigen::internal::set_is_malloc_allowed(true);
		
		
		fftw_free(input_array);
		fftw_free(output_array);
		
	

}
void derivk3D(Eigen::Tensor<std::complex<double>, 3>& vark, Eigen::Tensor<double, 3>& K, Eigen::Tensor<std::complex<double>, 3>& derivative){
	//Eigen::internal::set_is_malloc_allowed(false);
	derivative = vark * K * 1i; //.asDiagonal(); //correct
    //Eigen::internal::set_is_malloc_allowed(true);

}
void laplaciank3D(Eigen::Tensor<std::complex<double>, 3>& vark,Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<std::complex<double>, 3>& derivative){	
	derivative =  -1. * vark * ksqu;
}
void convolve3D(Eigen::Tensor<std::complex<double>, 3>& fk, Eigen::Tensor<std::complex<double>, 3>& gk, Eigen::Tensor<std::complex<double>, 3>& fgk){
	
	Eigen::Tensor<double, 3> f(nx,ny,nz);   
	//f.setZero();
	Eigen::Tensor<double, 3> g(nx,ny,nz);   
	//g.setZero();
	Eigen::Tensor<double, 3> fg(nx,ny,nz);    
	//fg.setZero();

	c2rfft3d(fk,f);
	c2rfft3d(gk,g);

	//multiply in real space
    //Eigen::internal::set_is_malloc_allowed(false);
	fg = f * g;
    //Eigen::internal::set_is_malloc_allowed(true);
	//take fft of product
	r2cfft3d(fg,fgk);

}



void calcCollFreqk_inertia3D( Eigen::Tensor<std::complex<double>, 3>& nek, Eigen::Tensor<std::complex<double>, 3>& Tik, Eigen::Tensor<std::complex<double>, 3>& Tek, double kb, double eps0, double mi, double me, double ri, double rn, double nn, double Oci, double Oce, double e,Eigen::Tensor<std::complex<double>, 3>& nuink, Eigen::Tensor<std::complex<double>, 3>&  nuiek, Eigen::Tensor<std::complex<double>, 3>&  nuiik, Eigen::Tensor<std::complex<double>, 3>& nuenk,Eigen::Tensor<std::complex<double>, 3>&  nueek,Eigen::Tensor<std::complex<double>, 3>&  nueik, Eigen::Tensor<std::complex<double>, 3>& isigPk, Eigen::Tensor<std::complex<double>, 3>& invnk, Eigen::Tensor<std::complex<double>, 3>&  hallIk,Eigen::Tensor<std::complex<double>, 3>& hallEk){	// Take inverse fft of nek, Tik, and Tek. 
	// The reason we have to convert all of this to real space is because we can't take square roots or reciprocals in Fourier space
    
	//Eigen::Matrix< double, (ny+1), (nx)> ne; 
	Eigen::Tensor<double, 3> ne(nx,ny,nz);   
	//ne.setZero();

	Eigen::Tensor<double, 3> Ti(nx,ny,nz);  
	//Ti.setZero();
	
	Eigen::Tensor<double, 3> Te(nx,ny,nz);   
	//Te.setZero();
	//clock_t start_time1 = clock();
    
	c2rfft3d(nek,ne);
	c2rfft3d(Tik,Ti);
	c2rfft3d(Tek,Te);
   

	

	// Set scalar doubles for variables that are needed in the loops
	double Vthi, Vthe, lambdaD, Lambda;
	
	// Initialize 2D physical arrays for the collision frequencies. Will take fft of them all at the end
	//Eigen::Matrix< double, (ny+1), (nx)> nuin; 
	Eigen::Tensor<double, 3> nuin(nx,ny,nz);  
	//nuin.setZero();
	Eigen::Tensor<double, 3> nuii(nx,ny,nz);  
	//nuii.setZero();
	Eigen::Tensor<double, 3> nuie(nx,ny,nz);   
	//nuie.setZero();
	Eigen::Tensor<double, 3> nuen(nx,ny,nz);  
	//nuen.setZero();
	Eigen::Tensor<double, 3> nuei(nx,ny,nz);   
	//nuei.setZero();
	Eigen::Tensor<double, 3> nuee(nx,ny,nz);   
	//nuee.setZero();
	Eigen::Tensor<double, 3> isigP(nx,ny,nz);   
	//isigP.setZero();
	Eigen::Tensor<double, 3> invn(nx,ny,nz);   
	//invn.setZero();

	Eigen::Tensor<double, 3> hallE(nx,ny,nz);   
	//hallE.setZero();
	Eigen::Tensor<double, 3> hallI(nx,ny,nz);  
	//hallI.setZero();
    
	Eigen::Tensor<double, 3> test(nx,ny,nz);  

	// Begin big loop to calculating everything.
	//#pragma omp parallel for collapse(3)
    
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 

				Vthi = sqrt( 2. * kb * Ti(k,i,j) / mi);
				Vthe = sqrt( 2. * kb * Te(k,i,j) / me);
				nuin(k,i,j) =  nn * Vthi* EIGEN_PI * (ri+rn) * (ri+rn);
				nuen(k,i,j) = nn * Vthe * EIGEN_PI * rn * rn;
				// Calculate "inverse Pedersen conductivity"
				isigP(k,i,j) = 1.0 / (e * ( nuin(k,i,j)/ Oci + nuen(k,i,j)/ Oce) );
				// Calculate Debye length
				lambdaD = sqrt( eps0 * kb * Te(k,i,j)/ (ne(k,i,j) * e * e) );					
				// Calculate plasma parameter
				Lambda = 12.0 * EIGEN_PI * ne(k,i,j) * lambdaD * lambdaD * lambdaD;
				// Calculate electron-electron collision frequency
				nuee(k,i,j) = ne(k,i,j)* e * e *e *e * log(Lambda/3.0) / ( 2. * EIGEN_PI * eps0 * eps0 * me * me * Vthe * Vthe * Vthe);
				//nuee(k,i,j) = ne(k,i,j)* pow(e,4.0) * log(Lambda/3.0) / ( 2. * EIGEN_PI * eps0 * eps0 * me * me * pow(Vthe,3.0));

				// Calculate ion-ion collision frequency
				nuii(k,i,j)= nuee(k,i,j) * sqrt(me/mi);
				// Calculate ion-electron collision frequency
				nuie(k,i,j) = nuee(k,i,j) * 0.5 * me / mi;
			
				// Calculate electron-ion collision frequency
				nuei(k,i,j) = nuee(k,i,j);
				// Calculate the inverse of the density
				// inverse of ne in Fourier space (which is needed for several terms in the temperature equation )	
				invn(k,i,j) = 1.0 / ne(k,i,j);

				//*************************	Hall parameters: TEST forb the inertial function ****************************
				hallE(k,i,j) = nuen(k,i,j)/Oce;
				hallI(k,i,j) = nuin(k,i,j)/Oci;
				//test(k,i,j)= nuen(k,i,j)/Oce;



			}
		}
	}
	
	
	// Take FFTs of everything now
	r2cfft3d(nuin, nuink);// nuink is NOT READ-only 
	r2cfft3d(nuie, nuiek);
	r2cfft3d(nuii, nuiik);
	r2cfft3d(nuen, nuenk);
	//std::cout<<nuenk<<endl;
	r2cfft3d(nuei, nueik);
	r2cfft3d(nuee, nueek);
	r2cfft3d(invn, invnk);
	r2cfft3d(isigP, isigPk);
	//***********test for the inertial function****************
	r2cfft3d(hallE, hallEk);
	r2cfft3d(hallI, hallIk);
	

}
void EcalcV_ExBk3D(Eigen::Tensor<std::complex<double>, 3>& dphidxk,Eigen::Tensor<std::complex<double>, 3>& dphidyk,Eigen::Tensor<std::complex<double>, 3>& dphidzk, double B[], double B2, Eigen::Tensor<std::complex<double>, 3>& vexbkx, Eigen::Tensor<std::complex<double>, 3>& vexbky,Eigen::Tensor<std::complex<double>, 3>& vexbkz){
	 //
	 //#pragma omp parallel
	 //{
		#pragma omp parallel for collapse(3)
		//#pragma omp parallel for schedule(static) collapse(3)
	 	for(int i = 0; i< nx; i++){
			for(int j = 0; j< ny; j++){
				for(int k = 0; k< nz; k++){ 
					vexbkx(k,i,j) = (dphidzk(k,i,j) * B[1] - dphidyk(k,i,j) * B[2]) /B2; //(-1.0 * dphidyk(k,i,j) * B[2] + dphidzk(k,i,j) * B[1]) /B2;
					vexbky(k,i,j) =  (dphidxk(k,i,j) * B[2] - dphidzk(k,i,j) * B[0]) /B2;
					vexbkz(k,i,j) =  (dphidyk(k,i,j) * B[0] - dphidxk(k,i,j) * B[1]) /B2;
				}
			}
		}
	 //}
	 

	
}
void Ecalc_diamag3D(Eigen::Tensor<std::complex<double>, 3>& dpdxk,Eigen::Tensor<std::complex<double>, 3>& dpdyk,Eigen::Tensor<std::complex<double>, 3>& dpdzk, double B[], double B2, double qa, Eigen::Tensor<std::complex<double>, 3>& nak, Eigen::Tensor<std::complex<double>, 3>& diamagxk,Eigen::Tensor<std::complex<double>, 3>& diamagyk, Eigen::Tensor<std::complex<double>, 3>& diamagzk){
	Eigen::Tensor<std::complex<double>, 3> predivx(nx,ny,nz);   
	//predivx.setZero();
	Eigen::Tensor<std::complex<double>, 3> predivy(nx,ny,nz);    
	//predivy.setZero();
	Eigen::Tensor<std::complex<double>, 3> predivz(nx,ny,nz);    
	//predivz.setZero();
	
	#pragma omp parallel for schedule(static) collapse(3)
    
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				predivx(k,i,j) = (dpdyk(k,i,j) * B[2] - dpdzk(k,i,j) * B[1])/(B2*qa*-1);
				predivy(k,i,j) = (dpdzk(k,i,j) * B[0] - dpdxk(k,i,j) * B[2])/(B2*qa*-1);
				predivz(k,i,j) = (dpdxk(k,i,j) * B[1] - dpdyk(k,i,j) * B[0])/(B2*qa*-1);


			}
		}
	}
    
	Eigen::Tensor<double, 3> diamagx(nx,ny,nz); 
	//diamagx.setZero();
	Eigen::Tensor<double, 3> diamagy(nx,ny,nz);   
	//diamagy.setZero();
	Eigen::Tensor<double, 3> diamagz(nx,ny,nz);   
	//diamagz.setZero();

	c2rfft3d(predivx,diamagx);
	c2rfft3d(predivy,diamagy);
	c2rfft3d(predivz,diamagz);

	Eigen::Tensor<double, 3> ne(nx,ny,nz);  
	//ne.setZero();

	c2rfft3d(nak,ne);
	//divide in real space
	Eigen::Tensor<double, 3> div1(nx,ny,nz);   
	//div1.setZero();
	Eigen::Tensor<double, 3> div2(nx,ny,nz);  
	//div2.setZero();
	Eigen::Tensor<double, 3> div3(nx,ny,nz);  
	//div3.setZero();
	Eigen::Tensor<double, 3> test(nx,ny,nz);  //adding this fixes temp seg fault error?
    
	div1 = diamagx / ne;
	div2 = diamagy / ne;
	div3 = diamagz / ne;
   
	
	r2cfft3d(div1,diamagxk);
	r2cfft3d(div2,diamagyk);
	r2cfft3d(div3,diamagzk);
	
	
	
}
double Ecalc_dt3D(double U[], Eigen::Tensor<double, 3> vexbx, Eigen::Tensor<double, 3> vexby, Eigen::Tensor<double, 3> vexbz,Eigen::Tensor<double, 3> diamagxi,Eigen::Tensor<double, 3> diamagyi,  Eigen::Tensor<double, 3> diamagzi, Eigen::Tensor<double, 3> diamagxe, Eigen::Tensor<double, 3> diamagye, Eigen::Tensor<double, 3> diamagze,double cfl, double kmax, double maxdt){
		double *absArr;
	    absArr = (double*) fftw_malloc(nx*ny *sizeof(double));
		//memset(absArr, 42, (nx*(ny+1))* sizeof(double));
		memset(absArr, 42, nx*ny* sizeof(double));

		//you can change it later so that U is eigen and uses abs
		for (int i=0; i < sizee; i++){
			if (U[i] < 0){
				absArr[i] = abs(U[i]);
	
			}
		}

		double vMaxArr[10];  
		Eigen::Tensor<double, 0> AbsMaxAsTensor = vexbx.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor1 = vexby.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor2 = vexbz.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor3 = diamagxi.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor4 = diamagyi.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor5 = diamagzi.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor6 = diamagxe.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor7 = diamagye.abs().maximum();
		Eigen::Tensor<double, 0> AbsMaxAsTensor8 = diamagze.abs().maximum();

		vMaxArr[0] = AbsMaxAsTensor(0);//vexbx.maxCoeff();//max2D(vexbx); //max of 2d arrays
		vMaxArr[1] = AbsMaxAsTensor1(0); //vexby.array().maxCoeff();//max2D(vexby);
		vMaxArr[2] = AbsMaxAsTensor2(0);
		vMaxArr[3] = AbsMaxAsTensor3(0); //diamagxi.array().maxCoeff();//max2D(diamagxi);
		vMaxArr[4] = AbsMaxAsTensor4(0); //diamagyi.array().maxCoeff();//max2D(diamagyi);
		vMaxArr[5] = AbsMaxAsTensor5(0);
		vMaxArr[6] = AbsMaxAsTensor6(0); //diamagxe.array().maxCoeff();//max2D(diamagxe);
		vMaxArr[7] = AbsMaxAsTensor7(0); //diamagye.array().maxCoeff();//max2D(diamagye);
		vMaxArr[8] = AbsMaxAsTensor8(0);
		vMaxArr[9] = Arr1DMax(U,3); // absolute of u, neutral wind //vMaxArr[6] = Arr1DMax((U), 3);
		//vMaxArr[6] = Arr1DMax(absArr,3);
		//std::cout<<vMaxArr[0]<<endl;
		
		double max = Arr1DMax(vMaxArr, 10);
		//std::cout<<max<<endl;
		//Eigen::internal::set_is_malloc_allowed(false);
		double dt = cfl / (max * kmax); // added 
        //Eigen::internal::set_is_malloc_allowed(true);
		// try using min 
		//dt = std::min(dt, maxdt);
		
		if (dt < maxdt){
		return dt;
		}else {
			return maxdt;
		}
		
		fftw_free(absArr);
	return 0;
}
void Ecalc_residualn3D(Eigen::Tensor<std::complex<double>, 3>& vexbxk, Eigen::Tensor<std::complex<double>, 3>& vexbyk, Eigen::Tensor<std::complex<double>, 3>& vexbzk, Eigen::Tensor<std::complex<double>, 3>& nink, Eigen::Tensor<std::complex<double>, 3>& residnoutk, Eigen::Tensor<double, 3>& kx, Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz){
	
	Eigen::Tensor<std::complex<double>, 3> dninxk(nx,ny,nz);   
	//dninxk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dninyk(nx,ny,nz);   
	//dninyk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dninzk(nx,ny,nz);   
	//dninzk.setZero();

	Eigen::Tensor<std::complex<double>, 3> mult1(nx,ny,nz);   
	//mult1.setZero();
	Eigen::Tensor<std::complex<double>, 3> mult2(nx,ny,nz);  
	//mult2.setZero();
	Eigen::Tensor<std::complex<double>, 3> mult3(nx,ny,nz);    
	//mult3.setZero();
    //Eigen::internal::set_is_malloc_allowed(false);
	derivk3D(nink, kx, dninxk);
	derivk3D(nink, ky, dninyk);
	derivk3D(nink, kz, dninzk);
    
	
	convolve3D(vexbxk, dninxk, mult1);
    //Eigen::internal::set_is_malloc_allowed(true);
	convolve3D(vexbyk, dninyk, mult2);
	convolve3D(vexbzk, dninzk, mult3);
	residnoutk = mult1 + mult2 + mult3;
    

	
}
void EcalcPotSourcek_inertia3D(Eigen::Tensor<std::complex<double>, 3>& ne0k,Eigen::Tensor<std::complex<double>, 3>&  nek,Eigen::Tensor<std::complex<double>, 3>& dndx0k, Eigen::Tensor<std::complex<double>, 3>& dndy0k,Eigen::Tensor<std::complex<double>, 3>& dndz0k,Eigen::Tensor<std::complex<double>, 3>& dndxk , Eigen::Tensor<std::complex<double>, 3>& dndyk, Eigen::Tensor<std::complex<double>, 3>& dndzk ,Eigen::Tensor<std::complex<double>, 3>& dphidx0k , Eigen::Tensor<std::complex<double>, 3>& dphidy0k, Eigen::Tensor<std::complex<double>, 3>& dphidz0k,Eigen::Tensor<std::complex<double>, 3>& dphidx1k,Eigen::Tensor<std::complex<double>, 3>& dphidy1k, Eigen::Tensor<std::complex<double>, 3>& dphidz1k,Eigen::Tensor<std::complex<double>, 3>& Pi1k, Eigen::Tensor<std::complex<double>, 3>& Pe1k, double uxB[], double e, double Cm, Eigen::Tensor<std::complex<double>, 3>& hallEk, Eigen::Tensor<std::complex<double>, 3>& hallIk , Eigen::Tensor<std::complex<double>, 3>& vexbkx0,Eigen::Tensor<std::complex<double>, 3>& vexbky0,  Eigen::Tensor<std::complex<double>, 3>& vexbkz0, Eigen::Tensor<std::complex<double>, 3>& vexbkx,Eigen::Tensor<std::complex<double>, 3>& vexbky, Eigen::Tensor<std::complex<double>, 3>& vexbkz,Eigen::Tensor<double, 3>& kx,Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz,Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<std::complex<double>, 3>& potSourcek_inertia){

	
 	// get Laplacians of ion pressure:
	Eigen::Tensor<std::complex<double>, 3> d2Pi1k(nx,ny,nz);  
	//d2Pi1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2Pe1k(nx,ny,nz);  
	//d2Pe1k.setZero();

	laplaciank3D(Pi1k, ksqu, d2Pi1k); 
	laplaciank3D(Pe1k, ksqu, d2Pe1k); 
	

 	// calculate ion and electron pressure terms:
	Eigen::Tensor<std::complex<double>, 3> pikTerm(nx,ny,nz); 
	//pikTerm.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> pekTerm(nx,ny,nz);  
	//pekTerm.setZero();
	#pragma omp parallel for schedule(static) collapse(3)
    
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				pikTerm(k,i,j) = -1.0*(Cm/e) * hallIk(k,i,j);
				pekTerm(k,i,j) = (Cm/e) * hallEk(k,i,j);
			}
		}
	}
	
	//std::cout << std::fixed << std::setprecision(4);
	//std::cout<<pikTerm<<endl;
	
 	// First convolution for ions to multiply by Laplacian of pressure
	convolve3D(d2Pi1k, pikTerm, pikTerm); // this is = convolution2D of hallIk with -ksqu.*Pik	
	//aapx(d2Pi1k, pikTerm, pikTerm); 
	
	//std::cout<<pikTerm<<endl;
 	// First convolution for electrons to multiply by Laplacian of pressure
	convolve3D(d2Pe1k, pekTerm, pekTerm);
	//aapx(d2Pe1k, pekTerm, pekTerm);
	//std::cout<<pekTerm<<endl;

	Eigen::Tensor<std::complex<double>, 3> ne1k(nx,ny,nz);  
	//ne1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dndx1k(nx,ny,nz); 
	//dndx1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> dndy1k(nx,ny,nz); 
	//dndy1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> dndz1k(nx,ny,nz); 
	//dndz1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dphidxk(nx,ny,nz); 
	//dphidxk.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dphidyk(nx,ny,nz); 
	//dphidyk.setZero();

	Eigen::Tensor<std::complex<double>, 3> dphidzk(nx,ny,nz); 
	//dphidzk.setZero();

	Eigen::Tensor<std::complex<double>, 3> vexbkx1(nx,ny,nz); 
	//vexbkx1.setZero();

	Eigen::Tensor<std::complex<double>, 3> vexbky1(nx,ny,nz); 
	//vexbky1.setZero();

	Eigen::Tensor<std::complex<double>, 3> vexbkz1(nx,ny,nz);  
	//vexbkz1.setZero();

	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				ne1k(k,i,j) = nek(k,i,j) - ne0k(k,i,j);
				//ne1k[j + nyk*i][k] = nek[j + nyk*i][k]-ne0k[j + nyk*i][k]; //not sure about [j + nyk*i][k] 
				dndx1k(k,i,j) = dndxk(k,i,j) - dndx0k(k,i,j);
				//dndx1k[j + nyk*i][k] = dndxk[j + nyk*i][k]-dndx0k[j + nyk*i][k];
				dndy1k(k,i,j) = dndyk(k,i,j) - dndy0k(k,i,j);
				//dndy1k[j + nyk*i][k] = dndyk[j + nyk*i][k]-dndy0k[j + nyk*i][k];
				dndz1k(k,i,j) = dndzk(k,i,j) - dndz0k(k,i,j);
				//dndz1k[j + nyk*i][k] = dndzk[j + nyk*i][k]-dndz0k[j + nyk*i][k];

				dphidxk(k,i,j) = dphidx0k(k,i,j) + dphidx1k(k,i,j);
				//dphidxk[j + nyk*i][k] = dphidx0k[j + nyk*i][k]+dphidx1k[j + nyk*i][k];
				dphidyk(k,i,j) = dphidy0k(k,i,j) + dphidy1k(k,i,j);
				//dphidyk[j + nyk*i][k] = dphidy0k[j + nyk*i][k]+dphidy1k[j + nyk*i][k];
				dphidzk(k,i,j) = dphidz0k(k,i,j) + dphidz1k(k,i,j);
				//dphidzk[j + nyk*i][k] = dphidz0k[j + nyk*i][k]+dphidz1k[j + nyk*i][k];

				vexbkx1(k,i,j) = vexbkx(k,i,j) - vexbkx0(k,i,j);
				//vexbkx1[j + nyk*i][k] = vexbkx[j + nyk*i][k]-vexbkx0[j + nyk*i][k];
				vexbky1(k,i,j) = vexbky(k,i,j) - vexbky0(k,i,j);
				//vexbky1[j + nyk*i][k] = vexbky[j + nyk*i][k]-vexbky0[j + nyk*i][k];
				vexbkz1(k,i,j) = vexbkz(k,i,j) - vexbkz0(k,i,j);
				//vexbkz1[j + nyk*i][k] = vexbkz[j + nyk*i][k]-vexbkz0[j + nyk*i][k];

			}
		}
	}
	//test MUST be changed
	

	
	Eigen::Tensor<std::complex<double>, 3> div_n_nabphi_x(nx,ny,nz);  
	//div_n_nabphi_x.setZero();

	Eigen::Tensor<std::complex<double>, 3> div_n_nabphi_y(nx,ny,nz); 
	//div_n_nabphi_y.setZero();

	Eigen::Tensor<std::complex<double>, 3> div_n_nabphi_z(nx,ny,nz);  
	//div_n_nabphi_z.setZero();
	
  // divergence (n nabla phi) term:
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_x1(nx,ny,nz);  
	//div_dummy_x1.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_y1(nx,ny,nz); 
	//div_dummy_y1.setZero();

	Eigen::Tensor<std::complex<double>, 3> div_dummy_z1(nx,ny,nz);  
	//div_dummy_z1.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_x0(nx,ny,nz);  
	//div_dummy_x0.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_y0(nx,ny,nz); 
	//div_dummy_y0.setZero();

	Eigen::Tensor<std::complex<double>, 3> div_dummy_z0(nx,ny,nz); 
	//div_dummy_z0.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummyTerm(nx,ny,nz);  
	//div_dummyTerm.setZero();

 	//**********************************************************
	Eigen::Tensor<std::complex<double>, 3> div_dummy_dx1(nx,ny,nz);  
	//div_dummy_dx1.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_dy1(nx,ny,nz);  
	//div_dummy_dy1.setZero();

	Eigen::Tensor<std::complex<double>, 3> div_dummy_dz1(nx,ny,nz); 
	//div_dummy_dz1.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_dx0(nx,ny,nz); 
	//div_dummy_dx0.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> div_dummy_dy0(nx,ny,nz);  
	//div_dummy_dy0.setZero();

	Eigen::Tensor<std::complex<double>, 3> div_dummy_dz0(nx,ny,nz);  
	//div_dummy_dz0.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> coeff(nx,ny,nz);  
	//coeff.setZero();
	
 //Do the convolution of nek and dphi1k , to find the sum of final div dummy term seperate them like this?
    convolve3D(nek, dphidx1k, div_dummy_x1);
	convolve3D(nek, dphidy1k, div_dummy_y1);
	convolve3D(nek, dphidz1k, div_dummy_z1); 

	
	convolve3D(ne1k, dphidx0k, div_dummy_x0);
	convolve3D(ne1k, dphidy0k, div_dummy_y0);
	convolve3D(ne1k, dphidz0k, div_dummy_z0); 


 //Take derivatives of dummy variable:
 	derivk3D(div_dummy_x1, kx, div_dummy_dx1);  
	derivk3D(div_dummy_y1, ky, div_dummy_dy1);
	derivk3D(div_dummy_z1, kz, div_dummy_dz1);
	
	derivk3D(div_dummy_x0, kx, div_dummy_dx0); 
	derivk3D(div_dummy_y0, ky, div_dummy_dy0);
	derivk3D(div_dummy_z0, kz, div_dummy_dz0);
	
	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				div_dummyTerm(k,i,j) = div_dummy_x1(k,i,j) + div_dummy_y1(k,i,j) + div_dummy_z1(k,i,j) + div_dummy_x0(k,i,j) + div_dummy_y0(k,i,j) + div_dummy_z0(k,i,j) ;
				//div_dummyTerm[j + nyk*i][k] = div_dummy_x1[j + nyk*i][k]+div_dummy_y1[j + nyk*i][k] + div_dummy_z1[j + nyk*i][k]+div_dummy_x0[j + nyk*i][k]+div_dummy_y0[j + nyk*i][k] + div_dummy_z0[j + nyk*i][k] ; 
				coeff(k,i,j) = -1.*(Cm)*(hallEk(k,i,j)+ hallIk(k,i,j));
				//coeff[j + nyk*i][k] = -1.*(Cm)*(hallEk[j + nyk*i][k]+ hallIk[j + nyk*i][k]); //not sure how to do this 

				div_n_nabphi_x(k,i,j) = div_dummy_dx1(k,i,j) + div_dummy_dx0(k,i,j); 
				//div_n_nabphi_x[j + nyk*i][k] = div_dummy_dx1[j + nyk*i][k]+div_dummy_dx0[j + nyk*i][k]; //change x =x+..
				div_n_nabphi_y(k,i,j) = div_dummy_dy1(k,i,j) + div_dummy_dy0(k,i,j);
				//div_n_nabphi_y[j + nyk*i][k] = div_dummy_dy1[j + nyk*i][k]+div_dummy_dy0[j + nyk*i][k];
				div_n_nabphi_z(k,i,j) = div_dummy_dz1(k,i,j) + div_dummy_dz0(k,i,j);
				//div_n_nabphi_z[j + nyk*i][k] = div_dummy_dz1[j + nyk*i][k]+div_dummy_dz0[j + nyk*i][k];
			}
		}
	}
	convolve3D(div_n_nabphi_x,coeff, div_n_nabphi_x); 
	convolve3D(div_n_nabphi_y, coeff, div_n_nabphi_y); 
	convolve3D(div_n_nabphi_z, coeff, div_n_nabphi_z);  
 //
 // uXb Term: do the same thing with coeff here but different sign 
 //  double uxBTerm; to find convov of this term then it must be complex
  //(think of uxBTerm as a place holder being updated at each iteration)

	Eigen::Tensor<std::complex<double>, 3> coeff1(nx,ny,nz);  
	//coeff1.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> uxBTerm(nx,ny,nz); 
	//uxBTerm.setZero();
	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				uxBTerm(k,i,j) = uxB[0] * dndx1k(k,i,j) + uxB[1] * dndy1k(k,i,j) + uxB[2] * dndz1k(k,i,j);
				//uxBTerm[j + nyk*i][k] = uxB[0]* dndx1k[j + nyk*i][k] + uxB[1]*dndy1k[j + nyk*i][k] + uxB[2]*dndz1k[j + nyk*i][k] ;  
				coeff1(k,i,j) = (Cm)*(hallEk(k,i,j)+ hallIk(k,i,j));
				//coeff1[j + nyk*i][k] = (Cm)*(hallEk[j + nyk*i][k]+ hallIk[j + nyk*i][k]); 

			}
		}
	}
	
	convolve3D(uxBTerm, coeff1, uxBTerm);
	//std::array<long,3> offsetA = {0,0,0};        //Starting point:(row,column,matrix)
	//std::array<long,3> extentA = {(ny),nx,1}; //Finish point:(row,column,matrix)
	//std::array<long,2> shapeA = {ny,(nx)};
    //std::cout << uxBTerm.slice(offsetA, extentA).reshape(shapeA) << std::endl;  //Extract slice and reshape it into a 10x10 matrix.

 // *****************************************************
	Eigen::Tensor<std::complex<double>, 3> v_nabnab_phiTerm(nx,ny,nz);  
	//v_nabnab_phiTerm.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_x(nx,ny,nz);  
	//nabnab_dummyTerm_n0v0phi1_x.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_y(nx,ny,nz);  
	//nabnab_dummyTerm_n0v0phi1_y.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_z(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_z.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_yz(nx,ny,nz);  
	//nabnab_dummyTerm_n0v0phi1_yz.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_x(nx,ny,nz);  
	//nabnab_dummyTerm_n1vphi_x.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_y(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_y.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_z(nx,ny,nz);  
	//nabnab_dummyTerm_n1vphi_z.setZero();
	//cross terms
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_yz(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_yz.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_xy(nx,ny,nz);  
	//nabnab_dummyTerm_n1vphi_xy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_xz(nx,ny,nz);  
	//nabnab_dummyTerm_n1vphi_xz.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_yx(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_yx.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_zx(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_zx.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_zy(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_zy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi(nx,ny,nz);  
	//nabnab_dummyTerm_n1vphi.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_x(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_x.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_y(nx,ny,nz);  
	//nabnab_dummyTerm_n0v1phi_y.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_z(nx,ny,nz);  
	//nabnab_dummyTerm_n0v1phi_z.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi(nx,ny,nz);  
	//nabnab_dummyTerm_n0v1phi.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2phidx1k(nx,ny,nz);  
	//d2phidx1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2phidy1k(nx,ny,nz); 
	//d2phidy1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidz1k(nx,ny,nz);  
	//d2phidz1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2phidxk(nx,ny,nz); 
	//d2phidxk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidyk(nx,ny,nz); 
	//d2phidyk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidzk(nx,ny,nz); 
	//d2phidzk.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2phidxyk(nx,ny,nz); 
	//d2phidxyk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidxzk(nx,ny,nz);  
	//d2phidxzk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidzyk(nx,ny,nz);  
	//d2phidzyk.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2phidyxk(nx,ny,nz); 
	//d2phidyxk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidyzk(nx,ny,nz); 
	//d2phidyzk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidzxk(nx,ny,nz); 
	//d2phidzxk.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidxy1k(nx,ny,nz); 
	//d2phidxy1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidzx1k(nx,ny,nz);  
	//d2phidzx1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidzy1k(nx,ny,nz); 
	//d2phidzy1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> d2phidxz1k(nx,ny,nz);  
	//d2phidxz1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidyx1k(nx,ny,nz); 
	//d2phidyx1k.setZero();

	Eigen::Tensor<std::complex<double>, 3> d2phidyz1k(nx,ny,nz); 
	//d2phidyz1k.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_xy(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_xy.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_yx(nx,ny,nz);  
	//nabnab_dummyTerm_n0v0phi1_yx.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_xz(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_xz.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_zy(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_zy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_zx(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_zx.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_xy(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_xy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_xz(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_xz.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_zy(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_zy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_yx(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_yx.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_yz(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_yz.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_zx(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_zx.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_yy(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_yy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v0phi1_zz(nx,ny,nz); 
	//nabnab_dummyTerm_n0v0phi1_zz.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_yy(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_yy.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n1vphi_zz(nx,ny,nz); 
	//nabnab_dummyTerm_n1vphi_zz.setZero();

	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_zz(nx,ny,nz);  
	//nabnab_dummyTerm_n0v1phi_zz.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> nabnab_dummyTerm_n0v1phi_yy(nx,ny,nz); 
	//nabnab_dummyTerm_n0v1phi_yy.setZero();
  // first term:

 //Derivatives:
 	derivk3D(dphidx1k, kx, d2phidx1k); 
	derivk3D(dphidy1k, ky, d2phidy1k);
	derivk3D(dphidz1k, kz, d2phidz1k);
	
	derivk3D(dphidz1k, ky, d2phidzy1k); 
	derivk3D(dphidz1k, kx, d2phidzx1k);
	//cross deriv of dphidxk
	derivk3D(dphidx1k, ky, d2phidxy1k); 
	derivk3D(dphidy1k, kx, d2phidyx1k);
	derivk3D(dphidx1k, kz, d2phidxz1k);
	derivk3D(dphidy1k, kz, d2phidyz1k);


 // convolutions:
 	convolve3D(vexbkx0, d2phidx1k, nabnab_dummyTerm_n0v0phi1_x);
	convolve3D(vexbky0, d2phidy1k, nabnab_dummyTerm_n0v0phi1_y);
	convolve3D(vexbkz0, d2phidyz1k, nabnab_dummyTerm_n0v0phi1_yz);

	//cross terms convo
	convolve3D(vexbky0, d2phidxy1k, nabnab_dummyTerm_n0v0phi1_xy);
	convolve3D(vexbkx0, d2phidyx1k, nabnab_dummyTerm_n0v0phi1_yx);
	convolve3D(vexbkz0, d2phidxz1k,nabnab_dummyTerm_n0v0phi1_xz); 

	//cross terms convo
	convolve3D(vexbky0, d2phidzy1k, nabnab_dummyTerm_n0v0phi1_zy);
	convolve3D(vexbkx0, d2phidzx1k, nabnab_dummyTerm_n0v0phi1_zx);
	convolve3D(vexbkz0, d2phidz1k,nabnab_dummyTerm_n0v0phi1_z); 
	
 	// Calculate other n1vphi term:
 	//Derivatives:
	derivk3D(dphidxk, kx, d2phidxk); 
	derivk3D(dphidyk, ky, d2phidyk);
	derivk3D(dphidzk, kz, d2phidzk);
	

	//cross deriv of dphidxk
	derivk3D(dphidxk, ky, d2phidxyk); 
	derivk3D(dphidxk, kz, d2phidxzk);
	derivk3D(dphidyk, kx, d2phidyxk);
	derivk3D(dphidyk, kz, d2phidyzk);
	derivk3D(dphidzk, ky, d2phidzyk); 
	derivk3D(dphidzk, kx, d2phidzxk);
	
 // convolutions:
	convolve3D(vexbkx, d2phidxk, nabnab_dummyTerm_n1vphi_x);
	convolve3D(vexbky, d2phidyk, nabnab_dummyTerm_n1vphi_y);
	convolve3D(vexbkz, d2phidzk, nabnab_dummyTerm_n1vphi_z);
	
	//cross terms convo
	convolve3D(vexbky, d2phidxyk, nabnab_dummyTerm_n1vphi_xy);
	convolve3D(vexbkx, d2phidyxk, nabnab_dummyTerm_n1vphi_yx);
	convolve3D(vexbkz, d2phidyzk, nabnab_dummyTerm_n1vphi_yz);
	convolve3D(vexbkz, d2phidxzk, nabnab_dummyTerm_n1vphi_xz);
	convolve3D(vexbky, d2phidzyk, nabnab_dummyTerm_n1vphi_zy);
	convolve3D(vexbkx, d2phidzxk, nabnab_dummyTerm_n1vphi_zx);
	

 //  Calculate other n0v1phi term:
 // convolutions:
	convolve3D(vexbkx1, d2phidxk, nabnab_dummyTerm_n0v1phi_x);
	convolve3D(vexbky1, d2phidyk, nabnab_dummyTerm_n0v1phi_y);

	// cross terms convo 
	convolve3D(vexbky1, d2phidxyk, nabnab_dummyTerm_n0v1phi_xy);
	convolve3D(vexbkx1, d2phidyxk, nabnab_dummyTerm_n0v1phi_yx);
	convolve3D(vexbkz1, d2phidxzk, nabnab_dummyTerm_n0v1phi_xz);
	convolve3D(vexbkz1, d2phidyzk, 	nabnab_dummyTerm_n0v1phi_yz);
	
	convolve3D(vexbkz1, d2phidzk,nabnab_dummyTerm_n0v1phi_z);
	convolve3D(vexbky1, d2phidzyk, nabnab_dummyTerm_n0v1phi_zy);
	convolve3D(vexbkx1, d2phidzxk, nabnab_dummyTerm_n0v1phi_zx);

	//aapx(vexbky1, d2phidxyk, nabnab_dummyTerm_n0v1phi_xy);
	//aapx(vexbkx1, d2phidyxk, nabnab_dummyTerm_n0v1phi_yx);

	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				nabnab_dummyTerm_n0v0phi1(k,i,j) = (nabnab_dummyTerm_n0v0phi1_x(k,i,j) ) + ( nabnab_dummyTerm_n0v0phi1_xy(k,i,j) ) + ( nabnab_dummyTerm_n0v0phi1_xz(k,i,j) );
				//nabnab_dummyTerm_n0v0phi1[j + nyk*i][k] =  (nabnab_dummyTerm_n0v0phi1_x[j + nyk*i][k])+ (nabnab_dummyTerm_n0v0phi1_xy[j + nyk*i][k]) + (nabnab_dummyTerm_n0v0phi1_xz[j + nyk*i][k]) ; 
				nabnab_dummyTerm_n0v0phi1_yy(k,i,j) = (nabnab_dummyTerm_n0v0phi1_yx(k,i,j) ) + ( nabnab_dummyTerm_n0v0phi1_y(k,i,j) ) + ( nabnab_dummyTerm_n0v0phi1_yz(k,i,j) );
				//nabnab_dummyTerm_n0v0phi1_yy[j + nyk*i][k] =  (nabnab_dummyTerm_n0v0phi1_yx[j + nyk*i][k])+ (nabnab_dummyTerm_n0v0phi1_y[j + nyk*i][k]) + (nabnab_dummyTerm_n0v0phi1_yz[j + nyk*i][k]) ; 
				nabnab_dummyTerm_n0v0phi1_zz(k,i,j) = (nabnab_dummyTerm_n0v0phi1_zx(k,i,j) ) + ( nabnab_dummyTerm_n0v0phi1_zy(k,i,j) ) + ( nabnab_dummyTerm_n0v0phi1_z(k,i,j) );
				//nabnab_dummyTerm_n0v0phi1_zz[j + nyk*i][k] = (nabnab_dummyTerm_n0v0phi1_zx[j + nyk*i][k])+ (nabnab_dummyTerm_n0v0phi1_zy[j + nyk*i][k]) + (nabnab_dummyTerm_n0v0phi1_z[j + nyk*i][k]) ; 

				nabnab_dummyTerm_n1vphi(k,i,j) = (nabnab_dummyTerm_n1vphi_x(k,i,j) ) + ( nabnab_dummyTerm_n1vphi_xy(k,i,j) ) + ( nabnab_dummyTerm_n1vphi_xz(k,i,j) );
				//nabnab_dummyTerm_n1vphi[j + nyk*i][k] =  (nabnab_dummyTerm_n1vphi_x[j + nyk*i][k])+ (nabnab_dummyTerm_n1vphi_xy[j + nyk*i][k]) + (nabnab_dummyTerm_n1vphi_xz[j + nyk*i][k]) ;
				nabnab_dummyTerm_n1vphi_yy(k,i,j) = (nabnab_dummyTerm_n1vphi_yx(k,i,j) ) + ( nabnab_dummyTerm_n1vphi_y(k,i,j) ) + ( nabnab_dummyTerm_n1vphi_yz(k,i,j) );
				//nabnab_dummyTerm_n1vphi_yy[j + nyk*i][k] =  (nabnab_dummyTerm_n1vphi_yx[j + nyk*i][k])+ (nabnab_dummyTerm_n1vphi_y[j + nyk*i][k]) + (nabnab_dummyTerm_n1vphi_yz[j + nyk*i][k]) ;
				nabnab_dummyTerm_n1vphi_zz(k,i,j) = (nabnab_dummyTerm_n1vphi_zx(k,i,j) ) + ( nabnab_dummyTerm_n1vphi_zy(k,i,j) ) + ( nabnab_dummyTerm_n1vphi_z(k,i,j) );
				//nabnab_dummyTerm_n1vphi_zz[j + nyk*i][k] =  (nabnab_dummyTerm_n1vphi_zx[j + nyk*i][k])+ (nabnab_dummyTerm_n1vphi_zy[j + nyk*i][k]) + (nabnab_dummyTerm_n1vphi_z[j + nyk*i][k]) ;

				nabnab_dummyTerm_n0v1phi(k,i,j) = (nabnab_dummyTerm_n0v1phi_x(k,i,j) ) + ( nabnab_dummyTerm_n0v1phi_xy(k,i,j) ) + ( nabnab_dummyTerm_n0v1phi_xz(k,i,j) );
				//nabnab_dummyTerm_n0v1phi[j + nyk*i][k] =  (nabnab_dummyTerm_n0v1phi_x[j + nyk*i][k])+  (nabnab_dummyTerm_n0v1phi_xy[j + nyk*i][k]) +  (nabnab_dummyTerm_n0v1phi_xz[j + nyk*i][k]);  
				nabnab_dummyTerm_n0v1phi_yy(k,i,j) = (nabnab_dummyTerm_n0v1phi_yx(k,i,j) ) + ( nabnab_dummyTerm_n0v1phi_y(k,i,j) ) + ( nabnab_dummyTerm_n0v1phi_yz(k,i,j) );
				//nabnab_dummyTerm_n0v1phi_yy[j + nyk*i][k] =  (nabnab_dummyTerm_n0v1phi_yx[j + nyk*i][k])+ (nabnab_dummyTerm_n0v1phi_y[j + nyk*i][k]) + (nabnab_dummyTerm_n0v1phi_yz[j + nyk*i][k]); 
				nabnab_dummyTerm_n0v1phi_zz(k,i,j) = (nabnab_dummyTerm_n0v1phi_zx(k,i,j) ) + ( nabnab_dummyTerm_n0v1phi_zy(k,i,j) ) + ( nabnab_dummyTerm_n0v1phi_z(k,i,j) );
				// nabnab_dummyTerm_n0v1phi_zz[j + nyk*i][k] =  (nabnab_dummyTerm_n0v1phi_zx[j + nyk*i][k])+ (nabnab_dummyTerm_n0v1phi_zy[j + nyk*i][k]) + (nabnab_dummyTerm_n0v1phi_z[j + nyk*i][k]);  
			}
		}
	}

 // Multiply by n0:
	convolve3D(nabnab_dummyTerm_n0v0phi1, ne0k, nabnab_dummyTerm_n0v0phi1);
	convolve3D(nabnab_dummyTerm_n0v0phi1_yy, ne0k, nabnab_dummyTerm_n0v0phi1_yy);
	convolve3D(nabnab_dummyTerm_n0v0phi1_zz, ne0k, nabnab_dummyTerm_n0v0phi1_zz);


 // Multiply by n1:
	convolve3D(nabnab_dummyTerm_n1vphi, ne1k, nabnab_dummyTerm_n1vphi);
	convolve3D(nabnab_dummyTerm_n1vphi_yy, ne1k, nabnab_dummyTerm_n1vphi_yy);
	convolve3D(nabnab_dummyTerm_n1vphi_zz, ne1k, nabnab_dummyTerm_n1vphi_zz);

	

 // Multiply by n0:
	convolve3D(ne0k,nabnab_dummyTerm_n0v1phi ,nabnab_dummyTerm_n0v1phi);
	convolve3D(nabnab_dummyTerm_n0v1phi_yy, ne0k, nabnab_dummyTerm_n0v1phi_yy);
	convolve3D(nabnab_dummyTerm_n0v1phi_zz, ne0k, nabnab_dummyTerm_n0v1phi_zz);


	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v0phi1_x(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v0phi1_x.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v0phi1_y(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v0phi1_y.setZero();

	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v0phi1_z(nx,ny,nz);  
	//dnabnab_dummyTerm_n0v0phi1_z.setZero();

	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n1vphi_z(nx,ny,nz);  
	//dnabnab_dummyTerm_n1vphi_z.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n1vphi_x(nx,ny,nz); 
	//dnabnab_dummyTerm_n1vphi_x.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n1vphi_y(nx,ny,nz);  
	//dnabnab_dummyTerm_n1vphi_y.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v1phi_x(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v1phi_x.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v1phi_y(nx,ny,nz);  
	//dnabnab_dummyTerm_n0v1phi_y.setZero();

	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v1phi_z(nx,ny,nz);  
	//dnabnab_dummyTerm_n0v1phi_z.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v0phi1_xy(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v0phi1_xy.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v0phi1_yx(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v0phi1_yx.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n1vphi_xy(nx,ny,nz); 
	//dnabnab_dummyTerm_n1vphi_xy.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n1vphi_yx(nx,ny,nz); 
	//dnabnab_dummyTerm_n1vphi_yx.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v1phi_xy(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v1phi_xy.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> dnabnab_dummyTerm_n0v1phi_yx(nx,ny,nz); 
	//dnabnab_dummyTerm_n0v1phi_yx.setZero();


  //*****************************
   // try nabnab_dummyTerm_n0v0phi1
   	derivk3D(nabnab_dummyTerm_n0v0phi1, kx, dnabnab_dummyTerm_n0v0phi1_x);
	derivk3D(nabnab_dummyTerm_n0v0phi1_yy, ky, dnabnab_dummyTerm_n0v0phi1_y);
	derivk3D(nabnab_dummyTerm_n0v0phi1_zz, kz, dnabnab_dummyTerm_n0v0phi1_z);
	
	//cross terms
	derivk3D(nabnab_dummyTerm_n1vphi, kx, dnabnab_dummyTerm_n1vphi_x);
	derivk3D(nabnab_dummyTerm_n1vphi_yy, ky, dnabnab_dummyTerm_n1vphi_y);
	derivk3D(nabnab_dummyTerm_n1vphi_zz, kz, dnabnab_dummyTerm_n1vphi_z);
	
	// cross terms
	derivk3D(nabnab_dummyTerm_n0v1phi, kx, dnabnab_dummyTerm_n0v1phi_x); // it was kx
	
	derivk3D(nabnab_dummyTerm_n0v1phi_yy, ky, dnabnab_dummyTerm_n0v1phi_y);
	derivk3D(nabnab_dummyTerm_n0v1phi_zz, ky, dnabnab_dummyTerm_n0v1phi_z);
	
	//Eigen::MatrixXcd dummyPot((ny+1),nxk);
    //dummyPot.setZero();
	#pragma omp parallel for schedule(static) collapse(3)
	for(int i = 0; i< nx; i++){
		for(int j = 0; j< ny; j++){
			for(int k = 0; k< nz; k++){ 
				v_nabnab_phiTerm(k,i,j) = - (dnabnab_dummyTerm_n0v0phi1_x(k,i,j) +  dnabnab_dummyTerm_n0v0phi1_y(k,i,j)  +  dnabnab_dummyTerm_n0v0phi1_z(k,i,j) 
				+ dnabnab_dummyTerm_n1vphi_x(k,i,j) +  dnabnab_dummyTerm_n1vphi_y(k,i,j)  +  dnabnab_dummyTerm_n1vphi_z(k,i,j) 
				+ dnabnab_dummyTerm_n0v1phi_x(k,i,j) +  dnabnab_dummyTerm_n0v1phi_y(k,i,j)  +  dnabnab_dummyTerm_n0v1phi_z(k,i,j) );
				
				potSourcek_inertia(k,i,j) = v_nabnab_phiTerm(k,i,j) + div_n_nabphi_x(k,i,j) + div_n_nabphi_y(k,i,j) +  pikTerm(k,i,j) + pekTerm(k,i,j) + uxBTerm(k,i,j);


			}
		}
	}




}
int potentialk3D(Eigen::Tensor<std::complex<double>, 3>& invnk, Eigen::Tensor<std::complex<double>, 3>& dndxk, Eigen::Tensor<std::complex<double>, 3>& dndyk, Eigen::Tensor<std::complex<double>, 3>& dndzk,Eigen::Tensor<std::complex<double>, 3>& phik, Eigen::Tensor<std::complex<double>, 3>& potSourcek,Eigen::Tensor<double, 3>& kx, Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz, Eigen::Tensor<double, 3>&  ninvksqu, double err_max, int max_iter){ 

	Eigen::Tensor<std::complex<double>, 3> dphidxk(nx,ny,nz);  
	dphidxk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dphidyk(nx,ny,nz);  
	dphidyk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dphidzk(nx,ny,nz);  
	dphidzk.setZero();
	
	
	Eigen::Tensor<std::complex<double>, 3> gradNgradPhi_x(nx,ny,nz);   
	gradNgradPhi_x.setZero(); 
	Eigen::Tensor<std::complex<double>, 3> gradNgradPhi_y(nx,ny,nz);   
	gradNgradPhi_y.setZero(); 
	Eigen::Tensor<std::complex<double>, 3> gradNgradPhi_z(nx,ny,nz);   
	gradNgradPhi_z.setZero();
	
	Eigen::Tensor<std::complex<double>, 3> RHS(nx,ny,nz);   
	RHS.setZero();
	Eigen::Tensor<std::complex<double>, 3> gradNgradPhi(nx,ny,nz); 
	gradNgradPhi.setZero();
	
	double phik_max, phik_max_old;
	double it_error;

	// Begin counter for the number of iterations it takes to calculate phi
	int count = 0;
	
	// Begin while loop

	do{			
		// Calculate phi derivatives
		derivk3D(phik, kx, dphidxk);
		derivk3D(phik, ky, dphidyk);
		derivk3D(phik, kz, dphidzk);
		

		convolve3D(dndxk, dphidxk, gradNgradPhi_x);
		convolve3D(dndyk, dphidyk, gradNgradPhi_y);
		convolve3D(dndzk, dphidzk, gradNgradPhi_z);
		
	
		#pragma omp parallel for collapse(3) //num_threads(12) 
		//#pragma omp parallel for reduction(+: gradNgradPhi[:nx])
		for(int i = 0; i< nx; i++){
			for(int j = 0; j< ny; j++){
				for(int k = 0; k< nz; k++){ 
		
					RHS(k,i,j) = potSourcek(k,i,j) - gradNgradPhi_x(k,i,j) - gradNgradPhi_y(k,i,j) - gradNgradPhi_z(k,i,j);
				}
			}
		}
	
		// Convolve RHS with invnk
		convolve3D(RHS, invnk, RHS);	
		
		// Calculate maximum of absolute value of previous phi
        
		Eigen::Tensor<double, 0> AbsMaxAsTensor = phik.abs().maximum();
		phik_max_old = AbsMaxAsTensor(0); //phik.abs().maximum();
		
		phik = RHS * ninvksqu;

		// Calculate maximum of absolute value of updated phi(new phi)
		Eigen::Tensor<double, 0> AbsMaxAsTensor1 = phik.abs().maximum();
		phik_max = AbsMaxAsTensor1(0); //max_absComp(phik); //by
	
		count = count + 1; //count += 1;
	
		it_error = fabs((phik_max-phik_max_old)/phik_max);	//err_max is the error we want to converge to	
		
		
	 }while( it_error > err_max && count  <= max_iter && phik_max > err_max); // and instead &&
	
	 return count;

}
void Ecalc_sourcen3D(Eigen::Tensor<double, 3>& ksqu, Eigen::Tensor<std::complex<double>, 3>& nk, double d,Eigen::Tensor<std::complex<double>, 3>& sourcenk){
 
	Eigen::Tensor<std::complex<double>, 3> lapnk(nx,ny,nz);   
	//lapnk.setZero();
	laplaciank3D(nk,ksqu,lapnk);	
	sourcenk = d * lapnk;  //for loop?
    

}
void ERK4(Eigen::Tensor<std::complex<double>, 3>& f, double dt, Eigen::Tensor<std::complex<double>, 3>& residual, Eigen::Tensor<std::complex<double>, 3>& source, int stage, Eigen::Tensor<std::complex<double>, 3>& fout){
 	double alpha[4] = {1./4, 1./3, 1./2, 1.};
 
		for(int i = 0; i< nx; i++){
			for(int j = 0; j< ny; j++){
				for(int k = 0; k< nz; k++){ 
					fout(k,i,j) = f(k,i,j) - (alpha[stage] * dt * (residual(k,i,j) - source(k,i,j) ) );
					//fout(j + (ny+1)*i) =  f(j + (ny+1)*i);// - (alpha[stage] * dt * (residual(j + (ny+1)*i) - source(j + (ny+1)*i)));				

				}
			}
		}
		
}
void Ecalc_residualt3D(Eigen::Tensor<std::complex<double>, 3>& voxk, Eigen::Tensor<std::complex<double>, 3>& voyk, Eigen::Tensor<std::complex<double>, 3>& vozk,Eigen::Tensor<std::complex<double>, 3>& tempink, Eigen::Tensor<std::complex<double>, 3>& tempoutk, Eigen::Tensor<double, 3>& kx, Eigen::Tensor<double, 3>& ky, Eigen::Tensor<double, 3>& kz){
	
	Eigen::Tensor<std::complex<double>, 3> dtempinxk(nx,ny,nz);  
	//dtempinxk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dtempinyk(nx,ny,nz);  
	//dtempinyk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dtempinzk(nx,ny,nz);   
	//dtempinzk.setZero();
	
	
	Eigen::Tensor<std::complex<double>, 3> dvoxk(nx,ny,nz);  
	//dvoxk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dvoyk(nx,ny,nz);  
	//dvoyk.setZero();
	Eigen::Tensor<std::complex<double>, 3> dvozk(nx,ny,nz);   
	//dvozk.setZero();

	
	Eigen::Tensor<std::complex<double>, 3> divvo(nx,ny,nz);   
	//divvo.setZero();
	Eigen::Tensor<std::complex<double>, 3> mult1(nx,ny,nz);   
	//mult1.setZero();
	Eigen::Tensor<std::complex<double>, 3> mult2(nx,ny,nz);   
	//mult2.setZero();
	Eigen::Tensor<std::complex<double>, 3> mult3(nx,ny,nz);  
	//mult3.setZero();
	Eigen::Tensor<std::complex<double>, 3> mult4(nx,ny,nz);  
	//mult4.setZero();
	

	Eigen::Tensor<std::complex<double>, 3> dninxk(nx,ny,nz);  
	//dninxk.setZero();
	
	derivk3D(tempink, kx, dtempinxk);
	derivk3D(tempink, ky, dtempinyk);
	derivk3D(tempink, kz, dtempinzk);
	derivk3D(voxk, kx, dvoxk);
	derivk3D(voyk, ky, dvoyk);
	derivk3D(vozk, kz, dvozk);
	
	//Arr3DArr3DAdd(dvoxk, dvoyk, divvo);
	divvo = dvoxk + dvoyk + dvozk;
	
	convolve3D(dtempinxk, voxk, mult1);
	convolve3D(dtempinyk, voyk, mult2);
	convolve3D(dtempinzk, vozk, mult4);
	convolve3D(tempink, divvo, mult3);
	//rscalArr3DMult(mult3, 2/3, mult3);
	Eigen::Tensor<std::complex<double>, 3> dummy(nx,ny,nz);   
	//dummy.setZero();
	dummy = 2/3 * mult3; 
	
	for(int i = 0; i< nx; i++){
			for(int j = 0; j< ny; j++){
				for(int k = 0; k< nz; k++){ 
					tempoutk(k,i,j) = mult1(k,i,j) + mult2(k,i,j) + dummy(k,i,j) + mult4(k,i,j) ;
					//tempoutk(k,i,j) = mult1(k,i,j) + mult2(k,i,j) + mult3(k,i,j) + mult4(k,i,j) ;
				}
			}
	}
	
}
