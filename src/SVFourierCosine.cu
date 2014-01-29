/*
 * SVFourierCosine.cpp
 *
 *  Created on: Aug 28, 2013
 *      Author: Sabbir
 *  Rewritten on: Jan 19, 2014
 *     Author: Matthew Dixon
 */

#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <vector>
#include <cuComplex.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include<stdlib.h>

#define pi 3.1415926535897932384626433832795
#define rpart(x)   (cuCreal(x))
#define ipart(x)   (cuCimag(x))
#define cmplx(x,y) (make_cuDoubleComplex(x,y))
float totaltime 	= 0.0;
float counter 		= 1.0;
double s0 		= 0.0;
double r0 		= 0.0; 
double q0		= 0.0;
size_t num_blocks 	= 0;
size_t block_size 	= 256;
int num_bytes		= 0;
bool bLoad 		= false;
//host global arrays
double *h_K,*h_T,*h_W,*h_OP;
char *h_Types;
//device variables
double* d_sums		= 0;
double *d_input_p0 	= 0;
double *d_input_K 	= 0;
double *d_input_T 	= 0;
char *d_input_Types 	= 0;
char *model		= 0;

__device__ cuDoubleComplex squareRoot(cuDoubleComplex x)

{
	cuDoubleComplex inc;
	cuDoubleComplex c = x;
	cuDoubleComplex r = c;
	
	for(int j=0; j < 10; j++)
	{
		inc = cuCadd(r,r);
		inc = cuCdiv((cuCsub(cuCmul(r,r), c)),inc);
		r = cuCsub(r,inc);
	}
	
	return r;
}

__device__ double carg(const cuDoubleComplex& z) {return (double)atan2(ipart(z), rpart(z));} // polar angle
__device__ double cabs(const cuDoubleComplex& z) {return (double)cuCabs(z);}

__device__ cuDoubleComplex cPow(const cuDoubleComplex& z, const int &n)
{
	return cmplx((pow(cabs(z), n)*cos(n*carg(z))), (pow(cabs(z), n)*sin(n*carg(z))));
}

__device__ cuDoubleComplex cPow(const cuDoubleComplex& z, const double &n)
{
	return cmplx((pow(cabs(z), n)*cos(n*carg(z))), (pow(cabs(z), n)*sin(n*carg(z))));
}
__device__ cuDoubleComplex my_complex_exp (cuDoubleComplex arg)
{
   cuDoubleComplex res;
   double s, c;
   double e = exp(arg.x);
   sincos(arg.y, &s, &c);
   res.x = c * e;
   res.y = s * e;
   return res;
}

__device__ cuDoubleComplex HestonCF(
        double u,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho
        ) 
{
	cuDoubleComplex j1={0.0,1.0};
        double a 			= lmbda*meanV;
        double b 			= lmbda;
        double sigma2 		        = sigma*sigma;
	cuDoubleComplex d   		= squareRoot(cuCadd(cPow((cuCsub(cuCmul(j1,cmplx(rho*sigma*u,0)),cmplx(b,0))),2),cuCmul((cuCadd(cmplx(u*u,0),cuCmul(j1,cmplx(u,0)))),cmplx(sigma2,0))));
	cuDoubleComplex g   		= cuCdiv(cuCsub(cmplx(b,0),cuCadd(cuCmul(j1,cmplx(rho*sigma*u,0)),d)), cuCadd(cuCsub(cmplx(b,0),cuCmul(j1,cmplx(rho*sigma*u,0))),d));
	cuDoubleComplex ret 		= my_complex_exp(cuCmul(j1, make_cuDoubleComplex(u*(r-q)*T,0)));	
	cuDoubleComplex temp2 		= cuCmul(cuCsub(cuCsub(cmplx(b,0),cuCmul(cmplx(rho*sigma*u,0),j1)),d),cmplx(T,0));   //  (b - rho*j1*sigma*u - d)*T
	double root_sqr 		= cabs(cuCdiv(cuCsub(cmplx(1.0,0), cuCmul(g,my_complex_exp(cuCmul(cuCsub(cmplx(0,0),d),cmplx(T,0))))),cuCsub(cmplx(1.0,0),g)));
	double logarithm		= log(root_sqr);
	double arg 			= carg(cuCdiv(cuCsub(cmplx(1.0,0), cuCmul(g,my_complex_exp(cuCmul(cuCsub(cmplx(0,0),d),cmplx(T,0))))),cuCsub(cmplx(1.0,0),g)));
	cuDoubleComplex temp3 		= cuCmul(cmplx(2.0,0),cuCadd(cmplx(logarithm,0),cmplx(0,arg)));  //  2.0*log((1.0-g*exp(-d*T))/(1.0-g));	
	cuDoubleComplex temp1 		= cuCsub(temp2,temp3);   // 	((b - rho*j1*sigma*u - d)*T - 2.0*log((1.0-g*exp(-d*T))/(1.0-g)))
	ret 				= cuCmul(ret,my_complex_exp(cuCmul(cmplx(a/sigma2,0),temp1)));
	temp1 				= cmplx(v0/sigma2,0);
	temp2 				= cuCsub(cuCsub(cmplx(b,0),cuCmul(cmplx(rho*sigma*u,0),j1)),d);    //  (b - rho*j1*sigma*u - d)
	temp3 				= cuCsub(cmplx(1.0,0),my_complex_exp(cuCmul(cuCsub(cmplx(0,0),d),cmplx(T,0))));  //   (1.0-exp(-d*T))
	cuDoubleComplex temp4 		= cuCsub(cmplx(1.0,0),cuCmul(g,my_complex_exp(cuCmul(cuCsub(cmplx(0,0),d),cmplx(T,0))))); //   (1.0-g*exp(-d*T))
	temp1 				= cuCmul(temp1, temp2);
	temp1 				= cuCmul (temp1,temp3);
	temp1 				= cuCdiv(temp1,temp4);
	temp1 				= my_complex_exp(temp1);
        return cuCmul(ret,temp1);
}

__device__ cuDoubleComplex BatesCF(
        double u,
        double T,
        double r,
        double q,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho,
        double a, 
        double b,
        double lmbda_prime
        ) 
{
        cuDoubleComplex HCF 	= HestonCF(u,T,r,q,sigma,lmbda,meanV,v0,rho);
	cuDoubleComplex j1	= {0.0,1.0};
        cuDoubleComplex temp1 	= cuCmul(cmplx(lmbda_prime*T*-a*u,0),j1);
        cuDoubleComplex uj 	= cuCmul(cmplx(u,0),j1);
        cuDoubleComplex temp2 	= cuCsub(uj,cmplx(1.0,0));
        cuDoubleComplex temp3  	= cuCmul(cmplx(0.5*b*b,0),uj);
       	temp3 			= cuCmul(temp3,temp2);
        temp3 			= cuCadd(cuCmul(uj,cmplx(log(1.0+a),0)),temp3);
       	temp3  		        = cuCsub(my_complex_exp(temp3),cmplx(1.0,0));
       	temp3 			= cuCmul(cmplx(lmbda_prime*T,0),temp3);
       	temp3 			= cuCmul(cmplx(lmbda_prime*T,0),temp3);
       	temp3 			= my_complex_exp(cuCadd(temp1,temp3));

        //emp = my_complex_exp(lmbda_prime*T*(-a*u*j + (exp(u*j*log(1.0+a)+0.5*b*b*u*j*(u*j-1.0))-1.0)));
        return (cuCmul(HCF,temp3));
        
}
__device__ cuDoubleComplex VGCF(
        double u,
        double T,
        double r,
        double q,
        double sigma,
        double theta,
        double nu
        ) 
{
	cuDoubleComplex j1	= {0.0,1.0};
        double sigma2 = sigma*sigma;
        double omega = (1.0/nu)*(log(1.0-theta*nu-sigma2*nu/2.0));
        double tmp = 1.0 + 0.5*sigma2*u*u*nu;
        //tmp <- 1.0 - j*theta*nu*u + 0.5*sigma*sigma*u*u*nu
        cuDoubleComplex temp= cuCsub(cmplx(tmp,0),cuCmul(cmplx(theta*nu*u,0),j1));
        temp = cPow(temp, T/nu); 
        temp = cuCdiv(my_complex_exp(cuCmul(cmplx((r+omega-q)*u*T,0),j1)),temp);
        //ret <- exp(j*u*(r + omega - q)*T   - T*log(tmp)/nu)
        //= exp(j*u*(r+omega-q)*T)/tmp**(T/nu)
        return temp;
}

__device__ cuDoubleComplex CGMYCF(
        double u,
        double T,
        double r,
        double q,
        double C,
        double G,
        double M, 
        double Y 
        ) 
{
	cuDoubleComplex j1	= {0.0,1.0};
        double tg = tgamma(-Y);
        double m = -C*tg*(pow(M-1.0,Y)-pow(M,Y)+pow((G+1.0),Y)-pow(G,Y));
        cuDoubleComplex uj = cuCmul(cmplx(u,0),j1);
        cuDoubleComplex tp= cuCsub(cPow(cuCsub(cmplx(M,0),uj),Y),cmplx(pow(M,Y),0));
        cuDoubleComplex tn= cuCsub(cPow(cuCadd(cmplx(G,0),uj),Y),cmplx(pow(G,Y),0));
        cuDoubleComplex temp = cuCmul(cmplx(C*T*tg,0),cuCadd(tp,tn)); 
        //tmp = C*T*gamma(-Y)*((M-j*u)**Y-M**Y+(G+j*u)**Y-G**Y)
        return(my_complex_exp(cuCadd(cuCmul(uj, cmplx((r-q+m)*T,0)), temp))); 
        // ret <- exp(j*u*(r-q+m)*T + tmp)
}
__device__ double xi(
	double k,
	double a,
	double b,
	double c,
	double d) 
{
    double ret = 1.0/(1.0+pow(k*pi/(b-a),2))*(cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a)*sin(k*pi*(d-a)/(b-a))*exp(d)-k*pi/(b-a)*sin(k*pi*(c-a)/(b-a))*exp(c));
    return ret;
}

__device__ double psi(
	double k,
	double a,
	double b,
	double c,
	double d)
{
    double ret = 0.0;
	
    if (k==0)
	ret = d-c;
    else
	ret = (sin(k*pi*(d-a)/(b-a))-sin(k*pi*(c-a)/(b-a)))*(b-a)/(k*pi);
       
    return ret;
}

__global__ void HestonCOS(
		double* p0, 	//parameters
		double S, 	//option underlying
		double r0,	//risk free inst. short rate
		double q0,	//annual dividend rate
		double* K, 	//option strike
		double* T, 	//option maturity
		char* Types, 	//option type {'C','P'}
	        bool bIncorporateNLContraint,
	        double* results)
{	
    extern __shared__ cuDoubleComplex sdata[];
    double kappa, theta, sigma, rho, v0;
    double U, unit = 0.5;
    kappa = p0[0]; theta =p0[1]; sigma = p0[2]; rho = p0[3]; v0=p0[4];
	
    if (bIncorporateNLContraint)
	kappa = (kappa + sigma*sigma)/(2.0*theta);
      
    double lmbda = kappa, meanV = theta;
    cuDoubleComplex j1 = {0.0, 1.0};
    double sigma2 = sigma*sigma;
    double lmbda2 = lmbda*lmbda;	
    double c1 = r0*T[blockIdx.x]+(1-exp(-lmbda*T[blockIdx.x]))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T[blockIdx.x];	
    double c2 = 1.0/(8.0*lmbda2*lmbda)*(sigma*T[blockIdx.x]*lmbda*exp(-lmbda*T[blockIdx.x])*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-exp(-lmbda*T[blockIdx.x]))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T[blockIdx.x]*(-4.0*lmbda*rho*sigma+sigma2+4.0*lmbda2)+sigma2*((meanV-2.0*v0)*exp(-2.0*lmbda*T[blockIdx.x])+meanV*(6.0*exp(-lmbda*T[blockIdx.x])-7.0)+2.0*v0)+8.0*lmbda2*(v0-meanV)*(1-exp(-lmbda*T[blockIdx.x])));
    double a = c1-12.0*sqrt(fabs(c2));
    double b = c1+12.0*sqrt(fabs(c2));
    double x = log(S/K[blockIdx.x]);
   
    int tx = threadIdx.x;	
	
    if (Types[blockIdx.x] == 'C')
	U = 2.0/(b-a)*(xi(tx,a,b,0,b) - psi(tx,a,b,0,b));
    else
	U = 2.0/(b-a)*(-xi(tx,a,b,a,0) + psi(tx,a,b,a,0));

    cuDoubleComplex HCF = HestonCF(tx*pi/(b-a),T[blockIdx.x],r0,q0,sigma,lmbda,meanV,v0,rho);
	
    if(threadIdx.x == 0)
    {
	sdata[tx] = cuCmul(HCF, cuCmul(cmplx(unit*U, 0),my_complex_exp(cuCdiv(cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)),cmplx((b-a),0)))));    // unit*HCF*exp(j1*double(k)*pi*(x-a)/(b-a))*U;		
    }
    else
    {
    	unit = 1.0;
	sdata[tx] = cuCmul(cuCmul(cuCmul(cmplx(unit,0),HCF),my_complex_exp(cuCdiv (cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)), cuCsub(cmplx(b,0),cmplx(a,0))))), cmplx(U,0));						
    }
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
		
		if(tx < offset)
		{
			sdata[tx] = cuCadd(sdata[tx], sdata[tx+ offset]);
		}
		__syncthreads();
    }
	
    if(threadIdx.x == 0)
    {
       results[blockIdx.x] = K[blockIdx.x]*exp(-r0*T[blockIdx.x])*rpart(sdata[0]);   	   
    }
}

__global__ void BatesCOS(
		double* p0, 	//parameters
		double S, 	//option underlying
		double r0,	//risk free inst. short rate
		double q0,	//annual dividend rate
		double* K, 	//option strike
		double* T, 	//option maturity
		char* Types, 	//option type {'C','P'}
	        bool bIncorporateNLContraint,
	        double* results)
{	
    extern __shared__ cuDoubleComplex sdata[];
    double kappa, theta, sigma, rho, v0, a_prime, b_prime, lmbda_prime;
    double U, unit = 0.5;
    kappa = p0[0]; theta =p0[1]; sigma = p0[2]; rho = p0[3]; v0=p0[4]; a_prime =p0[5]; b_prime = p0[6]; lmbda_prime = p0[7];
	
    if (bIncorporateNLContraint)
	kappa = (kappa + sigma*sigma)/(2.0*theta);
      
    double lmbda = kappa, meanV = theta;
    cuDoubleComplex j1 = {0.0, 1.0};
    double sigma2 = sigma*sigma;
    double lmbda2 = lmbda*lmbda;	
    double c1 = r0*T[blockIdx.x]+(1-exp(-lmbda*T[blockIdx.x]))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T[blockIdx.x];	
    double c2 = 1.0/(8.0*lmbda2*lmbda)*(sigma*T[blockIdx.x]*lmbda*exp(-lmbda*T[blockIdx.x])*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-exp(-lmbda*T[blockIdx.x]))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T[blockIdx.x]*(-4.0*lmbda*rho*sigma+sigma2+4.0*lmbda2)+sigma2*((meanV-2.0*v0)*exp(-2.0*lmbda*T[blockIdx.x])+meanV*(6.0*exp(-lmbda*T[blockIdx.x])-7.0)+2.0*v0)+8.0*lmbda2*(v0-meanV)*(1-exp(-lmbda*T[blockIdx.x])));
    double a = c1-12.0*sqrt(fabs(c2));
    double b = c1+12.0*sqrt(fabs(c2));
    double x = log(S/K[blockIdx.x]);
   
    int tx = threadIdx.x;	
	
    if (Types[blockIdx.x] == 'C')
	U = 2.0/(b-a)*(xi(tx,a,b,0,b) - psi(tx,a,b,0,b));
    else
	U = 2.0/(b-a)*(-xi(tx,a,b,a,0) + psi(tx,a,b,a,0));

    cuDoubleComplex BCF = BatesCF(tx*pi/(b-a),T[blockIdx.x],r0,q0,sigma,lmbda,meanV,v0,rho,a_prime, b_prime, lmbda_prime);
	
    if(threadIdx.x == 0)
    {
	sdata[tx] = cuCmul(BCF, cuCmul(cmplx(unit*U, 0),my_complex_exp(cuCdiv(cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)),cmplx((b-a),0)))));    // unit*HCF*exp(j1*double(k)*pi*(x-a)/(b-a))*U;		
    }
    else
    {
    	unit = 1.0;
	sdata[tx] = cuCmul(cuCmul(cuCmul(cmplx(unit,0),BCF),my_complex_exp(cuCdiv (cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)), cuCsub(cmplx(b,0),cmplx(a,0))))), cmplx(U,0));						
    }
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
		
		if(tx < offset)
		{
			sdata[tx] = cuCadd(sdata[tx], sdata[tx+ offset]);
		}
		__syncthreads();
    }
	
    if(threadIdx.x == 0)
    {
       results[blockIdx.x] = K[blockIdx.x]*exp(-r0*T[blockIdx.x])*rpart(sdata[0]);   	   
    }
}

__global__ void VGCOS(
		double* p0, 	//parameters
		double S, 	//option underlying
		double r0,	//risk free inst. short rate
		double q0,	//annual dividend rate
		double* K, 	//option strike
		double* T, 	//option maturity
		char* Types, 	//option type {'C','P'}
	        bool bIncorporateNLContraint,
	        double* results)
{	
    extern __shared__ cuDoubleComplex sdata[];
    double sigma, theta, nu; 
    double U, unit = 0.5;
    sigma = p0[0]; theta =p0[1]; nu = p0[2];
	
    cuDoubleComplex j1 = {0.0, 1.0};
    double sigma2 = sigma*sigma;
    double theta2 = theta*theta; 
    double nu2    = nu*nu; 
    double c1 = (r0+theta)*T[blockIdx.x];
    double c2 = (sigma2 + nu*theta2)*T[blockIdx.x];
    double c4 = 3.0*(sigma2*sigma2*nu +2.0*theta2*theta2*nu2*nu+4.0*sigma2*theta2*nu2)*T[blockIdx.x];
    double a = c1-10.0*sqrt(c2 + sqrt(c4));
    double b = c1+10.0*sqrt(c2 + sqrt(c4));
    double x = log(S/K[blockIdx.x]);
   
    int tx = threadIdx.x;	
	
    if (Types[blockIdx.x] == 'C')
	U = 2.0/(b-a)*(xi(tx,a,b,0,b) - psi(tx,a,b,0,b));
    else
	U = 2.0/(b-a)*(-xi(tx,a,b,a,0) + psi(tx,a,b,a,0));

    cuDoubleComplex VGCF_ = VGCF(tx*pi/(b-a),T[blockIdx.x],r0,q0,sigma,theta,nu);
	
    if(threadIdx.x == 0)
    {
	sdata[tx] = cuCmul(VGCF_, cuCmul(cmplx(unit*U, 0),my_complex_exp(cuCdiv(cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)),cmplx((b-a),0)))));    // unit*HCF*exp(j1*double(k)*pi*(x-a)/(b-a))*U;		
    }
    else
    {
    	unit = 1.0;
	sdata[tx] = cuCmul(cuCmul(cuCmul(cmplx(unit,0),VGCF_),my_complex_exp(cuCdiv (cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)), cuCsub(cmplx(b,0),cmplx(a,0))))), cmplx(U,0));						
    }
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
		
		if(tx < offset)
		{
			sdata[tx] = cuCadd(sdata[tx], sdata[tx+ offset]);
		}
		__syncthreads();
    }
	
    if(threadIdx.x == 0)
    {
       results[blockIdx.x] = K[blockIdx.x]*exp(-r0*T[blockIdx.x])*rpart(sdata[0]);   	   
    }
}
__global__ void CGMYCOS(
		double* p0, 	//parameters
		double S, 	//option underlying
		double r0,	//risk free inst. short rate
		double q0,	//annual dividend rate
		double* K, 	//option strike
		double* T, 	//option maturity
		char* Types, 	//option type {'C','P'}
	        bool bIncorporateNLContraint,
	        double* results)
{	
    extern __shared__ cuDoubleComplex sdata[];
    double sigma, theta, nu,Y; 
    double U, unit = 0.5;
    sigma = p0[0]; theta =p0[1]; nu = p0[2]; Y = p0[3];
	
    cuDoubleComplex j1 = {0.0, 1.0};
    double C = 1.0/nu;
    double sigma2 = sigma*sigma;
    double theta2 = theta*theta; 
    double G = theta/sigma2 + sqrt(theta2/(sigma2*sigma2) + 2.0/(nu*sigma2));
    double M = theta/sigma2 + sqrt(theta2/(sigma2*sigma2) + 2.0/(nu*sigma2));

    double c1 = r0*T[blockIdx.x] + C*T[blockIdx.x]*tgamma(1.0-Y)*(pow(M,Y-1) - pow(G,Y-1));
    double c2 = sigma2*T[blockIdx.x] + C*T[blockIdx.x]*tgamma(2.0-Y)*(pow(M,Y-2.0) + pow(G,Y-2.0));
    double c4 = C*T[blockIdx.x]*tgamma(4.0-Y)*(pow(M,Y-4) + pow(G,Y-4));
    double a = c1-10.0*sqrt(c2 + sqrt(c4));
    double b = c1+10.0*sqrt(c2 + sqrt(c4));
    double x = log(S/K[blockIdx.x]);
   
    int tx = threadIdx.x;	
	
    if (Types[blockIdx.x] == 'C')
	U = 2.0/(b-a)*(xi(tx,a,b,0,b) - psi(tx,a,b,0,b));
    else
	U = 2.0/(b-a)*(-xi(tx,a,b,a,0) + psi(tx,a,b,a,0));

    cuDoubleComplex CGMYCF_ = CGMYCF(tx*pi/(b-a),T[blockIdx.x],r0,q0,C,G,M,Y);
	
    if(threadIdx.x == 0)
    {
	sdata[tx] = cuCmul(CGMYCF_, cuCmul(cmplx(unit*U, 0),my_complex_exp(cuCdiv(cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)),cmplx((b-a),0)))));    // unit*HCF*exp(j1*double(k)*pi*(x-a)/(b-a))*U;		
    }
    else
    {
    	unit = 1.0;
	sdata[tx] = cuCmul(cuCmul(cuCmul(cmplx(unit,0),CGMYCF_),my_complex_exp(cuCdiv (cuCmul(j1,cmplx(pi*(x-a)*double(tx), 0)), cuCsub(cmplx(b,0),cmplx(a,0))))), cmplx(U,0));						
    }
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
		
		if(tx < offset)
		{
			sdata[tx] = cuCadd(sdata[tx], sdata[tx+ offset]);
		}
		__syncthreads();
    }
	
    if(threadIdx.x == 0)
    {
       results[blockIdx.x] = K[blockIdx.x]*exp(-r0*T[blockIdx.x])*rpart(sdata[0]);   	   
    }
}
#ifdef __cplusplus
  extern "C"
#endif

void copy_data (double* K, double* T, double* W, double* OP, char* Types,double s0_, double r0_, double q0_,int length){

    h_K =K; h_T=T; h_W=W; h_OP=OP;h_Types=Types;s0=s0_;r0=r0_;q0=q0_;
    num_blocks = length;  
    num_bytes = sizeof(double)*num_blocks;

    cudaMalloc((void**)&d_sums, num_bytes);
    cudaMalloc((void**)&d_input_K, num_bytes);
    cudaMalloc((void**)&d_input_T, num_bytes);
    cudaMalloc((void**)&d_input_Types, num_bytes);

    cudaMemcpy(d_input_K, h_K, num_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_input_T, h_T, num_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_input_Types, h_Types, sizeof(char) * num_blocks, cudaMemcpyHostToDevice);
	
    bLoad = true;
}

#ifdef __cplusplus
  extern "C"
#endif

void dealloc_data (){

    cudaFree(d_sums);
    cudaFree(d_input_K);
    cudaFree(d_input_T);
    cudaFree(d_input_Types);

    bLoad = false;
}

#ifdef __cplusplus
  extern "C"
#endif
void set_model(char* model_){
   model = model_;
}

#ifdef __cplusplus
  extern "C"
#endif
void set_block_size(int block_size_){
   block_size = (size_t)block_size_;
}

#ifdef __cplusplus
  extern "C"
#endif

double error_func (double *p0, int length)
{
    std::cout.precision(16);
    if (!bLoad){
       std::cerr<<"Error: call copy_data first to load option data"<<std::endl; 
       return 0.0;
    }
    //float elapsedTime;
    //cudaEvent_t start, stop;
    double rmse = 0.0; 		
    bool bIncorporateNLConstraint = true;
    double* h_result = new double[num_blocks];
    cudaMalloc((void**)&d_input_p0, sizeof(double) * length);
    cudaMemcpy(d_input_p0, &p0[0], sizeof(double) * length, cudaMemcpyHostToDevice);
	
    /*for(int z=0; z<5;z++)
	{
	   printf("p0[%d] = %f\n",z,p0[z]);
    }*/
	
    //cudaEventCreate(&start);
    //cudaEventRecord(start,0);
	
    if (strcmp(model,"Heston")==0)
        HestonCOS <<<num_blocks, block_size, (block_size*2) * sizeof(double)>>>(d_input_p0, s0, r0, q0, d_input_K, d_input_T, d_input_Types, bIncorporateNLConstraint, d_sums);
    else if (strcmp(model,"Bates")==0)
        BatesCOS <<<num_blocks, block_size, (block_size*2) * sizeof(double)>>>(d_input_p0, s0, r0, q0, d_input_K, d_input_T, d_input_Types, bIncorporateNLConstraint, d_sums);
    else if (strcmp(model,"VG")==0)
        VGCOS <<<num_blocks, block_size, (block_size*2) * sizeof(double)>>>(d_input_p0, s0, r0, q0, d_input_K, d_input_T, d_input_Types, bIncorporateNLConstraint, d_sums);
    else if (strcmp(model,"CGMY")==0)
        CGMYCOS <<<num_blocks, block_size, (block_size*2) * sizeof(double)>>>(d_input_p0, s0, r0, q0, d_input_K, d_input_T, d_input_Types, bIncorporateNLConstraint, d_sums);

    cudaMemcpy(h_result, d_sums, num_bytes, cudaMemcpyDeviceToHost);
	
    for(int index = 0 ; index < num_blocks ; index++)
    {
	rmse += pow(h_W[index]*(h_result[index] - h_OP[index]),2);
    }
	
    rmse = sqrt(rmse/num_blocks);
    //totaltime += elapsedTime;
    //std::cout<<"TotalTime= "<<totaltime<<std::endl;
    //std::cout<<"Iteration= "<<counter<<std::endl;
    counter++;
	
    delete [] h_result;
    cudaFree(d_input_p0);
	
    return rmse;
}
