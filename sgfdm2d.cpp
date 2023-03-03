/*
 * sgfdm2d.cpp
 * 2D elastic wave propagation exercise
 * J.D. De Basabe (jonas@cicese.mx)
 */

#include "macros.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// Macros to compute the x and z coordinates of grid nodes
#define XX(i) (XMIN + i*hx)
#define ZZ(i) (ZMIN + i*hz)
#define TT(i) (i*ht)

#define ZTOP(x) (0.5*(ZMIN+ZMAX))

//#define ZTOP(x) (ZMIN)

// P-wave velocity
double alpha(double x, double z)
{
	if( z<ZTOP(x) )
		return 1500.0;
	return 2000.0;
}

// S-wave velocity
double beta(double x, double z)
{
	if( z<ZTOP(x) )
		return 0.0;
	return 1000.0;
}

double rho(double x, double z)
{
	if( z<ZTOP(x) )
		return 1000.0;
	return 2000.0;
}

// Ricker wavelet
double Ricker(double t)
{
	double t0 = 0.78/PKFREQ;
	return -2.0*sqr(M_PI*PKFREQ)*(1.0-2.0*sqr(M_PI*PKFREQ*(t-t0)))
		*exp( -sqr(M_PI*PKFREQ*(t-t0)) );
}

// Gaussian distribution
double Gauss(double t, double mn, double sd)
{
	return exp( -sqr( (t-mn)/sd ) )/(sd*sqrt(2.0*M_PI));
}

// Source
double Source(double x, double z, double t)
{
	return Gauss(x,XSRC,XSW)*Gauss(z,ZSRC,ZSW)*Ricker(t);
}

// Wave propagation function
void sgfd2d(int nx, int nz, int nt)
{
	double hx, hz, ht, px, pz, l2m, l, m, b;
	int ix, iz, it;
	char fnsnapshot[256];
	if( nx<=1 || nz<=1 )
	{ // check for errors in input parameters
		printf("\nERROR in input parameters, nx=%i, nz=%i\n",nx,nz);
		return;
	}
	// Declare matrices for the pressure field
	dmat vx(nx,nz,0.0), vz(nx,nz,0.0), sxx(nx,nz,0.0), szz(nx,nz,0.0), sxz(nx,nz,0.0);

	printf("\n\tStaggered-Grid Finite-Diference 2D Elastic Wave Propagation\n");
	// Compute increments
	hx = (double) (XMAX - XMIN)/(nx-1.0);
	hz = (double) (ZMAX - ZMIN)/(nz-1.0);
	printf("Finite Difference Mesh = %i x %i\n",nx,nz);
	printf("Increments: hx = %g, hz = %g\n", hx, hz);
	if( nt<=0 )
	{// compute nt and ht if not provided
		ht= (hx<hz) ? 0.7*hx/VMAX : 0.7*hz/VMAX;  // get the smallest increments
		/*
        if( hx<hz )
            ht= 0.7*hx/VMAX;
        else
            ht= 0.7*hz/VMAX;
         */
		nt= TMAX/ht;
		if( nt*ht<TMAX )// make sure the time domain is fully covered
			++nt;
	}
	else// compute ht
		ht= (double) TMAX/nt;
	px= ht/hx;
	pz= ht/hz;
	printf("Time increment = %g, NT = %i\n", ht, nt);
	printf("Beginning time-stepping loop...\n");
	// Time-stepping loop - Not parallel
	for( it=1; it<=nt; ++it )
	{
		// Update velocity
#pragma omp parallel for private(ix, iz, b)
		for( ix=1; ix<(nx-1); ++ix )
			for( iz=1; iz<(nz-1); ++iz )
			{
				b  = 1.0/rho(XX(ix),ZZ(iz));

                vx[ix][iz]+= b*(px*(sxx[ix+1][iz] - sxx[ix][iz]) + pz*(sxz[ix][iz] - sxz[ix][iz-1]) )
                  + ht*NSX*Source(XX(ix),ZZ(iz),TT(it));
                vz[ix][iz]+= b*(pz*(szz[ix][iz+1] - szz[ix][iz]) + px*(sxz[ix][iz] - sxz[ix-1][iz]) )
                  + ht*NSZ*Source(XX(ix),ZZ(iz),TT(it));
			}

		// Update stress
#pragma omp parallel for private(ix, iz, b, m, l2m, l)
		for( ix=1; ix<(nx-1); ++ix )
			for( iz=1; iz<(nz-1); ++iz )
			{
				b  = 1.0/rho(XX(ix),ZZ(iz));
				m  = sqr(beta(XX(ix),ZZ(iz)))/b;
				l2m= sqr(alpha(XX(ix),ZZ(iz)))/b;
				l  = l2m - 2.0*m;

                sxx[ix][iz]+= l2m*px*(vx[ix][iz] - vx[ix-1][iz]) + l*pz*(vz[ix][iz] - vz[ix][iz-1]);
                szz[ix][iz]+= l*px*(vx[ix][iz] - vx[ix-1][iz]) + l2m*pz*(vz[ix][iz] - vz[ix][iz-1]);
                sxz[ix][iz]+= m*( px*(vz[ix+1][iz] - vz[ix][iz]) + pz*(vx[ix][iz+1] - vx[ix][iz]) );
			}

		// Save snapshots
		if( it%SAVEINT == 0 || it==nt)
		{
			sprintf(fnsnapshot,"./OUTPUT/ux%.8i.bin",it);
			vx.Save(fnsnapshot);
			printf("t = %8.4g, BIN file: %s\n", it*ht, fnsnapshot);
			sprintf(fnsnapshot,"./OUTPUT/uz%.8i.bin",it);
			vz.Save(fnsnapshot);
			printf("t = %8.4g, BIN file: %s\n", it*ht, fnsnapshot);
		}
	}
	printf("End of time-stepping loop.\n");
}
