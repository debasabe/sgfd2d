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
#define TT(i) (i*ht)

#define ZTOP(x) (ZMIN+50.0*(sin(2.0*M_PI*x/XMAX)+1.0))
// P-wave velocity
double alpha(double x, double z)
{
	return 1.0;
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
double Source(double x, double t)
{
	return Gauss(x,XSRC,XSW)*Ricker(t);
}

// Wave propagation function
void fdm1d(int nx, int nt)
{
	double hx, ht, p, c;
	int ix, it;
	char fnsnapshot[256];
	if( nx<=1 )
	{ // check for errors in input parameters
		printf("\nERROR in input parameters, nx=%i\n",nx);
		return;
	}
	// Declare matrices for the pressure field
	dvec u(nx,0.0);

	printf("\n\tFinite-Diference 1D Acoustic Wave Propagation\n");
	// Compute increments
	hx = (double) (XMAX - XMIN)/(nx-1.0);
	printf("Finite Difference Mesh: nx = %i, hx=%g\n",nx,hx);
	if( nt<=0 )
	{// compute nt and ht if not provided
		ht= hx/VMAX;
		nt= TMAX/ht;
		if( nt*ht<TMAX )// make sure the time domain is fully covered
			++nt;
	}
	else// compute ht
		ht= (double) TMAX/nt;
	px= ht/hx;
	printf("Time increment = %g, NT = %i\n", ht, nt);
	printf("Beginning time-stepping loop...\n");
	// Time-stepping loop
	for( it=1; it<=nt; ++it )
	{
		for( ix=1; ix<(nx-1); ++ix )
        {
        }

		// Save snapshots
		if( it%SAVEINT == 0 || it==nt)
		{
		}
	}
	printf("End of time-stepping loop.\n");
}
