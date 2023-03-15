/*
 * fdm1d.cpp
 * 1D acoustic wave propagation exercise
 * J.D. De Basabe (jonas@cicese.mx)
 */

#include "macros.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// Macros to compute the x and z coordinates of grid nodes
#define XX(i) (XMIN + i*hx)
#define TT(i) (i*ht)

// P-wave velocity
double alpha(double x)
{
    double v=VMAX;
    if( x<XMAX/5.0)
	{
        v=VMAX;
	}
	return v;
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
int fdm1d(int nx, int nt)
{
	double hx, ht, p;
	int ix, it;
    char image[100], fnsnapshot[]= "OUTPUT/out.bin";
	if( nx<=1 )
	{ // check for errors in input parameters
		printf("\nERROR in input parameters, nx=%i\n",nx);
		return 0;
	}
	// Declare matrices for the pressure field
	dmat u;

	printf("\n\tFinite-Diference 1D Acoustic Wave Propagation\n");
	// Compute increments
	hx = (double) (XMAX - XMIN)/(nx-1.0);
	printf("Finite Difference Mesh: nx = %i, hx=%g\n",nx,hx);
	if( nt<=0 )
	{// compute nt and ht if not provided
		ht= hx/VMAX; // Stability condition
		nt= TMAX/ht; // Careful with roundup
		if( nt*ht<TMAX )// make sure the time domain is fully covered
			++nt; // nt = nt + 1
	}
	else// compute ht
		ht= (double) TMAX/nt;
	p= ht/hx; // not exactly as defined in class
    u.Init(nx,nt+1,0.0); // u[m][n], m=0...nx-1; n=0...nt
	printf("Time increment = %g, NT = %i\n", ht, nt);
	printf("Beginning time-stepping loop...\n");
	// Time-stepping loop
	for( it=1; it<nt; ++it )
	{	// finite-differences stencil
		for( ix=1; ix<(nx-1); ++ix )
        {   // vel = alpha(ix*hx);
            u[ix][it+1]= 2.0*u[ix][it] - u[ix][it-1]
                + sqr(p*alpha(XX(ix)))*(u[ix+1][it] - 2.0*u[ix][it] + u[ix-1][it])
                + ht*Source(XX(ix),it*ht);
        }
	}

	if( !u.Save(fnsnapshot) )
    {
        printf("Error saving results to file %s\n",fnsnapshot);
        return 0;
    }
	printf("End of time-stepping loop.\n");
    sprintf(image, "ximage < %s n1=%i d1=%g d2=%g cmap=rgb1 title='c=%g' &",
            fnsnapshot,nt+1, ht, hx, VMAX);
    printf("%s\n", image);
    return system(image);
}

// argc = numero de parameters
// argv = parametros como lista de palabras
int main(int argc, char **argv)
{   // argv[0] = comando completo, e.g. 'fdm1d 101 600'
    int nx=0, nt=0;
    switch( argc )
    {
        case 2:
            sscanf(argv[1],"%i",&nx);
            break;
        case 3:
            sscanf(argv[1],"%i",&nx);
            sscanf(argv[2],"%i",&nt);
            break;
        default:
            nx= 201;
            nt= 0;
            break;
    }
    return fdm1d(nx,nt);
}

