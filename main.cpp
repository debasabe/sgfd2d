/*
 * wave.cpp
 * 2D acoustic wave propagation exercise
 * J.D. De Basabe (jonas@cicese.mx)
 */

#include <iostream>
#include <stdio.h>
#include <omp.h>


void sgfd2d(int nx, int nz, int nt);

int main(int argc, char **argv)
{
	int nx=0, nz=0, nt=0;
    double tstart = omp_get_wtime();
	switch( argc )
	{
		case 2:
			sscanf(argv[1],"%i",&nx);
			nz= nx;
			break;
		case 3:
			sscanf(argv[1],"%i",&nx);
			sscanf(argv[2],"%i",&nz);
			break;
		case 4:
			sscanf(argv[1],"%i",&nx);
			sscanf(argv[2],"%i",&nz);
			sscanf(argv[3],"%i",&nt);
			break;
		default:
			nx= nz= 200;
			nt= 0;
	}
	sgfd2d(nx,nz,nt);
    printf("\nEllapsed time = %g\n",omp_get_wtime() - tstart);
	return 0;
}
