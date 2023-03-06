/*
 * main.cpp
 * 2D wave propagation with FDM
 * J.D. De Basabe (jonas@cicese.mx)
 */

#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <sys/stat.h>
#include <errno.h>

// Prototype of finite-difference function
void sgfd2d(int nx, int nz, int nt);

// Entry point, gets parameters form the command line
// argc is the number of parameters plus one
// argv is a vector of strings of the parameters
int main(int argc, char **argv)
{
	int nx=0, nz=0, nt=0,            // FD discretization
        err=0;                       // Error code
    char szName[]= {"./OUTPUT"};     // Output folder
    double tstart = omp_get_wtime(); // Get initial time

	switch( argc ) // Check how many parameters were set
	{
		case 2: // One parameter sets nx and assumes nz is equal
			sscanf(argv[1],"%i",&nx);
			nz= nx;
			break;

		case 3: // Two parameters, reads nx and nz
			sscanf(argv[1],"%i",&nx);
			sscanf(argv[2],"%i",&nz);
			break;

		case 4: // Three parameters, reads nx, nz and nt
			sscanf(argv[1],"%i",&nx);
			sscanf(argv[2],"%i",&nz);
			sscanf(argv[3],"%i",&nt);
			break;

		default: // No parameters, set default values for the FD grid
			nx= nz= 200;
			nt= 0;
	}

	// Create output folder
	err= mkdir(szName,0777); // Modificar para correr en Win
    if( err && errno!=EEXIST ) // Check for errors creating folder
    {    printf("ERROR creating output dir %s, error code = %i\n", szName, errno);
	}
    // Call the main FD method
	sgfd2d(nx,nz,nt);

    // Report the time
    printf("\nEllapsed time = %g\n",omp_get_wtime() - tstart);
	return 0;
}
