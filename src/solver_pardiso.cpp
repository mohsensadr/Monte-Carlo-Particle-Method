/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on for unsymmetric linear systems                               */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Institute of Computational Science                 */
/*      Universita della Svizzera italiana, Lugano, Switzerland.        */
/*      Email: olaf.schenk@usi.ch                                       */
/* -------------------------------------------------------------------- */
//#include "cpp_headers.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
// #include "headers.h"
/* PARDISO prototype. */
// #include "llapack"

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <mkl_lapacke.h>
*/
#include "cpp_headers.h"
//int solving_phi(int *ia, int *ja, double *a, int n, double *b, double *x, int num_thr)
int solver_pardiso(int *ia, int *ja, double *a, int n, double *b, double *x)
{
  /* Matrix data. */
	//int n = 5;//MKL_INT n = 5;
	//int ia[ 6] = { 1, 4, 6, 9, 12, 14 };//	MKL_INT ia[ 6]
	//int ja[13] = { 1, 2, 4,
	//	1, 2,
	//	3, 4, 5,
	//	1, 3, 4,
	//	2, 5 };//MKL_INT ja[13]
	//double a[18] = { 1.0, -1.0, -3.0,
	//	-2.0, 5.0,
	//	4.0, 6.0, 4.0,
	//	-4.0, 2.0, 7.0,
	//	8.0, -5.0 };
	int mtype = 11; /* Real unsymmetric matrix */ //MKL_INT
	/* RHS and solution vectors. */
	//double b[5], x[5];
	int nrhs = 1; /* Number of right hand sides. */ //MKL_INT nrhs
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	int iparm[64];//MKL_INT iparm[64]
	int maxfct, mnum, phase, error, msglvl;//	MKL_INT maxfct
	/* Auxiliary variables. */
	int i;//MKL_INT i
	double ddum; /* Double dummy */
	int idum; /* Integer dummy. */ //MKL_INT idum
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 0; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	//printf("\nReordering completed ... ");
	//printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	//printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	//printf("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	/* Set right hand side to one. */
	//for (i = 0; i < n; i++) {
	//	b[i] = 1;
	//}
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, b, x, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	//printf("\nSolve completed ... ");
	//printf("\nThe solution of the system is: ");
	//for (i = 0; i < n; i++) {
	//	printf("\n x [%d] = % f", i, x[i] );
	//}
	//printf ("\n");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	return error;
}
