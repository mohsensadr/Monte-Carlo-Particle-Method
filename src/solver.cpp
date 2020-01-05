#include "cpp_headers.h"

void solver( )
{
  double a[] = {1.0,2.0,3.0,4.0}; //NO need for column-major mode
  double b[] = {5.0, 11.0}; //NO need for column-major mode
  int n = 2;
  int nrhs = 1;
  int lda = n;
  int ipiv[n];
  int ldb = n;
  int info;
  dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

  int i;
  for(i=0; i<n; i++){
    printf("x[i]=%lf\n\n\n\n",b[i]);
  }
  //printf("info = %ld\n", info);
}
