#include<stdio.h>
#include<stdlib.h>



void printMatrix(double**,int,int);
double** transpose(double**,int,int);
double** multiply(double**,double**,int,int,int,int);
double** multiplyV(double**,double**,int,int);
double** invert(double**,int);


int main(int args, char** argv){
  FILE* train = NULL;
 
  // CHANGE BACK TO 3 FOR FINAL SUBMIT
  if(args < 3){
    printf("error\n");
    exit(0);
  }
  
  char* trainName = argv[1];
  
  int i = 0;
  int j = 0;
  // read training data from file
  train = fopen(trainName, "r");
  int k = 0;
  int n = 0;
  // number of columns/rows 
  fscanf(train, "%d\n%d\n", &k, &n);
  k = k + 1;
  
  // create training matrix
  double** trainMat = (double**)malloc(n * sizeof(double*));
  for(i = 0; i<n; i++){
    trainMat[i] = (double*)malloc(k * sizeof(double));
  }


  // create house price matrix [n X 1], and populate with last entry of each row
  double** Y = (double**)malloc(n * sizeof(double*));
  for(i=0; i<n; i++){
    Y[i] = (double*)malloc(1 * sizeof(double));
  }
  // populate training matrix
  for(i = 0; i<n; i++){
    for(j = 0; j<k+1; j++){
      // last entry -- must read new line character
      if(j == k){
	fscanf(train, "%lf\n", &Y[i][0]);
      }else if(j == 0){
	trainMat[i][j] = 1.0;
      }else{
	fscanf(train, "%lf,", &trainMat[i][j]);
      }
    }
  }
  
 
  double** t = (double**)malloc(n * sizeof(double*));
  t = transpose(trainMat,n,k);

  double** p = (double**)malloc(n * sizeof(double*));
  p = multiply(t,trainMat, k, n, n, k);

 
  double** r = (double**)malloc(k * sizeof(double*));
  r = invert(p, k);

  // should be 5x7
  double** pseudo = (double**)malloc(k * sizeof(double*));
  for(i=0;i<k;i++){
    pseudo[i] = (double*)malloc(n * sizeof(double));
  }
  
  pseudo = multiply(r, t, k, k, k, n);

  double** W = (double**)malloc(n * sizeof(double*));
  W = multiply(pseudo, Y, k, n, n, 1);

  
  //free matrices
  for(i=0;i<n;i++){
    free(trainMat[i]);
  }
  free(trainMat);
  free(t);
  free(p);




  // test file operations
  FILE* test = NULL;
  char* testName = argv[2];
  test = fopen(testName, "r");
  
  int m;
  fscanf(test, "%d\n", &m);

  double** testMat = (double**)malloc(m * sizeof(double*));
  for(i=0; i<m; i++){
    testMat[i] = (double*)malloc((k-1) * sizeof(double));
  }

  for(i=0; i<m; i++){
    for(j=0; j<k-1; j++){
      if(j == (k-1)){
	fscanf(test, "%lf\n", &testMat[i][j]);
      }else{
	fscanf(test, "%lf,", &testMat[i][j]);
      }
    }
  }
   

  double** P = (double**)malloc(m * sizeof(double**));
  for(i=0; i<m; i++){
    P[i] = (double*)malloc(1 * sizeof(double));
  }

  double price;
  for(i=0;i<m;i++){
    for(j=0;j<(k-1);j++){
      price += (testMat[i][j] * W[j+1][0]);
    }
    price += W[0][0];
    P[i][0] = price;
    price = 0.0;
  }

  printMatrix(P, m, 1);


  return 0;
}


// WORKING
// arguments are original matrix, number of rows, number of columns
double ** transpose(double ** mat, int n, int k){
  int i = 0;
  int j = 0;

  // resultant transposed matrix, opposite dimensions of the original matrix
  double** tmat = (double**)malloc(k * sizeof(double*));
  for(i = 0; i<k; i++){
    tmat[i] = (double*)malloc(n * sizeof(double));
  }
  // transpose matrix
  for(i = 0; i<n; i++){
    for(j = 0; j<k; j++){
      tmat[j][i] = mat[i][j];
    }
  }
  return tmat;
}

// n and k specify dimensions of resultant matrix
double** multiply(double** m1, double** m2, int n1, int k1, int n2,  int k2){
  int i = 0;
  int j = 0;
  int l = 0;
  double sum = 0.0;

  double** r = (double**)malloc(n1 * sizeof(double*));
  for(i = 0; i < n1; i++){
    r[i] = (double*)malloc(k2 * sizeof(double));
  }
 
  for(i=0; i < n1; i++){
    for(j=0; j < k2; j++){
      for(l=0; l < k1; l++){
	sum += ((m1[i][l] * m2[l][j]));
      }
      r[i][j] = sum;
      sum = 0.0;
    }
    sum = 0.0;
  }
  return r;
}


double** multiplyV(double** m, double** v, int n, int k){
  int i = 0;
  int l = 0;
  double sum = 0.0;

  //n:5 k:1

  double** f = (double**)malloc(n * sizeof(double*));
  for(i=0; i<n; i++){
    f[i] = (double*)malloc(1 * sizeof(double));
  }

  for(i=0; i < n; i++){
      for(l=0; l < n; l++){
	sum += ((m[i][l] * v[l][0]));
      }
      f[i][0] = sum;
      sum = 0;
  }

  return f;

}


double** invert(double** m, int n){
  int i = 0;
  int j = 0;
  int k = 0;
  double temp;

  // create identity matrix I
  double** I = (double**)malloc(n * sizeof(double*));
  for(i=0; i<n; i++){
    I[i] = (double*)malloc(n * sizeof(double));
  }

  // fill identity matrix
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      if(i == j){
	I[i][j] = 1;
      }else{
	I[i][j] = 0;
      }
    }
  }

  // create augmented matrix A -> W=2n, H=n
  double** A = (double**)malloc(n * sizeof(double**));
  for(i=0; i<n; i++){
    A[i] = (double*)malloc((2*n) * sizeof(double));
  }
  // concatenate m and I to A
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      A[i][j] = m[i][j];
    }
    for(k = n; k<(2*n); k++){
      A[i][k] = I[i][k-n];
    }
  }
 
  for(i=0; i<n; i++){
    if(A[i][i] != 1.0){
      temp = A[i][i];
      for(j=0; j<(2*n); j++){
	A[i][j] /= temp;
      }
    }
    
    for(j=(i+1); j<n; j++){
      if(A[j][i] != 0){
	temp = A[j][i];
	for(k=0; k<(2*n); k++){
	  A[j][k] += (-A[i][k] * temp); // ex: A[1][0] = A[1][0] + (-A[0][0] * A[1][0]);
	}
      }
    }
  }

  // reduced row echelon form
  for(i=(n-1); i>=0; i--){
    for(j=(i-1); j>=0; j--){
      if(A[j][i] != 0){
	temp = A[j][i];
	for(k=0; k<(2*n); k++){
	  A[j][k] += (-temp * A[i][k]);
	} 
      }
    }
  }


  double** inv = (double**)malloc(n * sizeof(double*));
  for(i=0; i<n; i++){
    inv[i] = (double*)malloc(n* sizeof(double));
  }

  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      inv[i][j] = A[i][j+n];
    }
  }

  return inv;
}




// WORKING
void printMatrix(double** matrix, int m, int k){
  int i = 0;
  int j = 0;

  for(i=0;i<m;i++){
    for(j=0;j<k;j++){
      printf("%0.0lf ", matrix[i][j]);
    }
    printf("\n");
  }
}
