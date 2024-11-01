#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 512

float Mat[N][N], MatDD[N][N];
float V1[N], V2[N], V3[N];
void InitData(){

  int i,j;
  srand(334411);

  for( i = 0; i < N; i++ ){
    for( j = 0; j < N; j++ ){
      Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
      if ( (abs(i - j) <= 3) && (i != j)){
        MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
      }else{
        if ( i == j ){
          MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
        }else{
          MatDD[i][j] = 0.0;
        }
      }
    }
  }
  for( i = 0; i < N; i++ ){
    V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
    V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
    V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
  }
}

void PrintVect( float vect[N], int from, int numel ){
  for (int i = from; i < (from + numel); i ++){
    printf("%f", vect[i]);
    if ( i < (from + numel - 1) ) {
      printf(" ");
    }else{
      printf("\n");
    }
  }
}

void PrintRow( float mat[N][N], int row, int from, int numel ) {
  for (int  i=from; i < (from + numel); i++){
    printf("%f", mat[row][i]);
    if (  i < (from + numel - 1) ) {
      printf(" ");
    }else{
      printf("\n");
    }
  }
}

void  MultEscalar( float vect[N], float vectres[N], float alfa ){
  for (int i=0; i<N; i++){
    vectres[i] = vect[i] * alfa;
  }
}

float Scalar( float vect1[N], float vect2[N] ){
  float r = 0;
  for (int i=0; i<N; i++){
    r = r + ( vect1[i] * vect2[i] );
  }
  return r;
}

float Magnitude( float vect[N] ){
  float magn = 0;
  for (int i=0; i<N; i++){
    magn = magn + vect[i] * vect[i];
  }
  magn = sqrt(magn);
  return magn;
}

int Ortogonal( float vect1[N], float vect2[N] ){
  float producto;
  int r = 0;
  producto = Scalar( vect1, vect2 );
  if (producto == 0) {
  r = 1;
  }
  return r;
}

void  Projection( float vect1[N], float vect2[N], float vectres[N] ){
  float scalar_mult, magn_v2, div_scalar_magn;
  scalar_mult = Scalar(vect1,  vect2);
  magn_v2 = Magnitude(vect2);
  div_scalar_magn = scalar_mult / magn_v2;
  MultEscalar(vect2, vectres, div_scalar_magn);
}

float Infininorm( float M[N][N] ){
  float acumulador, inf_max;
  for (int i=0; i<N; i++){
    acumulador = 0;
    for (int j=0; j<N; j++){
      acumulador = acumulador + fabs(M[i][j]);
    }
    if ( acumulador > inf_max ){
      inf_max = acumulador;
    }
  }
  return inf_max;
}


float Onenorm( float M[N][N] ){
  float acumulador, u_max;
  for (int j=0; j<N; j++){
    acumulador = 0;
    for (int i=0; i<N; i++){
      acumulador = acumulador + fabs(M[i][j]);
    }
    if ( acumulador > u_max ){
      u_max = acumulador;
    }
  }
  return u_max;
}

float NormFrobenius( float M[N][N] ){
  float square_accumulator, result;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      square_accumulator = square_accumulator + pow(M[i][j], 2);
    }
  }
  result = sqrt(square_accumulator);
  return result;
}

int DiagonalDom( float M[N][N] ){
  float diagonal, matriu_exclosa;
  int b_diagonalDom = 0;
//  printf("INSIDE FUNC\n");
  for (int i=0; i<N; i++){
//    printf("INSIDE i LOOP\n");
    for (int j=0; j<N; j++){
//      printf("i=%d  j=%d\n", i, j);
      if (i==j){
       diagonal = diagonal + fabs( M[i][j] );
//       printf("%f  %f\n", diagonal, fabs( M[i][j] ));
      }else{
        matriu_exclosa = matriu_exclosa + fabs( M[i][j] );
      }
    }
  }
//  printf("%d  %f  %f\n", b_diagonalDom, diagonal, matriu_exclosa);
  if (diagonal >= matriu_exclosa){
    b_diagonalDom = 1;
  }
//  printf("%d\n", b_diagonalDom);
  return b_diagonalDom;
}

int main(){
  InitData();
  int a;
  a = DiagonalDom( MatDD );
  printf("%d\n", a);
}

