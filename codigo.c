#include <stdio.h>
#include <stdlib.h>

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

int main(){
  InitData();
  float res_vect[N];
  float escalar = 3.0;
  PrintVect(V1, 0, 10);
  PrintVect(V1, 256, 10);
  PrintVect(V2, 0, 10);
  PrintVect(V2, 256, 10);
  PrintVect(V3, 0, 10);
  PrintVect(V3, 256, 10);
  PrintRow(Mat, 0, 0, 10);
  PrintRow(Mat, 100, 0, 10);
  MultEscalar(V1, res_vect, escalar);
  escalar = Scalar(V1, V2);
  printf("%f\n", escalar); //Hauria de donar zero
}

