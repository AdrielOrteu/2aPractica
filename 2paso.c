
void PrintRow( float mat[N][N], int row, int from, int numel ) {
  for (i=from; i < (from + numel); i++){
    prinf("%f", mat[row][i];)
    if i < (from + numel - 1) {
      printf(" ");
    }else{
      printf("\n");
    }
  }
}
