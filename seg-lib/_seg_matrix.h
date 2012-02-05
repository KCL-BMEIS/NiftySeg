#pragma once

#include <iostream>
#include <math.h>

using namespace std;

template <class D> class matrix{

  int size; //numb of columns
  int size2; //numb of rows
  D* data;
  void allocate()   {
    delete[] data;
    data = new D [size*size2];
  };
  matrix() {
    this->size = 0;
    this->size2 = 0;
    this->data = NULL;
  };

public:
  matrix(int tmpsize) {
    this->size2 = 0;
    this->size = 0;
    this->data = NULL;
    if (tmpsize <= 0) tmpsize = 5;
    this->size = tmpsize;
    this->size2 = tmpsize;
    this->data = NULL;
    allocate();
  };
  matrix(int tmpsize, int tmpsize2)  { // row, column
    this->size2 = 0;
    this->size = 0;
    this->data = NULL;
    if (tmpsize <= 0){ cout << "ERROR: Size <= 0" << endl;
        return;
      }
    if (tmpsize2 <= 0){ cout << "ERROR: Size <= 0" << endl;
        return;
      }
    this->size2 = tmpsize;
    this->size = tmpsize2;
    this->data = NULL;
    allocate();
  };
  ~matrix() { delete[] data; };
  void comparetoidentity()  {
    if(size==size2){
        int worstdiagonal = 0;
        D maxunitydeviation = 0.0;
        D currentunitydeviation;
        for ( int i = 0; i < size; i++ )  {
            currentunitydeviation = data[i*size+i] - 1.;
            if ( currentunitydeviation < 0.0) currentunitydeviation *= -1.;
            if ( currentunitydeviation > maxunitydeviation )  {
                maxunitydeviation = currentunitydeviation;
                worstdiagonal = i;
              }
          }
        int worstoffdiagonalrow = 0;
        int worstoffdiagonalcolumn = 0;
        D maxzerodeviation = 0.0;
        D currentzerodeviation ;
        for ( int i = 0; i < size; i++ )  {
            for ( int j = 0; j < size; j++ )  {
                if ( i == j ) continue;  // we look only at non-diagonal terms
                currentzerodeviation = data[i*size+j];
                if ( currentzerodeviation < 0.0) currentzerodeviation *= -1.0;
                if ( currentzerodeviation > maxzerodeviation )  {
                    maxzerodeviation = currentzerodeviation;
                    worstoffdiagonalrow = i;
                    worstoffdiagonalcolumn = j;
                  }

              }
          }
        cout << "Max diagonal Residual: "
             << maxunitydeviation << " at row/column " << worstdiagonal << endl;
        cout << "Max off-diagonal Residual: "
             << maxzerodeviation << " at row = " << worstoffdiagonalrow
             << ", column = " << worstoffdiagonalcolumn << endl;
      }
    else{cout << "ERROR: Compare to identity - Matrix is not Square"<< endl;}
  }
  void settoproduct(matrix& left, matrix& right)  {
    int lsize = left.getsize();
    int lsize2 = left.getsize2();
    int rsize = right.getsize();
    int rsize2 = right.getsize2();

    if (lsize!=rsize2){cout << "ERROR: Can't multiply. Wrong sizes/" << endl;
        return;
      }
    size2=lsize2; // numb rows
    size=rsize; // numb cols
    allocate();
    D leftvalue=0;
    D rightvalue=0;
    for ( int i = 0; i < size2; i++ ) // for each row
      for ( int j = 0; j < size; j++ )  { // for each column
          D sum = 0.0;
          bool success;
          for (int c = 0; c < lsize; c++)  {
              left.getvalue(i,c,leftvalue,success);
              right.getvalue(c,j,rightvalue,success);
              sum += leftvalue * rightvalue;
            }
          setvalue(i,j,sum);
        }
  }

  void dumpmatrix()  {

    for (int i=0; i < size2; i++) // curr row
      {
        for (int j=0; j<size; j++){ // curr column
            cout << data[i*size+j] << " ";}
        cout << endl;
      }
  };
  void copymatrix(matrix&  source)  {
    size = source.getsize();
    size2 = source.getsize2();
    allocate();
    D value=0;
    for ( int i = 0; i < size2; i++ )// curr row
      for ( int j = 0; j < size; j++ )  {// curr column
          bool success;
          source.getvalue(i,j,value,success); // curr_row, curr_column
          data[i*size+j] = value;
        }
  };

  void setsize(int newsize) {
    if (newsize > 0 )
      {
        size = newsize ;
        size2 = newsize;
        allocate();
      }
  };

  void setsize(int newsize, int newsize2) {
    if (newsize > 0 && newsize2 > 0)
      {
        size = newsize ;
        size2 = newsize2;
        allocate();
      }
  };

  int getsize() { return size; }

  int getsize2() { return size2; }

  void getvalue(int row, int column, D& returnvalue, bool& success)   {
    if ( (row>=size2) || (column>=size)
         || (row<0) || (column<0) )
      {  success = false;
        return;    }
    returnvalue = data[ row * size + column ];
    success = true;
  };
  bool setvalue(int row, int column, D newvalue)  {
    if ( (row >= size2) || (column >= size)
         || (row<0) || (column<0) ) return false;
    data[ row * size + column ] = newvalue;
    return true;
  };
  void invert()  {
    if (size!=size2){
        cout << "Matrix in not square" << endl;
        return;
      }
    if (size <= 1) return;
    for (int i=1; i < size; i++) {
        data[i] /= data[0]; // normalize row 0
      }
    for (int i=1; i < size; i++)  {
        for (int j=i; j < size; j++)  { // do a column of L
            D sum = 0.0;
            for (int k = 0; k < i; k++){
                sum += data[j*size+k] * data[k*size+i];
              }
            data[j*size+i] -= sum;
          }
        if (i == size-1) continue;
        for (int j=i+1; j < size; j++)  {  // do a row of U
            D sum = 0.0;
            for (int k = 0; k < i; k++){
                sum += data[i*size+k]*data[k*size+j];
              }
            data[i*size+j] =(data[i*size+j]-sum) / data[i*size+i];
          }
      }
    for ( int i = 0; i < size; i++ )  // invert L
      for ( int j = i; j < size; j++ )  {
          D x = 1.0;
          if ( i != j ) {
              x = 0.0;
              for ( int k = i; k < j; k++ ){
                  x -= data[j*size+k]*data[k*size+i];
                }
            }
          data[j*size+i] = x / data[j*size+j];
        }
    for ( int i = 0; i < size; i++ )   // invert U
      for ( int j = i; j < size; j++ )  {
          if ( i == j ){continue;}
          D sum = 0.0;
          for ( int k = i; k < j; k++ ){
              sum += data[k*size+j]*( (i==k) ? 1.0 : data[i*size+k] );
            }
          data[i*size+j] = -sum;
        }
    for ( int i = 0; i < size; i++ )   // final inversion
      for ( int j = 0; j < size; j++ )  {
          D sum = 0.0;
          for ( int k = ((i>j)?i:j); k < size; k++ ){
              sum += ((j==k)?1.0:data[j*size+k])*data[k*size+i];
            }
          data[j*size+i] = sum;
        }
  };
  double determinant(){

    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;

    if (size < 1) {

      } else if (size == 1) { /* Shouldn't get used */
        det = data[0];
      } else if (size == 2) {
        det = data[0+size*0] * data[1+size*1] - data[1+size*0] * data[0+size*1];
      } else {
        det = 0;
        for (j1=0;j1<size;j1++) {
            m = (double **) malloc((size-1)*sizeof(double *));
            for (i=0;i<size-1;i++)
              m[i] = (double *) malloc((size-1)*sizeof(double));
            for (i=1;i<size;i++) {
                j2 = 0;
                for (j=0;j<size;j++) {
                    if (j == j1)
                      continue;
                    m[i-1][j2] = data[i+size*j];
                    j2++;
                  }
              }
            det += pow(-1.0,1.0+j1+1.0) * data[0+size*j1] * Determinant_lib(m,size-1);
            for (i=0;i<size-1;i++)
              free(m[i]);
            free(m);
          }
      }
    return det;
  }

  double Determinant_lib(double **a,int n)
  {
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;

    if (n < 1) { /* Error */

      } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
      } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
      } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = (double **) malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
              m[i] = (double *) malloc((n-1)*sizeof(double));
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                    if (j == j1)
                      continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                  }
              }
            det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant_lib(m,n-1);
            for (i=0;i<n-1;i++)
              free(m[i]);
            free(m);
          }
      }
    return(det);
  }
};

