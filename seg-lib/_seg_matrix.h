#pragma once

#include <iostream>
#include <cmath>

using namespace std;

#include "NiftySegWinExportHeader.h"


/// @brief My own very simple implementation of a matrix class.
///
/// Did this mostly for fun. This class will probably be replaced with Eigen, since Eigen is better and it has been recently added as a dependency.
template <class matrixDataType> class seg_Matrix
{

    int size; ///< @brief Number of columns
    int size2; ///< @brief Number of rows
    matrixDataType* data; ///< @brief Data pointer, stored as an array

    /// @brief Delete the current data and allocate a matrix of length size2*size
    void allocate()
    {
        if(data!=NULL)
        {
            delete [] data;
        }
        data = new matrixDataType [size*size2];
    }

    /// @brief Constructure. Initialise all to 0/NULL
    seg_Matrix()
    {
        this->size = 0;
        this->size2 = 0;
        this->data = NULL;
    }

public:

    /// @brief Initialise a square matrix of size squareSize
    seg_Matrix(int squareSize)
    {
        this->size2 = 0;
        this->size = 0;
        this->data = NULL;
        if (squareSize <= 2)
        {
            squareSize = 2;
        }
        this->size = squareSize;
        this->size2 = squareSize;
        this->data = NULL;
        allocate();
    }

    /// @brief Initialise a matrix of size _size1*_size2
    seg_Matrix(int _size1, int _size2)
    {
        this->size2 = 0;
        this->size = 0;
        this->data = NULL;
        if (_size1 <= 0)
        {
            cout << "ERROR: Size <= 0" << endl;
            return;
        }
        if (_size2 <= 0)
        {
            cout << "ERROR: Size <= 0" << endl;
            return;
        }
        this->size2 = _size1;
        this->size = _size2;
        this->data = NULL;
        allocate();
    }

     /// @brief Destructor. Delete data if not NULL
    ~seg_Matrix()
    {
        if(data!=NULL)
        {
        delete[] data;
        }
    }
    /// @brief Prints out the diagonal and off-diagonal resuduals of the matrix whn compared to identity. Provides insights about the singularity. Mostly used for debugging.
    void comparetoidentity()
    {
        if(size==size2)
        {
            int worstdiagonal = 0;
            matrixDataType maxunitydeviation = 0.0;
            matrixDataType currentunitydeviation;
            for ( int i = 0; i < size; i++ )
            {
                currentunitydeviation = data[i*size+i] - 1.;
                if ( currentunitydeviation < 0.0) currentunitydeviation *= -1.;
                if ( currentunitydeviation > maxunitydeviation )
                {
                    maxunitydeviation = currentunitydeviation;
                    worstdiagonal = i;
                }
            }
            int worstoffdiagonalrow = 0;
            int worstoffdiagonalcolumn = 0;
            matrixDataType maxzerodeviation = 0.0;
            matrixDataType currentzerodeviation ;
            for ( int i = 0; i < size; i++ )
            {
                for ( int j = 0; j < size; j++ )
                {
                    if ( i == j ) continue;  // we look only at non-diagonal terms
                    currentzerodeviation = data[i*size+j];
                    if ( currentzerodeviation < 0.0) currentzerodeviation *= -1.0;
                    if ( currentzerodeviation > maxzerodeviation )
                    {
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
        else
        {
            cout << "ERROR: Compare to identity - Matrix is not Square"<< endl;
        }
    }

     /// @brief Set the current matrix to the product of left and right
    void settoproduct(seg_Matrix& left, seg_Matrix& right)
    {
        int lsize = left.getSizeColumn();
        int lsize2 = left.getSizeRow();
        int rsize = right.getSizeColumn();
        int rsize2 = right.getSizeRow();

        if (lsize!=rsize2)
        {
            cout << "ERROR: Can't multiply. Wrong sizes/" << endl;
            return;
        }
        size2=lsize2; // numb rows
        size=rsize; // numb cols
        allocate();
        matrixDataType leftvalue=0;
        matrixDataType rightvalue=0;
        for ( int i = 0; i < size2; i++ ) // for each row
            for ( int j = 0; j < size; j++ )    // for each column
            {
                matrixDataType sum = 0.0;
                bool success;
                for (int c = 0; c < lsize; c++)
                {
                    left.getvalue(i,c,leftvalue,success);
                    right.getvalue(c,j,rightvalue,success);
                    sum += leftvalue * rightvalue;
                }
                setvalue(i,j,sum);
            }
    }

     /// @brief Print the matrix
    void dumpmatrix()
    {

        for (int i=0; i < size2; i++) // curr row
        {
            for (int j=0; j<size; j++)  // curr column
            {
                cout << data[i*size+j] << " ";
            }
            cout << endl;
        }
    }

    /// @brief Make "this" a copy of the input source
    void copymatrix(seg_Matrix&  source)
    {
        size = source.getSizeColumn();
        size2 = source.getSizeRow();
        allocate();
        matrixDataType value=0;
        for ( int i = 0; i < size2; i++ )// curr row
            for ( int j = 0; j < size; j++ )   // curr column
            {
                bool success;
                source.getvalue(i,j,value,success); // curr_row, curr_column
                data[i*size+j] = value;
            }
    }

    /// @brief Set the matrix as a diagonal matrix of size diagonalSize. Realocates the data afterwards.
    void setsize(int diagonalSize)
    {
        if (diagonalSize > 0 )
        {
            size = diagonalSize ;
            size2 = diagonalSize;
            allocate();
        }
    }

    /// @brief Set the matrix to the size _size1*_size2. Realocates the data afterwards.
    void setsize(int _size1, int _size2)
    {
        if (_size1 > 0 && _size2 > 0)
        {
            size = _size1 ;
            size2 = _size2;
            allocate();
        }
    }

    /// @brief Get the number of columns
    int getSizeColumn()
    {
        return size;
    }

    /// @brief Get the number of rows
    int getSizeRow()
    {
        return size2;
    }

    /// @brief Get the value at a specific location (row,column) and output it as returnvalue
    void getvalue(int row, int column, matrixDataType& returnvalue, bool& success)
    {
        if ( (row>=size2) || (column>=size)
             || (row<0) || (column<0) )
        {
            success = false;
            return;
        }
        returnvalue = data[ row * size + column ];
        success = true;
    }

    /// @brief Set the value at a specific location (row,column) to newValue
    bool setvalue(int row, int column, matrixDataType newvalue)
    {
        if ( (row >= size2) || (column >= size)
             || (row<0) || (column<0) ) return false;
        data[ row * size + column ] = newvalue;
        return true;
    }

    /// @brief Invert the matrix using a LU decomposition
    void invert()
    {
        if (size!=size2)
        {
            cout << "Matrix in not square" << endl;
            return;
        }
        if (size <= 1) return;
        for (int i=1; i < size; i++)
        {
            data[i] /= data[0]; // normalize row 0
        }
        for (int i=1; i < size; i++)
        {
            for (int j=i; j < size; j++)    // do a column of L
            {
                matrixDataType sum = 0.0;
                for (int k = 0; k < i; k++)
                {
                    sum += data[j*size+k] * data[k*size+i];
                }
                data[j*size+i] -= sum;
            }
            if (i == size-1) continue;
            for (int j=i+1; j < size; j++)     // do a row of U
            {
                matrixDataType sum = 0.0;
                for (int k = 0; k < i; k++)
                {
                    sum += data[i*size+k]*data[k*size+j];
                }
                data[i*size+j] =(data[i*size+j]-sum) / data[i*size+i];
            }
        }
        for ( int i = 0; i < size; i++ )  // invert L
            for ( int j = i; j < size; j++ )
            {
                matrixDataType x = 1.0;
                if ( i != j )
                {
                    x = 0.0;
                    for ( int k = i; k < j; k++ )
                    {
                        x -= data[j*size+k]*data[k*size+i];
                    }
                }
                data[j*size+i] = x / data[j*size+j];
            }
        for ( int i = 0; i < size; i++ )   // invert U
            for ( int j = i; j < size; j++ )
            {
                if ( i == j )
                {
                    continue;
                }
                matrixDataType sum = 0.0;
                for ( int k = i; k < j; k++ )
                {
                    sum += data[k*size+j]*( (i==k) ? 1.0 : data[i*size+k] );
                }
                data[i*size+j] = -sum;
            }
        for ( int i = 0; i < size; i++ )   // final inversion
            for ( int j = 0; j < size; j++ )
            {
                matrixDataType sum = 0.0;
                for ( int k = ((i>j)?i:j); k < size; k++ )
                {
                    sum += ((j==k)?1.0:data[j*size+k])*data[k*size+i];
                }
                data[j*size+i] = sum;
            }
    }

    /// @brief Estimate the matrix determinant using a recursive method
    double determinant()
    {

        int i,j,j1,j2;
        double det = 0;
        double **m = NULL;

        if (size < 1)
        {

        }
        else if (size == 1)     /* Shouldn't get used */
        {
            det = data[0];
        }
        else if (size == 2)
        {
            det = data[0+size*0] * data[1+size*1] - data[1+size*0] * data[0+size*1];
        }
        else
        {
            det = 0;
            for (j1=0; j1<size; j1++)
            {
                m = (double **) malloc((size-1)*sizeof(double *));
                for (i=0; i<size-1; i++)
                    m[i] = (double *) malloc((size-1)*sizeof(double));
                for (i=1; i<size; i++)
                {
                    j2 = 0;
                    for (j=0; j<size; j++)
                    {
                        if (j == j1)
                            continue;
                        m[i-1][j2] = data[i+size*j];
                        j2++;
                    }
                }
                det += pow(-1.0,1.0+j1+1.0) * data[0+size*j1] * Determinant_lib(m,size-1);
                for (i=0; i<size-1; i++)
                    free(m[i]);
                free(m);
            }
        }
        return det;
    }

    /// @brief The actual recusive function used to estimate the matrix determinant. \sa determinant
    double Determinant_lib(double **a,int n)
    {
        int i,j,j1,j2;
        double det = 0;
        double **m = NULL;

        if (n < 1)   /* Error */
        {

        }
        else if (n == 1)     /* Shouldn't get used */
        {
            det = a[0][0];
        }
        else if (n == 2)
        {
            det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
        }
        else
        {
            det = 0;
            for (j1=0; j1<n; j1++)
            {
                m = (double **) malloc((n-1)*sizeof(double *));
                for (i=0; i<n-1; i++)
                    m[i] = (double *) malloc((n-1)*sizeof(double));
                for (i=1; i<n; i++)
                {
                    j2 = 0;
                    for (j=0; j<n; j++)
                    {
                        if (j == j1)
                            continue;
                        m[i-1][j2] = a[i][j];
                        j2++;
                    }
                }
                det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant_lib(m,n-1);
                for (i=0; i<n-1; i++)
                    free(m[i]);
                free(m);
            }
        }
        return(det);
    }
};
