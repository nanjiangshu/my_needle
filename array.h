/*
 * =====================================================================================
 * 
 *       Filename:  array.h
 *    Description:  array implementation
 *        Version:  1.0
 *        Created:  05/02/2006 05:12:34 PM CEST
 *        Modified: 2008-02-13 15:23:55 Wednesday Week 06 
 *       Revision:  none
 *       Compiler:  gcc
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */

#ifndef HAS_ARRAY_H
#define HAS_ARRAY_H

#include <iostream>
//Array1D
template <class T>  class Array1D // can be accessed both by .array1D and []
{
    protected:	
        unsigned int size;
        //	T* array1D;
    public:	
        T* array1D; // can also accessed by array1D directory, quicker
        Array1D(unsigned int size); 
        ~Array1D(void); 
        T& operator [] (unsigned int x); 
        unsigned int GetSize(void);
        void Init(T initValue);
};

// Array1D
template <class T> Array1D<T>::Array1D(unsigned int size)/*{{{*/
{
	array1D = new T[size];
	this->size = size ;
}/*}}}*/
template <class T> Array1D<T>::~Array1D(void) /*{{{*/
{
    delete [] array1D ;
    array1D = NULL;
    this->size = 0;
}/*}}}*/
template <class T> void Array1D<T>::Init(T initValue) /*{{{*/
{
    unsigned int i;
    for(i = 0; i < size; i ++)
    {
        array1D[i] = initValue;
    }
}/*}}}*/
template <class T> unsigned int Array1D<T>::GetSize(void)  { return  size;} 
template <class T> T& Array1D<T>::operator [] (unsigned int i) { return array1D[i]; }


//array 2D
// old-style, Create2DArray and Delete2DArray pair
template <class T> T **Create2DArray(T **array, unsigned int xSize, unsigned int ySize);
template <class T> void Delete2DArray(T **array, unsigned int xSize);

/*****************************************************************************
 * Array1D, Array2D, Array3D, new style, will be deconstructed automatically
 * when the procedure is over
 ****************************************************************************/
template <class T> class Array2D //accessed by .array2D
{
    public:
        T** array2D;
        Array2D(unsigned int rowSize,unsigned int colSize);
        ~Array2D(void);
        unsigned int GetRowSize(void);
        unsigned int GetColSize(void);
        void Init(T initValue = T(0));
    private:	
        unsigned int rowSize;
        unsigned int colSize;
};

//**********************************************************
//Create 2-dimensional array, which can be accessed by array[row][col]
//**********************************************************
template <class T> T **Create2DArray(T **array, unsigned int xSize, unsigned int ySize)/*{{{*/
{
    unsigned int i ;
    array = new T*[xSize] ;
    for ( i = 0 ; i < xSize ; i ++ )
        array[i] = new T[ySize] ;
    return array ;
}/*}}}*/
template <class T> void Delete2DArray(T **array, unsigned int xSize)/*{{{*/
{
    unsigned int i ; 
    for ( i = 0 ; i < xSize ; i ++ )
        delete [] array[i] ;
    delete [] array ;
}/*}}}*/

// Array2D Class, 
template <class T> Array2D<T>::Array2D(unsigned int rowSize, unsigned int colSize)
{
	this->rowSize = rowSize;
	this->colSize = colSize;
	unsigned int i ;
	array2D = new T*[rowSize] ;
	for ( i = 0 ; i < rowSize ; i ++ )
		array2D[i] = new T[colSize] ;
}
// destroy the 2D array
template <class T> Array2D<T>::~Array2D(void)/*{{{*/
{
	unsigned int i ; 
	for ( i = 0 ; i < rowSize ; i ++ )
		delete [] array2D[i] ;
	delete [] array2D ;
    array2D = NULL;
    this->rowSize = 0;
    this->colSize = 0;
}/*}}}*/
template <class T>  void Array2D<T>::Init(T initValue) /*{{{*/
{ 
    unsigned int i,j;
    for (i = 0 ; i < rowSize; i ++)
    {
        for(j = 0; j < colSize; j ++)
        {
            array2D[i][j] = initValue;
        }
    }
} /*}}}*/
template <class T>  unsigned int Array2D<T>::GetRowSize(void) { return  rowSize;} 
template <class T>  unsigned int Array2D<T>::GetColSize(void) { return  colSize;} 


//Array3D
// old-style, Create3DArray and Delete3DArray pair
template <class T> T*** Create3DArray(T ***array, unsigned int xSize, unsigned int ySize, unsigned int zSize);
template <class T> void Delete3DArray(T ***array, unsigned int xSize, unsigned int ySize);

template <class T> class Array3D //accessed by .array3D
{
    public:
        T*** array3D;
        Array3D(unsigned int xSize,unsigned int ySize, unsigned int zSize);
        ~Array3D(void);
        unsigned int GetXSize(void);
        unsigned int GetYSize(void);
        unsigned int GetZSize(void);
        void Init(T initValue = T(0));
    private:	
        unsigned int xSize;
        unsigned int ySize;
        unsigned int zSize;
};

//**********************************************************
//Create 3-dimensional array, which can be accessed by array[x][y][z]
//**********************************************************
template <class T> T*** Create3DArray(T ***array, unsigned int xSize, unsigned int ySize, unsigned int zSize)/*{{{*/
{
    unsigned int i ;
    array = new T**[xSize] ;
    for ( i = 0 ; i < xSize ; i ++ )
        array[i] = Create2DArray(array[i], ySize, zSize);
    return array ;
}/*}}}*/

/*****************************************************************************
 * Delete3DArray()
 * free the memory allocated to 3d array
 ****************************************************************************/
template <class T> void Delete3DArray(T ***array, unsigned int xSize, unsigned int ySize)/*{{{*/
{
    unsigned int i ; 
    for ( i = 0 ; i < xSize ; i ++ )
        Delete2DArray(array[i], ySize)  ;
    delete [] array ;
}/*}}}*/

// Array3D
template <class T> Array3D<T>::Array3D(unsigned int xSize, unsigned int ySize, unsigned int zSize)/*{{{*/
{
	this->xSize = xSize;
	this->ySize = ySize;
	this->zSize = zSize;
	unsigned int i,j ;
	array3D = new T**[xSize] ;
	for ( i = 0 ; i < xSize ; i ++ )
    {
        array3D[i] = new T*[ySize] ;
        for( j = 0 ; j < ySize ; j ++ )
        { array3D[i][j] = new T[zSize]; }
    }
}/*}}}*/
template <class T> Array3D<T>::~Array3D(void)/*{{{*/
{
	unsigned int i,j ; 
	for ( i = 0 ; i < xSize ; i ++ )
    {
        for( j = 0 ; j < ySize ; j ++ )
        { delete [] array3D[i][j] ; }
        delete [] array3D[i];
    }
	delete [] array3D ;
    array3D = NULL;
    this->xSize = 0;
    this->ySize = 0;
    this->zSize = 0;
}/*}}}*/
template <class T>  void Array3D<T>::Init(T initValue) /*{{{*/
{ 
	unsigned int i,j,k ; 
	for ( i = 0 ; i < xSize ; i ++ )
    {
        for( j = 0 ; j < ySize ; j ++ )
        { 
            for(k = 0 ; k < zSize; k++)
            { array3D[i][j][k] = initValue; }
        }
    }
} /*}}}*/
template <class T>  unsigned int Array3D<T>::GetXSize(void) { return  xSize;} 
template <class T>  unsigned int Array3D<T>::GetYSize(void) { return  ySize;} 
template <class T>  unsigned int Array3D<T>::GetZSize(void) { return  zSize;} 


//array2d_sub
template <class T> class Array2D_Sub // can be accessed by [][]
                                     // using [] operator, but much slower
{
protected:
	unsigned int rowSize;
	unsigned int colSize;
	T* array;
public:
	class Row
	{
		Array2D_Sub& array2D;
		unsigned int const row;
	public:
		Row(Array2D_Sub& _array2D,unsigned int _row):
		  array2D (_array2D),row(_row) {}
		T& operator [] (unsigned int col) const
		{ return array2D.Select(row,col);}
	};
	Array2D_Sub(unsigned int rowSize,unsigned int colSize)
	{
		this->rowSize = rowSize;
		this->colSize = colSize;
		array = new T[rowSize*colSize];
	}
	~Array2D_Sub()
	{
		delete [] array;
		this->rowSize = 0 ;
		this->colSize = 0 ;
	}
	T& Select(unsigned int row, unsigned int col)
	{
		return array[row*colSize+col];
	}
	Row operator [] (unsigned int row) {return Row(*this,row);}
	unsigned int GetRowSize(void)	{return this->rowSize; }
	unsigned int GetColSize(void)    {return this->colSize; }
};

#endif /*HAS_ARRAY_H*/
