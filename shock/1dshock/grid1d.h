#include<mm_malloc.h>
#include<iostream>


template <int numfield, typename dataType> 
class grid1D 
 {
  public:
   grid1D(int size1=1, int stride1=0)
  {
   m1 = size1;
   actualGridSize = m1;
   nB1 = stride1;
   n1 = size1+2*stride1;
   nE1 = size1+ stride1-1;
   sizeGrid = n1;
   numField = numfield;


   numElem= sizeGrid*numField ;

   data = (dataType*) _mm_malloc(numElem*sizeof(dataType), 4096);
   
   numField_offset = sizeGrid*numField;
   }
   

  ~grid1D()
  { 
   delete data; 
  }
  
  int size() const
  {
   return numElem;
  }
     
  int getIndex (int i1, int dv)
  const{ return  ((i1)*numField + dv);//CHECK
   }
  
  dataType  Node(const int i1,  int dv) const  { return data[ getIndex(i1,dv)];} 
  dataType& Node(const int i1,  int dv)        { return data[ getIndex(i1,dv)];} 
   
  

  void initializeIntegers()
  {
   for (int i=nB1; i<=nE1; i++)
   for (int dv=0; dv<numfield; dv++)
   {
    Node(i, dv) = getIndex(i,dv);
    // Cell(i,j,k, dv) = getIndex(i,j,k,dv);
   }
  }
   
  template <int N, int bL, typename dataType1>
  friend std::ostream& operator<< (std::ostream& os,const grid1D<N, dataType1> &thisVec);
 
   
     
   
  public:
  int m1;                // Size of physical  grid
  int nB1;             // Number of ghost nodes at ends in  direction i and the beginning of the physical grid
  int nE1;             // Last element of physical grid nEi = mi + nBi-1
  int n1;                // Size of computational  grid in  direction i, ni = mi+ 2*nBi
  int numElem;                   // total number of element in data array sizeGrid*numField
  int sizeGrid;                  // total number of computational grid points n1*n2*n3
  int numField,actualGridSize; 
  int numField_offset; 
  dataType *data; 
 };
 
	
