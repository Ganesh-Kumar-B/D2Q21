#include<iostream>

#define SOLID true
#define FLUID false

template<typename T>
class Grid_N_C_2D{
    public:
int nNodes , nCell,dv,total_mem_node,total_mem_cell, noghost;

int n_x_node,n_y_node,
    n_x,n_y,
    n_x_cell,n_y_cell, //Plus the 2*ghost
    d_v;

T *data_node,*data_cell;

Grid_N_C_2D(int nx, int ny,int no_ghost, int dv){
    
    n_x = nx;
    n_y = ny;

    noghost = no_ghost;
  


    n_x_node = nx + 2*no_ghost;    //2 ghost nodes in each direction
    n_y_node = ny + 2*no_ghost;
   
    d_v = dv;                      //No of population
    

    n_x_cell = nx + 2*no_ghost;   // total no of cells 
    n_y_cell = ny + 2*no_ghost;


    nNodes = (n_x_node)*(n_y_node); 
    nCell  = (n_x_cell)*(n_y_cell);    
   

    total_mem_node = nNodes*dv;
    total_mem_cell = nCell*dv; 
   
    data_node = new T[total_mem_node]{0};
    data_cell = new T[total_mem_cell]{0};
}


~Grid_N_C_2D(){
    delete [] data_node;
    delete [] data_cell;
}
T Node(const int i,const int j,const int dv) const{return data_node[(i*n_y_node +j)*d_v + dv ];} //start 0 from ghost node
T& Node(const int i, const int j,const  int dv) {return data_node[(i*n_y_node +j)*d_v + dv ];}   // start 0 from ghost node

T Cell(const int i,const int j,const int dv) const{return data_cell[(i*n_y_cell +j)*d_v + dv ];} //start 0 from ghost cell
T& Cell(const int i, const int j,const  int dv) {return data_cell[(i*n_y_cell +j)*d_v + dv ];}   //start 0 from ghost cell




};