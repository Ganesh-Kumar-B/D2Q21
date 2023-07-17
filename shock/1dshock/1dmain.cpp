#include <fstream>
#include <iostream>
using namespace std;
#include "D1Q3.h"
#include "grid3D.h"


/*Not using this yet*/
template <int numfield, typename dataType>
void copyPopulationFROMgrid(lbModelD1Q3<dataType> &lbModel,
         gridBCC3D<numfield,dataType> &gridLB, int i, int j , int k) {

  for (int dv = 0 ; dv < NUM_DV; dv++)
    lbModel.fEq[dv] = gridLB.Node(i,j,k,dv);  
}


template <int numfield, typename dataType>
void computeMomentsFromGrid(lbModelD1Q3<dataType> &lbModel, gridBCC3D< numfield,dataType> &gridLB, 
  int i, int j , int k, double &rho, double &u1) {
  rho = 0;
  u1 = 0;

  for(int dv=0;dv<NUM_DV;dv++)
  {
  rho += gridLB.Node(i,j,k,dv);
  u1  += gridLB.Node(i,j,k,dv) * lbModel.linkVel(dv,DX1);
  }
      u1 = u1/rho;

// std::cout << "rho from grid "<< rho << std::endl;
// std::cout << "ux  from grid "<< u1  << std::endl;


}

template <int numfield, typename dataType>
void collide(lbModelD1Q3<dataType> &lbModel,
             gridBCC3D< numfield,dataType> &gridLB,dataType beta) {
  
  dataType twoBeta = 2.0*beta ;
  dataType rho,u1;

  dataType c=1;//REMEMBER used c here

/*In initialisation, our grid gets feqs using rhoINITIAL and uINITIAL values  for all nodes*/
  for(int k=gridLB.nB3; k<=gridLB.nE3; k++)
   for(int j=gridLB.nB2; j<=gridLB.nE2; j++)
    for(int i=gridLB.nB1; i<=gridLB.nE1; i++)
    { 

   computeMomentsFromGrid(lbModel,gridLB,i,j,k,rho,u1);
   getEqISO_O2(lbModel,rho,u1,c);
     for(int dv=0; dv<numfield; dv++)
      gridLB.Node(i,j,k,dv) += twoBeta*(lbModel.fEq[dv] - gridLB.Node(i,j,k,dv));
    }
}



template <int numfield, typename datatype>
void advection(gridBCC3D <numfield,datatype> &gridLB)
{
  //forward Sweep
  for (int i=1;i<=gridLB.m1;i++)
    {
     gridLB.Node(i,0,0,dv_M1)=gridLB.Node(i+1,0,0,dv_M1);
    }
 //Backward Sweep
  for (int i=gridLB.m1;i>=1;i--)
  {
  	gridLB.Node(i,0,0,dv_P1)=gridLB.Node(i-1,0,0,dv_P1);
  }
}



template <int numfield, typename dataType>
void initialConditions(lbModelD1Q3<dataType> &lbModel, gridBCC3D< numfield,dataType> &gridLB,
  dataType rhominus,dataType rhoplus, dataType U0)
 {



               
    for(int i=gridLB.nB1; i<=gridLB.nE1/2; i++) {
   dataType c=1;//REMEMBER used c here
   getEqISO_O2(lbModel,rhominus,U0,c);
   // computeMoments(lbModel,rhominus,u1); 
     for(int dv=0;dv<3;dv++)
      gridLB.Node(i,0,0,dv) = lbModel.fEq[dv];
  }

    for(int i=(gridLB.nE1/2)+1; i<=gridLB.nE1; i++) {
   dataType c=1;//REMEMBER used c here
   getEqISO_O2(lbModel,rhoplus,U0,c);
   // computeMoments(lbModel,rhoplus,u1); 
     for(int dv=0;dv<3;dv++)
      gridLB.Node(i,0,0,dv) = lbModel.fEq[dv];
  }

}




template<int numfield, typename datatype>
void prepareBoundary(gridBCC3D<numfield,datatype> &gridLB)
{
  for(int dv=0; dv<numfield; dv++)
  {
  gridLB.Node(0,0,0,dv)=gridLB.Node(1,0,0,dv);
  gridLB.Node(gridLB.m1+1,0,0,dv)=gridLB.Node(gridLB.m1,0,0,dv);
  }
}



template <int numfield, typename dataType>
void printpopulation(lbModelD1Q3<dataType> &lbModel, gridBCC3D< numfield,dataType> &gridLB){
for (int i=0;i<=gridLB.m1+1;i++){
{  for(int dv=0;dv<3;dv++){
cout <<"gridLB.Node("<< i<<",0,0,"<< dv<<")    "<< gridLB.Node(i,0,0,dv)<< endl;
}
// << endl;}
}
}
}





template <int numfield, typename dataType>
void writer(lbModelD1Q3<dataType> &lbModel, gridBCC3D< numfield,dataType> &gridLB){


// if(t%20==0){
// double rho,ux;
// ofstream write("data.dat");

// for (int i=1;i<=gridLB.m1;i++){


// computeMomentsFromGrid(lbModel, gridLB,  i, 0 , 0, rho, ux);

// write  << i<< " "<< rho << " "<< ux << endl;
// }
// write.close();

// }

  double rho,ux;
ofstream write("datavelocity400.dat");

for (int i=1;i<=gridLB.m1;i++){


computeMomentsFromGrid(lbModel, gridLB,  i, 0 , 0, rho, ux);

write  << i<< " "<< rho << " "<< ux << endl;
}
write.close();
}


int main(){

  double t=0;
  double tf=400;



lbModelD1Q3 <double> lbModel(1.0);





double rhominus=1.5;
double rhoplus=0.75;
double ux=0;

double c=1;

double beta=1;

// getEqISO_O2(lbModel,rho,ux,c);

// computeMoments(lbModel,rho,ux);


gridBCC3D<3,double> lbgrid(800,1,1,1,0,0);

initialConditions(lbModel,lbgrid,rhominus,rhoplus,ux);

while(t<=tf){
collide(lbModel,lbgrid,beta);

prepareBoundary(lbgrid);

advection(lbgrid);
t=t+1;
}

writer(lbModel,lbgrid);

// printpopulation(lbModel,lbgrid);









return 0;




}