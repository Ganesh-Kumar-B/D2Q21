#include <cmath>

#ifndef _ANSU_D1Q3_
#define _ANSU_D1Q3_
#define NUM_DV 3

enum dirTag{
	DX1,
	};
enum dvTag{
    dv_ZERO,
    dv_P1, //Positive Side
	  dv_M1, //Negative Side

};

template<typename dataType>
class lbModelD1Q3
{
  public:
    lbModelD1Q3(dataType spacing){
    c = spacing;
    for(int dv=0;dv<NUM_DV;dv++){
         this->linkVel(dv, DX1)	=0.0;
    }
	
	theta0 = 1.0/3;
  theta0inv = 1.0/theta0;
	
  weight[dv_ZERO] =  4.0/6.0;
  weight[dv_M1]= weight[dv_P1]= 1.0/6.0;
	
	 this->linkVel(dv_P1, DX1) =    1.0*c;
	 this->linkVel(dv_M1, DX1) =   -1.0*c;
}

  dataType& linkVel(int dv, int dir){
    return linkData[dv*1 + dir];
  }
  const dataType linkVel(int dv, int dir) const{
    return linkData[dv*1  + dir];
  }

  dataType c;
  dataType theta0;
  dataType theta0inv;
	dataType weight[NUM_DV];
  dataType linkData[NUM_DV]; 
  dataType fEq[NUM_DV];
	dataType dot[NUM_DV];
} ;


 //Base Implementation
template<typename T1>
void getEqISO_O2(lbModelD1Q3<T1> &D1Q3, T1 rho, T1 ux,T1 c)
{

T1 ubar=ux/c;//u bar is always less than 1

T1 b=(2*ubar + sqrt(1+3*ubar*ubar))/(1-ubar); //Taking b which will always be Positive for ubar<1
// T1 b=(2*ubar - sqrt(1+3*ubar*ubar))/(1-ubar);
// T1 a=6*rho*ubar*(b/(b*b-1));
T1 a=6*rho/(4+b+(1/b));



D1Q3.fEq[dv_ZERO] = D1Q3.weight[dv_ZERO]*a;
D1Q3.fEq[dv_P1]   = D1Q3.weight[dv_P1]  *(a*b);
D1Q3.fEq[dv_M1]   = D1Q3.weight[dv_M1]  *(a/b);
}

template <typename T1>
void computeMoments(lbModelD1Q3<T1> &D1Q3, T1 &rho, T1 &u1)
{
	rho=0.0;
	u1=0.0;

		rho= D1Q3.fEq[dv_ZERO] + D1Q3.fEq[dv_P1] + D1Q3.fEq[dv_M1];
		u1 = D1Q3.fEq[dv_ZERO]*D1Q3.linkVel(dv_ZERO,DX1) + D1Q3.fEq[dv_P1]*D1Q3.linkVel(dv_P1,DX1) + 
				 D1Q3.fEq[dv_M1]*D1Q3.linkVel(dv_M1,DX1);
 u1=u1/rho;

// std::cout << "rho "<< rho << std::endl;
// std::cout << "ux  "<< u1  << std::endl;

}












#endif
