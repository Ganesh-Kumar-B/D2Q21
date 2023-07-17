#ifndef _lbmD2Q21_H_
#define _lbmD2Q21_H_

enum velocityDir{   dV_ZERO_ZERO, dV_P1_ZERO, dV_ZERO_P1, dV_M1_ZERO, dV_ZERO_M1,
                    dV_P2_ZERO, dV_ZERO_P2, dV_M2_ZERO, dV_ZERO_M2,
                    dV_PH1_PH1, dV_MH1_PH1, dV_MH1_MH1, dV_PH1_MH1,
                    dV_P1_P1,   dV_M1_P1,   dV_M1_M1,   dV_P1_M1, 
                    dV_PH3_PH3, dV_MH3_PH3, dV_MH3_MH3, dV_PH3_MH3};  //MH1 and PH1 represent one-half velocities, i.e. (1/2)c

int oppdV[21] =  {  dV_ZERO_ZERO,  dV_M1_ZERO, dV_ZERO_M1, dV_P1_ZERO, dV_ZERO_P1,
                    dV_M2_ZERO, dV_ZERO_M2, dV_P2_ZERO, dV_ZERO_P2,
                     dV_MH1_MH1, dV_PH1_MH1, dV_PH1_PH1, dV_MH1_PH1,
                     dV_M1_M1,  dV_P1_M1,   dV_P1_P1,  dV_M1_P1,
                     dV_MH3_MH3, dV_PH3_MH3, dV_PH3_PH3, dV_MH3_PH3};                                                                                  

template<typename T >
struct lbmD2Q21
{
    int dvN;

     T W[21]; 
     T Cx[21]; 
     T Cy[21]; 
     T Cz[21]; 
     T CxF[21];
     T CxC[21];
     T CyF[21];
     T CyC[21];
    T theta0, thetaInverse ;
    lbmD2Q21(T c); //constructor
    
};



//constructor
template<typename T>
lbmD2Q21<T>::lbmD2Q21(T c1)
    {
        dvN = 21;
        theta0 = 0.39046855785225293;
        thetaInverse = 1.0/theta0;

//ZERO
W[dV_ZERO_ZERO]= 0.2541886612233226;	 Cx[dV_ZERO_ZERO]=0.0;	         Cy[dV_ZERO_ZERO]=0.0;	      

//SC-1
W[dV_P1_ZERO]  =0.0842214730545326;	             Cx[dV_P1_ZERO]  = c1;	         Cy[dV_P1_ZERO]  =0.0;	      
W[dV_ZERO_P1]  =0.0842214730545326;	             Cx[dV_ZERO_P1]  =0.0;	         Cy[dV_ZERO_P1]  = c1;	      
W[dV_M1_ZERO]  =0.0842214730545326;	             Cx[dV_M1_ZERO]  =-c1;	         Cy[dV_M1_ZERO]  =0.0;	      
W[dV_ZERO_M1]  =0.0842214730545326;	             Cx[dV_ZERO_M1]  =0.0;	         Cy[dV_ZERO_M1]  =-c1;	      
            
//SC-2
W[dV_P2_ZERO]  =0.00426526385104285;	         Cx[dV_P2_ZERO]  = 2*c1;	     Cy[dV_P2_ZERO]  = 0.0;	      
W[dV_ZERO_P2]  =0.00426526385104285;	         Cx[dV_ZERO_P2]  =0.0;	         Cy[dV_ZERO_P2]  = 2*c1;	      
W[dV_M2_ZERO]  =0.00426526385104285;	         Cx[dV_M2_ZERO]  =-2*c1;	     Cy[dV_M2_ZERO]  = 0.0;	      
W[dV_ZERO_M2]  =0.00426526385104285;	         Cx[dV_ZERO_M2]  =0.0;	         Cy[dV_ZERO_M2]  = -2*c1;	

//BCC-1/2
W[dV_PH1_PH1] = 0.07057903278838625;           Cx[dV_PH1_PH1]= +c1 * 0.5;       Cy[dV_PH1_PH1]=c1 * 0.5;      
W[dV_MH1_PH1] = 0.07057903278838625;           Cx[dV_MH1_PH1]= -c1 * 0.5;       Cy[dV_MH1_PH1]=c1 * 0.5;
W[dV_MH1_MH1] = 0.07057903278838625;           Cx[dV_MH1_MH1]= -c1 * 0.5;       Cy[dV_MH1_MH1]=-c1 * 0.5;
W[dV_PH1_MH1] = 0.07057903278838625;           Cx[dV_PH1_MH1]= +c1 * 0.5;       Cy[dV_PH1_MH1]=-c1 * 0.5;


//BCC1
W[dV_P1_P1] = 0.02583182337108208;              Cx[dV_P1_P1]= +c1 ;             Cy[dV_P1_P1]=c1  ; 
W[dV_M1_P1] = 0.02583182337108208;              Cx[dV_M1_P1]= -c1 ;             Cy[dV_M1_P1]=c1  ;
W[dV_M1_M1] = 0.02583182337108208;              Cx[dV_M1_M1]= -c1 ;             Cy[dV_M1_M1]=-c1 ;
W[dV_P1_M1] = 0.02583182337108208;              Cx[dV_P1_M1]= +c1 ;             Cy[dV_P1_M1]=-c1 ;





//BCC-3/2
W[dV_PH3_PH3] = 0.001555241629125597;          Cx[dV_PH3_PH3]= +c1 * 1.5;       Cy[dV_PH3_PH3]=c1  * 1.5;      
W[dV_MH3_PH3] = 0.001555241629125597;          Cx[dV_MH3_PH3]= -c1 * 1.5;       Cy[dV_MH3_PH3]=c1  * 1.5;
W[dV_MH3_MH3] = 0.001555241629125597;          Cx[dV_MH3_MH3]= -c1 * 1.5;       Cy[dV_MH3_MH3]=-c1 * 1.5;
W[dV_PH3_MH3] = 0.001555241629125597;          Cx[dV_PH3_MH3]= +c1 * 1.5;       Cy[dV_PH3_MH3]=-c1 * 1.5;




        CxF[dV_ZERO_ZERO]=0.0;	              
        CxF[dV_P1_ZERO]  = c1;	              
        CxF[dV_ZERO_P1]  =0.0;	              
        CxF[dV_M1_ZERO]  =-c1;	              
        CxF[dV_ZERO_M1]  =0.0;	              
        CxF[dV_P2_ZERO]  = 2*c1;		      
        CxF[dV_ZERO_P2]  =0.0;	              
        CxF[dV_M2_ZERO]  =-2*c1;		      
        CxF[dV_ZERO_M2]  =0.0;	              
        CxF[dV_PH1_PH1]  = 0;                   
        CxF[dV_MH1_PH1]  = -1;          
        CxF[dV_MH1_MH1]  = -1;          
        CxF[dV_PH1_MH1]  = 0;
        CxF[dV_P1_P1]    = +c1 ;
        CxF[dV_M1_P1]    = -c1 ;
        CxF[dV_M1_M1]    = -c1 ;
        CxF[dV_P1_M1]    = +c1 ;          
        CxF[dV_PH3_PH3]  = 1;          
        CxF[dV_MH3_PH3]  = -2;          
        CxF[dV_MH3_MH3]  = -2;          
        CxF[dV_PH3_MH3]  = 1;          

        CxC[dV_ZERO_ZERO]=0.0;	    
        CxC[dV_P1_ZERO]  = c1;	    
        CxC[dV_ZERO_P1]  =0.0;	    
        CxC[dV_M1_ZERO]  =-c1;	    
        CxC[dV_ZERO_M1]  =0.0;	    
        CxC[dV_P2_ZERO]  = 2*c1;	    
        CxC[dV_ZERO_P2]  =0.0;	        
        CxC[dV_M2_ZERO]  =-2*c1;	    
        CxC[dV_ZERO_M2]  =0.0;	        
        CxC[dV_PH1_PH1]  = 1;   
        CxC[dV_MH1_PH1]  = 0;   
        CxC[dV_MH1_MH1]  = 0;   
        CxC[dV_PH1_MH1]  = 1;
        CxC[dV_P1_P1]    = +c1 ;
        CxC[dV_M1_P1]    = -c1 ;
        CxC[dV_M1_M1]    = -c1 ;
        CxC[dV_P1_M1]    = +c1 ;   
        CxC[dV_PH3_PH3]  = 2; 
        CxC[dV_MH3_PH3]  = -1; 
        CxC[dV_MH3_MH3]  = -1; 
        CxC[dV_PH3_MH3]  = 2; 




        CyF[dV_ZERO_ZERO]=0.0;	       
        CyF[dV_P1_ZERO]  =0.0;	       
        CyF[dV_ZERO_P1]  = c1;	       
        CyF[dV_M1_ZERO]  =0.0;	       
        CyF[dV_ZERO_M1]  =-c1;	       
        CyF[dV_P2_ZERO]  = 0.0;	       
        CyF[dV_ZERO_P2]  = 2*c1;	       
        CyF[dV_M2_ZERO]  = 0.0;	       
        CyF[dV_ZERO_M2]  = -2*c1;	   
        CyF[dV_PH1_PH1]  =0;      
        CyF[dV_MH1_PH1]  =0;      
        CyF[dV_MH1_MH1]  =-1;     
        CyF[dV_PH1_MH1]  =-1;
        CyF[dV_P1_P1]    =c1  ;
        CyF[dV_M1_P1]    =c1  ;
        CyF[dV_M1_M1]    =-c1 ;
        CyF[dV_P1_M1]    =-c1 ;     
        CyF[dV_PH3_PH3]  =1;     
        CyF[dV_MH3_PH3]  =1;     
        CyF[dV_MH3_MH3]  =-2;     
        CyF[dV_PH3_MH3]  =-2;     
        
        
        CyC[dV_ZERO_ZERO]=0.0;	      
        CyC[dV_P1_ZERO]  =0.0;	      
        CyC[dV_ZERO_P1]  = c1;	      
        CyC[dV_M1_ZERO]  =0.0;	      
        CyC[dV_ZERO_M1]  =-c1;	      
        CyC[dV_P2_ZERO]  = 0.0;	      
        CyC[dV_ZERO_P2]  = 2*c1;	  
        CyC[dV_M2_ZERO]  = 0.0;	      
        CyC[dV_ZERO_M2]  = -2*c1;	
        CyC[dV_PH1_PH1]  =1;      
        CyC[dV_MH1_PH1]  =1;
        CyC[dV_MH1_MH1]  =0;
        CyC[dV_PH1_MH1]  =0;
        CyC[dV_P1_P1]    =c1  ;
        CyC[dV_M1_P1]    =c1  ;
        CyC[dV_M1_M1]    =-c1 ;
        CyC[dV_P1_M1]    =-c1 ;   
        CyC[dV_PH3_PH3]  =2;   
        CyC[dV_MH3_PH3]  =2;
        CyC[dV_MH3_MH3]  =-1;
        CyC[dV_PH3_MH3]  =-1;






    }



#endif

