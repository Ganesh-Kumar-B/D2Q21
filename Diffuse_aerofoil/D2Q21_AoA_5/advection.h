#ifndef  ADVECTION_H
#define  ADVECTION_H

#include "lbmD2Q21.h"
// enum velocityDir{   dV_ZERO_ZERO, dV_P1_ZERO, dV_ZERO_P1, dV_M1_ZERO, dV_ZERO_M1,  //SC1
//                     dV_P2_ZERO, dV_ZERO_P2, dV_M2_ZERO, dV_ZERO_M2,                //SC2
//                     dV_PH1_PH1, dV_MH1_PH1, dV_MH1_MH1, dV_PH1_MH1,                //SC1/2
//                     dV_PH3_PH3, dV_MH3_PH3, dV_MH3_MH3, dV_PH3_MH3};               //SC3/2





// template<typename T>
// void advection(Grid_N_C_2D<T> &grid){
    
//         double a,b;
//         for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//             a = grid.Cell(i,grid.noghost -1,dV_MH1_PH1);
//             // std::cout<<a<<" "<<i<<std::endl;
//         for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost); j++){ 
            
//             grid.Node(i,j,dV_ZERO_M1) = grid.Node(i  ,j+1, dV_ZERO_M1 ); //for SC1 -Node
//             grid.Node(i,j,dV_M1_ZERO) = grid.Node(i+1,j  , dV_M1_ZERO );

//             grid.Cell(i,j,dV_ZERO_M1) = grid.Cell(i  ,j+1, dV_ZERO_M1 ); //for SC1 -cell
//             grid.Cell(i,j,dV_M1_ZERO) = grid.Cell(i+1,j  , dV_M1_ZERO );

//             ///half velocities
//             //std::cout<<grid.Node(i,j,dV_MH1_MH1)<<"before "<<i<<j <<std::endl;
            
//             grid.Node(i,j,dV_MH1_MH1) = grid.Cell(i  ,j  , dV_MH1_MH1);
//             grid.Node(i,j,dV_MH1_PH1) = a;     
//             a = grid.Cell(i,j,dV_MH1_PH1);       

//             grid.Cell(i,j,dV_MH1_MH1) = grid.Node(i+1,j+1, dV_MH1_MH1);                        
//             grid.Cell(i,j,dV_MH1_PH1) = grid.Node(i+1,j  , dV_MH1_PH1);

//             //3 velocities

//             grid.Node(i,j,dV_ZERO_M3) = grid.Node(i  ,j+3, dV_ZERO_M3);
//             grid.Node(i,j,dV_M3_ZERO) = grid.Node(i+3,j  , dV_M3_ZERO);

//             grid.Cell(i,j,dV_ZERO_M3) = grid.Cell(i  ,j+3, dV_ZERO_M3);
//             grid.Cell(i,j,dV_M3_ZERO) = grid.Cell(i+3,j  , dV_M3_ZERO);
//             //std::cout<<grid.Node(i,j,dV_MH1_MH1) <<" "<<i<<j<<std::endl;
//         } 
//     }

//         for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){
//             b = grid.Node(i,grid.n_y_node - grid.noghost,dV_PH1_MH1);
//         for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){

//             //SC1
//             grid.Node(i,j, dV_ZERO_P1) = grid.Node(i  ,j-1, dV_ZERO_P1);
//             grid.Node(i,j, dV_P1_ZERO) = grid.Node(i-1,j  , dV_P1_ZERO);
//             grid.Cell(i,j,dV_ZERO_P1) = grid.Cell(i  ,j-1, dV_ZERO_P1 );
//             grid.Cell(i,j,dV_P1_ZERO) = grid.Cell(i-1,j  , dV_P1_ZERO );
      
//             ///half velocities
//            // std::cout<<grid.Node(i,j,dV_PH1_PH1)<<"before "<<i<<j <<std::endl;
//             grid.Cell(i,j,dV_PH1_MH1) = b;
//             grid.Cell(i,j,dV_PH1_PH1) = grid.Node(i  ,j  , dV_PH1_PH1);
//             b=grid.Node(i,j,dV_PH1_MH1);
//             grid.Node(i,j,dV_PH1_MH1) = grid.Cell(i-1,j  , dV_PH1_MH1);       
//             grid.Node(i,j,dV_PH1_PH1) = grid.Cell(i-1,j-1, dV_PH1_PH1);



//             //std::cout<<grid.Cell(i,j,dV_PH1_PH1) <<" "<<i<<j<<std::endl;
//             //SC3
//             grid.Node(i,j, dV_ZERO_P3) = grid.Node(i  ,j-3, dV_ZERO_P3);
//             grid.Node(i,j, dV_P3_ZERO) = grid.Node(i-3,j  , dV_P3_ZERO);


//             grid.Cell(i,j,dV_ZERO_P3) = grid.Cell(i  ,j-3, dV_ZERO_P3 );
//             grid.Cell(i,j,dV_P3_ZERO) = grid.Cell(i-3,j  , dV_P3_ZERO );        
//         }}

// }

// template<typename T>
// void advection1(Grid_N_C_2D<T> &grid){


//         #pragma omp parallel 
//         {


//         #pragma omp for schedule(static)
//         for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//         for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 

//             grid.Node(i,j,dV_ZERO_M1) = grid.Node(i  ,j+1, dV_ZERO_M1 ); //for SC1 -Node
//             grid.Node(i,j,dV_M1_ZERO) = grid.Node(i+1,j  , dV_M1_ZERO );

//             grid.Node(i,j,dV_ZERO_M3) = grid.Node(i  ,j+3, dV_ZERO_M3);
//             grid.Node(i,j,dV_M3_ZERO) = grid.Node(i+3,j  , dV_M3_ZERO);

//             grid.Node(i,j  ,dV_MH1_MH1) = grid.Cell(i  ,j   , dV_MH1_MH1);
//             grid.Node(i,j  ,dV_MH1_PH1) = grid.Cell(i  ,j -1, dV_MH1_PH1);

//         }
//         }


//         #pragma omp for schedule(static)
//         for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
//         for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 

//             grid.Cell(i,j,dV_ZERO_M1) = grid.Cell(i  ,j+1, dV_ZERO_M1 ); //for SC1 -cell
//             grid.Cell(i,j,dV_M1_ZERO) = grid.Cell(i+1,j  , dV_M1_ZERO );

//             grid.Cell(i,j,dV_ZERO_M3) = grid.Cell(i  ,j+3, dV_ZERO_M3);
//             grid.Cell(i,j,dV_M3_ZERO) = grid.Cell(i+3,j  , dV_M3_ZERO);

//             grid.Cell(i,j,dV_MH1_MH1) = grid.Node(i+1,j+1, dV_MH1_MH1);                        
//             grid.Cell(i,j,dV_MH1_PH1) = grid.Node(i+1,j  , dV_MH1_PH1);
            

//         }
//         }

//         #pragma omp for schedule(static)
//         for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
//         for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            
//             grid.Cell(i,j,dV_ZERO_P1) = grid.Cell(i  ,j-1, dV_ZERO_P1 );
//             grid.Cell(i,j,dV_P1_ZERO) = grid.Cell(i-1,j  , dV_P1_ZERO );

//             grid.Cell(i,j,dV_ZERO_P3) = grid.Cell(i  ,j-3, dV_ZERO_P3 );
//             grid.Cell(i,j,dV_P3_ZERO) = grid.Cell(i-3,j  , dV_P3_ZERO );

//             grid.Cell(i,j,dV_PH1_MH1) = grid.Node(i,j+1,dV_PH1_MH1);
//             grid.Cell(i,j,dV_PH1_PH1) = grid.Node(i,j  ,dV_PH1_PH1);

//         }}


//         #pragma omp for schedule(static)
//         for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
//         for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            
//             grid.Node(i,j, dV_ZERO_P1) = grid.Node(i  ,j-1, dV_ZERO_P1);
//             grid.Node(i,j, dV_P1_ZERO) = grid.Node(i-1,j  , dV_P1_ZERO);

//             grid.Node(i,j, dV_ZERO_P3) = grid.Node(i  ,j-3, dV_ZERO_P3);
//             grid.Node(i,j, dV_P3_ZERO) = grid.Node(i-3,j  , dV_P3_ZERO);

//             grid.Node(i,j,dV_PH1_MH1) = grid.Cell(i-1,j  ,dV_PH1_MH1);
//             grid.Node(i,j,dV_PH1_PH1) = grid.Cell(i-1,j-1,dV_PH1_PH1);

//         }}


//         }

// }


template<typename T>
void advection_21(Grid_N_C_2D<T> &grid){


        for(int i = 0 + grid.noghost ; i < grid.n_x_node - (grid.noghost);i++){
        for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 

            grid.Node(i,j,dV_ZERO_M1) = grid.Node(i  ,j+1, dV_ZERO_M1 ); //for SC1 -Node
            grid.Node(i,j,dV_M1_ZERO) = grid.Node(i+1,j  , dV_M1_ZERO );



            grid.Node(i,j,dV_ZERO_M2) = grid.Node(i  ,j+2, dV_ZERO_M2);   //for SC2 Node
            grid.Node(i,j,dV_M2_ZERO) = grid.Node(i+2,j  , dV_M2_ZERO);   

            grid.Node(i,j  ,dV_MH1_MH1) = grid.Cell(i  ,j   , dV_MH1_MH1);   ///for SC1/2
            grid.Node(i,j  ,dV_MH1_PH1) = grid.Cell(i  ,j -1, dV_MH1_PH1);

            grid.Node(i,j  ,dV_MH3_MH3) = grid.Cell(i+1,j+1, dV_MH3_MH3);    //for SC3/2
            grid.Node(i,j  ,dV_MH3_PH3) = grid.Cell(i+1,j-2, dV_MH3_PH3);
          
            grid.Node(i,j, dV_M1_M1) = grid.Node(i+1,j+1, dV_M1_M1);
            grid.Node(i,j, dV_M1_P1) = grid.Node(i+1,j-1, dV_M1_P1);

        }
        
        for(int j = 0 + grid.noghost ; j < grid.n_y_node - (grid.noghost);j++){ 
          
            grid.Cell(i,j,dV_ZERO_M1) = grid.Cell(i  ,j+1, dV_ZERO_M1 );
            grid.Cell(i,j,dV_M1_ZERO) = grid.Cell(i+1,j  , dV_M1_ZERO );            
            
            grid.Cell(i,j,dV_ZERO_M2) = grid.Cell(i  ,j+2, dV_ZERO_M2);
            grid.Cell(i,j,dV_M2_ZERO) = grid.Cell(i+2,j  , dV_M2_ZERO);

            grid.Cell(i,j,dV_MH1_MH1) = grid.Node(i+1,j+1, dV_MH1_MH1);  
            grid.Cell(i,j,dV_MH1_PH1) = grid.Node(i+1,j  , dV_MH1_PH1);

            grid.Cell(i,j,dV_MH3_MH3) = grid.Node(i+2,j+2, dV_MH3_MH3);  
            grid.Cell(i,j,dV_MH3_PH3) = grid.Node(i+2,j-1, dV_MH3_PH3);

            grid.Cell(i,j, dV_M1_M1) = grid.Cell(i+1,j+1, dV_M1_M1);
            grid.Cell(i,j, dV_M1_P1) = grid.Cell(i+1,j-1, dV_M1_P1);            

            }
        }

        for(int i = grid.n_x_node - (grid.noghost +1); i> grid.noghost-1; i--){                                                                                        //  #pragma omp parallel for shared(b)
        for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            
            grid.Cell(i,j,dV_ZERO_P1) = grid.Cell(i  ,j-1, dV_ZERO_P1 );
            grid.Cell(i,j,dV_P1_ZERO) = grid.Cell(i-1,j  , dV_P1_ZERO );

            grid.Cell(i,j,dV_ZERO_P2) = grid.Cell(i  ,j-2, dV_ZERO_P2 );
            grid.Cell(i,j,dV_P2_ZERO) = grid.Cell(i-2,j  , dV_P2_ZERO );

            grid.Cell(i,j,dV_PH1_MH1) = grid.Node(i,j+1,dV_PH1_MH1);
            grid.Cell(i,j,dV_PH1_PH1) = grid.Node(i,j  ,dV_PH1_PH1);

            grid.Cell(i,j,dV_PH3_MH3) = grid.Node(i-1,j+2,dV_PH3_MH3);
            grid.Cell(i,j,dV_PH3_PH3) = grid.Node(i-1,j-1,dV_PH3_PH3);

            grid.Cell(i,j,dV_P1_M1) = grid.Cell(i-1, j+1, dV_P1_M1);
            grid.Cell(i,j,dV_P1_P1) = grid.Cell(i-1, j-1, dV_P1_P1);

        }
        
        for(int j = grid.n_y_node - (grid.noghost +1); j> grid.noghost-1; j--){
            
            grid.Node(i,j, dV_ZERO_P1) = grid.Node(i  ,j-1, dV_ZERO_P1);
            grid.Node(i,j, dV_P1_ZERO) = grid.Node(i-1,j  , dV_P1_ZERO);

            grid.Node(i,j, dV_ZERO_P2) = grid.Node(i  ,j-2, dV_ZERO_P2);
            grid.Node(i,j, dV_P2_ZERO) = grid.Node(i-2,j  , dV_P2_ZERO);

            grid.Node(i,j,dV_PH1_MH1) = grid.Cell(i-1,j  ,dV_PH1_MH1);
            grid.Node(i,j,dV_PH1_PH1) = grid.Cell(i-1,j-1,dV_PH1_PH1);

            grid.Node(i,j,dV_PH3_MH3) = grid.Cell(i-2,j+1,dV_PH3_MH3);
            grid.Node(i,j,dV_PH3_PH3) = grid.Cell(i-2,j-2,dV_PH3_PH3);

            grid.Node(i,j,dV_P1_M1) = grid.Node(i-1, j+1, dV_P1_M1);
            grid.Node(i,j,dV_P1_P1) = grid.Node(i-1, j-1, dV_P1_P1);
        }   
        }
}



template<typename T,typename T1>
void G_calc(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lb, Grid_N_C_2D<T> &G_calc, Grid_N_C_2D<T1> &marker){

double feq[21] = {0};
get_equi(feq,lb,0.0,0.0,1.0,lb.theta0);

     for(int i = 0 + grid.noghost; i < grid.n_x_node - (grid.noghost) ; i++){
        for(int j = /*1 +*/ grid.noghost;j < grid.n_y_node - (grid.noghost) ; j++){

            
            
            // For SC
            for(int dv = 1; dv<5; dv++){

                if(marker.Node(i,j,0) == FLUID){
                if(marker.Node(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID){

                    G_calc.Node(i,j,0)                                      += feq[oppdV[dv]] *1 + feq[oppdV[dv + 4]]* 2;    //for the closest node all the sc velocites will hit the wall
                    G_calc.Node(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],0)    += feq[oppdV[dv+4]]*2;                             // for the next closest the sc2 velocities will hit the walll
                
                }
                }



                if(marker.Cell(i,j,0) == FLUID)
                if(marker.Cell(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID){
                 
                    G_calc.Cell(i,j,0)                                      += feq[oppdV[dv  ]]*1 + feq[oppdV[dv + 4]] *2;
                    G_calc.Cell(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],0)    += feq[oppdV[dv+4]]*2;

                }
            }
            


            // for BCC
            for(int dv = 9; dv<13; dv ++){
                if(marker.Node(i,j,0) == FLUID)
                if(marker.Cell(i + (int)lb.CxF[dv], j + (int)lb.CyF[dv],0) == SOLID){
                    
                    G_calc.Node(i,j,0)                                                    += feq[oppdV[dv   ]]*0.5 + feq[oppdV[dv +4]]*1.0 + feq[oppdV[dv +8]]*1.5;
                    G_calc.Cell(i + (int)lb.CxF[oppdV[dv]], j + (int)lb.CyF[oppdV[dv]],0) += feq[oppdV[dv +4]]*1.0 + feq[oppdV[dv +8]]*1.5;
                    G_calc.Node(i - (int)lb.CxF[dv+4], j - (int)lb.CyF[dv+4],0)           += feq[oppdV[dv +8]]*1.5;

                }

                if(marker.Cell(i,j,0) == FLUID)
                if(marker.Node(i + (int)lb.CxC[dv], j + (int)lb.CyC[dv],0) == SOLID){
                    G_calc.Cell(i,j,0)                                                          += feq[oppdV[dv]]   *0.5 + feq[oppdV[dv +4]]*1.0 + feq[oppdV[dv +8]]*1.5;
                    G_calc.Node(i + (int)lb.CxC[oppdV[dv]], j + (int)lb.CyC[oppdV[dv]],0)       += feq[oppdV[dv +4]]*1.0 + feq[oppdV[dv +8]]*1.5;
                    G_calc.Cell(i - (int)lb.CxC[dv +4],j - (int)lb.CyC[dv+4],0)                 += feq[oppdV[dv +8]]*1.5;
                }
            }
        }
        }
}



//meant to be done before advection
template<typename T,typename T1>
void diffuse_B_21_Solid(Grid_N_C_2D<T> &grid, lbmD2Q21<T> &lb, Grid_N_C_2D<T> &G_calc, Grid_N_C_2D<T1> &marker,int left , int right, int bottom, int top){

double feq[21] = {0};
Grid_N_C_2D<double> G_num(grid.n_x,grid.n_y,2,21);

get_equi(feq,lb,0.0,0.0,1.0,lb.theta0);

//note calculating the scattering probability
 for(int i = left;i < right ; i++){
        for(int j = bottom;j < top ; j++){
            
            // For SC
            for(int dv = 1; dv<5; dv++){

                if(marker.Node(i,j,0) == FLUID){
                if(marker.Node(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID){

                    G_num.Node(i,j,0)                                      += grid.Node(i,j,dv) *1 + grid.Node(i,j,dv + 4)* 2;              //for the closest node all the sc velocites will hit the wall
                    G_num.Node(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],0)    += grid.Node(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],dv + 4)*2;      // for the next closest the sc2 velocities will hit the walll
                
                }
                }

                if(marker.Cell(i,j,0) == FLUID)
                if(marker.Cell(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID){
                 
                    G_num.Cell(i,j,0)                                      += grid.Cell(i,j,dv) *1 + grid.Cell(i,j,dv + 4)* 2;          
                    G_num.Cell(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],0)    += grid.Cell(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],dv + 4)*2;

                }
            }
            
            // for BCC
            for(int dv = 9; dv<13; dv ++){
                if(marker.Node(i,j,0) == FLUID)
                if(marker.Cell(i + (int)lb.CxF[dv], j + (int)lb.CyF[dv],0) == SOLID){
                    G_num.Node(i,j,0)                                                    += grid.Node(i,j,dv)*0.5 + grid.Node(i,j,dv +4)*1.0 + grid.Node(i,j,dv +8)*1.5;
                    G_num.Cell(i + (int)lb.CxF[oppdV[dv]], j + (int)lb.CyF[oppdV[dv]],0) += grid.Cell( i + (int)lb.CxF[oppdV[dv]], j + (int)lb.CyF[oppdV[dv]] ,dv +4)*1.0 + grid.Cell( i + (int)lb.CxF[oppdV[dv]], j + (int)lb.CyF[oppdV[dv]],dv +8)*1.5;
                    G_num.Node(i - (int)lb.CxF[dv+4], j - (int)lb.CyF[dv+4],0)           += grid.Node(i - (int)lb.CxF[dv+4],j - (int)lb.CyF[dv+4],dv +8)*1.5;

                }

                if(marker.Cell(i,j,0) == FLUID)
                if(marker.Node(i + (int)lb.CxC[dv], j + (int)lb.CyC[dv],0) == SOLID){
                    G_num.Cell(i,j,0)                                                     += grid.Cell(i,j,dv)*0.5   + grid.Cell(i,j,dv +4)*1.0 + grid.Cell(i,j,dv +8)*1.5;
                    G_num.Node(i + (int)lb.CxC[oppdV[dv]], j + (int)lb.CyC[oppdV[dv]],0)  += grid.Node(i + (int)lb.CxC[oppdV[dv]],j + (int)lb.CyC[oppdV[dv]],dv +4)*1.0    + grid.Node(i + (int)lb.CxC[oppdV[dv]],j + (int)lb.CyC[oppdV[dv]],dv +8)*1.5;
                    G_num.Cell(i - (int)lb.CxC[dv +4],j - (int)lb.CyC[dv+4],0)            += grid.Cell(i - (int)lb.CxC[dv +4],j - (int)lb.CyC[dv+4],dv +8)*1.5;
                }
            }
        }
        }

    for(int i = left;i < right ; i++){
        for(int j = bottom;j < top ; j++){

            if(G_calc.Cell(i,j,0) != 0){
                    G_num.Cell(i,j,0) = G_num.Cell(i,j,0)/G_calc.Cell(i,j,0);

                    }
            
            if(G_calc.Node(i,j,0) != 0)
                    G_num.Node(i,j,0) = G_num.Node(i,j,0)/G_calc.Node(i,j,0);

        }}





        //note calculation of the population from the diffuse bounce back
        for(int i = left; i<= right; i++){
        for(int j = bottom; j<= top; j++){
            
            if(marker.Node(i,j,0) == FLUID)
                for(int dv =1; dv <5; dv++){
                    if(marker.Node(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID){
                        grid.Node(i + (int)lb.Cx[dv]  ,j + (int)lb.Cy[dv]  , oppdV[dv  ]) = G_num.Node(i,j,0) * feq[oppdV[dv  ]];
                        grid.Node(i + (int)lb.Cx[dv+4],j + (int)lb.Cy[dv+4], oppdV[dv+4]) = G_num.Node(i,j,0) * feq[oppdV[dv+4]];

                        grid.Node(i - (int)lb.Cx[dv] + (int)lb.Cx[dv+4] ,j - (int)lb.Cy[dv] + (int)lb.Cy[dv+4], oppdV[dv+4]) = G_num.Node(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],0) * feq[oppdV[dv+4]];

                    }



                //checking with respect to Node
                for(int dv = 9; dv < 13; dv++){
                    if(marker.Cell(i + (int)lb.CxF[dv], j + (int)lb.CyF[dv], 0) == SOLID){
                        grid.Cell(i + (int)lb.CxF[dv]   , j + (int)lb.CyF[dv]   ,oppdV[dv  ])   = G_num.Node(i,j,0) * feq[oppdV[dv  ]];
                        grid.Node(i + (int)lb.Cx[dv +4] , j + (int)lb.Cy[dv+ 4] ,oppdV[dv+4])   = G_num.Node(i,j,0) * feq[oppdV[dv+4]];
                        grid.Cell(i + (int)lb.CxF[dv +8], j + (int)lb.CyF[dv+ 8],oppdV[dv+8])   = G_num.Node(i,j,0) * feq[oppdV[dv+8]];


                        grid.Cell(i + (int)lb.CxF[oppdV[dv]] + (int)lb.Cx[dv +4] , j + (int)lb.CyF[oppdV[dv]] + (int)lb.Cy[dv+ 4] , oppdV[dv+4]) = G_num.Cell(i + (int)lb.CxF[oppdV[dv]] , j + (int)lb.CyF[oppdV[dv]],0) * feq[oppdV[dv+4]];
                        grid.Node(i + (int)lb.CxF[oppdV[dv]] + (int)lb.CxC[dv +8], j + (int)lb.CyF[oppdV[dv]] + (int)lb.CyC[dv+ 8], oppdV[dv+8]) = G_num.Cell(i + (int)lb.CxF[oppdV[dv]] , j + (int)lb.CyF[oppdV[dv]],0) * feq[oppdV[dv+8]];


                        grid.Cell(i - (int)lb.CxF[dv+4]  + (int)lb.CxF[dv +8]    , j - (int)lb.CyF[dv+4]   + (int)lb.CyF[dv+ 8]   , oppdV[dv+8]) = G_num.Node(i - (int)lb.CxF[dv+4] , j - (int)lb.CyF[dv+4], 0) * feq[oppdV[dv+8]];

                    }
                }
            }

            if(marker.Cell(i,j,0) == FLUID)
            for(int dv =1; dv <5; dv++){
                if(marker.Cell(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID){
                    grid.Cell(i + (int)lb.Cx[dv]  ,j + (int)lb.Cy[dv]  , oppdV[dv  ]) = G_num.Cell(i,j,0) * feq[oppdV[dv  ]];
                    grid.Cell(i + (int)lb.Cx[dv+4],j + (int)lb.Cy[dv+4], oppdV[dv+4]) = G_num.Cell(i,j,0) * feq[oppdV[dv+4]];

                    grid.Cell(i - (int)lb.Cx[dv] + (int)lb.Cx[dv+4] ,j - (int)lb.Cy[dv] + (int)lb.Cy[dv+4], oppdV[dv+4]) = G_num.Cell(i - (int)lb.Cx[dv],j - (int)lb.Cy[dv],0) * feq[oppdV[dv+4]];

                }
            }
            for(int dv = 9; dv < 13; dv++){
              
                //checking with respect to Cell
                if(marker.Node(i + (int)lb.CxC[dv], j + (int)lb.CyC[dv], 0) == SOLID){

                    grid.Node(i + (int)lb.CxC[dv]   , j + (int)lb.CyC[dv]   ,oppdV[dv])   = G_num.Cell(i,j,0) * feq[oppdV[dv]]  ;
                    grid.Cell(i + (int)lb.Cx [dv +4], j + (int)lb.Cy[dv+ 4] ,oppdV[dv+4]) = G_num.Cell(i,j,0) * feq[oppdV[dv+4]];
                    grid.Node(i + (int)lb.CxC[dv +8], j + (int)lb.CyC[dv+ 8],oppdV[dv+8]) = G_num.Cell(i,j,0) * feq[oppdV[dv+8]];

                    grid.Node(i + (int)lb.CxC[oppdV[dv]] + (int)lb.Cx [dv +4]  , j + (int)lb.CyC[oppdV[dv]] + (int)lb.Cy[dv+ 4] , oppdV[dv+4]) = G_num.Node(i + (int)lb.CxC[oppdV[dv]] , j + (int)lb.CyC[oppdV[dv]],0) * feq[oppdV[dv+4]];
                    grid.Cell(i + (int)lb.CxC[oppdV[dv]] + (int)lb.CxF[dv +8]  , j + (int)lb.CyC[oppdV[dv]] + (int)lb.CyF[dv+ 8], oppdV[dv+8]) = G_num.Node(i + (int)lb.CxC[oppdV[dv]] , j + (int)lb.CyC[oppdV[dv]],0) * feq[oppdV[dv+8]];

                    grid.Node(i - (int)lb.CxC[dv +4]     + (int)lb.CxC[dv +8]  ,j - (int)lb.CyC[dv+4]       + (int)lb.CyC[dv+ 8], oppdV[dv+8]) = G_num.Cell(i - (int)lb.CxC[dv +4] , j - (int)lb.CyC[dv+4],0) *feq[oppdV[dv+8]];

                }
            }
        }
        }
}





//setting momentum into the solid use it before the advection
template<typename T, typename T1>
void mom_in(Grid_N_C_2D<T> &grid, lbmD2Q21<T> &lb, Grid_N_C_2D<T1> &marker, double &jx, double &jy, int left , int right, int top, int bottom){
    
    int count = 0;

    jx = 0;
    jy = 0;

    for(int i = left ; i<right;i++){          //bounce back starting i 
        for(int j = bottom ; j<top;j++){





//SEEING WITH RESPECT TO NODE
             if(marker.Node(i,j,0)==FLUID){    
            //note contribution from sc1 and Sc3

            for(int dv = 1; dv< 5; dv++){       

                if(marker.Node(i + lb.Cx[dv], j + lb.Cy[dv],0) == SOLID){
                    jx += grid.Node(i,j,dv) *lb.Cx[dv];   count+=1;
                    jy += grid.Node(i,j,dv) *lb.Cy[dv];   count+=1;

                    jx += grid.Node(i,j,dv+4) *lb.Cx[dv+4];count+=1;
                    jy += grid.Node(i,j,dv+4) *lb.Cy[dv+4];count+=1;

                }
            }

            for(int dv = 9; dv<13; dv++){
                if(marker.Cell(i+ (int)lb.CxF[dv], j + (int)lb.CyF[dv],0)== SOLID){

                    jx += grid.Node(i,j,dv) * lb.Cx[dv];count+=1;
                    jy += grid.Node(i,j,dv) * lb.Cy[dv];count+=1;
                    
                    jx += grid.Node(i,j,dv+4) * lb.Cx[dv+4];count+=1;
                    jy += grid.Node(i,j,dv+4) * lb.Cy[dv+4];count+=1;

                    jx += grid.Node(i,j,dv+8) * lb.Cx[dv+8];count+=1;
                    jy += grid.Node(i,j,dv+8) * lb.Cy[dv+8];count+=1;



                    jx += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], dv +4) * lb.Cx[dv+4];count+=1;
                    jy += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], dv +4) * lb.Cy[dv+4];count+=1;

                    jx += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], dv +8) * lb.Cx[dv+8];count+=1;
                    jy += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], dv +8) * lb.Cy[dv+8];count+=1;

                    jx += grid.Node(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], dv+8)* lb.Cx[dv+8];count+=1;
                    jy += grid.Node(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], dv+8)* lb.Cy[dv+8];count+=1;



                }
            }
             }


//seeing with respect to cell
            if(marker.Cell(i,j,0) == FLUID){    
            //note contribution from sc1 and Sc3

            for(int dv = 1; dv< 5; dv++){

                if(marker.Cell(i + lb.Cx[dv], j + lb.Cy[dv],0) == SOLID){
                    jx += grid.Cell(i,j,dv) *lb.Cx[dv];count+=1;
                    jy += grid.Cell(i,j,dv) *lb.Cy[dv];count+=1;

                    jx += grid.Cell(i,j,dv+4) *lb.Cx[dv+4];count+=1;
                    jy += grid.Cell(i,j,dv+4) *lb.Cy[dv+4];count+=1;

                }
            }

            for(int dv = 9; dv<13; dv++){
                if(marker.Node(i+ (int)lb.CxC[dv], j + (int)lb.CyC[dv],0)== SOLID){


                    //bcc 
                    //checking the first boundary cell
                    jx += grid.Cell(i,j,dv) * lb.Cx[dv];count+=1;
                    jy += grid.Cell(i,j,dv) * lb.Cy[dv];count+=1;
                    
                    jx += grid.Cell(i,j,dv+4) * lb.Cx[dv+4];count+=1;
                    jy += grid.Cell(i,j,dv+4) * lb.Cy[dv+4];count+=1;

                    jx += grid.Cell(i,j,dv+8) * lb.Cx[dv+8];count+=1;
                    jy += grid.Cell(i,j,dv+8) * lb.Cy[dv+8];count+=1;



                    //checking the next node
                    jx += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], dv +4) * lb.Cx[dv+4];count+=1;
                    jy += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], dv +4) * lb.Cy[dv+4];count+=1;

                    jx += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], dv +8) * lb.Cx[dv+8];count+=1;
                    jy += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], dv +8) * lb.Cy[dv+8];count+=1;

                    //checking the next cell
                    jx += grid.Cell(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], dv+8)* lb.Cx[dv+8];count+=1;
                    jy += grid.Cell(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], dv+8)* lb.Cy[dv+8];count+=1;

                }
            }
             }
        }}

std::cout<<"count "<<count<< std:: endl;



}



//setting Momentum out of the solid ----------------use it after the advection
template<typename T, typename T1>
void mom_out(Grid_N_C_2D<T> &grid, lbmD2Q21<T> &lb, Grid_N_C_2D<T1> &marker, double &jx, double &jy, int left , int right, int top, int bottom){


   int count = 0;
    jx = 0;
    jy = 0;

    for(int i = left ; i<right;i++){          //bounce back starting i 
        for(int j = bottom ; j<top;j++){





//SEEING WITH RESPECT TO NODE
            if(marker.Node(i,j,0)==FLUID){    
            //note contribution from sc1 and Sc3

            for(int dv = 1; dv< 5; dv++){

                if(marker.Node(i + lb.Cx[dv], j + lb.Cy[dv],0) == SOLID){
                    jx += grid.Node(i,j,oppdV[dv]) *lb.Cx[oppdV[dv]];count+=1;
                    jy += grid.Node(i,j,oppdV[dv]) *lb.Cy[oppdV[dv]];count+=1;

                    jx += grid.Node(i,j,oppdV[dv]+4) *lb.Cx[oppdV[dv]+4];count+=1;
                    jy += grid.Node(i,j,oppdV[dv]+4) *lb.Cy[oppdV[dv]+4];count+=1;

                }
            }

            for(int dv = 9; dv<13; dv++){
                if(marker.Cell(i+ (int)lb.CxF[dv], j + (int)lb.CyF[dv],0)== SOLID){

                    jx += grid.Node(i,j,oppdV[dv]) * lb.Cx[oppdV[dv]];count+=1;
                    jy += grid.Node(i,j,oppdV[dv]) * lb.Cy[oppdV[dv]];count+=1;
                    
                    jx += grid.Node(i,j,oppdV[dv+4]) * lb.Cx[oppdV[dv+4]];count+=1;
                    jy += grid.Node(i,j,oppdV[dv+4]) * lb.Cy[oppdV[dv+4]];count+=1;
                    
                    jx += grid.Node(i,j,oppdV[dv+8]) * lb.Cx[oppdV[dv+8]];count+=1;
                    jy += grid.Node(i,j,oppdV[dv+8]) * lb.Cy[oppdV[dv+8]];count+=1;



                    jx += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], oppdV[dv+4]) * lb.Cx[oppdV[dv+4]];count+=1;
                    jy += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], oppdV[dv+4]) * lb.Cy[oppdV[dv+4]];count+=1;

                    jx += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], oppdV[dv+8]) * lb.Cx[oppdV[dv+8]];count+=1;
                    jy += grid.Cell(i + (int) lb.CxF[oppdV[dv]], j + (int) lb.CyF[oppdV[dv]], oppdV[dv+8]) * lb.Cy[oppdV[dv+8]];count+=1;

                    jx += grid.Node(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], oppdV[dv+8])* lb.Cx[oppdV[dv+8]];count+=1;
                    jy += grid.Node(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], oppdV[dv+8])* lb.Cy[oppdV[dv+8]];count+=1;


                }
            }
             }


//seeing with respect to cell
            if(marker.Cell(i,j,0) == FLUID){    
            //note contribution from sc1 and Sc3

            for(int dv = 1; dv< 5; dv++){

                if(marker.Cell(i + lb.Cx[dv], j + lb.Cy[dv],0) == SOLID){
                    jx += grid.Cell(i,j,oppdV[dv]) *lb.Cx[oppdV[dv]];count+=1;
                    jy += grid.Cell(i,j,oppdV[dv]) *lb.Cy[oppdV[dv]];count+=1;

                    jx += grid.Cell(i,j,oppdV[dv]+4) *lb.Cx[oppdV[dv]+4];count+=1;
                    jy += grid.Cell(i,j,oppdV[dv]+4) *lb.Cy[oppdV[dv]+4];count+=1;

                }
            }

            for(int dv = 9; dv<13; dv++){
                if(marker.Node(i+ (int)lb.CxC[dv], j + (int)lb.CyC[dv],0)== SOLID){



                    //bcc 
                    //checking the first boundary cell
                    jx += grid.Cell(i,j,oppdV[dv]) * lb.Cx[oppdV[dv]];count+=1;
                    jy += grid.Cell(i,j,oppdV[dv]) * lb.Cy[oppdV[dv]];count+=1;
                    
                    jx += grid.Cell(i,j,oppdV[dv]+4) * lb.Cx[oppdV[dv]+4];count+=1;
                    jy += grid.Cell(i,j,oppdV[dv]+4) * lb.Cy[oppdV[dv]+4];count+=1;

                    jx += grid.Cell(i,j,oppdV[dv]+8) * lb.Cx[oppdV[dv]+8];count+=1;
                    jy += grid.Cell(i,j,oppdV[dv]+8) * lb.Cy[oppdV[dv]+8];count+=1;



                    //checking the next node
                    jx += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], oppdV[dv +4]) * lb.Cx[oppdV[dv+4]];count+=1;
                    jy += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], oppdV[dv +4]) * lb.Cy[oppdV[dv+4]];count+=1;

                    jx += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], oppdV[dv +8]) * lb.Cx[oppdV[dv+8]];count+=1;
                    jy += grid.Node(i + (int) lb.CxC[oppdV[dv]], j + (int) lb.CyC[oppdV[dv]], oppdV[dv +8]) * lb.Cy[oppdV[dv+8]];count+=1;

                    //checking the next cell
                    jx += grid.Cell(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], oppdV[dv+8])* lb.Cx[oppdV[dv+8]];count+=1;
                    jy += grid.Cell(i - (int)2*lb.Cx[dv], j - (int)2*lb.Cy[dv], oppdV[dv+8])* lb.Cy[oppdV[dv+8]];count+=1;

                }
            }
             }
        }}

        std::cout<<"count "<<count<< std:: endl;

}



// setting made to run after advection
template<typename T, typename T1>
void BB_21(Grid_N_C_2D<T> &grid, lbmD2Q21<T> &lb, Grid_N_C_2D<T1> &marker, int left , int right, int top, int bottom){

    for(int i = left ; i<right;i++){                                                                                                                  //bounce back starting i 
        for(int j = bottom ; j<top; j++){

            if(marker.Node(i,j,0) == FLUID){

                //note we are on the node
                for(int dv = 1; dv< 5; dv++){

                    if(marker.Node(i + (int)lb.Cx[dv], j + (int)lb.Cy[dv],0) == SOLID){

                        grid.Node(i,j,oppdV[dv   ]) = grid.Node(i + lb.Cx[dv]   , j + lb.Cy[dv]   , dv  );
                        grid.Node(i,j,oppdV[dv +4]) = grid.Node(i + lb.Cx[dv +4], j + lb.Cy[dv +4], dv+4);

                        grid.Node(i + lb.Cx[oppdV[dv]], j + lb.Cy[oppdV[dv]], oppdV[dv+4]) = grid.Node(i + lb.Cx[dv], j+ lb.Cy[dv], dv+4);

                    }
                
                

                for(int dv = 9; dv<13; dv++){
                    if(marker.Cell(i + (int) lb.CxF[dv], j + (int) lb.CyF[dv], 0 )  == SOLID){

                        grid.Node(i , j, oppdV[dv  ]) = grid.Cell(i + (int) lb.CxF[dv  ], j + (int) lb.CyF[dv  ], dv  );
                        grid.Node(i , j, oppdV[dv+4]) = grid.Node(i + (int) lb.CxF[dv+4], j + (int) lb.CyF[dv+4], dv+4);
                        grid.Node(i , j, oppdV[dv+8]) = grid.Cell(i + (int) lb.CxF[dv+8], j + (int) lb.CyF[dv+8], dv+8);

                        grid.Cell(i + lb.CxF[oppdV[dv]], j + lb.CyF[oppdV[dv]], oppdV[dv+4]) = grid.Cell(i + (int) lb.CxF[dv  ], j + (int) lb.CyF[dv  ], dv+4);
                        grid.Cell(i + lb.CxF[oppdV[dv]], j + lb.CyF[oppdV[dv]], oppdV[dv+8]) = grid.Node(i + (int) lb.CxF[dv+4], j + (int) lb.CyF[dv+4], dv+8);

                        grid.Node(i + lb.CxF[oppdV[dv+4]], j + lb.CyF[oppdV[dv+4]], oppdV[dv+8]) = grid.Cell(i + (int) lb.CxF[dv  ], j + (int) lb.CyF[dv], dv+8);

                    }
                }
            }


        }
        
        if(marker.Cell(i,j,0) == FLUID){

            //note we are on the node
            for(int dv = 1; dv< 5; dv++){

                if(marker.Cell(i + (int)lb.Cx[dv], j + (int)lb.Cy[dv],0) == SOLID){

                    grid.Cell(i,j,oppdV[dv   ]) = grid.Cell(i + lb.Cx[dv]   , j + lb.Cy[dv]   , dv  );
                    grid.Cell(i,j,oppdV[dv +4]) = grid.Cell(i + lb.Cx[dv +4], j + lb.Cy[dv +4], dv+4);

                    grid.Cell(i + lb.Cx[oppdV[dv]], j + lb.Cy[oppdV[dv]], oppdV[dv+4]) = grid.Cell(i + lb.Cx[dv], j+ lb.Cy[dv], dv+4);

                }
            
            

            for(int dv = 9; dv<13; dv++){
                if(marker.Node(i + (int) lb.CxC[dv], j + (int) lb.CyC[dv], 0 )  == SOLID){

                    grid.Cell(i , j, oppdV[dv  ]) = grid.Node(i + (int) lb.CxC[dv  ], j + (int) lb.CyC[dv  ], dv  );
                    grid.Cell(i , j, oppdV[dv+4]) = grid.Cell(i + (int) lb.CxC[dv+4], j + (int) lb.CyC[dv+4], dv+4);
                    grid.Cell(i , j, oppdV[dv+8]) = grid.Node(i + (int) lb.CxC[dv+8], j + (int) lb.CyC[dv+8], dv+8);

                    grid.Node(i + lb.CxC[oppdV[dv]], j + lb.CyC[oppdV[dv]], oppdV[dv+4]) = grid.Node(i + (int) lb.CxC[dv  ], j + (int) lb.CyC[dv  ], dv+4);
                    grid.Node(i + lb.CxC[oppdV[dv]], j + lb.CyC[oppdV[dv]], oppdV[dv+8]) = grid.Cell(i + (int) lb.CxC[dv+4], j + (int) lb.CyC[dv+4], dv+8);

                    grid.Cell(i + lb.CxC[oppdV[dv+4]], j + lb.CyC[oppdV[dv+4]], oppdV[dv+8]) = grid.Node(i + (int) lb.CxC[dv  ], j + (int) lb.CyC[dv], dv+8);


                }
            }
        }


        }




        }
  }


}









/*as in the paper*/
template<typename T>
void diffuse_B_21_new(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lb, double u_top,double u_bottom, double theta_top, double theta_bottom
                            ){


double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
int j;
double feq[21] = {0};



//refactor -------------------------TOP WALL-------------------------------/// 
get_equi(feq,lb,u_top,0.0,rho,theta_top);
for(int i = grid.noghost; i<grid.n_x_node-grid.noghost;i++){
 nx =  0; ny= -1;
 
 //outgoing populations Cell or Node for last and second last
 int outC_last[] = {2,6,9,10,13,14,17,18};
 int outC_2nd_last[] = {6,17,18};
 int outN_last[] = {2,6,13,14,17,18};
 int outN_2nd_last[] = {6};

double G_num_cell = 0, G_den_cell = 0, G_num_node = 0, G_den_node = 0,G_node = 0,G_cell = 0;

//calculating G

    for(int k = 0; k<sizeof(outC_last)/sizeof(outC_last[0]); k++){
        G_num_cell += grid.Cell(i,grid.n_y_cell - grid.noghost -1, outC_last[k]) * abs(lb.Cy[outC_last[k]]);
        G_den_cell += feq[outC_last[k]] *abs(lb.Cy[outC_last[k]]);       
    }
        
    for(int k = 0; k<sizeof(outC_2nd_last)/sizeof(outC_2nd_last[0]); k++){
        G_num_cell += grid.Cell(i,grid.n_y_cell - grid.noghost -2, outC_2nd_last[k]) * abs(lb.Cy[outC_2nd_last[k]]);
        G_den_cell += feq[outC_2nd_last[k]] *abs(lb.Cy[outC_2nd_last[k]]) ;       
    }

        G_cell = G_num_cell/G_den_cell;

    for(int k = 0; k<sizeof(outN_last)/sizeof(outN_last[0]); k++){
        G_num_node += grid.Node(i,grid.n_y_node - grid.noghost -1, outN_last[k]) * abs(lb.Cy[outN_last[k]]);
        G_den_node += feq[outN_last[k]] *abs(lb.Cy[outN_last[k]]) ;       
    }

    for(int k = 0; k<sizeof(outN_2nd_last)/sizeof(outN_2nd_last[0]); k++){
        G_num_node += grid.Node(i,grid.n_y_node - grid.noghost -2, outN_2nd_last[k]) * abs(lb.Cy[outN_2nd_last[k]]);    
        G_den_node += feq[outN_2nd_last[k]] *abs(lb.Cy[outN_2nd_last[k]]) ;       
    }

        G_node = G_num_node / G_den_node;

    // pushing to the last node and cell
    int j = grid.n_y_node - grid.noghost -1;

    for(int k = 0; k < 2; k++) 
        grid.Cell(i + (int)lb.Cx[outC_last[k]], j +  (int)lb.Cy[outC_last[k]], oppdV[outC_last[k]]) = G_cell * feq[oppdV[outC_last[k]]];

    for(int k = 2; k < 8;k++)
        grid.Node(i + (int)lb.CxC[outC_last[k]], j + (int)lb.CyC[outC_last[k]], oppdV[outC_last[k]]) = G_cell * feq[oppdV[outC_last[k]]];

    for(int k = 0; k < 2;k++)
        grid.Node(i + (int)lb.Cx[outN_last[k]], j + (int)lb.Cy[outN_last[k]] , oppdV[outN_last[k]]) = G_node * feq[oppdV[outN_last[k]]];

    for(int k =2; k < 6; k++)
        grid.Cell(i + (int)lb.CxF[outN_last[k]], j + (int)lb.CyF[outN_last[k]],oppdV[outN_last[k]]) = G_node * feq[oppdV[outN_last[k]]];

    //pushing to the last cell
    j = grid.n_y_node -grid.noghost - 2;

    for(int k = 0; k<1; k++)
        grid.Cell(i + (int)lb.Cx[outC_2nd_last[k]], j +  (int)lb.Cy[outC_2nd_last[k]], oppdV[outC_2nd_last[k]]) = G_cell * feq[oppdV[outC_2nd_last[k]]];

    for(int k = 1; k<3; k++)
        grid.Node(i + (int)lb.CxC[outC_2nd_last[k]], j + (int)lb.CyC[outC_2nd_last[k]], oppdV[outC_2nd_last[k]]) = G_cell * feq[oppdV[outC_2nd_last[k]]];

    for(int k = 0; k<1; k++)
        grid.Node(i + (int)lb.Cx[outN_2nd_last[k]], j +  (int)lb.Cy[outN_2nd_last[k]], oppdV[outN_2nd_last[k]]) = G_node * feq[oppdV[outN_2nd_last[k]]];

    for(int k = 0; k<0; k++)
        grid.Cell(i + (int)lb.CxF[outN_2nd_last[k]], j +  (int)lb.CyF[outN_2nd_last[k]], oppdV[outN_2nd_last[k]]) = G_node * feq[oppdV[outN_2nd_last[k]]];



}
    


//refactor -----------bottom-------------------------------//

get_equi(feq,lb,u_bottom,0.0,rho,theta_bottom);
for(int i = grid.noghost; i<grid.n_x_node-grid.noghost;i++){
 nx =  0; ny= 1;
 
 //outgoing populations Cell or Node for last and second last
 int outC_last[] = {4,8,15,16,19,20};
 int outC_2nd_last[] = {8};
 int outN_last[] = {4,8,11,12,15,16,19,20};
 int outN_2nd_last[] = {8,19,20};

double G_num_cell = 0, G_den_cell = 0, G_num_node = 0, G_den_node = 0,G_node = 0,G_cell = 0;

//calculating G

    for(int k = 0; k<sizeof(outC_last)/sizeof(outC_last[0]); k++){
        G_num_cell += grid.Cell(i,grid.noghost, outC_last[k]) * abs(lb.Cy[outC_last[k]]);
        G_den_cell += feq[outC_last[k]] *abs(lb.Cy[outC_last[k]]) ;       
    }
        
     for(int k = 0; k<sizeof(outC_2nd_last)/sizeof(outC_2nd_last[0]); k++){
        G_num_cell += grid.Cell(i,grid.noghost +1, outC_2nd_last[k]) * abs(lb.Cy[outC_2nd_last[k]]);
        G_den_cell += feq[outC_2nd_last[k]] *abs(lb.Cy[outC_2nd_last[k]]) ;       
    }

        G_cell = G_num_cell/G_den_cell;

     for(int k = 0; k<sizeof(outN_last)/sizeof(outN_last[0]); k++){
        G_num_node += grid.Node(i,grid.noghost, outN_last[k]) * abs(lb.Cy[outN_last[k]]);
        G_den_node += feq[outN_last[k]] *abs(lb.Cy[outN_last[k]]) ;       
    }

    for(int k = 0; k<sizeof(outN_2nd_last)/sizeof(outN_2nd_last[0]); k++){
            G_num_node += grid.Node(i,grid.noghost +1, outN_2nd_last[k]) * abs(lb.Cy[outN_2nd_last[k]]);    
            G_den_node += feq[outN_2nd_last[k]] *abs(lb.Cy[outN_2nd_last[k]]) ;       
        }

        G_node = G_num_node / G_den_node;

// reassigning the populations
    // pushing to the last node and cell
    int j = grid.noghost;

    for(int k=0; k< 2; k++) 
        grid.Cell(i + (int)lb.Cx[outC_last[k]], j +  (int)lb.Cy[outC_last[k]], oppdV[outC_last[k]]) = G_cell * feq[oppdV[outC_last[k]]];

    for(int k = 2; k<6;k++)
        grid.Node(i + (int)lb.CxC[outC_last[k]], j + (int)lb.CyC[outC_last[k]], oppdV[outC_last[k]]) = G_cell * feq[oppdV[outC_last[k]]];

    
    for(int k = 0; k < 2;k++)
        grid.Node(i + (int)lb.Cx[outN_last[k]], j + (int)lb.Cy[outN_last[k]] , oppdV[outN_last[k]]) = G_node * feq[oppdV[outN_last[k]]];

    for(int k =2; k<8; k++)
        grid.Cell(i + (int)lb.CxF[outN_last[k]], j + (int)lb.CyF[outN_last[k]],oppdV[outN_last[k]]) = G_node * feq[oppdV[outN_last[k]]];

    //pushing to the last cell
    j = grid.noghost +1;

    for(int k = 0; k<1; k++)
        grid.Cell(i + (int)lb.Cx[outC_2nd_last[k]], j +  (int)lb.Cy[outC_2nd_last[k]], oppdV[outC_2nd_last[k]]) = G_cell * feq[oppdV[outC_2nd_last[k]]];

    for(int k = 0; k<0; k++)
        grid.Node(i + (int)lb.CxC[outC_2nd_last[k]], j + (int)lb.CyC[outC_2nd_last[k]], oppdV[outC_2nd_last[k]]) = G_cell * feq[oppdV[outC_2nd_last[k]]];

    for(int k = 0; k<1; k++)
        grid.Node(i + (int)lb.Cx[outN_2nd_last[k]], j +  (int)lb.Cy[outN_2nd_last[k]], oppdV[outN_2nd_last[k]]) = G_node * feq[oppdV[outN_2nd_last[k]]];

    for(int k = 1; k<3; k++)
        grid.Cell(i + (int)lb.CxF[outN_2nd_last[k]], j +  (int)lb.CyF[outN_2nd_last[k]], oppdV[outN_2nd_last[k]]) = G_node * feq[oppdV[outN_2nd_last[k]]];

}
}

template<typename T,typename T1>
void equi_correction(Grid_N_C_2D<T> &grid,Grid_N_C_2D<T1> &marker,lbmD2Q21<T> &lb, int left , int right, int bottom, int top){

double feq[21] = {0};
get_equi(feq, lb,0.0,0.0,1.0,lb.theta0);
for(int i = left ; i < right;i++)
    for(int j = bottom ; j < top; j++){
        
        if(marker.Node(i,j,0) == SOLID)
        for(int dv = 0; dv<21; dv++)
            grid.Node(i,j,dv) = feq[dv];

        if(marker.Cell(i,j,0) == SOLID)
        for(int dv = 0; dv<21; dv++)
            grid.Cell(i,j,dv) = feq[dv];
    }

}






template<typename T>
void stationary_correction(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21){
    //refactor ---------TOP -----------------//
    //note   stationary popuulation correction
    ///   last------------------//
    for(int i = grid.noghost; i<grid.n_x_node-grid.noghost; i++){
    for(int j = grid.n_y_node - grid.noghost -5 ; j<= grid.n_y_node - grid.noghost -1; j++ ){
    double sum = 0;
    for(int k = 1; k<grid.d_v; k++)
        sum += grid.Node(i,j,k);

    grid.Node(i,j,dV_ZERO_ZERO) = 1 - sum;
    
    sum = 0;
    for(int k = 1; k<grid.d_v; k++)
        sum += grid.Cell(i,j,k);

    grid.Cell(i,j,dV_ZERO_ZERO) = 1 - sum;
    }
    }
   
   
    //refactor  ------BOTTOM-------------//
    //note   stationary popuulation correction
    ///   last------------------//
    for(int i = grid.noghost ; i < grid.n_x_node-grid.noghost; i++){
    for(int j = grid.noghost ; j<= grid.noghost+5; j++ ){
    double sum = 0;
    for(int k = 1; k < grid.d_v; k++)
        sum += grid.Node(i,j,k);

    grid.Node(i,j,dV_ZERO_ZERO) = 1 - sum;
    
    sum = 0;
    for(int k = 1; k<grid.d_v; k++)
        sum += grid.Cell(i,j,k);

    grid.Cell(i,j,dV_ZERO_ZERO) = 1 - sum;
    }
    }
}










template<typename T>
void diffuse_B_21(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21, double u_top,double u_bottom, double theta_bottom, double theta_top){
double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
int j;
double feq[21] = {0};


////--------------------------------------------------top wall---------------------------------------------////
//---------------------------------------------------------------------------------------------------------///
for(int i = grid.noghost+2; i<grid.n_x_node-grid.noghost-2;i++){
    nx =  0; ny= -1;
    ux =u_top; uy = 0;
    j = grid.n_y_node - 1;

    get_equi(feq,lbD2Q21,ux,uy,rho, theta_top);

    //for the top most node

G = (grid.Node(i,j -2, dV_ZERO_P2) + grid.Cell(i+1,j -2, dV_MH3_PH3)  + grid.Cell(i -2,j -2,dV_PH3_PH3)) /
    (feq[dV_ZERO_M2] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);
    
    grid.Node(i,j,dV_MH3_MH3) = feq[dV_MH3_MH3] *G;
    grid.Node(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;
    grid.Node(i,j,dV_PH3_MH3) = feq[dV_PH3_MH3] *G;
    // for the cell
    G = (grid.Cell(i,j -2, dV_ZERO_P2))/ (feq[dV_ZERO_M2]);
    grid.Cell(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;



//next bottom node
    j = grid.n_y_node -2;
    G =(grid.Node(i ,j -1,dV_ZERO_P1) + grid.Node(i  ,j-2, dV_ZERO_P2)+
        grid.Cell(i ,j-1, dV_MH1_PH1) + grid.Cell(i-1,j-1, dV_PH1_PH1) +
        grid.Node(i -1, j-1,dV_P1_P1) + grid.Node(i+1, j-1, dV_M1_P1)+
        grid.Cell(i +1, j -2, dV_MH3_PH3) +grid.Cell(i-2,j-2,dV_PH3_PH3))/
        (feq[dV_ZERO_M1] + feq[dV_ZERO_M2] + feq[dV_PH1_MH1] +feq[dV_MH1_MH1]+feq[dV_M1_M1] + feq[dV_P1_M1]+ feq[dV_PH3_MH3] + feq[dV_MH3_MH3]) ; 
        // std::cout<<G<<std::endl;
    for(int dv = 0; dv<21; dv++){
        if(lbD2Q21.Cx[dv]*nx + lbD2Q21.Cy[dv]*ny > 0){
            grid.Node(i,j,dv) = feq[dv] * G;
        }
    }

// for the cell
G = (   grid.Cell(i ,j-1, dV_ZERO_P1) + grid.Cell(i,j -1, dV_ZERO_P2) +
        grid.Cell(i -1, j-1,dV_P1_P1) + grid.Cell(i+1, j-1, dV_M1_P1) +
        grid.Node(i+2,j-1,dV_MH3_PH3) + grid.Node(i-1,j-1,dV_PH3_PH3))/
        (feq[dV_ZERO_M1] + feq[dV_ZERO_M2] +feq[dV_M1_M1] + feq[dV_P1_M1]+feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);

    grid.Cell(i,j,dV_ZERO_M1) = feq[dV_ZERO_M1] *G;
    grid.Cell(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;
    grid.Cell(i,j,dV_P1_M1)   = feq[dV_P1_M1]   *G;
    grid.Cell(i,j,dV_M1_M1)   = feq[dV_M1_M1]   *G;
    grid.Cell(i,j,dV_PH3_MH3) = feq[dV_PH3_MH3] *G;
    grid.Cell(i,j,dV_MH3_MH3) = feq[dV_MH3_MH3] *G;

}



//-------------------------------------------------------------bottom wall-------------------------------------//
//-------------------------------------------------------------bottom wall-------------------------------------//
for(int i = grid.noghost +2 ; i<grid.n_x_node-grid.noghost-2;i++){

    get_equi(feq,lbD2Q21,ux,uy,rho,theta_bottom);
    nx = 0;  ny = 1;
    ux = u_bottom; uy = 0;
    j = 0;

    //0th node
    G = (grid.Node(i,j+2,dV_ZERO_M2))/(feq[dV_ZERO_P2]);
    grid.Node(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] *G;

    //Cell
    G = (grid.Cell(i, j +2, dV_ZERO_M2) + grid.Node(i+2, j+2, dV_MH3_MH3) + grid.Node(i-1, j +2, dV_PH3_MH3))/
        (feq[dV_ZERO_P2] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);

    grid.Cell(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] * G;
    grid.Cell(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] * G;
    grid.Cell(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] * G;

    j =1;
    //1st
    G = (   grid.Node(i  ,j +1, dV_ZERO_M1) + grid.Node(i,j +2, dV_ZERO_M2) +
            grid.Node(i+1, j+1, dV_M1_M1  ) + grid.Node(i-1,j+1,dV_P1_M1) +
            grid.Cell(i +1, j+1,dV_MH3_MH3) + grid.Cell(i-2,j+1, dV_PH3_MH3))/
        (feq[dV_ZERO_P1] + feq[dV_ZERO_P2] +feq[dV_P1_P1] + feq[dV_M1_P1]+ feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);

    grid.Node(i,j,dV_ZERO_P1) = feq[dV_ZERO_P1] *G;
    grid.Node(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] *G;
    grid.Node(i,j,dV_P1_P1)   = feq[dV_P1_P1]   *G;
    grid.Node(i,j,dV_M1_P1)   = feq[dV_M1_P1]   *G;
    grid.Node(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] *G;
    grid.Node(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] *G;

    G =(grid.Cell(i, j +1, dV_ZERO_M1) +grid.Cell(i,j +2, dV_ZERO_M2) + 
        grid.Node(i+1,j+1,dV_MH1_MH1) + grid.Node(i,j+1, dV_PH1_MH1) + 
        grid.Cell(i+1, j+1, dV_M1_M1)  +grid.Cell(i-1,j+1,dV_P1_M1) +
        grid.Node(i+2, j +2, dV_MH3_MH3) +grid.Node(i-1, j+2, dV_PH3_MH3))/
        (feq[dV_ZERO_P1] + feq[dV_ZERO_P2] + feq[dV_PH1_PH1] + feq[dV_MH1_PH1] +feq[dV_P1_P1] + feq[dV_M1_P1]+ feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);

for(int dv = 0; dv<21; dv++){
        if(lbD2Q21.Cx[dv]*nx + lbD2Q21.Cy[dv]*ny > 0){
            grid.Cell(i,j,dv) = feq[dv] * G;
        }
    }
}
}



// template<typename T>
// void BB_solid(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21){



//     for(int i = grid.noghost +2 ; i<grid.n_x_node-grid.noghost-2;i++){          //bounce back starting i 
//         for(int i = grid.noghost +2 ; i<grid.n_x_node-grid.noghost-2;i++){      //bounce back stanting and ending j

//             if(lbgrid.marker_node(i,j)==0){ //////checking whether the node is fluid        
//             ///----------------------enters for the fluid nodes
//                             //////////////////////---------------------------------SC1, SC3, bounceback
//                     if(lbgrid.marker_node(i-1,j) == 1){

                                                


//             }
//             }

//         }
//     }



// }










template<typename T>
void slip_wall_bb(Grid_N_C_2D<T> &grid,lbmD2Q21<T> d2q21 ){

int j;

/*! \brief  velocities*/
int wcv[21]= {0,3,4,1,2,7,8,5,6,12,11,10,9,16,15,14,13,20,19,18,17}; //wall corresponding velocities

int edgeshift  =0;

for(int i= grid.noghost + edgeshift ; i< grid.n_x_node - grid.noghost -edgeshift ; i++){

    //====================================TOP===============================//
    j = grid.n_y_node -grid.noghost -1;

    //for sc1 and sc2
    for(int dv = 1; dv < 9; dv++){  //only go tilll the sc2+ve velocity
        //will shift both the single speed and multispeed 
        
        grid.Cell(i , j +1, oppdV[dv]) = grid.Cell(i  ,j  ,dv);   //for the last cell
        grid.Cell(i , j +2, oppdV[dv]) = grid.Cell(i  ,j-1,dv);   //for the second last cell only the sc2 velocity

        grid.Node(i , j +1, oppdV[dv]) = grid.Node(i  ,j  ,dv);
        grid.Node(i , j +2, oppdV[dv]) = grid.Node(i  ,j-1,dv);

    }
    }




                           
 for(int i= grid.noghost + edgeshift ; i< grid.n_x_node - grid.noghost  -edgeshift ; i++){
        //for the BCC
        for(int dv = 9; dv < 11 ; dv ++){

            grid.Node(i + (int)d2q21.CxC[dv], j + (int) d2q21.CyC[dv],  wcv[dv])      = grid.Cell(i,j,dv); //BCC 0.5

            grid.Cell(i , j +1  ,wcv[dv+4])   = grid.Cell(i ,j  ,dv+4);  //BCC 1.0
            grid.Node(i , j +2  ,wcv[dv+4])   = grid.Node(i ,j  ,dv+4);  //BCC 1.0
            

            grid.Node(i + (int)d2q21.CxC[dv]    , j + (int) d2q21.CyC[dv]   ,   wcv[dv +8]) = grid.Cell(i,j,dv+8);
        }

            grid.Node(i +1, j+2, dV_PH3_MH3) = grid.Cell(i, j -1, dV_PH3_PH3);
            grid.Cell(i   , j+1, dV_PH3_MH3) = grid.Node(i  ,j  ,dV_PH3_PH3);

            grid.Node(i   , j +2, dV_MH3_MH3) = grid.Cell(i, j -1, dV_MH3_PH3);
            grid.Cell(i-1 , j +1, dV_MH3_MH3) = grid.Node(i, j   , dV_MH3_PH3);

    }



    //====================  bottom  =====================================================//////////
    j = grid.noghost;
 for(int i= grid.noghost + edgeshift ; i< grid.n_x_node - grid.noghost -edgeshift; i++){

        for(int dv = 0; dv < 9; dv++){  

        grid.Node(i , j -1, oppdV[dv]) = grid.Node(i  ,j  ,dv);   //for the last cell
        grid.Node(i , j -2, oppdV[dv]) = grid.Node(i  ,j+1,dv);   //fot the second last cell only the sc2 velocity 

        grid.Cell(i , j -1, oppdV[dv]) = grid.Cell(i  ,j  ,dv);
        grid.Cell(i , j -2, oppdV[dv]) = grid.Cell(i  ,j+1,dv);
    
    }}
    


    j = grid.noghost;
 for(int i= grid.noghost + edgeshift; i< grid.n_x_node - grid.noghost -edgeshift ; i++){
          for(int dv = 11; dv < 13 ; dv ++){

            grid.Cell(i + (int)d2q21.CxF[dv], j + (int) d2q21.CyF[dv],  wcv[dv])   = grid.Node(i,j,dv);

            grid.Node(i , j -1  ,wcv[dv+4])   = grid.Node(i ,j  ,dv+4);  
            grid.Cell(i , j -2  ,wcv[dv+4])   = grid.Cell(i ,j  ,dv+4);

            grid.Cell(i + (int)d2q21.CxF[dv], j + (int) d2q21.CyF[dv],  wcv[dv+8])   = grid.Node(i,j,dv+8);
        }

            grid.Cell(i  ,j -2, dV_PH3_PH3) = grid.Node(i , j+1, dV_PH3_MH3);
            grid.Cell(i-1,j -2, dV_MH3_PH3) = grid.Node(i , j+1, dV_MH3_MH3);

            grid.Node(i +1,j - 1, dV_PH3_PH3) = grid.Cell(i,j, dV_PH3_MH3);
            grid.Node(i   ,j - 1, dV_MH3_PH3) = grid.Cell(i,j, dV_MH3_MH3);
    }
}







template<typename T>
void periodic_21(Grid_N_C_2D<T> &grid){

        //periodic
for(int j = 0; j<grid.n_y_node; j++)
{       
        //-------------------------------------------right to left -----------//
        //SC1
        grid.Node(grid.noghost -1,j,dV_P1_ZERO) = grid.Node( grid.n_x_node - (grid.noghost +1),j,dV_P1_ZERO); ///copying to first adjacent ghost node from last node physical domain
        grid.Cell(grid.noghost -1,j,dV_P1_ZERO) = grid.Cell( grid.n_x_node - (grid.noghost +1),j,dV_P1_ZERO);
        
        //SC2
        grid.Node(0,j,dV_P2_ZERO) = grid.Node((grid.n_x_cell - grid.noghost) -2,j,dV_P2_ZERO);         
        grid.Cell(0,j,dV_P2_ZERO) = grid.Cell((grid.n_x_cell - grid.noghost) -2,j,dV_P2_ZERO);

        grid.Node(1,j,dV_P2_ZERO) = grid.Node((grid.n_x_cell - grid.noghost) -1,j,dV_P2_ZERO);
        grid.Cell(1,j,dV_P2_ZERO) = grid.Cell((grid.n_x_cell - grid.noghost) -1,j,dV_P2_ZERO);

        //BCC1/2
        grid.Cell(grid.noghost -1,j,dV_PH1_PH1) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_PH1_PH1); 
        grid.Cell(grid.noghost -1,j,dV_PH1_MH1) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_PH1_MH1);

        //BCC1
        grid.Node(grid.noghost -1, j, dV_P1_P1) = grid.Node(grid.n_x_node - (grid.noghost +1), j, dV_P1_P1);
        grid.Node(grid.noghost -1, j, dV_P1_M1) = grid.Node(grid.n_x_node - (grid.noghost +1), j, dV_P1_M1);

        grid.Cell(grid.noghost -1, j, dV_P1_P1) = grid.Cell(grid.n_x_node - (grid.noghost +1), j, dV_P1_P1);
        grid.Cell(grid.noghost -1, j, dV_P1_M1) = grid.Cell(grid.n_x_node - (grid.noghost +1), j, dV_P1_M1);



        //BCC3/2
        grid.Cell(1,j,dV_PH3_PH3) = grid.Cell((grid.n_x_cell-grid.noghost) - 1,j, dV_PH3_PH3);
        grid.Cell(1,j,dV_PH3_MH3) = grid.Cell((grid.n_x_cell-grid.noghost) - 1,j, dV_PH3_MH3);

        grid.Cell(0,j,dV_PH3_PH3) = grid.Cell((grid.n_x_cell-grid.noghost) - 2,j, dV_PH3_PH3);
        grid.Cell(0,j,dV_PH3_MH3) = grid.Cell((grid.n_x_cell-grid.noghost) - 2,j, dV_PH3_MH3);

        grid.Node(1,j,dV_PH3_PH3) = grid.Node((grid.n_x_cell-grid.noghost) - 1,j, dV_PH3_PH3);
        grid.Node(1,j,dV_PH3_MH3) = grid.Node((grid.n_x_cell-grid.noghost) - 1,j, dV_PH3_MH3);

        //------------------------------left to right---------------------//
        //SC1
        grid.Node(grid.n_x_cell - grid.noghost,j,dV_M1_ZERO) = grid.Node(grid.noghost,j,dV_M1_ZERO);
        grid.Cell(grid.n_x_cell - grid.noghost,j,dV_M1_ZERO) = grid.Cell(grid.noghost,j,dV_M1_ZERO);

        //Sc2
        grid.Node(grid.n_x_cell - (grid.noghost)  ,j,dV_M2_ZERO)   = grid.Node(grid.noghost  ,j,dV_M2_ZERO);
        grid.Cell(grid.n_x_cell - (grid.noghost)  ,j,dV_M2_ZERO)   = grid.Cell(grid.noghost  ,j,dV_M2_ZERO);

        grid.Node((grid.n_x_cell - grid.noghost) +1,j,dV_M2_ZERO)   = grid.Node(grid.noghost+1,j,dV_M2_ZERO);
        grid.Cell((grid.n_x_cell - grid.noghost) +1,j,dV_M2_ZERO)   = grid.Cell(grid.noghost+1,j,dV_M2_ZERO);

        //BCC-1/2
        grid.Node(grid.n_x_cell - (grid.noghost),j,dV_MH1_PH1) = grid.Node(grid.noghost,j,dV_MH1_PH1);
        grid.Node(grid.n_x_cell - (grid.noghost),j,dV_MH1_MH1) = grid.Node(grid.noghost,j,dV_MH1_MH1);

        //BCC1
        grid.Node(grid.n_x_node - grid.noghost, j ,dV_M1_M1) = grid.Node(grid.noghost, j , dV_M1_M1);
        grid.Node(grid.n_x_node - grid.noghost, j ,dV_M1_P1) = grid.Node(grid.noghost, j , dV_M1_P1);

        grid.Cell(grid.n_x_node - grid.noghost, j ,dV_M1_M1) = grid.Cell(grid.noghost, j , dV_M1_M1);
        grid.Cell(grid.n_x_node - grid.noghost, j ,dV_M1_P1) = grid.Cell(grid.noghost, j , dV_M1_P1);


        //BCC 3/2
        grid.Node(grid.n_x_node - (grid.noghost),j,dV_MH3_MH3) = grid.Node(grid.noghost  , j, dV_MH3_MH3);
        grid.Node(grid.n_x_node - (grid.noghost),j,dV_MH3_PH3) = grid.Node(grid.noghost  , j, dV_MH3_PH3);


        grid.Node(grid.n_x_node - (grid.noghost) +1,j,dV_MH3_MH3) = grid.Node(grid.noghost +1 , j, dV_MH3_MH3);
        grid.Node(grid.n_x_node - (grid.noghost) +1,j,dV_MH3_PH3) = grid.Node(grid.noghost +1 , j, dV_MH3_PH3);

        grid.Cell(grid.n_x_cell - (grid.noghost),j,dV_MH3_MH3) = grid.Cell(grid.noghost  , j, dV_MH3_MH3);
        grid.Cell(grid.n_x_cell - (grid.noghost),j,dV_MH3_PH3) = grid.Cell(grid.noghost  , j, dV_MH3_PH3);

       



}


for(int i = 0; i<grid.n_x_node; i++)
{       
        //-------------------------------------------top to bottom -----------//
        //SC1
        grid.Node(i,grid.noghost -1,dV_ZERO_P1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1); ///copying to first adjacent ghost node from last node physical domain
        grid.Cell(i,grid.noghost -1,dV_ZERO_P1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1);
        
        //SC2
        grid.Node(i,0,dV_ZERO_P2) = grid.Node(i,(grid.n_y_cell - grid.noghost) -2,dV_ZERO_P2);         
        grid.Cell(i,0,dV_ZERO_P2) = grid.Cell(i,(grid.n_y_cell - grid.noghost) -2,dV_ZERO_P2);

        grid.Node(i,1,dV_ZERO_P2) = grid.Node(i,(grid.n_y_cell - grid.noghost) -1,dV_ZERO_P2);
        grid.Cell(i,1,dV_ZERO_P2) = grid.Cell(i,(grid.n_y_cell - grid.noghost) -1,dV_ZERO_P2);

        //SC1/2
        grid.Cell(i,grid.noghost -1,dV_PH1_PH1) = grid.Cell(i,grid.n_y_cell - (grid.noghost +1),dV_PH1_PH1); 
        grid.Cell(i,grid.noghost -1,dV_MH1_PH1) = grid.Cell(i,grid.n_y_cell - (grid.noghost +1),dV_MH1_PH1);


        //BCC1
        grid.Node(i, grid.noghost -1, dV_P1_P1) = grid.Node(i ,grid.n_y_node - (grid.noghost +1), dV_P1_P1);
        grid.Node(i, grid.noghost -1, dV_M1_P1) = grid.Node(i ,grid.n_y_node - (grid.noghost +1), dV_M1_P1);
        
        grid.Cell(i, grid.noghost -1, dV_P1_P1) = grid.Cell(i ,grid.n_y_node - (grid.noghost +1), dV_P1_P1);
        grid.Cell(i, grid.noghost -1, dV_M1_P1) = grid.Cell(i ,grid.n_y_node - (grid.noghost +1), dV_M1_P1);



        //SCi,3
        grid.Cell(i,1,dV_PH3_PH3) = grid.Cell(i,(grid.n_y_cell-grid.noghost) - 1, dV_PH3_PH3);
        grid.Cell(i,1,dV_MH3_PH3) = grid.Cell(i,(grid.n_y_cell-grid.noghost) - 1, dV_MH3_PH3);

        grid.Cell(i,0,dV_PH3_PH3) = grid.Cell(i,(grid.n_y_cell-grid.noghost) - 2, dV_PH3_PH3);
        grid.Cell(i,0,dV_MH3_PH3) = grid.Cell(i,(grid.n_y_cell-grid.noghost) - 2, dV_MH3_PH3);

        grid.Node(i,1,dV_PH3_PH3) = grid.Node(i,(grid.n_y_cell-grid.noghost) - 1, dV_PH3_PH3);
        grid.Node(i,1,dV_MH3_PH3) = grid.Node(i,(grid.n_y_cell-grid.noghost) - 1, dV_MH3_PH3);

        //------------------------------bottom top---------------------//
        //SC1
        grid.Node(i,grid.n_y_cell - grid.noghost,dV_ZERO_M1) = grid.Node(i,grid.noghost,dV_ZERO_M1);
        grid.Cell(i,grid.n_y_cell - grid.noghost,dV_ZERO_M1) = grid.Cell(i,grid.noghost,dV_ZERO_M1);

        //Sc2
        grid.Node(i,grid.n_y_cell - (grid.noghost)  ,dV_ZERO_M2)   = grid.Node(i,grid.noghost  ,dV_ZERO_M2);
        grid.Cell(i,grid.n_y_cell - (grid.noghost)  ,dV_ZERO_M2)   = grid.Cell(i,grid.noghost  ,dV_ZERO_M2);

        grid.Node(i,(grid.n_y_cell - grid.noghost) +1,dV_ZERO_M2)   = grid.Node(i,grid.noghost+1,dV_ZERO_M2);
        grid.Cell(i,(grid.n_y_cell - grid.noghost) +1,dV_ZERO_M2)   = grid.Cell(i,grid.noghost+1,dV_ZERO_M2);

        //BCC-1/2
        grid.Node(i,grid.n_y_cell - (grid.noghost),dV_PH1_MH1) = grid.Node(i,grid.noghost,dV_PH1_MH1);
        grid.Node(i,grid.n_y_cell - (grid.noghost),dV_MH1_MH1) = grid.Node(i,grid.noghost,dV_MH1_MH1);



        //BCC 1
        grid.Node(i , grid.n_y_node - grid.noghost,dV_M1_M1) = grid.Node(i , grid.noghost, dV_M1_M1);
        grid.Node(i , grid.n_y_node - grid.noghost,dV_P1_M1) = grid.Node(i , grid.noghost, dV_P1_M1);
    
        grid.Cell(i , grid.n_y_node - grid.noghost,dV_M1_M1) = grid.Cell(i , grid.noghost, dV_M1_M1);
        grid.Cell(i , grid.n_y_node - grid.noghost,dV_P1_M1) = grid.Cell(i , grid.noghost, dV_P1_M1);



        //BCC 3/2
        grid.Node(i,grid.n_y_node - (grid.noghost),dV_MH3_MH3) = grid.Node(i,grid.noghost , dV_MH3_MH3);
        grid.Node(i,grid.n_y_node - (grid.noghost),dV_PH3_MH3) = grid.Node(i,grid.noghost , dV_PH3_MH3);


        grid.Node(i,grid.n_y_node - (grid.noghost) +1,dV_MH3_MH3) = grid.Node(i,grid.noghost +1 , dV_MH3_MH3);
        grid.Node(i,grid.n_y_node - (grid.noghost) +1,dV_PH3_MH3) = grid.Node(i,grid.noghost +1 , dV_PH3_MH3);

        grid.Cell(i,grid.n_y_cell - (grid.noghost),dV_MH3_MH3) = grid.Cell(i,grid.noghost, dV_MH3_MH3);
        grid.Cell(i,grid.n_y_cell - (grid.noghost),dV_PH3_MH3) = grid.Cell(i,grid.noghost, dV_PH3_MH3);

       
}

}









template<typename T, typename T1>
void inlet(Grid_N_C_2D<T> &grid,lbmD2Q21<T1> &lb, double ux, double uy,double rho){

double feq[21] = {0};

get_equi(feq,lb,ux,uy,rho,lb.theta0);

    for(int i = grid.noghost; i < grid.noghost+2  ;i++){
        for(int j = 0; j<grid.n_y_node;j++){

            for(int dv = 0; dv <grid.d_v; dv++){
                grid.Node(i,j,dv) = feq[dv];
            }
                        
            for(int dv = 0; dv <grid.d_v; dv++){
                grid.Cell(i,j,dv) = feq[dv];
            }
            
        }
        }
}



template<typename T, typename T1>
void outlet(Grid_N_C_2D<T> &grid,lbmD2Q21<T1> &lb, double ux, double uy,double rho){

double feq[21] = {0};

get_equi(feq,lb,ux,uy,rho,lb.theta0);

    for(int i = grid.n_x_node - grid.noghost -4; i < grid.n_x_node -grid.noghost  ;i++){
        for(int j = grid.noghost; j<grid.n_y_node - grid.noghost;j++){



            ///============for equilibrizing the outer wall also
            // for(int dv = 0; dv <grid.d_v; dv++){
            //     grid.Node(i,j,dv) = feq_Node[dv];
            // }
                        

            // for(int dv = 0; dv <grid.d_v; dv++){
            //     grid.Cell(i,j,dv) = feq_Cell[dv];
            // }
            

            //============for copying the populations in the outer wall
            for(int dv = 0; dv< 21; dv++){
                grid.Node(i,j,dv) = grid.Node(i-1, j,dv);

            }
            for(int dv = 0; dv< 21; dv++){
               grid.Cell(i,j,dv) = grid.Cell(i-1, j,dv);
            }


        }
        }
}


// for(int i = left; i<= right; i++ ){
//     for(int j = bottom; j<= top; j++){
        
//         //Node 
//         double G_num = 0; //G_numerator

//         if(marker.Node(i,j,0) == FLUID){
//             for(int dv = 1; dv<5; dv++)
//                 if(marker.Node(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID)
//                     G_num += grid.Node(i,j,dv) *1.0 + grid.Node(i,j, dv+4) * 2.0;
    
//             for(int dv = 9; dv<13; dv ++)
//                 if(marker.Cell(i + (int)lb.CxF[dv], j + (int)lb.CyF[dv],0) == SOLID)
//                     G_num += grid.Node(i,j,dv) *0.5 + grid.Node(i,j, dv +4) * 1.0 + grid.Node(i,j,dv +8) * 1.5;

//         }
//         G_calc.Node(i,j,0) = G_num/G_calc.Node(i,j,0);
        


//         //Cell
//         G_num = 0;

//         if(marker.Cell(i,j,0) == FLUID){
//         for(int dv = 1; dv<5; dv++)
//             if(marker.Cell(i + (int)lb.Cx[dv],j + (int)lb.Cy[dv], 0) == SOLID)
//                 G_num += grid.Cell(i,j,dv) *1.0+ grid.Cell(i,j, dv+4) *2.0;

//         for(int dv = 9; dv<13; dv++)
//             if(marker.Node(i + (int)lb.CxC[dv], j + (int)lb.CyC[dv],0) == SOLID)
//                 G_num += grid.Cell(i,j,dv) *0.5 + grid.Cell(i,j, dv +4) *1.0 + grid.Cell(i,j,dv +8) *1.5 ;
//         }
//         G_calc.Cell(i,j,0) = G_num/G_calc.Cell(i,j,0);
        
//     }
// }



template<typename T, typename T1>
void u_free_stream(Grid_N_C_2D<T>  &grid, lbmD2Q21<T1> &lb, double &u_free, double &rho_free){

    u_free = 0;
    int count =0;
    rho_free = 0;
    for(int i = 10; i< 11; i++){
        for(int j =  grid.noghost; j < grid.n_y_node - grid.noghost; j++){

                for(int dv = 0; dv<grid.d_v; dv++){
                 
                    u_free += grid.Node(i,j,dv) * lb.Cx[dv];
                    rho_free += grid.Node(i,j,dv);
                }

            count+=1;
        }
    }

    u_free /= count;
    rho_free /= count;

}










;
#endif