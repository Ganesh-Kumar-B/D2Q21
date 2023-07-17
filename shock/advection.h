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



//// ====================================================FIRST ORDER==================================================///////      
////====================================================FIRST ORDER==================================================///////
////====================================================FIRST ORDER==================================================///////      
////====================================================FIRST ORDER==================================================///////
////====================================================FIRST ORDER==================================================///////
// template<typename T>
// void diffuse_B_21_with_1st_order_wall_correction(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21, double u0){

// double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
// int j;u0= u0;
// double feq[21] = {0};


// ////--------------------------------------------------top wall---------------------------------------------////
// //---------------------------------------------------------------------------------------------------------///
// for(int i = grid.noghost+2; i<grid.n_x_node-grid.noghost-2;i++){
//     nx =  0; ny= -1;
//     ux = u0; uy = 0;
//     j = grid.n_y_node - 1;

//     get_equi(feq,lbD2Q21,ux,uy,rho);

//     //for the top most node

// G = (grid.Node(i,j -2, dV_ZERO_P2) + grid.Cell(i+1,j -2, dV_MH3_PH3)  + grid.Cell(i -2,j -2,dV_PH3_PH3)) /
//     (feq[dV_ZERO_M2] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);
    
//     grid.Node(i,j,dV_MH3_MH3) = feq[dV_MH3_MH3] *G;
//     grid.Node(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;
//     grid.Node(i,j,dV_PH3_MH3) = feq[dV_PH3_MH3] *G;
//     // for the cell
//     G = (grid.Cell(i,j -2, dV_ZERO_P2))/ (feq[dV_ZERO_M2]);
//     grid.Cell(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;



// //next bottom node
//     j = grid.n_y_node -2;
    
//     G = (grid.Node(i,j -1,dV_ZERO_P1) + grid.Node(i,j -2, dV_ZERO_P2)+ grid.Cell(i,j-1, dV_MH1_PH1) + grid.Cell(i-1,j-1, dV_PH1_PH1) + grid.Cell(i +1, j -2, dV_MH3_PH3) +grid.Cell(i-2,j-2,dV_PH3_PH3))/
//         (feq[dV_ZERO_M1] + feq[dV_ZERO_M2] + feq[dV_PH1_MH1] +feq[dV_MH1_MH1] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]) ; 
//         // std::cout<<G<<std::endl;
//     for(int dv = 0; dv<21; dv++){
//         if(lbD2Q21.Cx[dv]*nx + lbD2Q21.Cy[dv]*ny > 0){
//             grid.Node(i,j,dv) = feq[dv] * G;
//         }
//     }
// // for the cell
// G = (grid.Cell(i ,j-1, dV_ZERO_P1) + grid.Cell(i,j -1, dV_ZERO_P2) +  grid.Node(i+2,j-1,dV_MH3_PH3) + grid.Node(i-1,j-1,dV_PH3_PH3))/
//     (feq[dV_ZERO_M1] + feq[dV_ZERO_M2] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);

//     grid.Cell(i,j,dV_ZERO_M1) = feq[dV_ZERO_M1] *G;
//     grid.Cell(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;
//     grid.Cell(i,j,dV_PH3_MH3) = feq[dV_PH3_MH3] *G;
//     grid.Cell(i,j,dV_MH3_MH3) = feq[dV_MH3_MH3] *G;

// }



// //-------------------------------------------------------------bottom wall-------------------------------------//
// //-------------------------------------------------------------bottom wall-------------------------------------//
// for(int i = grid.noghost +2 ; i<grid.n_x_node-grid.noghost-2;i++){

//     get_equi(feq,lbD2Q21,ux,uy,rho);
//     nx = 0;  ny = 1;
//     ux = 0; uy = 0;
//     j = 0;

//     //0th node
//     G = (grid.Node(i,j+2,dV_ZERO_M2))/(feq[dV_ZERO_P2]);
//     grid.Node(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] *G;

//     //Cell
//     G = (grid.Cell(i, j +2, dV_ZERO_M2) + grid.Node(i+2, j+2, dV_MH3_MH3) + grid.Node(i-1, j +2, dV_PH3_MH3))/
//         (feq[dV_ZERO_P2] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);

//     grid.Cell(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] * G;
//     grid.Cell(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] * G;
//     grid.Cell(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] * G;


//     j =1;
//     //1st
//     G = (grid.Node(i, j +1, dV_ZERO_M1) +grid.Node(i,j +2, dV_ZERO_M2) + grid.Cell(i +1, j+1,dV_MH3_MH3) + grid.Cell(i-2,j+1, dV_PH3_MH3))/
//         (feq[dV_ZERO_P1] + feq[dV_ZERO_P2] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);

//     grid.Node(i,j,dV_ZERO_P1) = feq[dV_ZERO_P1] *G;
//     grid.Node(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] *G;
//     grid.Node(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] *G;
//     grid.Node(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] *G;

//     G = (grid.Cell(i, j +1, dV_ZERO_M1) +grid.Cell(i,j +2, dV_ZERO_M2) + grid.Node(i+1,j+1,dV_MH1_MH1) + grid.Node(i,j+1, dV_PH1_MH1) + grid.Node(i+2, j +2, dV_MH3_MH3) +grid.Node(i-1, j+2, dV_PH3_MH3))/
//         (feq[dV_ZERO_P1] + feq[dV_ZERO_P2] + feq[dV_PH1_PH1] + feq[dV_MH1_PH1] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);

//     grid.Cell(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] * G;
//     grid.Cell(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] * G;
//     grid.Cell(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] * G;
//     grid.Cell(i,j,dV_ZERO_P1) = feq[dV_ZERO_P1] * G;
//     grid.Cell(i,j,dV_PH1_PH1) = feq[dV_PH1_PH1] * G;
//     grid.Cell(i,j,dV_MH1_PH1) = feq[dV_MH1_PH1] * G;
// }
//     ////--------------------------------->>>>>>>>>>>>>>>>>>>>> first order corrections<<<<<<<<<<<<<<<<<<<<<<<--------------------------///////
 
//     const int wall_vel[21] = {0,3,4,1,2,7,8,5,6,10,9,12,11,14,13,16,15}; //when you call this at ghost node and backtrace it you can go to the origin node
//                                                                     //wall corresponding velocity

//     // const int dv_sc[] = {}


//     double del_1 = del +0.5;    // big one use it for top node andn bottom cell
//     double del_2 = del ;    // small one use it for top cell and bottom node    

//     //////////===============for top velocities
//     for(int i = 0; i < grid.n_x_node - grid.noghost -2; i++){

//         for(int j =grid.n_y_node-grid.noghost; j<grid.n_y_node-grid.noghost +2;j++){

//             for(int dv = 1; dv< 9; dv++){   ///runs only for the SC
//                 grid.Node(i,j,dv) = -(1.0/2*(del_1))*(grid.Node(i + (int) lbD2Q21.Cx[wall_vel[dv]],j + (int) lbD2Q21.cy[wall_vel[dv]], oppdV[wall_corresponding_velocity[dv]]) - grid.Node(i,j,dv)) 
//                                        +(grid.Node(i + (int) lbD2Q21.Cx[wall_vel[dv]],j + (int) lbD2Q21.cy[wall_vel[dv]], oppdV[wall_corresponding_velocity[dv]]));
            
//                 grid.Cell(i,j,dv) = -(1.0/2*(del_2))*(grid.Cell(i + (int) lbD2Q21.Cx[wall_vel[dv]],j + (int) lbD2Q21.cy[wall_vel[dv]], oppdV[wall_corresponding_velocity[dv]]) - grid.Cell(i,j,dv)) 
//                                        +(grid.Cell(i + (int) lbD2Q21.Cx[wall_vel[dv]],j + (int) lbD2Q21.cy[wall_vel[dv]], oppdV[wall_corresponding_velocity[dv]]));
            
            
//             }

//             for(int dv = 9; dv< 21; dv++){ //runs for the bcc    
//                 grid.Node(i,j,dv) = -(1.0/del_1) *(grid.Cell(i +(int)))
            
//             }

//         }
//     }
// }








// template<typename T>
// void diffuse_B_21_out(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21,double u0){ ///////with respect to te velocity


// double G_top_node[2] = {},G_top_cell[2] = {},G_bot_node[2] = {}, G_bot_cell[2] = {};        ///0,1,2,3,4;---while call subract n_Y_node
// double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
// int j;u0= u0;
// double feq[21] = {0};
// int dv_top_node0[6] ={4,8,12,11,16,15},dv_top_node1[3] = {8,15,16},         
//     dv_top_cell0[4] ={4,8,15,16}      ,dv_top_cell1[1] = {8};
// int dv_bot_cell0[6] ={2,6,10,14,9,13} ,dv_bot_cell1[3] = {6,13,14},
//     dv_bot_node0[4] = {2,6,13,14}     ,dv_bot_node1[1] = {6};

// ux = u0; uy = 0;                                                                   //top velocity declaration
// get_equi(feq,lbD2Q21,ux,uy,rho);                                                //defining the equilibrium for the top




// //G declarations

// G_top_node[0] = 1.0/(feq[dV_ZERO_M1] + feq[dV_ZERO_M2] + feq[dV_PH1_MH1] +feq[dV_MH1_MH1] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);
// G_top_node[1] = 1.0/(feq[dV_ZERO_M2] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);

// G_top_cell[0] = 1.0/(feq[dV_ZERO_M1] + feq[dV_ZERO_M2] + feq[dV_PH3_MH3] + feq[dV_MH3_MH3]);
// G_top_cell[1] = 1.0/(feq[dV_ZERO_M2]);

// ux = u0; uy = 0;                    ///bottom velocity declarations
// get_equi(feq,lbD2Q21,ux,uy,rho)  ;                                                          //define the equilibrium for the bottom one


// G_bot_node[0] = 1/(feq[dV_ZERO_P1] + feq[dV_ZERO_P2] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);
// G_bot_node[1] = 1/(feq[dV_ZERO_P2]);

// G_bot_cell[0] = 1/(feq[dV_ZERO_P1] + feq[dV_ZERO_P2] + feq[dV_PH1_PH1] +feq[dV_MH1_PH1] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);
// G_bot_cell[1] = 1/(feq[dV_ZERO_P2] + feq[dV_PH3_PH3] + feq[dV_MH3_PH3]);



// for(int i = grid.noghost+2; i<grid.n_x_node- grid.noghost - 2; i++){
//         j = grid.n_x_node - grid.noghost;
//         for(int dv = 0; dv< grid.d_v ; dv++){
//             grid.Node(i,j,dv) = 0.0;
//             grid.Cell(i,j,dv) =0.0;
            
//         }
//         j = grid.n_x_node - grid.noghost +1;
//         for(int dv = 0; dv< grid.d_v ; dv++){
//             grid.Node(i,j,dv) = 0.0;
//             grid.Cell(i,j,dv) =0.0;
            
//         }
//         j = 0;
//           for(int dv = 0; dv< grid.d_v ; dv++){
//             grid.Node(i,j,dv) = 0.0;
//             grid.Cell(i,j,dv) =0.0;
            
//         }
//         j = 1;
//                 j = 0;
//           for(int dv = 0; dv< grid.d_v ; dv++){
//             grid.Node(i,j,dv) = 0.0;
//             grid.Cell(i,j,dv) =0.0;
            
//         }


// }



// // ================================top wall===============================================//
// // ///node
// // ////dv = dvzerop1
        
//         for(int i = grid.noghost+2; i<grid.n_x_node-grid.noghost-2;i++){
//         j = grid.n_y_node -grid.noghost -1;                                         //-------->>>>>>for the boundary node
//         //=======node
//         //sc1,3 
//         for(int dvc = 0; dvc<6; dvc++){
//             grid.Node(i , j+1, dv_top_node0[dvc]) += feq[dv_top_node0[dvc]]*grid.Node(i,j  ,dV_ZERO_P1) *G_top_node[0];
//             grid.Node(i , j+1, dv_top_node0[dvc]) += feq[dv_top_node0[dvc]]*grid.Node(i,j-1,dV_ZERO_P2) *G_top_node[0];
//         }
        
//         for(int dvc = 0; dvc<3; dvc++){
//             grid.Node(i , j+2, dv_top_node1[dvc]) += feq[dv_top_node1[dvc]]*grid.Node(i,j,dV_ZERO_P2) *G_top_node[1];
//         }

//         //=========cell
//         //sc1
//         j = grid.n_y_node - grid.noghost -1;
//         for(int dvc = 0; dvc<4; dvc++){
//             grid.Cell(i , j+1, dv_top_cell0[dvc]) += feq[dv_top_cell0[dvc]] *grid.Cell(i,j,dV_ZERO_P1) * G_top_cell[0];
//             grid.Cell(i , j+1, dv_top_cell0[dvc]) += feq[dv_top_cell0[dvc]]*grid.Cell(i,j-1,dV_ZERO_P2) *G_top_cell[0];
//         }
//         for(int dvc = 0; dvc<1; dvc++){
//             grid.Cell(i , j+2, dv_top_cell1[dvc]) += feq[dv_top_cell1[dvc]]*grid.Cell(i,j,dV_ZERO_P2) *G_top_cell[1];
//         }
//         }


//         //dv = ph2ph2
//         for(int i = grid.noghost; i<grid.n_x_node-grid.noghost-2-2;i++){   ///should run till -5 only
//             j = grid.n_y_node -grid.noghost -1;   
            
//             for(int dvc = 0; dvc<3; dvc++){
//                 grid.Node(i +2, j+2, dv_top_node1[dvc]) += feq[dv_top_node1[dvc]] *grid.Cell(i,j,dV_PH3_PH3) *G_top_node[1];
//             }

//             for(int dvc = 0; dvc<6 ;dvc++){
//                 grid.Node(i+2,j+2, dv_top_node0[dvc]) += feq[dv_top_node0[dvc]] * grid.Cell(i,j,dV_PH3_PH3) *G_top_node[0];
//             }
        
//         }

//         for(int i = grid.noghost +1 ; i<grid.n_x_node-grid.noghost-3;i++){
//             for(int dvc = 0; dvc<4; dvc++){
//             grid.Cell(i+1, j+1, dv_top_cell0[dvc]) += feq[dv_top_cell0[dvc]] *grid.Node(i,j,dV_PH3_PH3) * G_top_cell[0];
//             }
//         } 
        
        
//         j = grid.n_y_node -grid.noghost -2;
//         for(int i = grid.noghost +2; i<grid.n_x_node-grid.noghost-3;i++){
//             for(int dvc = 0; dvc<6 ;dvc++){
//                 grid.Node(i+1,j+1, dv_top_node0[dvc]) += feq[dv_top_node0[dvc]] * grid.Cell(i,j,dV_PH1_PH1) *G_top_node[0];
//             }
//         }

//         //dv = mh3ph3

//         j = grid.n_y_node -grid.noghost -1;
//         for(int i = grid.noghost+3; i<grid.n_x_node-grid.noghost-1;i++){
//              for(int dvc = 0; dvc<3; dvc++){
//                 grid.Node(i -1, j+2, dv_top_node1[dvc]) += feq[dv_top_node1[dvc]] *grid.Cell(i,j,dV_MH3_PH3) *G_top_node[1];
//              }
//         }

//         j = grid.n_y_node -grid.noghost -2;
//         for(int i = grid.noghost+3; i<grid.n_x_node-grid.noghost-1;i++){
//             for(int dvc = 0; dvc<6; dvc++){
//                 grid.Node(i -1, j+2, dv_top_node0[dvc]) += feq[dv_top_node0[dvc]] *grid.Cell(i,j,dV_MH3_PH3) *G_top_node[0];
//              }
//         }

//         j = grid.n_y_node -grid.noghost -1;
//         for(int i = grid.noghost+4; i<grid.n_x_node-grid.noghost-1;i++){
//             for(int dvc = 0; dvc<4; dvc++){
//                 grid.Cell(i +1, j -2,dv_top_cell0[dvc]) += feq[dv_top_cell0[dvc]] *grid.Node(i,j,dV_MH3_PH3) * G_top_cell[0];
//             }
//         }

//         for(int i = grid.noghost +2; i<grid.n_x_node-grid.noghost-4;i++){
//             for(int dvc = 0; dvc<6; dvc++){
//                 grid.Node(i , j+1, dv_top_node0[dvc]) += feq[dv_top_node0[dvc]] *grid.Cell(i,j,dV_MH1_PH1) *G_top_node[0];
//              }     
//         }
// // 

// //================================================================bottom wall======================================//
// //================================================================bottom wall======================================//

//      for(int i = grid.noghost+2; i<grid.n_x_node-grid.noghost-2;i++){
//         j = grid.noghost;       
    
//         //=======Cell
//         //sc1,3 
//         for(int dvc = 0; dvc<6; dvc++){
//             grid.Cell(i , j-1, dv_bot_cell0[dvc]) += feq[dv_bot_cell0[dvc]]*grid.Cell(i,j  ,dV_ZERO_M1) *G_bot_cell[0];
//             grid.Cell(i , j-1, dv_bot_cell0[dvc]) += feq[dv_bot_cell0[dvc]]*grid.Cell(i,j+1,dV_ZERO_M2) *G_bot_cell[0];
//         }
        
//         for(int dvc = 0; dvc<3; dvc++){
//             grid.Cell(i , j-2, dv_bot_cell1[dvc]) += feq[dv_bot_cell1[dvc]]*grid.Cell(i,j,dV_ZERO_M2) *G_bot_cell[1];
//         }

//         //=========cell
//         //sc1
//         j =grid.noghost;
//         for(int dvc = 0; dvc<4; dvc++){
//             grid.Node(i , j-1, dv_bot_node0[dvc]) += feq[dv_bot_node0[dvc]] *grid.Node(i,j  ,dV_ZERO_M1) * G_bot_node[0];
//             grid.Node(i , j-1, dv_bot_node0[dvc]) += feq[dv_bot_node0[dvc]]* grid.Node(i,j+1,dV_ZERO_M2) *G_bot_node[0];
//         }
//         for(int dvc = 0; dvc<1; dvc++){
//             grid.Node(i , j-2, dv_bot_node1[dvc]) += feq[dv_bot_node1[dvc]]*grid.Node(i,j,dV_ZERO_M2) *G_bot_node[1];
//         }
//         }

//         // dv = ph2mh2
//         j = grid.noghost;
//         for(int i = grid.noghost+1; i<grid.n_x_node-grid.noghost-2-1;i++){
//             j = grid.noghost;
//             for(int dvc=0; dvc<3;dvc++){
//                     grid.Cell(i+1,j-2,dv_bot_cell1[dvc]) += feq[dv_bot_cell1[dvc]]*grid.Node(i,j,dV_PH3_MH3) * G_bot_node[1];

//             }
//         }
//         j = grid.noghost;
//         for(int i = grid.noghost; i<grid.n_x_node - grid.noghost-4; i++){
            
//                 for(int dvc=0; dvc<4;dvc++){
//                     grid.Node(i+2,j-1, dv_bot_node0[dvc]) += feq[dv_bot_node0[dvc]] * grid.Cell(i,j,dV_PH3_PH3) * G_bot_node[0] ;
//                 }
            
//         }

//         j = grid.noghost +1;
//         for(int i = grid.noghost+1; i < grid.n_x_node - grid.noghost-3; i++){

//             for(int dvc = 0; dvc<6; dvc++ ){

//                 grid.Cell(i+1, j-2, dv_top_cell0[dvc]) += feq[dv_top_cell0[dvc]]*grid.Node(i,j,dV_PH3_MH3)*G_bot_cell[0];

//             }
//         }

//         //dv = ph1mh1
//         j = grid.noghost;
//         for(int i= grid.noghost+2; i<grid.n_x_node-2; i++){
            
//                 for(int dvc = 0; dvc<6;dvc++){

//                     grid.Cell(i,j-1, dv_bot_node0[dvc]) += feq[dv_bot_node0[dvc]] *  grid.Node(i,j,dV_PH1_PH1) * G_bot_node[0]; 
//                 }
            
//         }

//         //dvmh3mh3
//         j = grid.noghost;
//         for(int i = grid.noghost +4 ; i<grid.n_x_node - grid.noghost; i++){
//             for(int dvc = 0; dvc<3; dvc++){
//                 grid.Cell(i -2, j -2, dv_bot_cell1[dvc]) += feq[dv_bot_cell1[dvc]] * grid.Node(i,j,dV_MH3_MH3) * G_bot_cell[1];
            
//             }
//         }

//         for(int i = grid.noghost +3; i< grid.n_x_node - grid.noghost -1;i++){
//             for(int dvc = 0; dvc<4; dvc++){
//                 grid.Node(i -1 , j-1, dv_bot_node0[dvc]) += feq[dv_bot_node0[dvc]] * grid.Cell(i,j,dV_MH3_MH3) * G_bot_node[0];
//             }
//         }

//         j = grid.noghost +1;
//         for(int i = grid.noghost +4; i< grid.n_x_node - grid.noghost ; i++){
//             for(int dvc =0; dvc< 6; dvc++){
//                 grid.Cell(i -2 , j -2, dv_bot_cell0[dvc]) += feq[dv_bot_cell0[dvc]] * grid.Node(i,j, dV_MH3_MH3) * G_bot_cell[0];
//             }
//         }

//         //dv mh1mh1
//         j = grid.noghost;
//         for(int i = grid.noghost +3 ; i< grid.n_x_node - grid.noghost -1; i++){

//             for(int dvc = 0; dvc<6;dvc++){
//                 grid.Cell(i,j,dv_bot_cell0[dvc]) += feq[dv_bot_cell0[dvc]] *grid.Node(i,j,dV_MH1_MH1) * G_bot_cell[0];
//             }
//         }

// }


// template<typename T>
// void diffuse_B_21_t(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21, double u0){  //diffuse considering at different times different populatin get hit in the wall



// /////considering the top first ghost node is the wall
// for(int i = 0; i<grid.n_x_node - grid.noghost; i++){

    
//     //at time t = 1/3 for populatios landing on the first node
// int j = grid.n_y_node - grid.noghost;

//     G = ()

// }




// }





//////------------with the norm multiplied/// but the norm given in the krithivasan's paper is wrong check from manjusha thesis
// template<typename T>
// void diffuse_B_21(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lbD2Q21, double u0){

// double nx; double ny,ux ,uy = 0,rho = 1.0,G = 0; //top
// int j;u0= u0;
// double feq[21] = {0};


// ////--------------------------------------------------top wall---------------------------------------------////
// //---------------------------------------------------------------------------------------------------------///
// for(int i = grid.noghost+2; i<grid.n_x_node-grid.noghost-2;i++){
//     nx =  0; ny= -1;
//     ux = u0; uy = 0;
//     j = grid.n_y_node - 1;

//     get_equi(feq,lbD2Q21,ux,uy,rho);

//     //for the top most node

// G = (grid.Node(i,j -2, dV_ZERO_P2)*2.0 + grid.Cell(i+1,j -2, dV_MH3_PH3) *1.5 + grid.Cell(i -2,j -2,dV_PH3_PH3)*1.5) /
//     (feq[dV_ZERO_M2]*2.0 + feq[dV_PH3_MH3]*1.5 + feq[dV_MH3_MH3]*1.5);
    
//     grid.Node(i,j,dV_MH3_MH3) = feq[dV_MH3_MH3] *G;
//     grid.Node(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;
//     grid.Node(i,j,dV_PH3_MH3) = feq[dV_PH3_MH3] *G;
//     // for the cell
//     G = (grid.Cell(i,j -2, dV_ZERO_P2))/ (feq[dV_ZERO_M2]);
//     grid.Cell(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;

// //next bottom node
//     j = grid.n_y_node -2;
    
//     G = (grid.Node(i,j -1,dV_ZERO_P1) + grid.Node(i,j -2, dV_ZERO_P2)* 2.0+ grid.Cell(i,j-1, dV_MH1_PH1) *0.5+ grid.Cell(i-1,j-1, dV_PH1_PH1)*0.5 + grid.Cell(i +1, j -2, dV_MH3_PH3)*1.5 +grid.Cell(i-2,j-2,dV_PH3_PH3)*1.5)/
//         (feq[dV_ZERO_M1] + feq[dV_ZERO_M2]*2.0 + feq[dV_PH1_MH1]*0.5 +feq[dV_MH1_MH1]*0.5 + feq[dV_PH3_MH3]*1.5 + feq[dV_MH3_MH3]*1.5) ; 
//         // std::cout<<G<<std::endl;
//     for(int dv = 0; dv<21; dv++){
//         if(lbD2Q21.Cx[dv]*nx + lbD2Q21.Cy[dv]*ny > 0){

//             grid.Node(i,j,dv) = feq[dv] * G;

//         }

//     }
// // for the cell
// G = (grid.Cell(i ,j-1, dV_ZERO_P1) + grid.Cell(i,j -1, dV_ZERO_P2)*2.0 +  grid.Node(i+2,j-1,dV_MH3_PH3) *1.5+ grid.Node(i-1,j-1,dV_PH3_PH3)*1.5)/
//     (feq[dV_ZERO_M1]*1.0 + feq[dV_ZERO_M2]*2.0 + feq[dV_PH3_MH3]*1.5 + feq[dV_MH3_MH3]*1.5);

//     grid.Cell(i,j,dV_ZERO_M1) = feq[dV_ZERO_M1] *G;
//     grid.Cell(i,j,dV_ZERO_M2) = feq[dV_ZERO_M2] *G;
//     grid.Cell(i,j,dV_PH3_MH3) = feq[dV_PH3_MH3] *G;
//     grid.Cell(i,j,dV_MH3_MH3) = feq[dV_MH3_MH3] *G;

// }

// //-------------------------------------------------------------bottom wall-------------------------------------//
// //-------------------------------------------------------------bottom wall-------------------------------------//
// for(int i = grid.noghost +2 ; i<grid.n_x_node-grid.noghost-2;i++){

//     get_equi(feq,lbD2Q21,ux,uy,rho);
//     nx = 0;  ny = 1;
//     ux = 0; uy = 0;
//     j = 0;

//     //0th node
//     G = (grid.Node(i,j+2,dV_ZERO_M2))/(feq[dV_ZERO_P2]);
//     grid.Node(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] *G;

//     //Cell
//     G = (grid.Cell(i, j +2, dV_ZERO_M2)*2.0 + grid.Node(i+2, j+2, dV_MH3_MH3)*1.5 + grid.Node(i-1, j +2, dV_PH3_MH3)*1.5)/
//         (feq[dV_ZERO_P2]*2.0 + feq[dV_PH3_PH3]*1.5 + feq[dV_MH3_PH3]*1.5);

//     grid.Cell(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] * G;
//     grid.Cell(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] * G;
//     grid.Cell(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] * G;


//     j =1;
//     //1st
//     G = (grid.Node(i, j +1, dV_ZERO_M1) +grid.Node(i,j +2, dV_ZERO_M2)*2.0 + grid.Cell(i +1, j+1,dV_MH3_MH3) *1.5+ grid.Cell(i-2,j+1, dV_PH3_MH3)*1.5)/
//         (feq[dV_ZERO_P1] + feq[dV_ZERO_P2]*2.0 + feq[dV_PH3_PH3]*1.5 + feq[dV_MH3_PH3]*1.5);

//     grid.Node(i,j,dV_ZERO_P1) = feq[dV_ZERO_P1] *G;
//     grid.Node(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] *G;
//     grid.Node(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] *G;
//     grid.Node(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] *G;

//     G = (grid.Cell(i, j +1, dV_ZERO_M1) +grid.Cell(i,j +2, dV_ZERO_M2)*2.0 + grid.Node(i+1,j+1,dV_MH1_MH1)*0.5 + grid.Node(i,j+1, dV_PH1_MH1)*0.5 + grid.Node(i+2, j +2, dV_MH3_MH3)*1.5 +grid.Node(i-1, j+2, dV_PH3_MH3)*1.5)/
//         (feq[dV_ZERO_P1] + feq[dV_ZERO_P2]*2.0 + feq[dV_PH1_PH1]*0.5 + feq[dV_MH1_PH1] *0.5+ feq[dV_PH3_PH3]*1.5 + feq[dV_MH3_PH3]*1.5);

//     grid.Cell(i,j,dV_ZERO_P2) = feq[dV_ZERO_P2] * G;
//     grid.Cell(i,j,dV_PH3_PH3) = feq[dV_PH3_PH3] * G;
//     grid.Cell(i,j,dV_MH3_PH3) = feq[dV_MH3_PH3] * G;
//     grid.Cell(i,j,dV_ZERO_P1) = feq[dV_ZERO_P1] * G;
//     grid.Cell(i,j,dV_PH1_PH1) = feq[dV_PH1_PH1] * G;
//     grid.Cell(i,j,dV_MH1_PH1) = feq[dV_MH1_PH1] * G;
// }

// }

template<typename T>
void slip_wall_bb(Grid_N_C_2D<T> &grid,lbmD2Q21<T> d2q21 ){

int j;

/*! \brief  velocities*/
int wcv[21]= {0,3,4,1,2,7,8,5,6,12,11,10,9,16,15,14,13} ; 
 for(int i= grid.noghost +2 ; i< grid.n_x_node - grid.noghost -1 -2 ; i++){


    //====================================TOP===============================//
    j = grid.n_y_node -grid.noghost -1;
    //for sc1 and sc2
    for(int dv = 0; dv < 7; dv++){  //only go tilll the sc2+ve velocity

        grid.Cell(i , j +1, oppdV[dv]) = grid.Cell(i  ,j  ,dv);   //for the last cell
        grid.Cell(i , j +2, oppdV[dv]) = grid.Cell(i  ,j-1,dv);   //fot the second last cell only the sc2 velocity 

        grid.Node(i , j +1, oppdV[dv]) = grid.Node(i  ,j  ,dv);
        grid.Node(i , j +2, oppdV[dv]) = grid.Node(i  ,j-1,dv);
    }
    }
    



                           
 for(int i= grid.noghost +2 ; i< grid.n_x_node - grid.noghost -1 -2 ; i++){
     //for the BCC
    for(int dv = 9; dv < 11 ; dv ++){

            grid.Node(i + (int)ceil(d2q21.Cx[dv]), j + (int) ceil(d2q21.Cy[dv]),  wcv[dv])      = grid.Cell(i,j,dv  );
            grid.Cell(i , j +1  ,wcv[dv+4])   = grid.Cell(i ,j  ,dv+4);  
            grid.Node(i , j +2  ,wcv[dv+4])   = grid.Node(i ,j  ,dv+4);
        }
            grid.Node(i +1  , j +2 ,dV_PH3_MH3)    = grid.Cell(i   ,j-1,   dV_PH3_PH3);
            grid.Node(i     , j +2 ,dV_MH3_MH3)    = grid.Cell(i   ,j-1,   dV_MH3_PH3) ;
    }



    /////////////////====================bottom =====================================================//////////
    j = grid.noghost;
 for(int i= grid.noghost +2 ; i< grid.n_x_node - grid.noghost -1 -2 ; i++){

        for(int dv = 0; dv < 7; dv++){  //only go tilll the sc2+ve velocity

        grid.Node(i , j -1, oppdV[dv]) = grid.Node(i  ,j  ,dv);   //for the last cell
        grid.Node(i , j -2, oppdV[dv]) = grid.Node(i  ,j+1,dv);   //fot the second last cell only the sc2 velocity 

        grid.Cell(i , j -1, oppdV[dv]) = grid.Cell(i  ,j  ,dv);
        grid.Cell(i , j -2, oppdV[dv]) = grid.Cell(i  ,j+1,dv);
    
    }}
    

    j = grid.noghost;
 for(int i= grid.noghost +2 ; i< grid.n_x_node - grid.noghost -1 -2 ; i++){
          for(int dv = 11; dv < 13 ; dv ++){

            grid.Cell(i + (int)floor(d2q21.Cx[dv]), j + (int) floor(d2q21.Cy[dv]),  wcv[dv])      = grid.Node(i,j,dv  );
            grid.Node(i , j -1  ,wcv[dv+4])   = grid.Node(i ,j  ,dv+4);  
            grid.Cell(i , j -2  ,wcv[dv+4])   = grid.Cell(i ,j  ,dv+4);
        }

            grid.Cell(i  ,j -2 ,dV_PH3_PH3)      = grid.Node(i   ,j+1,   dV_PH3_MH3);        
            grid.Cell(i-1,j    ,dV_MH3_PH3)      = grid.Node(i   ,j+1,   dV_MH3_MH3);
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

// template<typename T>
// void periodic(Grid_N_C_2D<T> &grid){

  
//     // //---------------------------------top bottom---------------------------------------//can it be applied separately?????
//     // for(int i = 0; i<grid.n_x_node; i++){
//     //     ///-------top to bottom----///
//     //     //SC1
//     //     grid.Node(i,grid.noghost -1,dV_ZERO_P1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1); ///copying to first adjacent ghost node from last node physical domain
//     //     grid.Cell(i,grid.noghost -1,dV_ZERO_P1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1); 

//     //     //SC3
//     //     grid.Node(i,0,dV_ZERO_P3) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);         ///*****not generalized
//     //     grid.Cell(i,0,dV_ZERO_P3) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);

//     //     grid.Node(i,1,dV_ZERO_P3) = grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);
//     //     grid.Cell(i,1,dV_ZERO_P3) = grid.Cell(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);

//     //     grid.Node(i,2,dV_ZERO_P3) = grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);
//     //     grid.Cell(i,2,dV_ZERO_P3) = grid.Cell(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);

//     //     //SC-1/2
//     //     grid.Cell(i,grid.noghost -1,dV_PH1_PH1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_PH1_PH1); //top cell to bottom
//     //     grid.Cell(i,grid.noghost -1,dV_MH1_PH1) = grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_MH1_PH1);

//     //     ///-------bottom to top------///
//     //     //SC -1
//     //     grid.Node(i,grid.n_y_node - grid.noghost,dV_ZERO_M1) = grid.Node(i, grid.noghost,dV_ZERO_M1);
//     //     grid.Cell(i,grid.n_y_node - grid.noghost,dV_ZERO_M1) = grid.Cell(i, grid.noghost,dV_ZERO_M1);

//     //     //SC3
//     //     grid.Node(i,grid.n_y_node - (grid.noghost),dV_ZERO_M3)   = grid.Node(i,grid.noghost,dV_ZERO_M3);
//     //     grid.Cell(i,grid.n_y_node - (grid.noghost),dV_ZERO_M3)   = grid.Cell(i,grid.noghost,dV_ZERO_M3);

//     //     grid.Node(i,grid.n_y_node - (grid.noghost-1),dV_ZERO_M3) = grid.Node(i,grid.noghost+1,dV_ZERO_M3);
//     //     grid.Cell(i,grid.n_y_node - (grid.noghost-1),dV_ZERO_M3) = grid.Cell(i,grid.noghost+1,dV_ZERO_M3);
        
//     //     grid.Node(i,grid.n_y_node - (grid.noghost-2),dV_ZERO_M3) = grid.Node(i,grid.noghost+2,dV_ZERO_M3);
//     //     grid.Cell(i,grid.n_y_node - (grid.noghost-2),dV_ZERO_M3) = grid.Cell(i,grid.noghost+2,dV_ZERO_M3);


//     //     //SC-1/2
//     //     grid.Node(i,grid.n_y_node - (grid.noghost),dV_PH1_MH1) = grid.Node(i,grid.noghost,dV_PH1_MH1);
//     //     grid.Node(i,grid.n_y_node - (grid.noghost),dV_MH1_MH1) = grid.Node(i,grid.noghost,dV_MH1_MH1);

//     // }

//     //--------------------------------------------------left right-----------------------------------------------//

//     for(int j =0; j<grid.n_y_node  ;j++){
//         //--------right to left-------//
//         grid.Node(grid.noghost -1,j,dV_P1_ZERO) = grid.Node( grid.n_x_node - (grid.noghost +1),j,dV_P1_ZERO); ///copying to first adjacent ghost node from last node physical domain
//         grid.Cell(grid.noghost -1,j,dV_P1_ZERO) = grid.Cell( grid.n_x_node - (grid.noghost +1),j,dV_P1_ZERO);

//         //SC3
//         grid.Node(0,j,dV_P3_ZERO) = grid.Node(grid.n_x_cell - (grid.noghost +3),j,dV_P3_ZERO);         
//         grid.Cell(0,j,dV_P3_ZERO) = grid.Cell(grid.n_x_cell - (grid.noghost +3),j,dV_P3_ZERO);

//         grid.Node(1,j,dV_P3_ZERO) = grid.Node(grid.n_x_cell - (grid.noghost +2),j,dV_P3_ZERO);
//         grid.Cell(1,j,dV_P3_ZERO) = grid.Cell(grid.n_x_cell - (grid.noghost +2),j,dV_P3_ZERO);

//         grid.Node(2,j,dV_P3_ZERO) = grid.Node(grid.n_x_cell - (grid.noghost +1),j,dV_P3_ZERO);
//         grid.Cell(2,j,dV_P3_ZERO) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_P3_ZERO);

//         //BCC-1/2
//         grid.Cell(grid.noghost -1,j,dV_PH1_PH1) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_PH1_PH1); 
//         grid.Cell(grid.noghost -1,j,dV_PH1_MH1) = grid.Cell(grid.n_x_cell - (grid.noghost +1),j,dV_PH1_MH1);

//         ///-------left to right-----///
//         //SC -1
//         grid.Node(grid.n_x_cell - grid.noghost,j,dV_M1_ZERO) = grid.Node(grid.noghost,j,dV_M1_ZERO);
//         grid.Cell(grid.n_x_cell - grid.noghost,j,dV_M1_ZERO) = grid.Cell(grid.noghost,j,dV_M1_ZERO);

//         //SC3
//         grid.Node(grid.n_x_cell - (grid.noghost)  ,j,dV_M3_ZERO)   = grid.Node(grid.noghost  ,j,dV_M3_ZERO);
//         grid.Cell(grid.n_x_cell - (grid.noghost)  ,j,dV_M3_ZERO)   = grid.Cell(grid.noghost  ,j,dV_M3_ZERO);

//         grid.Node(grid.n_x_cell - (grid.noghost-1),j,dV_M3_ZERO)   = grid.Node(grid.noghost+1,j,dV_M3_ZERO);
//         grid.Cell(grid.n_x_cell - (grid.noghost-1),j,dV_M3_ZERO)   = grid.Cell(grid.noghost+1,j,dV_M3_ZERO);

//         grid.Node(grid.n_x_cell - (grid.noghost-2),j,dV_M3_ZERO)   = grid.Node(grid.noghost+2,j,dV_M3_ZERO);
//         grid.Cell(grid.n_x_cell - (grid.noghost-2),j,dV_M3_ZERO)   = grid.Cell(grid.noghost+2,j,dV_M3_ZERO);


//         //BCC-1/2
//         grid.Node(grid.n_x_cell - (grid.noghost),j,dV_MH1_PH1) = grid.Node(grid.noghost,j,dV_MH1_PH1);
//         grid.Node(grid.n_x_cell - (grid.noghost),j,dV_MH1_MH1) = grid.Node(grid.noghost,j,dV_MH1_MH1);

//     }

// }


//as


// template<typename T>
// void bounce_back_prep(Grid_N_C_2D<T> &grid){
// //________________________________________________bounceback_________________________________________________________//
   
//     for(int i =1; i<grid.n_x_node-1 ; i++){
        
//         //------SC1
//         //top________(grid.ny_node - grid.noghost) --->> 1st ghost node
//         grid.Node(i, grid.n_y_node - (grid.noghost),dV_ZERO_P1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1);
//         //taking the last cell as the boundary
//         grid.Cell(i, grid.n_y_node - (grid.noghost),dV_ZERO_M1) = grid.Cell(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P1);

//         //bottom
//         grid.Node(i,grid.noghost -1,dV_ZERO_M1) = grid.Node(i,grid.noghost,dV_ZERO_M1);
//         grid.Cell(i,grid.noghost -1,dV_ZERO_M1) = grid.Cell(i,grid.noghost,dV_ZERO_M1);


//         //----BCC1/2
//         //top
//         grid.Node(i, grid.n_y_node - (grid.noghost),dV_PH1_PH1) = grid.Cell(i-1, grid.n_y_node - (grid.noghost+1),dV_PH1_PH1);
//         grid.Node(i, grid.n_y_node - (grid.noghost),dV_MH1_PH1) = grid.Cell(i  , grid.n_y_node - (grid.noghost+1),dV_MH1_PH1);
//         //bottom
//         grid.Cell(i,grid.noghost -1,dV_MH1_MH1) = grid.Node(i+1,grid.noghost,dV_MH1_MH1);
//         grid.Cell(i,grid.noghost -1,dV_PH1_MH1) = grid.Node(i,grid.noghost,dV_PH1_MH1);

//         //---SC3
//         //top
//         grid.Node(i,  grid.n_y_node - (grid.noghost)  ,dV_ZERO_P3) =  grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);
//         grid.Node(i,  grid.n_y_node - (grid.noghost-1),dV_ZERO_P3) =  grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);
//         grid.Node(i,  grid.n_y_node - (grid.noghost-2),dV_ZERO_P3) =  grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);
                   
//         grid.Cell(i,  grid.n_y_cell - (grid.noghost)  ,dV_ZERO_P3) =  grid.Cell(i, grid.n_y_cell - (grid.noghost +3),dV_ZERO_P3);
//         grid.Cell(i,  grid.n_y_cell - (grid.noghost-1),dV_ZERO_P3) =  grid.Cell(i, grid.n_y_cell - (grid.noghost +2),dV_ZERO_P3);
//         grid.Cell(i,  grid.n_y_cell - (grid.noghost-2),dV_ZERO_P3) =  grid.Cell(i, grid.n_y_cell - (grid.noghost +1),dV_ZERO_P3);

//         //bottom
//         grid.Node(i, grid.noghost - 3,dV_ZERO_M3) =  grid.Node(i, grid.noghost   ,dV_ZERO_M3);
//         grid.Node(i, grid.noghost - 2,dV_ZERO_M3) =  grid.Node(i, grid.noghost +1,dV_ZERO_M3);
//         grid.Node(i, grid.noghost - 1,dV_ZERO_M3) =  grid.Node(i, grid.noghost +2,dV_ZERO_M3);

//         grid.Cell(i, grid.noghost - 3,dV_ZERO_M3) =  grid.Cell(i, grid.noghost   ,dV_ZERO_M3);
//         grid.Cell(i, grid.noghost - 2,dV_ZERO_M3) =  grid.Cell(i, grid.noghost +1,dV_ZERO_M3);
//         grid.Cell(i, grid.noghost - 1,dV_ZERO_M3) =  grid.Cell(i, grid.noghost +2,dV_ZERO_M3);        
 
//     }
    

// }



// template<typename T>
// void bounce_back_wall(Grid_N_C_2D<T> &grid){




//     for(int i =1; i<grid.n_x_node-1 ; i++){

//     //SC1
//     //top
//     grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_M1) = grid.Node(i, grid.n_y_node - (grid.noghost),dV_ZERO_P1);
//     grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_M1) = grid.Cell(i, grid.n_y_node - (grid.noghost),dV_ZERO_P1);
//     //bottom
//     grid.Node(i,grid.noghost,dV_ZERO_P1) = grid.Node(i,grid.noghost -1,dV_ZERO_M1);
//     grid.Cell(i,grid.noghost,dV_ZERO_P1) = grid.Cell(i,grid.noghost -1,dV_ZERO_M1);
//     //BCC1/2
//     //top
//     grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_MH1_MH1) = grid.Node(i+1, grid.n_y_node - (grid.noghost),dV_PH1_PH1);
//     grid.Cell(i, grid.n_y_node - (grid.noghost +1),dV_PH1_MH1) = grid.Node(i  , grid.n_y_node - (grid.noghost),dV_MH1_PH1);
//     //bottom
//     grid.Node(i,grid.noghost,dV_PH1_PH1) = grid.Cell(i-1,grid.noghost -1,dV_MH1_MH1);
//     grid.Node(i,grid.noghost,dV_MH1_PH1) = grid.Cell(i,grid.noghost -1,dV_PH1_MH1);

//     ///SC3
//     //top
//     grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost)  ,dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost-1),dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost-2),dV_ZERO_P3);

//     grid.Cell(i, grid.n_y_cell - (grid.noghost +1),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost)  ,dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost +2),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost-1),dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost +3),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost-2),dV_ZERO_P3);

//     //bottom
//     grid.Node(i, grid.noghost  ,dV_ZERO_P3) = grid.Node(i, grid.noghost - 1,dV_ZERO_M3);
//     grid.Node(i, grid.noghost+1,dV_ZERO_P3) = grid.Node(i, grid.noghost - 2,dV_ZERO_M3);
//     grid.Node(i, grid.noghost+2,dV_ZERO_P3) = grid.Node(i, grid.noghost - 3,dV_ZERO_M3);
    
//     grid.Cell(i, grid.noghost  ,dV_ZERO_P3) = grid.Cell(i, grid.noghost - 1,dV_ZERO_M3);
//     grid.Cell(i, grid.noghost+1,dV_ZERO_P3) = grid.Cell(i, grid.noghost - 2,dV_ZERO_M3);
//     grid.Cell(i, grid.noghost+2,dV_ZERO_P3) = grid.Cell(i, grid.noghost - 3,dV_ZERO_M3);

//     }
// }


// //////////////--------------------------------------------/working bounce back--------------------------------------------------------//////////
// template<typename T>
// void bounce_back(Grid_N_C_2D<T> &grid){
    

//     for ( int i = 1 ;i< grid.n_x_node -1; i++){

//     //TOP
//     //SC1
//      grid.Node(i, grid.n_y_node - (grid.noghost),dV_ZERO_M1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P1);
//      grid.Cell(i, grid.n_y_cell - (grid.noghost),dV_ZERO_M1) = grid.Cell(i, grid.n_y_cell - (grid.noghost+2),dV_ZERO_P1);


//     //SC3
//     grid.Node(i, grid.n_y_node - (grid.noghost  ),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost-1),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost +2),dV_ZERO_P3);
//     grid.Node(i, grid.n_y_node - (grid.noghost-2),dV_ZERO_M3) = grid.Node(i, grid.n_y_node - (grid.noghost +3),dV_ZERO_P3);

//     grid.Cell(i, grid.n_y_cell - (grid.noghost  ),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost +2),dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost-1),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost +3),dV_ZERO_P3);
//     grid.Cell(i, grid.n_y_cell - (grid.noghost-2),dV_ZERO_M3) = grid.Cell(i, grid.n_y_cell - (grid.noghost +4),dV_ZERO_P3);


//     ///SC1/2
// // //noslip bounce back
//     grid.Node(i + 1, grid.n_y_node - (grid.noghost), dV_MH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_PH1_PH1);     //////top right and top left corners are bouncing back in to the domain
//     grid.Node(i - 1, grid.n_y_node - (grid.noghost), dV_PH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_MH1_PH1);

//         //slip bounceback
//     // grid.Node(i , grid.n_y_node - (grid.noghost), dV_PH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_PH1_PH1);     
//     // grid.Node(i , grid.n_y_node - (grid.noghost), dV_MH1_MH1) = grid.Node(i, grid.n_y_node - (grid.noghost +1),dV_MH1_PH1);



//     //Bottom
//     //SC1
   
//     grid.Node(i, grid.noghost - 1,dV_ZERO_P1) = grid.Node(i,grid.noghost +1, dV_ZERO_M1);
//     grid.Cell(i, grid.noghost - 1,dV_ZERO_P1) = grid.Cell(i,grid.noghost   , dV_ZERO_M1);



//     //SC3
//     grid.Node(i, grid.noghost - 1,dV_ZERO_P3) =  grid.Node(i, grid.noghost+1, dV_ZERO_M3);
//     grid.Node(i, grid.noghost - 2,dV_ZERO_P3) =  grid.Node(i, grid.noghost+2, dV_ZERO_M3);
//     grid.Node(i, grid.noghost - 3,dV_ZERO_P3) =  grid.Node(i, grid.noghost+3, dV_ZERO_M3);

//     grid.Cell(i, grid.noghost - 1,dV_ZERO_P3) =  grid.Cell(i, grid.noghost   , dV_ZERO_M3);
//     grid.Cell(i, grid.noghost - 2,dV_ZERO_P3) =  grid.Cell(i, grid.noghost +1, dV_ZERO_M3);
//     grid.Cell(i, grid.noghost - 3,dV_ZERO_P3) =  grid.Cell(i, grid.noghost +2, dV_ZERO_M3);

//     //SC1/2
//     // //noslip bounceback
//     grid.Cell(i + 1, grid.noghost-1, dV_MH1_PH1) =  grid.Cell(i, grid.noghost, dV_PH1_MH1);            //bottom right and left also bounces back                    
//     grid.Cell(i - 1, grid.noghost-1, dV_PH1_PH1) =  grid.Cell(i, grid.noghost, dV_MH1_MH1);

//     // slip bounce back
//     // grid.Cell(i , grid.noghost-1, dV_PH1_PH1) =  grid.Cell(i, grid.noghost, dV_PH1_MH1);
//     // grid.Cell(i , grid.noghost-1, dV_MH1_PH1) =  grid.Cell(i, grid.noghost, dV_MH1_MH1);

// }

// }







template<typename T, typename T1>
void inlet_shock(Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T1> &lbD3Q13, double ux, double uy,double rho){

    for(int i = 0; i < lbgrid.noghost  ;i++){
        for(int j = 0; j<lbgrid.n_y_node - lbgrid.noghost+2;j++){

            for(int dv = 0; dv <lbgrid.d_v; dv++){
                lbgrid.Node(i,j,dv) = lbgrid.Node(lbgrid.noghost,j ,dv);
            }
            for(int dv = 0; dv <lbgrid.d_v; dv++){
                lbgrid.Cell(i,j,dv) = lbgrid.Cell(lbgrid.noghost,j ,dv);
            }   
            
        }
        }
}



template<typename T, typename T1>
void outlet_shock(Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T1> &lbD3Q13, double ux, double uy,double rho){

    for(int i = lbgrid.n_x_node - lbgrid.noghost ; i < lbgrid.n_x_node - lbgrid.noghost +2  ;i++){
        for(int j = 0; j<lbgrid.n_y_node - lbgrid.noghost+2;j++){



            ///============for equilibrizing the outer wall also
            for(int dv = 0; dv <lbgrid.d_v; dv++){
                lbgrid.Node(i,j,dv) = lbgrid.Node(lbgrid.n_x_node - lbgrid.noghost -1,j,dv);
            }
                        
            for(int dv = 0; dv <lbgrid.d_v; dv++){
                lbgrid.Cell(i,j,dv) = lbgrid.Cell(lbgrid.n_x_cell - lbgrid.noghost -1,j,dv);
            }
            


        }
        }
}

















;
#endif