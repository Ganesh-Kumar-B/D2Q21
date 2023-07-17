#ifndef FORCE_H
#define FORCE_H

#include<iostream>

#include"lbmD2Q21.h"
//--------------------------------------------------------DRTT------------------------------------////

template<typename T>
void J_domain( Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T> &lb, double &Jx, double &Jy , int left, int right , int top, int bottom){//if we take nodes as boundary

Jx = 0; Jy = 0;



///////////for no
for(int i = left; i<=right; i++){
    for (int j = bottom; j<= top  ; j++){

// if(lbgrid.marker_node(i,j) == 0){

        for(int dv = 0 ; dv< lbgrid.d_v; dv++){

            Jx += lbgrid.Node(i,j, dv) *lb.Cx[dv];
            Jy += lbgrid.Node(i,j, dv) *lb.Cy[dv]; 
        // }
}
    }
}

for(int i = left; i<= right - 1; i++){
    for(int j = bottom  ; j<= top -1; j++){

// if(lbgrid.marker_cell(i,j) == 0){

        for(int dv = 0; dv< lbgrid.d_v; dv++){

            Jx += lbgrid.Cell(i,j,dv) * lb.Cx[dv];
            Jy += lbgrid.Cell(i,j,dv) * lb.Cy[dv];

        }
// }
// 
    }
}
} 

// //--------------------------------------momentumn in and out of the system-----------------------------//
template<typename T, typename T1>
void Force_calc( Grid_N_C_2D<T> &grid,lbmD2Q21<T1> &lb,int left, int right, int top, int bottom,
                    double &Fx_in, double &Fy_in, double &fx_out, double &fy_out ){

Fx_in = 0; Fy_in = 0, fx_out = 0; fy_out = 0;



int count = 0;
int nx=0, ny=0;

int offset = 5; // this offset is just starting of the loop  *not the domain ends



//SETTING  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>INCOMING POPULATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//SC SCSCSCSCSCSSCS

 for(int i = left - offset; i<right + offset; i++)
    if(i <= left -1 || i>= right +1)
    for(int j = bottom - offset; j< top + offset; j ++)
        for(int dv = 0; dv< 9; dv++)
            if(i + (int)lb.Cx[dv] >= left && i + (int)lb.Cx[dv] <= right && j +(int) lb.Cy[dv] >= bottom && j + (int) lb.Cy[dv] <= top){count +=1;   // checkfor the node

                    Fx_in += grid.Node(i,j,dv) *lb.Cx[dv];
                    Fy_in += grid.Node(i,j,dv) *lb.Cy[dv];

        }

for(int j = bottom - offset; j< top + offset; j ++)
    if(j <= bottom -1 || j>= top +1)
     for(int i = left - offset; i<right + offset; i++)
        for(int dv = 0; dv< 9; dv++)
            if(i + (int)lb.Cx[dv] >= left && i + (int)lb.Cx[dv] <= right && j +(int) lb.Cy[dv] >= bottom && j + (int) lb.Cy[dv] <= top){ count +=1;    // checkfor the node

                    Fx_in += grid.Node(i,j,dv) *lb.Cx[dv];
                    Fy_in += grid.Node(i,j,dv) *lb.Cy[dv];

        }



//     for the cells
 for(int i = left - offset; i<right + offset; i++)
    if(i <= left -1 || i>= right )
    for(int j = bottom - offset; j< top + offset; j ++)
        for(int dv = 0; dv< 9; dv++)
            if(i + (int)lb.Cx[dv] >= left && i + (int)lb.Cx[dv] <= right-1 && j +(int) lb.Cy[dv] >= bottom && j + (int) lb.Cy[dv] <= top -1){ count +=1;   // checkfor the Cell

                    Fx_in += grid.Cell(i,j,dv) *lb.Cx[dv];
                    Fy_in += grid.Cell(i,j,dv) *lb.Cy[dv];

        }

 for(int j = bottom - offset; j< top + offset; j ++)
    if(j <= bottom -1 || j>= top )
     for(int i = left - offset; i<right + offset; i++)
        for(int dv = 0; dv< 9; dv++)
            if(i + (int)lb.Cx[dv] >= left && i + (int)lb.Cx[dv] <= right-1 && j +(int) lb.Cy[dv] >= bottom && j + (int) lb.Cy[dv] <= top -1){ count +=1;  // checkfor the cell

                    Fx_in += grid.Cell(i,j,dv) *lb.Cx[dv];
                    Fy_in += grid.Cell(i,j,dv) *lb.Cy[dv];

        }
 


//....................//setting BCC BCC BCC BCCBCC BCCBCC BCCBCC BCCBCC BCCBCC BCC
//Bcc populations that come from a Node and land in the cell in the C.V.
 for(int i = left - offset; i<right + offset; i++)
    if(i <= left -1 || i>= right +1)
    for(int j = bottom - offset; j< top + offset; j ++)
        for(int dv = 9; dv< 21; dv++)
            if(i + (int)lb.CxF[dv] >= left && i + (int)lb.CxF[dv] <= right-1 && j +(int) lb.CyF[dv] >= bottom && j + (int) lb.CyF[dv] <= top -1){count +=1;   // checkfor the cell
                // if it landed on one one of the cells in the C.V from the cell

                Fx_in += grid.Node(i,j,dv) *lb.Cx[dv];
                Fy_in += grid.Node(i,j,dv) *lb.Cy[dv];

            }








 for(int j = bottom - offset; j< top + offset; j ++)
    if(j <= bottom -1 || j>= top +1 )
     for(int i = left ; i<= right; i++)
        for(int dv = 9; dv< 21; dv++)
            if(i + (int)lb.CxF[dv] >= left && i + (int)lb.CxF[dv] <= right-1 && j +(int) lb.CyF[dv] >= bottom && j + (int) lb.CyF[dv] <= top -1){count +=1;   // checkfor the cell
                                                    
                    Fx_in += grid.Node(i,j,dv) *lb.Cx[dv];
                    Fy_in += grid.Node(i,j,dv) *lb.Cy[dv];

        }


///  fixme  BCC Cells
 for(int i = left - offset; i<right + offset; i++)
    if(i <= left -1 || i>= right)
    for(int j = bottom - offset; j< top + offset; j ++)
        for(int dv = 9; dv< 21; dv++)
        if(i + (int)lb.CxC[dv] >= left && i + (int)lb.CxC[dv] <= right && j +(int) lb.CyC[dv] >= bottom && j + (int) lb.CyC[dv] <= top){count +=1; 


                Fx_in += grid.Cell(i,j,dv) *lb.Cx[dv];
                Fy_in += grid.Cell(i,j,dv) *lb.Cy[dv];
        }

 for(int j = bottom - offset; j< top + offset; j ++)
    if(j <= bottom -1 || j>= top )
     for(int i = left ; i<= right -1; i++)
        for(int dv = 9; dv< 21; dv++)
        if(i + (int)lb.CxC[dv] >= left && i + (int)lb.CxC[dv] <= right && j +(int) lb.CyC[dv] >= bottom && j + (int) lb.CyC[dv] <= top){count +=1; 

            Fx_in += grid.Cell(i,j,dv) *lb.Cx[dv];
            Fy_in += grid.Cell(i,j,dv) *lb.Cy[dv];
        }

// std::cout<<"before"<<count<<std::endl;
count = 0;


//setting    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>   OUTGOING POPULATIONS<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<///
 

//SC

// .....for the NODES
for(int i = left; i<=right; i++)
    for (int j = bottom; j<= top  ; j++)
        for(int dv = 0; dv<9; dv++)
            if(!(i + (int)lb.Cx[dv] >= left && i + (int)lb.Cx[dv] <= right && j +(int) lb.Cy[dv] >= bottom && j + (int) lb.Cy[dv] <= top)){ count +=1;  // checkfor the node

                fx_out += grid.Node(i,j,dv) * lb.Cx[dv];
                fy_out += grid.Node(i,j,dv) * lb.Cy[dv];

            }

//......for the CELLS
for(int i = left; i<=right-1; i++)
    for (int j = bottom; j<= top-1  ; j++)
        for(int dv = 0; dv<9; dv++)
            if(!(i + (int)lb.Cx[dv] >= left && i + (int)lb.Cx[dv] <= right-1 && j +(int) lb.Cy[dv] >= bottom && j + (int) lb.Cy[dv] <= top -1)){ count +=1;  // checkfor the cell

                fx_out += grid.Cell(i,j,dv) * lb.Cx[dv];
                fy_out += grid.Cell(i,j,dv) * lb.Cy[dv];

            }       
        
//...........BCC
//from the NODES
for(int i = left; i<=right; i++)
    for (int j = bottom; j<= top  ; j++)
        for(int dv = 9; dv<21; dv++)
            if(!(i + (int)lb.CxF[dv] >= left && i + (int)lb.CxF[dv] <= right-1 && j +(int) lb.CyF[dv] >= bottom && j + (int) lb.CyF[dv] <= top -1)){count +=1; 

                fx_out += grid.Node(i,j,dv) * lb.Cx[dv];
                fy_out += grid.Node(i,j,dv) * lb.Cy[dv];    

            }



for(int i = left; i<=right-1; i++)
    for (int j = bottom; j<= top-1  ; j++)
        for(int dv = 9; dv<21; dv++)
            if(!(i + (int)lb.CxC[dv] >= left && i + (int)lb.CxC[dv] <= right && j +(int) lb.CyC[dv] >= bottom && j + (int) lb.CyC[dv] <= top)){count +=1; 

//std::cout<<i<<","<<j<<","<<dv<<std::endl;
                fx_out += grid.Cell(i,j,dv) * lb.Cx[dv];
                fy_out += grid.Cell(i,j,dv) * lb.Cy[dv];


            }


//  std::cout<<"after"<<count<<std::endl;

 
 }




template<typename T, typename T1>
void Force_SI_cylinder( Grid_N_C_2D<T> &grid,Grid_N_C_2D<T1> &marker,lbmD2Q21<T> &lb,int left, int right, int top, int bottom, int ic, int jc,
                    double &FxSI, double &FySI,  double beta ){
    FxSI = 0; FySI = 0;
                        double ux,uy, rho = 1.0, feq[21] = {0};

    double theta = 0;

double nx = 0, ny = 0; double r = 0;


    for(int i = left; i<=right; i++){
        for(int j = bottom; j<=top; j++){

            int ib,jb;

            if(marker.Node(i,j,0) == FLUID){
                ib = 0; jb = 0;
                for(int dv = 0; dv<5; dv++)
                if(marker.Node(i +(int) lb.Cx[dv], j +(int) lb.Cy[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }

                for(int dv = 9; dv<13; dv++)
                if(marker.Cell(i +(int) lb.CxF[dv], j +(int) lb.CyF[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }



                double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;

                if(ib != 0 && jb!= 0){ 
                    // this we are working on the boundary points if ib == 0 and jb == 0 then we enter  in this condition; 
                    get_moments_node(grid,lb,ux,uy,rho,theta, ib,jb);
                    get_equi(feq, lb, ux, uy, rho,theta);

                    double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                    double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;
                    for(int dv = 0; dv<21; dv++){

                        tau_xx += (1 - beta)*(grid.Node(i,j,dv) - feq[dv]) *(lb.Cx[dv] * lb.Cx[dv] - 0.5 * lb.theta0);
                        tau_yy += (1 - beta)*(grid.Node(i,j,dv) - feq[dv]) *(lb.Cy[dv] * lb.Cy[dv] - 0.5 * lb.theta0);

                        tau_xy += (1 - beta)*(grid.Node(i,j,dv) - feq[dv]) *(lb.Cx[dv] * lb.Cy[dv]);

                    }

                    r = sqrt((ib - ic)*(ib - ic) +(jb - jc)*(jb - jc));

                    nx = (ib - ic)/r;
                    ny = (jb - jc)/r;

                    
                    sigma_xx = tau_xx - beta*rho*lb.theta0;  // the total stress tensor on the object added with - p = rho*theta
                    sigma_yy = tau_yy - beta*rho*lb.theta0;
                    sigma_xy = tau_xy;

                    FxSI =  FxSI + sigma_xx*nx + sigma_xy*ny;
                    FySI =  FySI + sigma_xy*nx + sigma_yy*ny;

                    


                }
        }

        ib = 0; jb = 0;

            if(marker.Cell(i,j,0) == FLUID){
                for(int dv = 0; dv<5; dv++)
                if(marker.Cell(i +(int) lb.Cx[dv], j +(int) lb.Cy[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }

                for(int dv = 9; dv<13; dv++)
                if(marker.Node(i +(int) lb.CxC[dv], j +(int) lb.CyC[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }



                double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;

                if(ib != 0 && jb!= 0){
                    get_moments_cell(grid,lb,ux,uy,rho,theta, ib,jb);
                    get_equi(feq, lb, ux, uy, rho,theta);

                    double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                    double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;
                    for(int dv = 0; dv<21; dv++){

                        tau_xx += (1 - beta)*(grid.Cell(i,j,dv) - feq[dv]) *(lb.Cx[dv] * lb.Cx[dv] - 0.5 * lb.theta0);
                        tau_yy += (1 - beta)*(grid.Cell(i,j,dv) - feq[dv]) *(lb.Cy[dv] * lb.Cy[dv] - 0.5 * lb.theta0);

                        tau_xy += (1 - beta)*(grid.Cell(i,j,dv) - feq[dv]) *(lb.Cx[dv] * lb.Cy[dv]);


                    }

                    r = sqrt((ib+0.5 - ic)*(ib+0.5 - ic) +(jb+0.5 - jc)*(jb+0.5 - jc));

                    nx = (ib+0.5 - ic)/r;
                    ny = (jb+0.5 - jc)/r;

                    
                    sigma_xx = tau_xx - beta*rho*lb.theta0;  // the total stress tensor on the object added with - p = rho*theta
                    sigma_yy = tau_yy - beta*rho*lb.theta0;
                    sigma_xy = tau_xy;

                    FxSI =  FxSI + sigma_xx*nx + sigma_xy*ny;
                    FySI =  FySI + sigma_xy*nx + sigma_yy*ny;
                    


                }
        }




        }
    }

    }

















template<typename T, typename T1>
void Force_SI_cylinder1( Grid_N_C_2D<T> &grid,Grid_N_C_2D<T1> &marker,lbmD2Q21<T> &lb,int left, int right, int top, int bottom, int ic, int jc,
                    double &FxSI, double &FySI,  double beta ){
    FxSI = 0; FySI = 0;
                        double ux,uy, rho = 1.0, feq[21] = {0};

    double theta = 0;

double nx = 0, ny = 0; double r = 0;


    for(int i = left; i<=right; i++){
        for(int j = bottom; j<=top; j++){

            int ib,jb;

            if(marker.Node(i,j,0) == FLUID){
                ib = 0; jb = 0;
                for(int dv = 0; dv<5; dv++)
                if(marker.Node(i +(int) lb.Cx[dv], j +(int) lb.Cy[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }

                for(int dv = 9; dv<13; dv++)
                if(marker.Cell(i +(int) lb.CxF[dv], j +(int) lb.CyF[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }



                double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;

                if(ib != 0 && jb!= 0){ 
                    // this we are working on the boundary points if ib == 0 and jb == 0 then we enter  in this condition; 
                    get_moments_node(grid,lb,ux,uy,rho,theta, ib,jb);
                    get_equi(feq, lb, ux, uy, rho,theta);

                    double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                    double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;
                    for(int dv = 0; dv<21; dv++){

                        tau_xx += (grid.Node(i,j,dv) * (lb.Cx[dv] - ux) *(lb.Cx[dv] - ux) );
                        tau_yy += (grid.Node(i,j,dv) * (lb.Cy[dv] - uy) *(lb.Cy[dv] - uy) );

                        tau_xy += (grid.Node(i,j,dv) * (lb.Cx[dv] - ux) *(lb.Cy[dv] - uy) );

                    }

                    r = sqrt((ib - ic)*(ib - ic) +(jb - jc)*(jb - jc));

                    nx = (ib - ic)/r;
                    ny = (jb - jc)/r;

                    
                    sigma_xx = - beta * rho*lb.theta0 - (1.0 - beta)*tau_xx;  // the total stress tensor on the object added with - p = rho*theta
                    sigma_yy = - beta * rho*lb.theta0 - (1.0 - beta)*tau_yy;
                    sigma_xy = - (1.0 - beta)*tau_xy ;

                    FxSI =  FxSI + sigma_xx*nx + sigma_xy*ny;
                    FySI =  FySI + sigma_xy*nx + sigma_yy*ny;

                    


                }
        }

        ib = 0; jb = 0;

            if(marker.Cell(i,j,0) == FLUID){
                for(int dv = 0; dv<5; dv++)
                if(marker.Cell(i +(int) lb.Cx[dv], j +(int) lb.Cy[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }

                for(int dv = 9; dv<13; dv++)
                if(marker.Node(i +(int) lb.CxC[dv], j +(int) lb.CyC[dv],0) == SOLID){
                    ib = i;
                    jb = j;  //this will imply that  we are on the boundary points
                }



                double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;

                if(ib != 0 && jb!= 0){
                    get_moments_cell(grid,lb,ux,uy,rho,theta, ib,jb);
                    get_equi(feq, lb, ux, uy, rho,theta);

                    double tau_xx = 0, tau_yy = 0, tau_xy = 0;
                    double sigma_xx = 0, sigma_yy = 0, sigma_xy = 0;
                    for(int dv = 0; dv<21; dv++){

                        tau_xx += (grid.Cell(i,j,dv) * (lb.Cx[dv] - ux) *(lb.Cx[dv] - ux) );
                        tau_yy += (grid.Cell(i,j,dv) * (lb.Cy[dv] - uy) *(lb.Cy[dv] - uy) );

                        tau_xy += (grid.Cell(i,j,dv) * (lb.Cx[dv] - ux) *(lb.Cy[dv] - uy) );

                    }

                    r = sqrt((ib+0.5 - ic)*(ib+0.5 - ic) +(jb+0.5 - jc)*(jb+0.5 - jc));

                    nx = (ib+0.5 - ic)/r;
                    ny = (jb+0.5 - jc)/r;

                    
                    sigma_xx = - beta * rho*lb.theta0 - (1.0 - beta)*tau_xx;  // the total stress tensor on the object added with - p = rho*theta
                    sigma_yy = - beta * rho*lb.theta0 - (1.0 - beta)*tau_yy;
                    sigma_xy = - (1.0 - beta)*tau_xy ;

                    FxSI =  FxSI + sigma_xx*nx + sigma_xy*ny;
                    FySI =  FySI + sigma_xy*nx + sigma_yy*ny;

                    


                }
        }




        }
    }

    }

























;

#endif FORCE_H