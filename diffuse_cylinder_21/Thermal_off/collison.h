#include<iostream>
#include<cmath>
#include<vector>
#include<math.h>
#include<fstream>
#include<algorithm>
#include <sstream>
#include<string>
#include "lbmD2Q21.h"
#include "GRID_2D.h"
#define PI 3.14159265

//



template<typename T, typename T1>
void collide(Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T1> &lbD3Q21,double beta,double tau ){

// double dtheta = (theta - lbD3Q21.theta0)*lbD3Q21.thetaInverse;
double feq_Node[21] = {0},feq_Cell[21]={0}, ux = 0, uy = 0,theta = 0, rho = 0, G[21]={0};
int a = 0,b = 1; //just for representation of node or cell 





/// first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost) ; i++){
        for(int j =  lbgrid.noghost;j < lbgrid.n_y_node - (lbgrid.noghost) ; j++){
            
            
            get_moments_node(lbgrid, lbD3Q21,ux,uy,rho,theta, i, j);                             //for the node
            get_equi(feq_Node,lbD3Q21, ux, uy, rho,theta);
        

            for (int dv = 0; dv< 21; dv++){

                lbgrid.Node(i,j,dv) =  lbgrid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - lbgrid.Node(i,j,dv));
            
            }


            get_moments_cell(lbgrid, lbD3Q21,ux,uy,rho,theta,i,j);
            get_equi(feq_Cell,lbD3Q21, ux, uy, rho,theta);                          //for the cell
       
            for (int dv = 0; dv<21; dv++){

                lbgrid.Cell(i,j,dv) =  lbgrid.Cell(i,j,dv) + 2.0* beta*(feq_Cell[dv] - lbgrid.Cell(i,j,dv));
          
            }



        }}


}





template<typename T>
inline void get_equi(double feq[21], lbmD2Q21<T> &lb, double ux, double uy, double rho, double theta ){

    // double dtheta = 0; //set it zero for isothermal case




    double thetainv=  lb.thetaInverse;

    double u2 = ux*ux + uy*uy;
    double a1=0, ci2 = 0 , ci4 = 0 ;
    double  first,second,third,feq0=0;

    for (int dv = 0; dv< 21; dv++){

        ci2 = lb.Cx[dv] * lb.Cx[dv] + lb.Cy[dv] * lb.Cy[dv];
        ci4 = ci2*ci2;
        feq0 = rho*lb.W[dv]*(1.0);
        // std::cout<<rho<<std::endl;
        
        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv])*thetainv;
        second = 0.5*(first * first);
        third = -0.5*u2*thetainv;
        feq[dv] = feq0*(1+ first + second + third);
        
    }
        

}



template<typename T,typename T1>
inline void get_moments_node(Grid_N_C_2D<T> &lbgrid, lbmD2Q21<T1> &lb,double &Ux, double &Uy,double &Rho, double &theta,  int X, int Y){ ///node or cell 0-Node 1- cell
    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;
    theta = 0.0;

    for(int dv = 0; dv <21; dv++){
        Ux  += lbgrid.Node(X,Y,dv)*lb.Cx[dv];
        Uy  += lbgrid.Node(X,Y,dv)*lb.Cy[dv];
        Rho += lbgrid.Node(X,Y,dv);
    }  
    
    Ux = Ux/Rho;
    Uy = Uy/Rho;
  
    for(int dv = 0; dv<21; dv++)
        theta += lbgrid.Node(X,Y,dv)*(lb.Cx[dv]*lb.Cx[dv] + lb.Cy[dv]*lb.Cy[dv]);    //=>> f.ci^2

        theta = 0.5*(theta - Rho*(Ux*Ux + Uy*Uy))/(Rho);        //(sum_fi.ci^2 - rho u^2)/(2* Rho)


}






template<typename T,typename T1>
inline void get_moments_cell(Grid_N_C_2D<T> &lbgrid, lbmD2Q21<T1> &lb,double &Ux, double &Uy,double &Rho, double &theta,  int X, int Y){
    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;
    theta = 0.0;



    for(int dv = 0; dv <21; dv++){
        Ux  += lbgrid.Cell(X,Y,dv)*lb.Cx[dv];
        Uy  += lbgrid.Cell(X,Y,dv)*lb.Cy[dv];
        Rho += lbgrid.Cell(X,Y,dv);
    }   
    

    Ux = Ux/Rho;
    Uy = Uy/Rho;

    for(int dv = 0; dv<21; dv++)
        theta += lbgrid.Cell(X,Y,dv)*(lb.Cx[dv]*lb.Cx[dv] + lb.Cy[dv]*lb.Cy[dv]);    //=>> f.ci^2

        theta = 0.5*(theta - Rho*(Ux*Ux + Uy*Uy))/(Rho);        //(sum_fi.ci^2 - rho u^2)/(2* Rho)

}

//setting   getting the 9 moments for the grad
template<typename T,typename T1>
void get_grad_mom_thermal_node(Grid_N_C_2D<T> &lbgrid, lbmD2Q21<T1> &lb,double &Ux, double &Uy,double &Rho,double &theta, double &pxx,double &pyy,double &pxy, int X, int Y){

    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;
    pxx = 0; pyy = 0;
    pxy = 0;



    for(int dv = 0; dv <21; dv++){
        Ux  += lbgrid.Node(X,Y,dv)*lb.Cx[dv];
        Uy  += lbgrid.Node(X,Y,dv)*lb.Cy[dv];
        Rho += lbgrid.Node(X,Y,dv);

        pxx += lbgrid.Node(X,Y,dv)*lb.Cx[dv] * lb.Cx[dv];
        pyy += lbgrid.Node(X,Y,dv)*lb.Cy[dv] * lb.Cy[dv];

        pxy += lbgrid.Node(X,Y,dv)*lb.Cx[dv] * lb.Cy[dv];
    }   
    

    Ux = Ux/Rho;
    Uy = Uy/Rho;

    for(int dv = 0; dv<21; dv++)
        theta += lbgrid.Node(X,Y,dv)*(lb.Cx[dv]*lb.Cx[dv] + lb.Cy[dv]*lb.Cy[dv]);    //=>> f.ci^2

        theta = 0.5*(theta - Rho*(Ux*Ux + Uy*Uy))/(Rho);        //(sum_fi.ci^2 - rho u^2)/(2* Rho)



}

template<typename T,typename T1>
void get_grad_mom_thermal_Cell(Grid_N_C_2D<T> &lbgrid, lbmD2Q21<T1> &lb,double &Ux, double &Uy,double &Rho,double &theta, double &pxx,double &pyy,double &pxy, int X, int Y){

    Ux  = 0.0;
    Uy  = 0.0;
    Rho = 0.0;
    pxx = 0; pyy = 0;
    pxy = 0;



    for(int dv = 0; dv <21; dv++){
        Ux  += lbgrid.Cell(X,Y,dv)*lb.Cx[dv];
        Uy  += lbgrid.Cell(X,Y,dv)*lb.Cy[dv];
        Rho += lbgrid.Cell(X,Y,dv);

        pxx += lbgrid.Cell(X,Y,dv)*lb.Cx[dv] * lb.Cx[dv];
        pyy += lbgrid.Cell(X,Y,dv)*lb.Cy[dv] * lb.Cy[dv];

        pxy += lbgrid.Cell(X,Y,dv)*lb.Cx[dv] * lb.Cy[dv];
    }   
    Ux = Ux/Rho;
    Uy = Uy/Rho;


    for(int dv = 0; dv<21; dv++)
        theta += lbgrid.Cell(X,Y,dv)*(lb.Cx[dv]*lb.Cx[dv] + lb.Cy[dv]*lb.Cy[dv]);    //=>> f.ci^2

        theta = 0.5*(theta - Rho*(Ux*Ux + Uy*Uy))/(Rho);        //(sum_fi.ci^2 - rho u^2)/(2* Rho)



}





//setting getting the grad population
template<typename T>
void get_grad_isothermal(double feq[21], lbmD2Q21<T> &lb, double ux, double uy, double rho, double pxx, double pyy, double pxy){


for (int dv = 0 ; dv< 21; dv++){

        feq[dv] = lb.W[dv]*(
                        rho +   
                        (ux*lb.Cx[dv] + uy * lb.Cy[dv])*lb.thetaInverse*rho   + 
                        0.5*lb.thetaInverse*lb.thetaInverse*((pxx - rho*lb.theta0)*(lb.Cx[dv] *lb.Cx[dv] - lb.theta0) +
                                                             (pyy - rho*lb.theta0)*(lb.Cy[dv] *lb.Cy[dv] - lb.theta0) + 
                                                             2.0*(pxy)*(lb.Cx[dv]* lb.Cy[dv])
                                                            )                       
                        );
}
}



template<typename T>
void grad_inlet(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lb, double u0){

double f[21] = {0};

double ux, uy, rho, pxx, pyy, pxy,theta;

for(int j = 0; j< grid.n_y_node - 0; j ++){

    get_grad_mom_thermal_node(grid, lb, ux, uy, rho,theta, pxx, pyy,pxy,grid.noghost, j);

    ux = u0, uy = 0.0;

    get_grad_isothermal(f,lb,ux,uy, rho, pxx, pyy, pxy);

    for(int dv = 0; dv< 21; dv++){

        grid.Node(0,j,dv) = f[dv];
        grid.Node(1,j,dv) = f[dv];
        
        grid.Cell(0,j,dv) = f[dv];
        grid.Cell(1,j,dv) = f[dv];

    }
}
}



template<typename T>
void grad_outlet(Grid_N_C_2D<T> &grid,lbmD2Q21<T> &lb){

double f[21] = {0};

double ux, uy, rho, pxx, pyy, pxy,theta;

for(int j = 0; j< grid.n_y_node -0 ; j ++){

    get_grad_mom_thermal_Cell(grid, lb, ux, uy, rho,theta, pxx, pyy,pxy, grid.n_x_node - grid.noghost -1 , j);
    get_grad_isothermal(f,lb,ux,uy, rho, pxx, pyy, pxy);

    for(int dv = 0; dv< 21; dv++){
        grid.Node(grid.n_x_node - grid.noghost   ,j,dv) = f[dv];
        grid.Node(grid.n_x_node - grid.noghost +1,j,dv) = f[dv];
        
        grid.Cell(grid.n_x_node - grid.noghost   ,j,dv) = f[dv];
        grid.Cell(grid.n_x_node - grid.noghost +1,j,dv) = f[dv];
    }
}
}










template<typename T, typename T1>
void initialization(Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T1> &lbD2Q21, double U0){
     std::ofstream myfile ("ux_uy.dat"); //for the kida initialization
  if (myfile.is_open())
  { 
    double Feq_node[21] = {0},Feq_cell[21] = {0},Rho = 1.0;
    double x,y,
           n_to_n_dist = (2*PI)/lbgrid.n_x;    ///distance between nodes 

    double ux_node=0,uy_node=0,ux_cell=0,uy_cell=0;
    

    
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost); i++){
        for(int j = 0 + lbgrid.noghost; j < lbgrid.n_y_node - (lbgrid.noghost); j++){

            get_equi(Feq_node,lbD2Q21,ux_node,uy_node,Rho,lbD2Q21.theta0);
            for (int dv = 0; dv<21; dv++){
                lbgrid.Node(i,j,dv) = Feq_node[dv];
            }


            get_equi(Feq_cell,lbD2Q21,ux_cell,uy_cell, Rho,lbD2Q21.theta0);
            for (int dv = 0; dv<21; dv++){
                 lbgrid.Cell(i,j,dv) = Feq_cell[dv];
            }
        }
        }
        }
        }



template<typename T>
void printMass(Grid_N_C_2D<T> &lbgrid){    
    double a = 0;
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost) ; i++){
        for(int j = 0 + lbgrid.noghost;j < lbgrid.n_y_node - (lbgrid.noghost) ; j++){
            for (int dv = 0; dv< lbgrid.d_v; dv++){
                a += lbgrid.Node(i,j,dv) + lbgrid.Cell(i,j,dv);
            }
        }}
std::cout<<"   "<<a<<std::endl;

}


void get_coordinates(std::vector<double> &x,std::vector<double> &y, const std::string& filename){

    int j = 0;
   
    std::string a1,a2;
    std::ifstream my_file;
    my_file.open(filename);
    


        while(!my_file.eof())
        {
             getline(my_file, a1, ',');
             
             x.push_back(stod(a1));  ///////////////check whether there is no blank line at the end
            
             getline(my_file, a2, '\n');
             y.push_back(stod(a2));       
            
             
            j += 1;
            
        }     }


void get_markers(std::vector<double> &x,std::vector<double> &y, std::vector<int> &marker, const std::string& filename){

    int j = 0;
   
    std::string a1,a2,a3;
    std::ifstream my_file;
    my_file.open(filename);

    
        while(!my_file.eof())
        {                   

            getline(my_file, a1, ',');
            x.push_back(stod(a1));  ///////////////check whether there is no blank line at the end
                           

            getline(my_file, a2, ',');
            y.push_back(stod(a2)); 

            getline(my_file, a3, '\n');
            marker.push_back(stoi(a3));       
        
            j += 1;
            
        }     }



template<typename T1>
void coordinates_writer(Grid_N_C_2D<T1> &marker,int ic, int jc, double D){
    std::ofstream Myfiles("nodes_coordinates.csv");
    std::ofstream Myfiles1("cells_coordinates.csv");

    double delx = 1/D;
    double x_node,y_node, x_cell, y_cell;

    for(int i = 0 + marker.noghost; i < marker.n_x_node - (marker.noghost) ; i++){
        for(int j = 0 + marker.noghost; j < marker.n_y_node - (marker.noghost) ; j++){

            x_node = (i - ic)*delx;
            y_node = (j - jc)*delx; 

            Myfiles<<x_node<<","<<y_node<<std::endl;

            x_cell = x_node + (delx/2);
            y_cell = y_node + (delx/2);

            Myfiles1<<x_cell<<","<<y_cell<<std::endl;

        }}

}


//takes the marker from the files
template<typename T1>
void Mark_for_solid(Grid_N_C_2D<T1> &marker, int ic , int jc, double L){

    double delx = 1/L;

    std::vector<double> X_node,Y_node;  std::vector<int> marker_node;
    std::vector<double> X_cell,Y_cell;  std::vector<int> marker_cell;

    get_markers(X_node,Y_node,marker_node,"nodes_marker.csv");
    get_markers(X_cell,Y_cell,marker_cell,"cells_marker.csv");

    int index = 0;

    for(int i = 0 + marker.noghost; i < marker.n_x_node - (marker.noghost) ; i++){
        for(int j = 0 + marker.noghost;j < marker.n_y_node - (marker.noghost) ; j++){


            if(marker_node[index] == 1){
                marker.Node(i,j,0) = 1;
            }

            if(marker_cell[index] ==1){
                marker.Cell(i,j,0) = 1;
            }

            index +=1;

        }}


    std::ofstream Myfiles("values.txt");

    for(int i = 0 + marker.noghost; i < marker.n_x_node - (marker.noghost) ; i++){
        for(int j = 0 + marker.noghost;j < marker.n_y_node - (marker.noghost) ; j++){

            if(marker.Node(i,j,0)==1)
                Myfiles<<(i - ic)*delx<<","<<(j - jc)*delx<<std::endl;

            if(marker.Cell(i,j,0) ==1)
                Myfiles<<(i - ic)*delx + delx* 0.5<<","<<(j - jc)*delx + delx* 0.5 <<std::endl;


        }}

}





///////---------------works for solid between x axis  of 0 and 1-----------///
template<typename T1>
void Mark_solid(Grid_N_C_2D<T1> &marker){
    std::ofstream Myfiles("values.txt");
    std::vector<double> X_U,Y_U,X_L,Y_L;
    get_coordinates(X_U,Y_U,"0012_U.csv");
    get_coordinates(X_L,Y_L,"0012_L.csv");
    double delx = 0.025;
    double y_inter;
    int pos_xr;
  

    for(int i = 0 + marker.noghost; i < marker.n_x_node - (marker.noghost) ; i++){
        for(int j = 0 + marker.noghost;j < marker.n_y_node - (marker.noghost) ; j++){
            
            
            
            ///marking the nodes
            double x_node, y_node;
            x_node = (i - 3) * delx - delx*marker.n_x*0.25;
            y_node = (j - 3) * delx - delx*marker.n_y*0.5;          /////////for aerofoil 1.5

        

            std::vector<double>::iterator lower;
            if (x_node >0.0 && x_node < 1.0 ){

                if (y_node >= 0.0){   //////////

                /////upper half plane
                lower = std::lower_bound(X_U.begin(), X_U.end(), x_node);
                
                pos_xr = (lower - X_U.begin());  ////position of next largest no
               
                y_inter= ((Y_U[pos_xr]-Y_U[pos_xr-1])/(X_U[pos_xr] - X_U[pos_xr-1]))*(x_node - X_U[pos_xr-1])+ Y_U[pos_xr-1];
                if (y_node<y_inter){
                    

                    //add the i,j to the solid list 
                    marker.Node(i,j,0) = SOLID;
                  Myfiles<<i<<","<<j<<std::endl;

                }}

                if (y_node<0.0){
                ////lower plane
                lower = std::lower_bound(X_L.begin(), X_L.end(), x_node);

                pos_xr = (lower - X_L.begin());  ////position of next largest no

                y_inter= ((Y_L[pos_xr]-Y_L[pos_xr-1])/(X_L[pos_xr] - X_L[pos_xr-1]))*(x_node - X_L[pos_xr-1]) + Y_L[pos_xr-1];
                            
                if (y_node>y_inter){
                    //add the i,j to the solid list 
                    marker.Node(i,j,0) = SOLID;
                  Myfiles<<i<<","<<j<<std::endl;
                }}
                
            }





            double x_cell, y_cell;
            x_cell = x_node + (delx/2);
            y_cell = y_node + (delx/2);
            //std::cout<<x_cell<<std::endl;
            if (x_cell >0.0 && x_cell < 1.0 ){
                // std::cout<<"yes"<<" "<<y_node<<std::endl;

                if (y_cell > 0.0){


                /////upper half plane
                lower = std::lower_bound(X_U.begin(), X_U.end(), x_cell);
                
                pos_xr = (lower - X_U.begin());  ////position of next largest no
               
                y_inter= ((Y_U[pos_xr]-Y_U[pos_xr-1])/(X_U[pos_xr] - X_U[pos_xr-1]))*(x_cell - X_U[pos_xr-1])+ Y_U[pos_xr-1];
                if (y_cell<y_inter){
                    
                    //add the i,j to the solid list 
                    marker.Cell(i,j,0) = SOLID;

                }}

                if (y_cell<0.0){
                ///////////----------------------------------------lower plane--------------------------------------------/////////////////
                lower = std::lower_bound(X_L.begin(), X_L.end(), x_cell);

                pos_xr = (lower - X_L.begin());  ////position of next largest no

                y_inter= ((Y_L[pos_xr]-Y_L[pos_xr-1])/(X_L[pos_xr] - X_L[pos_xr-1]))*(x_cell - X_L[pos_xr-1]) + Y_L[pos_xr-1];
                            
                if (y_cell>y_inter){
                    //add the i,j to the solid list 
                    marker.Cell(i,j,0) = SOLID;
                }}
            }           
}}}
















































//    }