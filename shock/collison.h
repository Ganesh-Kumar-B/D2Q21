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





template<typename T, typename T1>
void collide(Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T1> &lbD3Q21,double beta,double tau, double gx, double gy ){

// double dtheta = (theta - lbD3Q21.theta0)*lbD3Q21.thetaInverse;
double feq_Node[21] = {0},feq_Cell[21]={0}, ux = 0, uy = 0,theta = 0, rho = 0, G[21]={0};
int a = 0,b = 1; //just for representation of node or cell 





/// first the population of nodes are resetted and second  the population of the cells are resetted
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost) ; i++){
        for(int j = /*1 +*/ lbgrid.noghost;j < lbgrid.n_y_node - (lbgrid.noghost) ; j++){
            
            
            get_moments_node(lbgrid, lbD3Q21,  ux, uy,rho,theta, i, j);                             //for the node
            
            // std::cout<<theta<<" i "<<i<<" j "<<j<<std::endl;
            get_equi(feq_Node,lbD3Q21, ux, uy, rho,theta);
            // std::cout<<ux<<std::endl;
        

            for (int dv = 0; dv< 21; dv++){
                G[dv] = rho * lbD3Q21.W[dv] * lbD3Q21.thetaInverse * ( (lbD3Q21.Cx[dv])*gx /*+ (lbD3Q21.Cy[dv])*gy */) ;
                
                     
            lbgrid.Node(i,j,dv) =  lbgrid.Node(i,j,dv) + 2.0* beta*(feq_Node[dv] - lbgrid.Node(i,j,dv))+ (2.0*beta*tau*G[dv]);
         //std::cout<<(2.0*beta*tau*G[dv])<<std::endl;
            }}}



            for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost)   ; i++){
            for(int j = 0 + lbgrid.noghost; j < lbgrid.n_y_node - (lbgrid.noghost)  /*- 1*/; j++){
            
             ux = 0; uy = 0;



            get_moments_cell(lbgrid, lbD3Q21,ux,uy,rho,theta,i,j);
            get_equi(feq_Cell,lbD3Q21, ux, uy, rho,theta);                          //for the cell
       
            for (int dv = 0; dv<21; dv++){
                G[dv] = rho * lbD3Q21.W[dv] * lbD3Q21.thetaInverse * ( (lbD3Q21.Cx[dv])*gx /*+ (lbD3Q13.Cy[dv])*gy */) ;
                // std::cout<<G[dv]<<"  "<<dv<<std::endl;
            lbgrid.Cell(i,j,dv) =  lbgrid.Cell(i,j,dv) + 2.0* beta*(feq_Cell[dv] - lbgrid.Cell(i,j,dv))+ (2.0*beta*tau*G[dv]);
          
            }

        }//std::cout<<"d"<<std::endl;
    }
}





template<typename T>
void get_equi(double feq[21], lbmD2Q21<T> &lb, double ux, double uy, double rho, double theta ){

    double dtheta  = (theta - lb.theta0)/ lb.theta0;
    double thetainv=  1/theta;

    double u2 = ux*ux + uy*uy;
    double a1=0, ci2 = 0 , ci4 = 0 ;
    double  first,second,third,feq0=0;
    for (int dv = 0; dv< 21; dv++){

        ci2 = lb.Cx[dv] * lb.Cx[dv] + lb.Cy[dv] * lb.Cy[dv];
        ci4 = ci2*ci2;
        feq0 = rho*lb.W[dv]*(1 + 0.5 * dtheta *(ci2 * lb.thetaInverse -2) + dtheta*dtheta*(1 - ci2*lb.thetaInverse + 0.125 * ci4 * lb.thetaInverse* lb.thetaInverse));
        // std::cout<<rho<<std::endl;
        
        first  = (ux*lb.Cx[dv] + uy*lb.Cy[dv])*thetainv;
        second = 0.5*(first * first);
        third = -0.5*u2*thetainv;
        feq[dv] = feq0*(1+ first + second + third);
        


    }
         
}



template<typename T,typename T1>
void get_moments_node(Grid_N_C_2D<T> &lbgrid, lbmD2Q21<T1> &lb,double &Ux, double &Uy,double &Rho, double &theta,  int X, int Y){ ///node or cell 0-Node 1- cell
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
void get_moments_cell(Grid_N_C_2D<T> &lbgrid, lbmD2Q21<T1> &lb,double &Ux, double &Uy,double &Rho, double &theta,  int X, int Y){
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





template<typename T, typename T1>
void initialization_shock(Grid_N_C_2D<T> &lbgrid,lbmD2Q21<T1> &lbD2Q21, double U0){
     std::ofstream myfile ("ux_uy.dat"); //for the kida initialization
  if (myfile.is_open())
  { 
    double Feq_node[21] = {0},Feq_cell[21] = {0},Rho = 1.0;
    double x,y,theta1,theta2,
           n_to_n_dist = (2*PI)/lbgrid.n_x;    ///distance between nodes 



    double ux_node=0,uy_node=0,ux_cell=0,uy_cell=0;
    
    std::cout <<"half" <<(lbgrid.n_x_node - (lbgrid.noghost))<<std::endl;   

    Rho = 1;
    theta1 = 1* lbD2Q21.theta0;

    get_equi(Feq_node,lbD2Q21,ux_node,uy_node,Rho,theta1);


    for(int i = 0 + lbgrid.noghost; i < (lbgrid.n_x_node - (lbgrid.noghost))/2; i++){
        for(int j = 0 + lbgrid.noghost; j < (lbgrid.n_y_node - (lbgrid.noghost)); j++){


            for (int dv = 0; dv<21; dv++){
                lbgrid.Node(i,j,dv) = Feq_node[dv];
            }

            for (int dv = 0; dv<21; dv++){
                 lbgrid.Cell(i,j,dv) = Feq_node[dv];
            }
        }
        }
        

            
    Rho     = 0.5;
    theta2  = 1.0 * lbD2Q21.theta0;
    get_equi(Feq_node,lbD2Q21,ux_node,uy_node,Rho,theta2);


    for(int i = (lbgrid.n_x_node - (lbgrid.noghost))/2; i < (lbgrid.n_x_node - (lbgrid.noghost)); i++){
        for(int j = 0 + lbgrid.noghost; j < (lbgrid.n_y_node - (lbgrid.noghost)); j++){
            
            
        

            for (int dv = 0; dv<21; dv++){
                lbgrid.Node(i,j,dv) = Feq_node[dv];
            }


            for (int dv = 0; dv<21; dv++){
                 lbgrid.Cell(i,j,dv) = Feq_node[dv];

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
                // if (lbgrid.Node(i,j,dv) != 0.0){std::cout<<i<<" "<<j<<" "<<dv<<" "<<"Node"<<std::endl;}
                // if (lbgrid.Cell(i,j,dv) != 0.0){std::cout<<j<<" "<<j<<" "<<dv<<" "<< "Cell"<<std::endl;}
            }
        }}
std::cout<<"   "<<a<<std::endl;

}


void get_coordinates(std::vector<double> &x,std::vector<double> &y, const std::string& filename){

    int j = 0;
   
    std::string a1,a2;
    std::ifstream my_file;
    my_file.open(filename);
    

       // string line;
       //f getline(my_file, line);
        while(!my_file.eof())
        {
             getline(my_file, a1, ',');
             
             x.push_back(stod(a1));  ///////////////check whether there is no blank line at the end
            
             getline(my_file, a2, '\n');
             y.push_back(stod(a2));       
            
             
             // std::cout<<x[j]<<" "<<j<<std::endl;
            j += 1;
            
        }     }

///////---------------works for solid between x axis  of 0 and 1-----------///
template<typename T1>
void Mark_solid(std::vector<std::vector<int>>& solid_node,std::vector<std::vector<int>>& solid_cell,const Grid_N_C_2D<T1> &lbgrid){
  std::ofstream Myfiles("values.txt");
  std::vector<double> X_U,Y_U,X_L,Y_L;
    get_coordinates(X_U,Y_U,"cylinder_u.csv");
    get_coordinates(X_L,Y_L,"cylinder_l.csv");
     double delx = 0.025;
     double y_inter;
     int pos_xr;
  
    for(int i = 0 + lbgrid.noghost; i < lbgrid.n_x_node - (lbgrid.noghost) ; i++){
        for(int j = 0 + lbgrid.noghost;j < lbgrid.n_y_node - (lbgrid.noghost) ; j++){
            

 
           
///marking the nodes
            double x_node, y_node;
            x_node = (i - 3) * delx - 1;
            y_node = (j - 3) * delx - 2.5;          /////////for aerofoil 1.5
            
            std::vector<double>::iterator lower;
            if (x_node >0.0 && x_node < 1.0 ){

                if (y_node >= 0.0){   //////////

                /////upper half plane
                lower = std::lower_bound(X_U.begin(), X_U.end(), x_node);
                
                pos_xr = (lower - X_U.begin());  ////position of next largest no
               
                y_inter= ((Y_U[pos_xr]-Y_U[pos_xr-1])/(X_U[pos_xr] - X_U[pos_xr-1]))*(x_node - X_U[pos_xr-1])+ Y_U[pos_xr-1];
                // std::cout<<y_inter<<" "<<y_node<<std::endl;
                if (y_node<y_inter){
                    
                                    //std::cout<<"yes"<<" "<<y_node<<std::endl;

                    //add the i,j to the solid list 
                    solid_node.push_back({i,j});
                 //Myfiles<<i<<std::endl;
                  Myfiles<<i<<","<<j<<std::endl;
                //   Myfiles<<x_node<<","<<y_node<<std::endl;

                }}

                if (y_node<0.0){
                ////lower plane
                lower = std::lower_bound(X_L.begin(), X_L.end(), x_node);

                pos_xr = (lower - X_L.begin());  ////position of next largest no

                y_inter= ((Y_L[pos_xr]-Y_L[pos_xr-1])/(X_L[pos_xr] - X_L[pos_xr-1]))*(x_node - X_L[pos_xr-1]) + Y_L[pos_xr-1];
                            
                if (y_node>y_inter){
                    //add the i,j to the solid list 
                    solid_node.push_back({i,j});
                        // Myfiles<<x_node<<","<<y_node<<std::endl;
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
                 //std::cout<<y_inter<<" "<<y_cell<<std::endl;
                if (y_cell<y_inter){
                    
                    
                    //add the i,j to the solid list 
                    solid_cell.push_back({i,j});
                    //   Myfiles<<i<<","<<j<<std::endl;
                    //  Myfiles<<x_node<<","<<y_node<<std::endl;

                }}

                if (y_cell<0.0){
                ///////////----------------------------------------lower plane--------------------------------------------/////////////////
                lower = std::lower_bound(X_L.begin(), X_L.end(), x_cell);

                pos_xr = (lower - X_L.begin());  ////position of next largest no

                y_inter= ((Y_L[pos_xr]-Y_L[pos_xr-1])/(X_L[pos_xr] - X_L[pos_xr-1]))*(x_cell - X_L[pos_xr-1]) + Y_L[pos_xr-1];
                            
                if (y_cell>y_inter){
                    //add the i,j to the solid list 
                    solid_cell.push_back({i,j});
                        // Myfiles<<i<<","<<j<<std::endl;
                        //  Myfiles<<x_node<<","<<y_node<<std::endl;
                }}
            }           
}}}






// template<typename T>
// void Mark_bb_velocities(std::vector<std::vector<int>> &node_vel, std::vector<std::vector<int>> &cell_vel,
//                        std::vector<std::vector<int>> &solid_node, std::vector<std::vector<int>> &solid_cell,
//                          Grid_N_C_2D<T> &lbgrid ){
    
//     std::vector<int> line;
//     bool isPresent;


//     for(int i = 35+ lbgrid.noghost; i < 90  ; i++){    ///standar run 490 to -1490 for the 4412 aerofoil
//     for(int j = 70+ lbgrid.noghost; j < 130; j++){      /// 730 to -190


        
        
        
//         //////////---------------------------------------FOR NODES---------------------//////////////////
//         line = {i,j}; //////checking whether the node is fluid        
//         isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end(); 
//         if(!isPresent){ 
//             ///enters for the non solid nodes
            
           
           
//             line = {i + 1,j};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//            std::cout<<isPresent<<" "<<i<<std::endl;
//             if(isPresent){
//             //    std::cout<<"yes"<<std::endl;
//                 //--------sc1
//                 node_vel.push_back({i  ,j,dV_M1_ZERO});  //3 -
                
//                 //------sc3
//                 node_vel.push_back({i  ,j, 5} );  ///SC3 markers, near wall M3_Zero = 5, previous M3_zero = 6,previous 7;
//                 node_vel.push_back({i-1,j, 6} );
//                 node_vel.push_back({i-2,j, 7} );
            
//             }

//             line = {i -1,j};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if(isPresent){

//                 ///__________sc1
//                 node_vel.push_back({i  ,j,dV_P1_ZERO});

//                 //-----------SC3
//                 node_vel.push_back({i  ,j, 8  }); ///sc3 markers, near wal p3zero = 8, next right = 9, next right = 10
//                 node_vel.push_back({i+1,j, 9  });
//                 node_vel.push_back({i+2,j, 10 });
//             }

//             line = {i, j+1};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if(isPresent){

//                 //--------------------SC1
//                 node_vel.push_back({i,j  ,dV_ZERO_M1});
                
//                 //--------------------SC3
//                 node_vel.push_back({i,j  , 11} );  ///
//                 node_vel.push_back({i,j-1, 12} );
//                 node_vel.push_back({i,j-2, 13} );
//             }

//             line = {i ,j-1};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if (isPresent){
//                 //-----------------sc1
//                 node_vel.push_back({i,j  ,dV_ZERO_P1});


//                 //---------------sc3
//                 node_vel.push_back({i,j  , 14 });
//                 node_vel.push_back({i,j+1, 15 });
//                 node_vel.push_back({i,j+2, 16 });
//             }

//             ///////////////--------------------------SC1/2----------------------------------///////////
//             line = {i,j};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(),line) != solid_cell.end();
//             if(isPresent){
//                 //---SC1/2
//                 node_vel.push_back({i, j , 17});
                
//             }
//             line = {i -1 , j-1};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(),line) != solid_cell.end();
//             if(isPresent){
//                 //----SC1/2
//                 node_vel.push_back({i,j,18});

//             }
//             line = {i-1 , j};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(),line) != solid_cell.end();
//             if(isPresent){
//                 node_vel.push_back({i,j, 19});
//             }

//             line = {i,j-1};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(),line) != solid_cell.end();
//             if(isPresent){
//                 node_vel.push_back({i,j,20});

//             }




       

            


        
       
       
//         }

//         ////////////--------------------------------------FOR CELLS------------------------//////////////////
         
//         line = {i,j}; //////checking whether the cell is fluid        
//         isPresent = std::find(solid_cell.begin(), solid_cell.end(), line) != solid_cell.end(); 
//         if(!isPresent){ 
//             ///enters for the non solid cells
           
           
//             line = {i + 1,j};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(), line) != solid_cell.end();
//             if(isPresent){
//                 //--------sc1
//                 cell_vel.push_back({i  ,j,dV_M1_ZERO});  //3 -
                
//                 //------sc3
//                 cell_vel.push_back({i  ,j, 5} );  ///SC3 markers, near wall M3_Zero = 5, previous M3_zero = 6,previous 7;
//                 cell_vel.push_back({i-1,j, 6} );
//                 cell_vel.push_back({i-2,j, 7} );
            
//             }

//             line = {i -1,j};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(), line) != solid_cell.end();
//             if(isPresent){

//                 ///__________sc1
//                 cell_vel.push_back({i  ,j,dV_P1_ZERO});

//                 //-----------SC3
//                 cell_vel.push_back({i  ,j, 8  }); ///sc3 markers, near wal p3zero = 8, next right = 9, next right = 10
//                 cell_vel.push_back({i+1,j, 9  });
//                 cell_vel.push_back({i+2,j, 10 });
//             }

//             line = {i, j+1};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(), line) != solid_cell.end();
//             if(isPresent){

//                 //--------------------SC1----------------------------//
//                 cell_vel.push_back({i,j  ,dV_ZERO_M1});
                
//                 //--------------------SC3----------------------------//
//                 cell_vel.push_back({i,j  , 11 });  ///
//                 cell_vel.push_back({i,j-1, 12 });
//                 cell_vel.push_back({i,j-2, 13 });
//             }

//             line = {i ,j-1};
//             isPresent = std::find(solid_cell.begin(), solid_cell.end(), line) != solid_cell.end();
//             if (isPresent){
//                 //-----------------sc1
//                 cell_vel.push_back({i,j  ,dV_ZERO_P1});


//                 //---------------sc3
//                 cell_vel.push_back({i,j  , 14} );
//                 cell_vel.push_back({i,j+1, 15} );
//                 cell_vel.push_back({i,j+2, 16} );
//             }

//             //////////-----------------------------SC1/2----------------------------------//////////////////





//             line = {i+1,j+1};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if (isPresent){
//                 cell_vel.push_back({i,j,17});
//             }

//             line ={i,j};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if (isPresent){
//                 cell_vel.push_back({i,j,18});
//             }
            
//             line = {i,j+1};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if (isPresent){
//                 cell_vel.push_back({i,j,19});
//             }

//             line = {i+1,j};
//             isPresent = std::find(solid_node.begin(), solid_node.end(), line) != solid_node.end();
//             if (isPresent){
//                 cell_vel.push_back({i,j,20});
//             }



//     }}


//     }






// }


// ///////node _vel has the components that needs to be changed
// template<typename T>
// void bb_vel(std::vector<std::vector<int>> &Node_vel, std::vector<std::vector<int>> &Cell_vel,
//                        std::vector<std::vector<int>> &solid_node, std::vector<std::vector<int>> &solid_cell,
//                          Grid_N_C_2D<T> &lbgrid)
// {

//         for(int i = 0; i<Node_vel.size();i++){
//             ///--------------------------------------sc1
//             if(Node_vel[i][2]==1) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], 1) = lbgrid.Node(Node_vel[i][0]-1 ,Node_vel[i][1]    ,3);       }

//             if(Node_vel[i][2]==2) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], 2) = lbgrid.Node(Node_vel[i][0]   , Node_vel[i][1]-1 ,4);       }

//             if(Node_vel[i][2]==3) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], 3) = lbgrid.Node(Node_vel[i][0]+1 ,Node_vel[i][1]    ,1);        }

//             if(Node_vel[i][2]==4) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], 4) = lbgrid.Node(Node_vel[i][0]   , Node_vel[i][1]+1 ,2);       }
           
//            ////////////////---------------------------sc3
//             if(Node_vel[i][2]==5) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_M3_ZERO ) = lbgrid.Node(Node_vel[i][0]+1 ,Node_vel[i][1], dV_P3_ZERO);  }
//             if(Node_vel[i][2]==6) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_M3_ZERO ) = lbgrid.Node(Node_vel[i][0]+3 ,Node_vel[i][1], dV_P3_ZERO);  }
//             if(Node_vel[i][2]==7) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_M3_ZERO ) = lbgrid.Node(Node_vel[i][0]+5 ,Node_vel[i][1], dV_P3_ZERO);        }

//             if(Node_vel[i][2]==8) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_P3_ZERO ) = lbgrid.Node(Node_vel[i][0]-1 ,Node_vel[i][1], dV_M3_ZERO);          }
//             if(Node_vel[i][2]==9) { lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_P3_ZERO ) = lbgrid.Node(Node_vel[i][0]-3 ,Node_vel[i][1], dV_M3_ZERO);      }
//             if(Node_vel[i][2]==10){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_P3_ZERO ) = lbgrid.Node(Node_vel[i][0]-5 ,Node_vel[i][1], dV_M3_ZERO);      }

//             if(Node_vel[i][2]==11){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_ZERO_M3 ) = lbgrid.Node(Node_vel[i][0]   ,Node_vel[i][1] +1, dV_ZERO_P3);         }
//             if(Node_vel[i][2]==12){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_ZERO_M3 ) = lbgrid.Node(Node_vel[i][0]   ,Node_vel[i][1] +3, dV_ZERO_P3);         }
//             if(Node_vel[i][2]==13){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_ZERO_M3 ) = lbgrid.Node(Node_vel[i][0]   ,Node_vel[i][1] +5, dV_ZERO_P3);         }

//             if(Node_vel[i][2]==14){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_ZERO_P3 ) = lbgrid.Node(Node_vel[i][0]   ,Node_vel[i][1] -1, dV_ZERO_M3);        }
//             if(Node_vel[i][2]==15){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_ZERO_P3 ) = lbgrid.Node(Node_vel[i][0]   ,Node_vel[i][1] -3, dV_ZERO_M3);         }
//             if(Node_vel[i][2]==16){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_ZERO_P3 ) = lbgrid.Node(Node_vel[i][0]   ,Node_vel[i][1] -5, dV_ZERO_M3);         }
           
           
//            //////////////////-------------------sc1/2
//             if(Node_vel[i][2]==17){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_MH1_MH1 ) = lbgrid.Cell(Node_vel[i][0]  ,Node_vel[i][1]  , dV_PH1_PH1);         }
//             if(Node_vel[i][2]==18){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_PH1_PH1 ) = lbgrid.Cell(Node_vel[i][0]-1,Node_vel[i][1]-1, dV_MH1_MH1);        }
//             if(Node_vel[i][2]==19){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_MH1_PH1 ) = lbgrid.Cell(Node_vel[i][0]  ,Node_vel[i][1]-1, dV_PH1_MH1);        }
//             if(Node_vel[i][2]==20){ lbgrid.Node(Node_vel[i][0] ,Node_vel[i][1], dV_PH1_MH1 ) = lbgrid.Cell(Node_vel[i][0]-1,Node_vel[i][1]  , dV_MH1_PH1);       }



//         }



//   for(int i = 0; i<Cell_vel.size();i++){
//             ///--------------------------------------sc1
//             if(Cell_vel[i][2]==1) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], 1) = lbgrid.Cell(Cell_vel[i][0]-1 ,Cell_vel[i][1]    ,3);       }

//             if(Cell_vel[i][2]==2) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], 2) = lbgrid.Cell(Cell_vel[i][0]   , Cell_vel[i][1]-1 ,4);       }

//             if(Cell_vel[i][2]==3) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], 3) = lbgrid.Cell(Cell_vel[i][0]+1 ,Cell_vel[i][1]    ,1);        }

//             if(Cell_vel[i][2]==4) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], 4) = lbgrid.Cell(Cell_vel[i][0]   , Cell_vel[i][1]+1 ,2);       }
           
//            ////////////////---------------------------sc3
//             if(Cell_vel[i][2]==5) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_M3_ZERO ) = lbgrid.Cell(Cell_vel[i][0]+1 ,Cell_vel[i][1], dV_P3_ZERO);  }
//             if(Cell_vel[i][2]==6) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_M3_ZERO ) = lbgrid.Cell(Cell_vel[i][0]+3 ,Cell_vel[i][1], dV_P3_ZERO);  }
//             if(Cell_vel[i][2]==7) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_M3_ZERO ) = lbgrid.Cell(Cell_vel[i][0]+5 ,Cell_vel[i][1], dV_P3_ZERO);        }

//             if(Cell_vel[i][2]==8) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_P3_ZERO ) = lbgrid.Cell(Cell_vel[i][0]-1 ,Cell_vel[i][1], dV_M3_ZERO);          }
//             if(Cell_vel[i][2]==9) { lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_P3_ZERO ) = lbgrid.Cell(Cell_vel[i][0]-3 ,Cell_vel[i][1], dV_M3_ZERO);      }
//             if(Cell_vel[i][2]==10){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_P3_ZERO ) = lbgrid.Cell(Cell_vel[i][0]-5 ,Cell_vel[i][1], dV_M3_ZERO);      }

//             if(Cell_vel[i][2]==11){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_ZERO_M3 ) = lbgrid.Cell(Cell_vel[i][0]   ,Cell_vel[i][1] +1, dV_ZERO_P3);         }
//             if(Cell_vel[i][2]==12){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_ZERO_M3 ) = lbgrid.Cell(Cell_vel[i][0]   ,Cell_vel[i][1] +3, dV_ZERO_P3);         }
//             if(Cell_vel[i][2]==13){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_ZERO_M3 ) = lbgrid.Cell(Cell_vel[i][0]   ,Cell_vel[i][1] +5, dV_ZERO_P3);         }

//             if(Cell_vel[i][2]==14){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_ZERO_P3 ) = lbgrid.Cell(Cell_vel[i][0]   ,Cell_vel[i][1] -1, dV_ZERO_M3);        }
//             if(Cell_vel[i][2]==15){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_ZERO_P3 ) = lbgrid.Cell(Cell_vel[i][0]   ,Cell_vel[i][1] -3, dV_ZERO_M3);         }
//             if(Cell_vel[i][2]==16){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_ZERO_P3 ) = lbgrid.Cell(Cell_vel[i][0]   ,Cell_vel[i][1] -5, dV_ZERO_M3);         }



//             //-----------SC1/2-------------------------//
//             if(Cell_vel[i][2]==17){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_MH1_MH1 ) = lbgrid.Node(Cell_vel[i][0]+1 , Cell_vel[i][1]+1   ,  dV_PH1_PH1 );}
//             if(Cell_vel[i][2]==18){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_PH1_PH1 ) = lbgrid.Node(Cell_vel[i][0]   , Cell_vel[i][1]     ,  dV_MH1_MH1 );}
//             if(Cell_vel[i][2]==19){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_MH1_PH1 ) = lbgrid.Node(Cell_vel[i][0]+1 , Cell_vel[i][1]     ,  dV_PH1_MH1 );}
//             if(Cell_vel[i][2]==20){ lbgrid.Cell(Cell_vel[i][0] ,Cell_vel[i][1], dV_PH1_MH1 ) = lbgrid.Node(Cell_vel[i][0]   , Cell_vel[i][1] +1  ,  dV_MH1_PH1 );}





//         }













































//    }