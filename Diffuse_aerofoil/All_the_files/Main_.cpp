
#include<iostream>
#include<fstream>
#include<cmath>
#include"collison.h"
#include"advection.h"
#include<iomanip>
#include<sstream>
#include"force.h"



// template<typename T, typename T1>
// void printVtk(lbmD2Q13<T1> &lbModel,  Grid_N_C_2D<T> &gridLB,  int step)
// {
//   T u1,u2,u3, rho;
//   std::ofstream file;
//   char fileName[250];
//   sprintf(fileName,"./result/velocity_%d.vtk",step) ;
//   file.open(fileName) ;

//      //   vtk file header
//   file<<"# vtk DataFile Version 3.0"<<std::endl<<"Velocity"<<std::endl;
//   file<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
//   file<<"DIMENSIONS "<<gridLB.n_x*2<<" "<<gridLB.n_y*2<<" "<<1<<std::endl;
//   file<<"POINTS "<< gridLB.n_x*gridLB.n_y*2 <<" double"<<std::endl;

//   for(int i = 0 + gridLB.noghost; i < gridLB.n_x_node - (gridLB.noghost); i++) 
// 	for (int j = 0 + gridLB.noghost; j < gridLB.n_y_node - (gridLB.noghost); j++)
// 		{
//          file<<i-3<<" "<< j-3<<" "<<0<< "\n" ;         //for Nodes
//        file<<i+0.5-3<<" "<<j + 0.5-3<<" "<<0<<"\n";       //cell
//         }
  
//   file<<"POINT_DATA "<< gridLB.n_x*gridLB.n_y*2 <<std::endl ;
//   file<<"VECTORS"<<" "<<"velocity"<<" "<<"double"<<std::endl ;

// for(int i = 0 + gridLB.noghost; i < gridLB.n_x_node - (gridLB.noghost); i++) 
// 	for (int j = 0 + gridLB.noghost; j < gridLB.n_y_node - (gridLB.noghost); j++){ 

//     get_moments(gridLB,lbModel,u1, u2, rho, i,j,0);
//     if(fabs(u1) < 0.000000000001){u1 = 0; }
//     if(fabs(u2) < 0.000000000001){u2 = 0; }
// 	  file << u1 << " " <<u2<<" "<<0.0<< std::endl ;
   
//      get_moments(gridLB,lbModel,u1,u2, rho, i,j,1);
//       if(fabs(u1) < 0.000000001){u1 = 0;}
//     if(fabs(u2) < 0.000000001){u2 = 0;}
//      file << u1 << " " << u2  <<" "<<0.0<< std::endl ;
//     }
//   file.close();
// }



template<typename T, typename T1>
void printdata(lbmD2Q21<T1> &lbModel,  Grid_N_C_2D<T> &gridLB,  int step, double u0)
{ 
  double u_inv,nx_inv, ny_inv;
  u_inv = 1/ u0;
  
  nx_inv = 1/(double)gridLB.n_x;
  ny_inv = 1/(double)gridLB.n_y;

  std::vector<int> line;bool isPresent;
  T u1,u2,u3,u4,um, rho1,rho2,theta1,theta2,del=0.05;
  std::ofstream file;
  char fileName[250];
  
  sprintf(fileName,"./Re_1000_Ma_0.075/velocity_%d.txt",step) ;
  
  file.open(fileName) ;
  file<<"x,y,ux,uy,rho,theta"<<std::endl;
  for(int i = 0 + gridLB.noghost; i < gridLB.n_x_node - (gridLB.noghost); i++) 
	for (int j = 0 + gridLB.noghost; j < gridLB.n_y_node - (gridLB.noghost); j++)
		{{

    get_moments_node(gridLB,lbModel,u1, u2, rho1,theta1, i,j);

	  

    get_moments_cell(gridLB,lbModel,u3,u4, rho2,theta2, i,j);


      file<<(double)(i-gridLB.noghost)*ny_inv  <<"," <<(double)(j-gridLB.noghost) *ny_inv   <<","<<u1*u_inv<<","<<u2*u_inv<<","<<rho1<<","<<theta1<<std::endl;
      file<<(double)(i+0.5-gridLB.noghost)*ny_inv <<","<<(double)(j + 0.5-gridLB.noghost) *ny_inv <<","<<u3* u_inv<<","<<u4 * u_inv<<","<<rho2<<","<<theta2<<std::endl;



}
}
}










int main()
{
    std::vector<std::vector<int>> Solid_Nodes,Solid_Cell;         /////////////try to use typedef next time
    std::vector<std::vector<int>> bb_node,bb_cell;
    int Nx = 800;int Ny =300;
    Grid_N_C_2D<double> grid(Nx,Ny,2,21);
    Grid_N_C_2D<double> G_num(Nx,Ny,2,1);
    Grid_N_C_2D<bool> marker(Nx,Ny,2,1);

    Grid_N_C_2D<double> G_deno(Nx,Ny,2,21);



    lbmD2Q21<double> d2q21(1.0);
    // double beta = 0.53;
    // double tau =0.5*(1/beta-1);
    // double Re=100.0;
    double cs = sqrt(d2q21.theta0);
    // std::cout<<cs<<std::endl;
    // double u0 = 0.05   *cs;
    // double gx = 8.0*u0*u0/(40*Re);
    // double gy = 0.0*u0*u0/(40*Re);
    
    double Re = 1000;
    std::cout<<"Re "<<Re<<std::endl;
    double L = 40;
    double Kn = 0.001;
    double Ma = 0.04;
    std::cout<<"Ma "<< Ma <<std::endl;
    double u0 = Ma * cs;
    std::cout<<"u0 = "<<u0<<std::endl;

    double Kin_Vis = u0*L/Re;
    double tau = Kin_Vis/(cs*cs);
    std::cout<<"Kinematic Viscosity "<<Kin_Vis<<std::endl;
    std::cout<<"tau "<<tau<<std::endl;

    double beta = 1.0/(2.0*tau + 1);
    std::cout<<beta<<std::endl;
    double gx = 0.0*u0*u0/(L*Re);
    double gy = 0.0*u0*u0/(L*Re);


    double jx_in, jy_in, jx_out, jy_out, u_free;


    double theta_bottom = d2q21.theta0 /*+ 0.3 * d2q21.theta0*/; 
    double theta_top = d2q21.theta0; 



    std::cout<<"theta_bottom "<<theta_bottom<<std::endl;
    std::cout<<"theta_top "<<theta_top<<std::endl;
    double a;

    double Rho_avg_fluid,Fx_dt,Fy_dt,Fx,Fy,ux_stream,
        Jx  ,Jy  ,Jx_dt  ,Jy_dt  ,Fx_in  ,Fy_in, Fx_boundary,Fy_boundary,mass  ;


    coordinates_writer(marker);  exit(0); //note for marking the aerofoil



    Mark_for_solid(marker); 
    

    G_calc(grid,d2q21,G_num,marker);

    int left = 225, right = 300, bottom =120 , top = 190; //2

    int step_large = 1000;
    
    //------------------------------Main code--------------------------//
    initialization(grid,d2q21,0.0);
    printMass(grid);
    printdata(d2q21,grid,0,u0);





    int sim_time = 10*Nx/u0; double time;



      std::cout<<"Simulation time "<< sim_time<<std::endl;
      for(int t = 1; t <= sim_time;t++){

        if(t%step_large == 0)                                        //DRTT
        J_domain(grid,d2q21, Jx, Jy, left, right, top, bottom);  //DRTT



        collide(grid,d2q21,beta,tau,gx,gy);
        periodic_21(grid);
        // bounce_back(grid);





        // slip_wall_bb(grid,d2q21);


        // diffuse_B_21_new(grid,d2q21,0.0,u0,theta_top, theta_bottom);
        // diffuse_B_21(grid,d2q21,0.0,u0,theta_bottom,theta_top);      //with respect to incoming populations
        // diffuse_B_21_out(grid,d2q21,u0);   //with respect to outgoing populations
        


        if(t%step_large == 0)                                                                       //MEA
        mom_in(grid,d2q21,marker, jx_in, jy_in, left, right, top, bottom);                          //MEA 


        if(t%step_large == 0)                                                                          //DRTT           
        Force_calc(grid,d2q21,left,right, top, bottom, Fx_in, Fy_in, Fx_boundary, Fy_boundary);  //DRTT


        diffuse_B_21_Solid(grid,d2q21,G_num,marker,left, right, bottom, top);        //diffuse bounce on a solid


        // grad_inlet(grid,d2q21,u0);
        // grad_outlet(grid,d2q21);

        advection_21(grid);
        // stationary_correction(grid,d2q21);

        // BB_21(grid,d2q21,marker, left, right, top, bottom);                 //bounce back for the solid 


        if(t%step_large == 0)                                                                               //MEA
        mom_out(grid,d2q21,marker, jx_out, jy_out, left, right, top, bottom);                               //MEA



        equi_correction(grid,marker,d2q21,left,right,bottom,top);



        //   // bb_vel(bb_node,bb_cell,Solid_Nodes,Solid_Cell,grid);
        // //  // bounce_back_wall(grid);



        if(t%step_large== 0){                                                                  //DRTT                             
            J_domain (grid, d2q21, Jx_dt, Jy_dt, left, right, top, bottom);                    //DRTT                                                                                 
                                                                                               //DRTT      
            Fx = Jx - Jx_dt + Fx_in -Fx_boundary;                                             //DRTT                                                        
            Fy = Jy - Jy_dt + Fy_in -Fy_boundary;                                              //DRTT                                                       
                                                                                               //DRTT                  
                                                                      //DRTT                                          
            u_free_stream(grid,d2q21, u_free);                                              //DRTT      
            // std::cout<< "Fx = "<<Jx - Jx_dt<< "Fy "<< Jy - Jy_dt<<std::endl;                    //DRTT                                                                                 
            // std::cout<< "Fx = "<<Fx_in -Fx_boundary<< "Fy "<< Fy_in -Fy_boundary<<std::endl;    //DRTT                                                                                                 
            std::cout<<" DRTT Cd = "<<Fx *(2/(u0*u0*L))<<" Cl = "<<( Fy)*(2/(u0*u0*L))<<std::endl;       //DRTT                                                                                              
        
        } 








        inlet(grid,   d2q21, u0,0,1.0);
        outlet(grid,  d2q21, u0,0,1.0);
    
      time = t*Kin_Vis/(L*L);
      //  std::cout<<t<<std::endl;
      //  if(time >=0.0){
        
        if(t %step_large == 0){
        std::cout<<" time "<<t<<" "<<time;

        u_free_stream(grid,d2q21, u_free);
        std::cout<<u_free<<std::endl;         
        
        std::cout<<"ufree_stream"<<u_free<<std::endl;                                                                                        //MEA
        std::cout<<" MEcd = "<< (jx_in - jx_out) *(2.0/(u0 * u0 *L))<<" cl = "<<(jy_in - jy_out)*(2.0/(u0 * u0 *L))<<std::endl;    //MEA

        printMass(grid);
        printdata(d2q21,grid,t,u0);
       }

    }
//------------------------------------------------------------------------//



//  for(int dv = 0; dv<17; dv++){
// for (int i = 0+grid.noghost; i<grid.n_x_node-grid.noghost;i++){
//                 for( int j = 0 +grid.noghost; j< grid.n_y_node -grid.noghost; j++){
                   
//                         if(grid.Node(i,j,dv) >0 || grid.Cell(i,j,dv) >0.0 ){
//                         std::cout<<grid.Node(i,j,dv)<<" Node "<<i<<" "<<j<<"-- "<<grid.Cell(i,j,dv)<<"cell  "<<i<<" "<<j<<" "<<dv<<std::endl;


//                         }
//                     }}}





std::cout<<std::endl;
    




   }
   


  //  initialization(grid,d3q19,0.0);
   
  //    grid.Cell(5,11,dV_ZERO_P3)=0.1;

  //  printVtk(d3q19,grid,0);

  //   for(int t = 1; t < 510;t++){std::cout<<t<<" ";
  //    collide(grid,d3q19,beta,tau,gx,gy);
        
       
  //       periodic(grid);
  //        bounce_back(grid);

  //       advection(grid);
  //      // bounce_back_wall(grid);
  //      printMass(grid);
  //      if(t%5 == 0){
  //     printVtk(d3q19,grid,t);
  //      }
  //   }


// for (int i =0; i<grid.n_x_node;i++){
//                 for( int j = 0; j< grid.n_y_node; j++){
//                     for(int dv = 0; dv<17; dv++){
//                         grid.Node(i,j,dv) = 0;
//                         grid.Cell(i,j,dv) = 0;
//                     }}}



// for (int i =0; i<grid.n_x_node;i++){
//                 for( int j = 0; j< grid.n_x_node; j++){
//                     for(int dv = 0; dv<13; dv++){
//                         if(grid.Node(i,j,dv) ==0.1 || grid.Cell(i,j,dv) == 0.1 ){
//                         std::cout<<grid.Node(i,j,dv)<<" Node "<<i<<" "<<j<<" "<<grid.Cell(i,j,dv)<<"cell  "<<i<<" "<<j<<" "<<dv<<std::endl;


//                         }
//                     }}}












    // for (int t = 0; t< 10; t++){
    //     std::cout<<t<<" ";
    //   // collide(grid,d3q19,beta,tau,gx,gy);
        
    //    periodic(grid);
    //    bounce_back_prep(grid);
    //    advection(grid);
    //    bounce_back_wall(grid);
    //     printMass(grid);
    //     if(t % 1 ==0){
    //     printVtk(d3q19,grid,t);
    //     }


    // }








//----------------------------------------------------------------writing to a file----------------------//////////////
// std::ofstream my_File("bb_vel_node.txt");

// for(int i = 0 ; i<bb_node.size();i++){
//     for(int j = 0; j<bb_node[i].size();j++){
//      if (j == 2){my_File<<bb_node[i][j]<<std::endl;}
//           else{my_File<<bb_node[i][j]<<",";}

// //     }}
// for(int i = 0 ; i<bb_node.size();i++){
//     for(int j = 0; j<bb_node[i].size();j++){
//      if (j == 2){my_File<<bb_node[i][j]<<std::endl;}
//           else{my_File<<bb_node[i][j]<<",";}

//     }}

//       std::ofstream my_File1("bb_vel_cell.txt");

// for(int i = 0 ; i<bb_cell.size();i++){
//     for(int j = 0; j<bb_cell[i].size();j++){
//      if (j == 2){my_File1<<bb_cell[i][j]<<std::endl;}
//           else{my_File1<<bb_cell[i][j]<<",";}

//     }}