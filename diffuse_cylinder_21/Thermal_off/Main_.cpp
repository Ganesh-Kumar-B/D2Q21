
#include<iostream>
#include<fstream>
#include<cmath>
#include"collison.h"
#include"advection.h"
#include<iomanip>
#include<sstream>
#include"force.h"
#include"print.h"



int main()
{
    
    int Ny = 250 ; int Nx = 5*Ny;   
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
    
    double Re = 10;
    std::cout<<"Re "<<Re<<std::endl;
    double L = 0.05*Ny;
    double Kn = 0.001;
    double Ma = 0.02;
    std::cout<<"Ma "<< Ma <<std::endl;
    double u0 = Ma * cs;
    std::cout<<"u0 = "<<u0<<std::endl;

    double Kin_Vis = u0*L/Re;
    double tau = Kin_Vis/(cs*cs);
    std::cout<<"Kinematic Viscosity "<<Kin_Vis<<std::endl;
    std::cout<<"tau "<<tau<<std::endl;

    double beta = 1.0/(2.0*tau + 1);
    std::cout<<beta<<std::endl;
    


    double jx_in, jy_in, jx_out, jy_out, u_free;


    double theta_bottom = d2q21.theta0 /*+ 0.3 * d2q21.theta0*/; 
    double theta_top = d2q21.theta0; 



    std::cout<<"theta_bottom "<<theta_bottom<<std::endl;
    std::cout<<"theta_top "<<theta_top<<std::endl;
    double a;

    double  Rho_avg_fluid,Fx_dt,Fy_dt,Fx,Fy,ux_stream,Fx_SI,Fy_SI,
            Jx  ,Jy  ,Jx_dt  ,Jy_dt  ,Fx_in  ,Fy_in, Fx_boundary,Fy_boundary,mass  ;



    int ic = Nx/3.3333; int jc = Ny/2;

    int shift = 5;
    int left    = ic - 0.5*L - shift,
        right   = ic + 0.5*L + shift,
        bottom  = jc - 0.5*L - shift,
        top     = jc + 0.5*L + shift;

    std::cout<<ic<<" "<<jc<<std::endl;
    
    // coordinates_writer(marker,ic,jc,L);  exit(0); //note for marking the aerofoil



    Mark_for_cylinder(marker,ic,jc,L); 
    // Mark_for_solid(marker, ic, jc, L); 
    




    G_calc(grid,d2q21,G_num,marker);

    


    
    //------------------------------Main code--------------------------//
    initialization(grid,d2q21,0.0);
    printMass(grid);
    printdata(d2q21,grid,0,u0);





    int sim_time = 5*Nx/u0; double time;
    int step_large = 2000;


      std::cout<<"Simulation time "<< sim_time<<std::endl;
      for(int t = 1; t <= sim_time;t++){

                  if(t%step_large == 0)                                        //DRTT
                  J_domain(grid,d2q21, Jx, Jy, left, right, top, bottom);  //DRTT



        collide(grid,d2q21,beta,tau);
        // periodic_21(grid);
        // bounce_back(grid);

        // periodic_x(grid,d2q21);
        periodic_y(grid,d2q21);





        // slip_wall_bb(grid,d2q21);


        // diffuse_B_21_new(grid,d2q21,0.0,u0,theta_top, theta_bottom);
        // diffuse_B_21(grid,d2q21,0.0,u0,theta_bottom,theta_top);      //with respect to incoming populations
        // diffuse_B_21_out(grid,d2q21,u0);   //with respect to outgoing populations
        


              if(t%step_large == 0)                                                                       //MEA
              mom_in(grid,d2q21,marker, jx_in, jy_in, left, right, top, bottom);                          //MEA 


              if(t%step_large == 0)                                                                          //DRTT           
              Force_calc(grid,d2q21,left,right, top, bottom, Fx_in, Fy_in, Fx_boundary, Fy_boundary);  //DRTT


        // diffuse_B_21_Solid(grid,d2q21,G_num,marker,left, right, bottom, top);        //diffuse bounce on a solid


        grad_inlet(grid,d2q21,u0);
        grad_outlet(grid,d2q21);

        advection_21(grid);
        // stationary_correction(grid,d2q21);



        BB_21(grid,d2q21,marker, left, right, top, bottom);                 //bounce back for the solid 

        if(t%step_large == 0)
        Force_SI_cylinder1(grid,marker,d2q21,left,right,top,bottom,ic,jc,Fx_SI,Fy_SI,beta);



                if(t%step_large == 0)                                                                               //MEA
                mom_out(grid,d2q21,marker, jx_out, jy_out, left, right, top, bottom);                               //MEA



        equi_correction(grid,marker,d2q21,left,right,bottom,top);





                  if(t%step_large== 0){                                                                 //DRTT                             
                      J_domain (grid, d2q21, Jx_dt, Jy_dt, left, right, top, bottom);                   //DRTT                                                                                 
                                                                                                        //DRTT      
                      Fx = Jx - Jx_dt + Fx_in -Fx_boundary;                                             //DRTT                                                        
                      Fy = Jy - Jy_dt + Fy_in -Fy_boundary;                                             //DRTT                                                       
                                                                                                        //DRTT                  
                                                                                                        //DRTT                                          
                      u_free_stream(grid,marker,d2q21, u_free,Rho_avg_fluid);                                                //DRTT
                      double cd_factor = (2.0 / (u_free*u_free*Rho_avg_fluid*L)) ;

                                     
                      std::cout<<" DRTT Cd = "<<Fx *(cd_factor)<<" Cl = "<<( Fy)*(cd_factor)<<std::endl;       //DRTT                                                                                              
                  
                  } 




              
                    time = t*Kin_Vis/(L*L);
                    
        
                      if(t %step_large == 0){
                      std::cout<<" time "<<t<<" "<<time;

                      u_free_stream(grid,marker,d2q21, u_free,Rho_avg_fluid);
                      double cd_factor = (2.0 / (u_free*u_free*Rho_avg_fluid*L)) ;



                      std::cout<<"ufree_stream"<<u_free<<"   rho_fluid  "<<Rho_avg_fluid<<std::endl;                                                                                        //MEA
                      std::cout<<" MEcd = "<< (jx_in - jx_out) *(cd_factor)<<" cl = "<<(jy_in - jy_out)*(cd_factor)<<std::endl;    //MEA
                      std::cout<<" SI  Cd = "<<Fx_SI *(cd_factor)<<" Cl = "<<( Fy_SI)*(cd_factor)<<std::endl;       //DRTT                                                                                              

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