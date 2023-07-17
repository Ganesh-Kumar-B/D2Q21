#pragma once

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
  
  sprintf(fileName,"./check/velocity_%d.txt",step) ;
  
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
;

