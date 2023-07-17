#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
#include <cmath>
typedef double myreal;

class riemann{
public:
// p2=x;//Pressure at left of the shock
myreal p1=100000*0.5;//Right side or Low Pressure or PRESSURE TO THE RIGHT OF THE SHOCK
myreal p4=100000;//Left side or High Pressure


myreal l=800;//length of shock tube in mm

myreal gamma=1.4;
myreal R=287;

myreal T1=348.4320 ;
myreal T4=348.4320;

// myreal T1=300*0.8;
// myreal T4=300;

myreal u1=0;
myreal u4=0;

myreal rho1=0.5;
myreal rho4=1;

myreal a1=sqrt(gamma*p1/rho1);
myreal a4=sqrt(gamma*p4/rho4);



/*derived parameters*/
myreal k1=(gamma-1)/(2*gamma); 
myreal k2=(gamma+1)/(2*gamma);
myreal k3=(gamma+1)/(gamma-1);

myreal k4=(2*gamma)/(gamma-1);


myreal operator()(myreal x){

myreal temp= u4 - u1 - ( (a1/gamma)*(x/p1-1)/(sqrt(k2*(x/p1-1)+1)) );


// myreal f=pow(p4/p1,k1)-pow(x/p1,k1)*pow(( 1 + ((gamma-1)*temp)/(2*a4) ),-1);



myreal f=pow(p4/p1,k1)*( 1 + ((gamma-1)*temp)/(2*a4) ) -pow(x/p1,k1);

return f;
}



};








int main(){


riemann f;
cout << "a in right region " << f.a1 << endl;
cout << "a in left region " << f.a4 << endl;




myreal x0=0.05*f.p4;

// cout << x0 << endl;

myreal x1=0.5*f.p4;
// cout << x1 << endl;

myreal x2;



myreal tol=1e-5;
myreal err=10;
myreal temp1;

int i=0;
while(err>tol){ 

//One more check:To see whether guesses gives values of functions of oppsite signs so 
//that a root lies between the 2 guesses
// if(f(x0)<0)
// {cout << "f(x0) is negative "<< endl;}
// if(f(x0)>0)
// {cout << "f(x0) is positive "<< endl;}
// if(f(x1)<0)
// {cout << "f(x1) is negative "<< endl;}
// if(f(x1)>0)
// {cout << "f(x1) is  positve "<< endl;}




x2= x1 -  f(x1)*( (x1-x0)/(f(x1)-f(x0)) );

// cout << " x2 "<< x2 << endl;
// cout << " x1 "<< x1 << endl;
// cout << " x0 "<< x0 << endl;
err=fabs(x2-x1);

temp1=x1;
// cout << "temp1 "<< temp1<< endl;
x1=x2;
x0=temp1;
// cout << " x0 "<< x0 << endl;
i=i+1;

}

cout << "No. OF ITERATIONS "<< "for tolerance " << tol << " is "<< i << endl;
cout << "The root or value of pressure P2  is "<< x2 << endl;


myreal p2=x2;

myreal u2;
u2= f.u1 + (f.a1/f.gamma)*(p2/f.p1 -1)/(sqrt(f.k2* (p2/f.p1 -1) + 1 ));

cout << "        Surface Contact Velocity: u2  "<< u2 << endl;


// myreal T2= f.T1*(f.p1/p2)*((1+f.k3*p2/f.p1)/(f.k3+p2/f.p1));

myreal T2= f.T1*(p2/f.p1)*((f.k3+p2/f.p1) / (1+f.k3*p2/f.p1));

cout << " Temperature in region 2 is  "<< T2 << endl;
myreal rho2=p2/(f.R*T2);
cout << " Density in region 2 is      "<< rho2 << endl;

myreal a2=f.a1*sqrt( p2*(f.k3 + p2/f.p1) / ( f.p1*( 1 + f.k3*p2/f.p1 ) ) );
cout << " a in region 2 is            "<< a2 << endl;


myreal S= f.u1 + f.a1*(sqrt(f.k2* (p2/f.p1 -1) + 1 ));
cout << "        Shock Speed: S "<< S << endl;


myreal t=0.4;//time in milliseconds

myreal x21=S*t + f.l/2;
cout << "x21 is  "<< x21 <<" at "<< t <<   " milliseconds"<<  endl;

myreal x22dash=u2*t +f.l/2;
cout << "x22dash is "<< x22dash <<" at "<< t << " milliseconds"<<  endl;



myreal p3,p2dash;
p3=p2;
p2dash=p2;
cout << p2dash << endl;
cout << f.p4 << endl;
cout << f.rho4  <<endl;
myreal rho2dash=f.rho4*pow(p2dash/f.p4,1/f.gamma);
cout << "Density at tail of expansion is "<< rho2dash << endl;

myreal a2dash=sqrt(f.gamma*p2dash/rho2dash);
cout << "Speed of Sound at tail of expansion is "<< a2dash << endl;

myreal x43=-f.a4*t +f.l/2;
cout << "x43(Expansion Head) is  "<< x43 <<" at "<< t <<   " milliseconds"<<  endl;

myreal x32dash= (u2-a2dash)*t + f.l/2;
cout << "x32dash(Expansion Tail) is  "<< x32dash <<" at "<< t <<   " milliseconds"<<  endl;
cout << endl;

cout <<setw(20)<< "Postion "<< x43 << " "<< x32dash << " "<< x22dash << " "<< x21 << endl; 
cout <<setw(20)<< "Pressure "<< f.p4 <<" "<<  p2<< " " << p2 << " "<< f.p1 << endl;
cout <<setw(20)<< "Velocity "<< f.u4 << " "<< u2 << " "<< u2 << " "<< f.u1 << endl;
cout <<setw(20)<< "Density "<< f.rho4<< " "<< rho2dash << " "<< rho2 << " "<< f.rho1 << endl;
cout <<setw(20)<< "Sound speed "<< f.a4 << " " << a2dash << " "<< a2 << " " << f.a1 << endl;


int n=100;//no. of points to be interpolated
myreal u[n+1];
myreal a[n+1];
myreal p[n+1];
myreal rho[n+1];


myreal dx=(x32dash-x43)/n;
cout << "dx "<< dx << endl;

u[0]=f.u4;
a[0]=f.a4;
p[0]=f.p4;
rho[0]=f.rho4;

ofstream myfile("density.dat");
myfile << 0 << " "<< f.rho4 << endl;
myfile << x43 << " "<< f.rho4 << endl;


ofstream myfile1("properties.dat");
for(int i=1;i<=n;i++){

myreal temp=x43+ dx*i-f.l/2;
// cout << "temp "<< temp<< endl;
u[i]= (2.0/(f.gamma+1))*( temp/t + ((f.gamma-1)/2.0)*f.u4 + f.a4 );
a[i]= u[i]-(temp/t);
p[i]=f.p4*pow(a[i]/f.a4,f.k4);

rho[i]=f.rho4*pow(p[i]/f.p4,1/f.gamma);


// cout <<"u:"<<  u[i]<< " "<< "a:"<< a[i]<< " "<< "p:"<< p[i]<< " "<<"rho:"<<  rho[i] << endl;
myreal x=temp+f.l/2;
myfile1 << x << " "<<  u[i]<< " "<< a[i]<< " "<< p[i]<< " "<<  rho[i] << endl;

myfile << x << " "<< rho[i]<< endl;
}






myfile << x22dash << " "<< rho2dash << endl;
myfile << x22dash << " "<< rho2 << endl;
myfile << x21 << " "<< rho2 << endl;
myfile << x21 << " "<< f.rho1 << endl;
myfile << f.l << " "<< f.rho1 << endl;

myfile.close();

myfile1.close();


return 0;


}