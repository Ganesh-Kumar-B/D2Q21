#include <iostream>
using namespace std;
#include <cmath>
typedef double myreal;


myreal f(myreal x){

myreal f=x*x*x - x -1;

return 	f;
}




int main(){





myreal x0,x1,x2;
cout << " Enter 1st guess "<< endl;
cin>> x0;
cout << " Enter 2nd guess "<< endl;
cin>> x1;

myreal tol=1e-3;
myreal err=10;
myreal temp1;

int i=0;
while(err>tol){ 
	
//NOT USED IN SECANT METHOD
// /*CHECK:Make x0 to be for smaller f*/
// if(f(x0)<f(x1)){
// x0=x0;
// x1=x1;	
// }
// else
// {myreal temp;
// temp=x0;	
// x0=x1;
// x1=temp;
// }

//One more check:To see whether guesses gives values of functions of oppsite signs so 
//that a root lies between the 2 guesses
if(f(x0)<0)
{cout << "f(x0) is negative "<< endl;}
if(f(x0)>0)
{cout << "f(x0) is positive "<< endl;}
if(f(x1)<0)
{cout << "f(x1) is negative "<< endl;}
if(f(x1)>0)
{cout << "f(x1) is  positve "<< endl;}




x2= x1 -  f(x1)*( (x1-x0)/(f(x1)-f(x0)) );

cout << " x2 "<< x2 << endl;
cout << " x1 "<< x1 << endl;
cout << " x0 "<< x0 << endl;
err=fabs(x2-x1);

temp1=x1;
cout << "temp1 "<< temp1<< endl;
x1=x2;
x0=temp1;
cout << " x0 "<< x0 << endl;
i=i+1;

}

cout << "No. OF ITERATIONS is "<< i << endl;
cout << "The root is "<< x2 << endl;



return 0;


}