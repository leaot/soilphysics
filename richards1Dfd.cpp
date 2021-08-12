#include<iostream>
#include<cmath>

// Calculations soil physics book chapter 3 
// In august 2021 
// static parameters for Warwick hydraulic functions see Zarba 1988 MIT Masters 
	#define	thetas  0.287
	#define	thetar  0.075
	#define	alpha  1.611e06
	#define	beta1  3.96

	#define	Ks 34 
	#define	A  1.175e06
	#define	beta2  4.74


using namespace std;

// define the functions for water retention and hydraulic conductivity
double theta(double h){
	double theta;
	theta = alpha * (thetas-thetar)/(alpha+ pow(abs(h), beta1)) + thetar;
	return theta;
}

double K(double h){
	double K;
	K = Ks * A/(A + (pow(abs(h), beta2)));
	return K;
}

double C(double h){
	double C;
	C = -alpha*(thetas-thetar) * pow(abs(h), beta1-1) * beta1/((alpha+ pow(abs(h), beta1))* (alpha+ pow(abs(h), beta1)));
	return C;
}


int main(){

       // double dz, dt, thetas, thetar, alpha, beta1, beta2, Ks, A; // parameters
       // double theta, K; // functions	
	double dz, dt;
// time and space discretization parameters
	dz = 1;
	dt = .00005;
int n, i;
double h[120][65];

/**
double r = 100;
	cout << theta(r) << "\t" << K(r) << "\t" << C(r) << "\n";

**/


for (n = 1; n < 120; n++){
	h[n][0] = -20.73;
}

for (i = 0; i <65; i++){
	h[0][i] = -61.5;
}

/**
for (n = 0; n < 120; n++){
cout  <<  n << "\t" << h[n][1] <<  "\t" << h[n][65] << "\n";
}
**/



for (n=0; n < 119; n++){
	for (i=1; i <=65; i++){

h[n+1][i] = h[n][i] + dt*1/(C(h[n][i])*dz)   * (  (K(h[n][i+1])+K(h[n][i]) )/2 * ((h[n][i+1]-h[n][i])/dz  -1 )
-
  (K(h[n][i])+K(h[n][i-1]) )/2 * ((h[n][i]-h[n][i-1])/dz  -1 ));

	}


}


for (i = 0; i <= 65; i++){
cout  <<  i << "\t" << h[10][i] << "\n";
}

return 0;
}

