// Solution to the heat equation 1D using an explicit method
// For a general reference see: 
// Applied numerical methods Carhanhan, Luther, Wilkes 1969 John Wiley
// Explicit solution p.436
// Do not use to make money, share teach and give credit
// tpl 25 sept 2021


#include<iostream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>

int main(){

// define the simulation parameters

double dx = 0.05;  
double dtau = 0.001; 
double TMAX = 0.95; 
double L = 1.00;  
double lambda = dtau/pow(dx,2);

int  nx = L/dx+1;
int  ntau = TMAX/dtau+1;

int i, j;

// define the grid 

double u[ntau][nx];  

// initialize the array with zeros

for(i=0; i < nx; i++){
for(j=0; j < ntau; j++){
	 u[j][i] = 0.0;
}}


std::cout << "-------------------------------------" <<  std::endl;

// define the initial conditions

for(i=0; i < nx; i++){
	u[0][i] = 0.0; //sin(M_PI*i*dx/L);// 0.5;
}

for(j=1; j < ntau; j++){
u[j][0] = 1.0; //100
u[j][nx-1] = 1.0; //100
}

// this is the explicit solution 
for(j=0; j < ntau-1; j++){
for(i=1; i < nx-1; i++){
		 u[j+1][i] = lambda*u[j][i-1] + (1-2*lambda)*u[j][i] + lambda*u[j][i+1]; 
} } 
	
char const* m[ntau][nx];

// color code for the terminal plot
for(j=0; j < ntau; j++){
for(i=0; i < nx; i++){
		 if  (u[j][i] > 0.9 and u[j][i] <=1){
		 m[j][i] = "\033[1;31m*\033[0m";}
		 else if  (u[j][i] > 0.8 and u[j][i] <=0.9){
		 m[j][i] = "\033[1;32m+\033[0m";}
		 else if  (u[j][i] > 0.7 and u[j][i] <=0.8){
		 m[j][i] = "\033[1;33m-\033[0m";}
		 else if  (u[j][i] > 0.6 and u[j][i] <=0.7){
		 m[j][i] = "\033[1;34m#\033[0m";}
		 else if  (u[j][i] > 0.5 and u[j][i] <=0.6){
		 m[j][i] = "\033[1;35m%\033[0m";}
		 else if  (u[j][i] > 0.4 and u[j][i] <=0.5){
		 m[j][i] = "\033[1;36m&\033[0m";}
		 else if  (u[j][i] > 0.3 and u[j][i] <=0.4){
		 m[j][i] = "\033[1;37m@\033[0m";}
		 else if  (u[j][i] > 0.2 and u[j][i] <=0.3){
		 m[j][i] = "=";}
		 else if  (u[j][i] > 0.1 and u[j][i] <=0.2){
		 m[j][i] = "!";}
		 else if  (u[j][i] >= 0.0 and u[j][i] <=0.1){
		 m[j][i] = ".";}
		 else{
		 m[i][j] = "err";}
} } 

// terminal plot


std::cout << " Explicit solution to the 1D heat equation in time   ";
std::cout << "\n";
std::cout << " _____________________________________________   ";
std::cout << "\n";


for(j=0; j < ntau; j++){ std::cout << "|  "; 
for(i=0; i < nx; i++){
		 std::cout << m[j][i] << " " ;
		 
}std::cout << "  |          t = " << std::setprecision(3) << j*dtau; std::cout << std::endl;} 


// prints the output data in terminal in gnuplot format
// uncomment to see the outoput data
// I have not included the code for exporting as *.dat, it is pretty simple though
/*
std::cout << "\n";
std::cout << " Output data in gnuplot format   ";
std::cout << "\n";
std::cout << " _____________________________________________   ";
std::cout << "\n";

// choose your interval of time by replacing "ntau"
for(j=0; j < ntau; j++){ 
for(i=0; i < nx; i++){
		 std::cout <<  std::setprecision(4) << u[j][i] << std::endl;
	}std::cout <<  std::endl;

} 
*/
 
return 0;

}
