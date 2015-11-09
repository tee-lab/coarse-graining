	#include <iostream>
	#include <math.h>
	#include <fstream>
	#include <string>
	#include <cstdlib>
	using namespace std;

#define N 1024
#define range 91


//-------------------mean----------------------------------------------//
float density(int* A){
	int i,j;
	float sum=0,count=0; 
	
	for(i=0;i<N;++i){
		for(j=0;j<N;++j){
			if (A[i*N+j] ==0 || A[i*N+j]==1){
				sum += A[i*N+j];
				count += 1;
			}
		}
	}
	float average = sum/count;
return average;
}

//=============================main function=================================//	 
int main () {
	int*A = new int[N*N];
	int*full = new int[N*N*range];
	float mean_cover, rainfall;
	
	  ifstream myfile;
	  myfile.open ("random_corresponding_to_tcp_p=0.28600:-0.00001:0.28510_q0.92_n1024.dat");
	  
	  if (myfile.is_open() ) {
		  for (int p=0 ;p<range; p++){
		  	
		  	for (int i =0; i<N; ++i){
		  		for (int j=0 ; j<N; ++j) {
					myfile >> full[p*N*N+i*N+j];
					A[i*N+j] = full[p*N*N+i*N+j];
				}
			}
			
			mean_cover = density(A);
			
			ofstream outdata;
			outdata.open("random_mean_correspondingTo_p=0.28600:-0.00001:0.28510.dat",ios::app);
			
			if (outdata.is_open()){
		 	outdata<< mean_cover<< endl;
			cout<<mean_cover<<endl;
			outdata.close();
			}
			
			else cout<< "oudata not open";	 	
	  		}
	  myfile.close();

	  }
	  else cout << "Unable to open file";

delete[] A;
delete[] full;	 
return 0;
}
