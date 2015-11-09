#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std;

const int na = -9999;
const int range = 126;	
const int N = 512;



//-------------------------initialize the matrix-----------------------//

void create_random_matrix (int* A, float x) {
	int i,j;
	for (i=0;i<N;i++) {
		for (j=0;j<N;j++){
			float number = rand() / (float)RAND_MAX;
			if (number < x) A[i*N+j] =1;
			else A[i*N+j] = 0;
		}
	}

}


//-----------------------mean cover-----------------------------------//

double density (double* A, int size) {
	
	int i,j;
	double count=0,sum=0;
	for (i=0;i<size;i++) {
		for (j=0;j<size;j++){	
			if (A[i*size+j] != na) {
				sum += A[i*size+j];
				count += 1;
			}
		}
	}
	
	double average = sum/count;
	//cout<<"mean="<<average<<endl;
return average;
}

//-----------------------coarse-graining-------------------------------//

void coarse_grain (int* A, double* A_cg, int n, int cgLength) {  // n is the size of coarse-grained matrix

	int i,j,k,l;
	
	for ( i=0; i<n; i++) {
		for ( j=0; j<n; j++) {
			
			int init_row = i*cgLength;
			int end_row = (i+1)*cgLength-1;
			int init_col = j*cgLength;
			int end_col = (j+1)*cgLength-1;
			
			double sum = 0, count=0;
			
			for ( k=init_row ; k<end_row+1; k++) { 
				for ( l=init_col; l<end_col+1; l++) {
					if (A[k*N+l]==0 || A[k*N+l]==1) { 
						sum += A[k*N+l];
						count += 1;
					}
				}
			}
			
			double reduced_mean = sum/count;
			A_cg[i*n+j] = reduced_mean;
			
		}
	}
}

//---------------------------------varaince--------------------------------------//

double spVar(double* A, int size, double mean) {

	double sum = 0,count = 0;
	for (int i=0; i< size; i++) {
		for (int j=0; j<size; j++) {
		
		if (A[i*size+j]!=na) {
			sum = sum + A[i*size+j]*A[i*size+j];
			count +=1;
		}
		}
	}
		
		
	double var = sum/count - mean*mean;
	//cout<<"var="<<var<<endl;
	return var;
}



//-------------------------------Correlation at lag1---------------------------------//

double spCor (double* A ,int size){

	double sum=0;
	int i,j;
	int down, right;
	
	for (i=0; i< size; i++) {
		for (j=0; j< size; j++) {
			
			down = (i+1)%size; 
			right = (j+1)%size; 
			
			sum = sum +  (A[i*size+j] * A[down*size+j]) + (A[i*size+j]* A[i*size+right]) ; 
			
		}
	}

	double corLag1 = sum/(double)(2*size*size);
	//cout<<"meansq="<<mean*mean<<"\t"<<"cor="<<corLag1<<endl;

return corLag1;

}


//--------------------------------------int to string converter----------------------------//
string IntToStr(int n){
    	stringstream result;
    	result << n;
    	return result.str();
}

						
//----------------------------------main function-----------------------------//

int main() {

	srand (time(NULL));
	float initial_density;
	double mean, correlation, variance; 
	int* mat = new int [N*N];
	int*full = new int[N*N*range];
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~reading data from a file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//This part reads the data from a file and saves it in the array "full"

	ifstream myfile;
	myfile.open ("tcp_rep2_snapshots_p=0.725:-0.001:0.6_q0_n512.dat");
		  
	if (myfile.is_open()){
		for (int x=0; x<range; x++){
			for (int i=0; i<N; i++){
		  		for (int j=0 ; j<N; j++) {
					myfile >> full[x*N*N+i*N+j];
				}
			}
		}
	}
	else cout << "Unable to open file";
	myfile.close();
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
	//Loop over different cg lengths
	
	for (int cgLength=1;cgLength<33;cgLength++) {
	
	//~~~~~~~~~~~~~ declare cg matrix~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
		int rem;
		if (N%cgLength ==0) rem = 0;
		else rem = N%cgLength ;
		int L = N-rem;
		int n = L/cgLength;
		
		double* mat_cg = new double[n*n];
		
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	
		// Loop over different 'p' values
		
		for (int x=0; x<range; x++){
		
			for (int i=0; i<N; i++){
		  		for (int j=0 ; j<N; j++) {
						mat[i*N+j] = full[x*N*N+i*N+j];
				}
			}

			coarse_grain (mat,mat_cg,n,cgLength);
		
			mean = density(mat_cg,n);
			correlation = spCor (mat_cg,n) - mean*mean;
			variance = spVar(mat_cg,n,mean);
		

			//~~~~~~~~~~~~~~~~~~~~~~~ save the correlation and variance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
			
			string filename1="tcp_q0_rep2_varianceCg" + IntToStr(cgLength)+ "_p=0.725:-0.001:0.6_n512.dat";
			string filename2="tcp_q0_rep2_correlationCg" + IntToStr(cgLength) + "_p=0.725:-0.001:0.6_n512.dat";
		    
		
			ofstream outdata1; 
			ofstream outdata2; 
			outdata1.open(filename1.c_str(),ios::app);
			outdata2.open(filename2.c_str(),ios::app);
			
			
			if ( outdata1.is_open()){
		 		
		 		outdata1<< variance<< endl;
		 		outdata2<<correlation<<endl;
				outdata1.close();
				outdata2.close();

			}
			else cout<< "oudata not open";
			
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
			
		  }
		  
		  delete [] mat_cg;
	}

	delete [] mat;
	delete [] full;
	 
	return 0;
	
}

