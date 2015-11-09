#include <iostream>
#include <math.h>
using std::cerr;
using namespace std;
#include <fstream>
using std::ofstream;
#include <cstdlib>       // for exit function
#include <ctime>


// parameter definitions
#define N 256    // system size
#define p_range 1000   // no. of p values in range 0-1      
#define T 6000      // No. of time steps
#define rep 1


//====== function to initialize the matrix ===============//
void create_random_matrix(int*U,float x) {  

	int i1, i2;
	
	for (i1 = 0; i1 < N; i1++) {
		for (i2 = 0; i2 < N; i2++) {
			float number = rand()/(float)RAND_MAX;
			if (number<x) U[i1*N+i2]=1;
			else U[i1*N+i2]=0;
		}
	}

}


//=======================average density=========================//
float density(int* A){
	int i,j;
	float sum=0; 
	for(i=0;i<N;++i){
		for(j=0;j<N;++j){
			sum+= A[i*N+j];
		}
	}
	float average= sum/(N*N);
return average;
}

//=========================average along rows======================//
void average_along_row(float* A,float*A1){
	int l,k;
	float sum=0;
	
	for(k=0;k<p_range;++k){
		sum=0;
		for(l=0;l<rep;++l){
			sum+=A[l*p_range+k];
		}
		A1[k]=sum/rep;
	}

}	

//====================select random site==============================//
void select_neighbor_of_site(int i ,int j , int*neighbor){
	int left,right,bottom,top ;
	int in = i ,jn = j;

	float test = rand()/(float)RAND_MAX ;
	
	if (i==0) top=N-1;
	else top = i-1;
	
	if (i==N-1) bottom = 0;
	else bottom = i+1 ;
	
	if (j==0) left = N-1;
	else left = j-1 ;
	
	if (j==N-1) right = 0;
	else right = j+1;
				
	if (test <= 0.25) in = top;
	else if ( test <= 0.5) in = bottom;
	else if ( test <=0.75) jn =left;
	else jn = right;
	
	neighbor[0] = in;
	neighbor[1] = jn;
	
}

//========select neighbor of pair=====================================//
int select_neighbor_of_pair(int in ,int jn, int i, int j){
	int left,right,top,bottom,leftn,rightn,topn,bottomn , neighbor_of_pair;
	
	if (i==0) top=N-1;
	else top = i-1;
	
	if (i==N-1) bottom = 0;
	else bottom = i+1 ;
	
	if (j==0) left = N-1;
	else left = j-1 ;
	
	if (j==N-1) right = 0;
	else right = j+1;
	
	if (in==0) topn=N-1;
	else topn = in-1;
	
	if (in==N-1) bottomn = 0;
	else bottomn = in+1 ;
	
	if (jn==0) leftn = N-1;
	else leftn = jn-1 ;
	
	if (jn==N-1) rightn = 0;
	else rightn = jn+1;
	
	int nn[6] ,c=0;
	
	if ((top*N +j) != (in*N+jn)) {
		nn[c]=top*N + j;
		c+=1;
	}
	if ((bottom*N + j) != (in*N+jn)) {
		nn[c]=bottom*N  + j;
		c+=1;
	}
	if ((i*N +right) != (in*N+jn)) {
		nn[c]= i*N + right;
		c+=1;
	}
	if ((i*N +left) != (in*N+jn)) {
		nn[c] = i*N + left;
		c+=1;
	}
	if ((topn*N +jn) != (i*N+j)) {
		nn[c]=topn*N + jn;
		c+=1;
	}
	if ((bottomn*N +jn) != (i*N+j)) {
		nn[c]=bottomn*N  + jn;
		c+=1;
	}
	if ((in*N +rightn) != (i*N+j)) {
		nn[c]= in*N + rightn;
		c+=1;
	}
	if ((in*N +leftn) != (i*N+j)) {
		nn[c] = in*N + leftn;
		c+=1;
	}
	
	float test =rand()/(float)RAND_MAX ;
	
	if (test <=(0.1666)) neighbor_of_pair= nn[0];
	else if ( test <= (2*0.1666)) neighbor_of_pair= nn[1];
	else if ( test <= (3*0.1666)) neighbor_of_pair= nn[2];
	else if ( test <= (4*0.1666)) neighbor_of_pair= nn[3];
	else if ( test <= (5*0.1666)) neighbor_of_pair= nn[4];
	else neighbor_of_pair = nn[5];
	

return neighbor_of_pair;

} 

//-------------------coarse grained variance-------------------//
float cg_var (int* A, int n, float mean) {				//n is the dimension of coarse graining i.e. dimension of the submatrix
	int i,j,k,l,count = 0;
	float reduced_mean[N*N/(n*n)] ;
	for (i=0; i<N-n+1; i+=n) {
		for (j=0; j<N-n+1; j+=n) {	
			float sum = 0;
			for (k=0; k<n; ++k) {
				for (l=0; l<n ; ++l) {
					sum+= (A[(i+k)*N+(j+l)]);
				}
			}
			
			reduced_mean[count] = sum/(float)(n*n);
			count+=1;
		}
	}
	float sum_sq=0;
	for (i=0; i<N*N/(n*n); i++) {
		sum_sq+= reduced_mean[i]*reduced_mean[i];
	}
	float mean_sq = sum_sq/(float)(N*N/(n*n)); 
	float variance = mean_sq - (mean*mean) ;
	//cout<<variance<<endl;
return variance;
}

//-------------------coarse grained skewness-------------------//
float cg_skew (int* A, int n, float mean, float variance) {				//n is the dimension of coarse graining i.e. dimension of the submatrix
	int i,j,k,l,count = 0;
	float reduced_mean[N*N/(n*n)] ;
	float skewness;
	for (i=0; i<N-n+1; i+=n) {
		for (j=0; j<N-n+1; j+=n) {	
			float sum = 0;
			for (k=0; k<n; ++k) {
				for (l=0; l<n ; ++l) {
					
					sum+= A[(i+k)*N+(j+l)] ;
				}
			}
			
			reduced_mean[count] = sum/(float)(n*n);
			count+=1;
		}
	}
	
	float sum_cube=0;
	for (i=0; i<N*N/n/n; i++) {
		sum_cube+= reduced_mean[i]*reduced_mean[i]*reduced_mean[i];
	}
	float mean_cube = sum_cube/(float)(N*N/(n*n));
	
	if (variance == 0) skewness = 0;
	else skewness = ( mean_cube - 3.0*mean*variance - mean*mean*mean )/ pow(variance,1.5) ;
	
	cout<<mean<<"\t"<<variance<<"\t"<<skewness<<endl;
return skewness;
}
			
		
////////////// main function //////////////////////////////////
int main(){
	
	srand(time(NULL));
	int x,l,t,i,j,z;
	float p[p_range],q=0.0;
	//float* den= new float[rep*p_range];
	float* av_skew= new float[p_range];
	float* skewness = new float[rep*p_range];
	int* neighbor = new int[2];
	
	for(i=0;i<p_range;++i) {
		p[i]= i/(float)p_range;
		//variance[i] = 0;
	}
	 
	int*A= new int[N*N];
	//int*A1= new int[N*N];
	
	for(x=600;x<700;++x){
		for(l=0;l<rep;++l){ 
			
			create_random_matrix(A,0.8);
			
			for(t=0;t<T;t++){
		
				for(z=0;z<N*N ; ++z){                // so that each site gets selected once on an average
				
					i = rand()%N;           // selecting one random site
					j = rand()%N;
					
					float test = rand()/(float)RAND_MAX;
					float test1 = rand()/(float)RAND_MAX;
					
					if (A[i*N+j]==1){     //if the site is occupied
						
						select_neighbor_of_site(i, j ,neighbor);    //look for a neighbor
						int in = neighbor[0] , jn = neighbor[1];
						
						if (A[in*N +jn]==0) {                     //if neighbor is empty
							if (test < p[x]) 
								A[in*N+jn]=1;                 //regular cp
							else A[i*N+j]=0;
						} 

						else {							  
							if (test < q){
								
								int neighbor_of_pair=select_neighbor_of_pair (in, jn, i, j);  //look for the neighbor of pair 
								A[neighbor_of_pair]=1;
							}
							else if (test1 < 1-p[x] )
								A[i*N+j]=0;	
						}
					}
			
					
				}	
				
			}
			float mean = density(A);
					
			skewness[l*p_range+x] = cg_skew (A,4,mean, cg_var(A,4,mean)); 
		}
		
		
	}
	
	average_along_row(skewness,av_skew);
	
	////////// saving data in a file///////////////////
	ofstream outdata; // outdata is like cin
	outdata.open("tcp_skewness_cg4_p=0.600:0.001:0.699_q0.dat"); // opens the file
	/*if( !outdata ) { // file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	} */

	for(i=600;i<700;++i){
		outdata<< skewness[i]<<"\n";
		//cout<< skewness[i]<<"\n";
	}
	outdata.close();
	
	
	delete[] A;
	//delete [] A1;
	//delete [] den;
	delete[] av_skew;
	delete[] neighbor;
	delete[] skewness;
	return 0;
}
