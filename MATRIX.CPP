#include "matrix.h"
extern Mesh *mesh1;
extern Mesh *mesh2;
extern Mesh *mesh3;
extern Vector3 vector;
Matrix::Matrix():element(0){}
Matrix::Matrix(Matrix3 &m3){
	element=matrix();
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			element[i][j]=m3.element[i][j];
		}
	}
}
Matrix::Matrix(Matrix4 &m4){
	element=matrix(0,4,0,4);
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			element[i][j]=m4.element[i][j];
		}
	}
}
void Matrix::nrerror(char error_text[]) { 
		 fprintf(stderr,"Numerical Recipes Run-Time Error:\n"); 
		 fprintf(stderr,"%s\n",error_text); 
		 fprintf(stderr,"Now Exiting to system\n"); 
		 exit(1); 
} 

double *Matrix::vectorm(int nl,int nh) { 
		 double *v; 
  
		 v = (double *) malloc((unsigned) (nh-nl+1) * sizeof(double)); 
		 if (!v) nrerror("allocation failure in vectorm()"); 
		 return v-nl; 
} 

double **Matrix::matrix(int nrl,int nrh,int ncl,int nch) { 
		 int i; 
		 element = (double **) malloc((unsigned) (nrh-nrl+1) * sizeof(double*)); 
		 if (!element) nrerror("row allocation failure in matrix(r,c)"); 
		 element -= nrl; 
		 for (i=nrl;i<=nrh;i++) { 
			 element[i] = (double *) malloc((unsigned) (nch-ncl+1)*sizeof(double)); 
			 if (!element[i]) nrerror("column allocation failure in matrix(r,c)"); 
			 element[i] -= ncl; 
		 } 
		 return element; 
} 

void Matrix::free_vector(double *v,int nl,int nh){ 
		 free((char *) (v+nl)); 
} 

void Matrix::free_matrix(double **m,int nrl,int nrh,int ncl,int nch){ 
		 int i; 

		 for (i=nrh;i>=nrl;i--) free((char *) (m[i] + ncl)); 
		 free((char *) (m + nrl)); 
} 
void Matrix::inverse(double **mat, int dim) { 
	 int i,j,*indx; 
	 double **y,d,*col; 

	 y = matrix(0,dim-1,0,dim-1); 
	 indx = (int *)malloc((unsigned)(dim*sizeof(int))); 
	 col = vectorm(0,dim-1); 
	 ludcmp(mat,dim,indx,&d); 
	 for (j=0;j<dim;j++) { 
		 for (i=0;i<dim;i++) col[i] = 0.0; 
		 col[j] = 1.0; 
		 lubksb(mat,dim,indx,col); 
		 for (i=0;i<dim;i++) y[i][j] = col[i]; 
	 } 

		   	 for (i=0;i<dim;i++) 
		   		for (j=0;j<dim;j++) 
		   		   mat[i][j] = y[i][j]; 
		   
	 free_matrix(y,0,dim-1,0,dim-1); 
	 free_vector(col,0,dim-1); 
	 free(indx); 
} 

void Matrix::ludcmp(double **a, int n, int *indx, double *d) { 
		 int i,imax,j,k; 
		 double   big,dum,sum,temp; 
		 double   *vv; 

		 vv = (double*)malloc((unsigned)(n*sizeof(double))); 
		 if (!vv)  
		   { 
			fprintf(stderr,"Error Allocating Vector Memory\n"); 
			exit(1); 
		   } 
		 *d = 1.0; 
		 for (i=0;i<n;i++) 
			{ 
			 big = 0.0; 
			 for (j=0;j<n;j++) 
			 { 
			  if ((temp=fabs(a[i][j])) > big) big = temp; 
			 } 
			 if (big == 0.0) 
			   { 
				fprintf(stderr,"Singular Matrix in Routine LUDCMP\n"); 
			 for (j=0;j<n;j++) printf(" %f ",a[i][j]); printf("/n"); 
			 exit(1); 
			   } 
			 vv[i] = 1.0/big; 
			} 
		 for (j=0;j<n;j++) 
			{ 
			 for (i=0;i<j;i++) 
			 { 
			  sum = a[i][j]; 
			  for (k=0;k<i;k++) sum -= a[i][k] * a[k][j]; 
			  a[i][j] = sum; 
			 } 
			 big = 0.0; 
			 for (i=j;i<n;i++) 
				{ 
			  sum = a[i][j]; 
			  for (k=0;k<j;k++) sum -= a[i][k] * a[k][j]; 
			  a[i][j] = sum; 
			  if ((dum=vv[i]*fabs(sum)) >= big) 
				{ 
				 big = dum; 
				 imax = i; 
				} 
			 } 
			 if (j != imax) 
			   { 
			 for (k=0;k<n;k++) 
				{ 
				 dum = a[imax][k]; 
				 a[imax][k] = a[j][k]; 
				 a[j][k] = dum; 
				} 
				*d = -(*d); 
			 vv[imax] = vv[j]; 
			   } 
			 indx[j] = imax; 
			 if (a[j][j] == 0.0) a[j][j] = TINY; 
			 if (j != n-1) 
			   { 
			 dum = 1.0 / a[j][j]; 
			 for (i=j+1;i<n;i++) a[i][j] *= dum; 
			   } 
			} 
		 free(vv); 
		} 

void Matrix::lubksb(double **a, int n, int *indx, double *b) { 
		 int i,ip,j,ii=-1; 
		 double   sum; 

		 for (i=0;i<n;i++) 
			{ 
			 ip = indx[i]; 
			 sum = b[ip]; 
			 b[ip] = b[i]; 
			 if (ii>=0) 
			   for (j=ii;j<i;j++) sum -= a[i][j] * b[j]; 
			 else if (sum) ii = i; 
			 b[i] = sum; 
			} 
		 for (i=n-1;i>=0;i--) 
			{ 
			 sum = b[i]; 
			 for (j=i+1;j<n;j++) sum -= a[i][j] * b[j]; 
			 b[i] = sum / a[i][i]; 
			} 
		} 
void Matrix::print(double **m,int nrl,int nrh,int ncl,int nch){
		for(int i=nrl;i<nrh;i++){
			for(int j=ncl;j<nch;j++){
				printf("%4.4e\t",m[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
Matrix3 Matrix::AddMatrix3(Matrix3 &M1,Matrix3 &M2){
	Matrix3 matrix3;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			matrix3.element[i][j]=M1.element[i][j]+M2.element[i][j];
		}
	}
		return matrix3;
}
Matrix3 Matrix::Add(Matrix3 &M1,float fValue){
	Matrix3 matrix3;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			matrix3.element[i][j]=M1.element[i][j]+fValue;
		}
	}
	return matrix3;
}
Matrix3 Matrix::Mat3Div(Matrix3 &M,float n){
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			M.element[i][j]/=n;
		}
	}
	return M;
}
Matrix3 Matrix::Sqrt(Matrix3 &M){
	Matrix3 M1;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			M1.element[i][j]=sqrt(M.element[i][j]);
		}
	}
	return M1;
}
float Matrix::Mat3lemda(Matrix3 &M){
	Vector3 v(1.0f,1.0f,1.0f);
	static int iteration=0;
	float lemda,lemda1;
	lemda=ComMat3lemda(M,&v,iteration);
	lemda1=ComMat3lemda(M,&v,iteration);
	if(lemda<0||lemda1<0){
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	while(fabs(lemda-lemda1)>(0.001*fabs(lemda))){
		lemda=lemda1;
		lemda1=ComMat3lemda(M,&v,iteration);
		if (iteration>10000){
			printf("warning:iteration may not be convergent\n");
			return (float)((fabs(lemda1)<TINY)?TINY:lemda1);
		}
	}
	return lemda1;
}
float Matrix::ComMat3lemda(Matrix3 &M,Vector3 *v,int iteration){
	iteration++;
	*v=vector.Mat3MulV3(M,*v);
	float MaxLemda=vector.Max(*v);
	*v=*v/MaxLemda;
	return MaxLemda;
}
