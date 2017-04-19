/* Divide and Conquer Matrix multiplication */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

//function to calculate execution time in milliseconds
double calc_time(struct timespec start, struct timespec stop) 
{
	double t;
	t = (stop.tv_sec - start.tv_sec) * 1000; 
	t += (stop.tv_nsec - start.tv_nsec) * 0.000001; 
	return t;
}

//function to free allocated memory
void freeThat(float ** a, float ** b, float ** c, float ** d, int n)
{
	for(int i=1; i<=n; ++i)
	{
		free(a[i]);
		free(b[i]);
		free(c[i]);
		free(d[i]);
	}
	free(a);
	free(b);
	free(c);
	free(d);
}

//allocates memory to matrix variables
float ** initMatrix(int n)
{
	float ** m = (float **) calloc(n, sizeof(float *));
	for(int i=1; i<=n; ++i)
	{
		m[i] = (float *) calloc(n, sizeof(float));
	}
	return m;
}

//function to read in input matrices
float ** initMatrixV(int n)
{
	float ** m = (float **) calloc(n, sizeof(float *));
	for(int i=1; i<=n; ++i)
	{
		m[i] = (float *) calloc(n, sizeof(float));
	}
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			scanf("%f", &m[i][j]);
		}
	}
	return m;
}

//function to add two matrices and return the sum matrix
float ** addMatrix(float ** a, float ** b, int n)
{
	float ** c = initMatrix(n);
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			c[i][j] = a[i][j] + b[i][j];
		}
	}
	return c;
}

//recursive function to multiply two matrices using Divide and Conquer approach	
float ** multiply(float ** a, float ** b, int n)
{
	int N = n/2;
	float ** c = initMatrix(n);
	if(n==1)
	{
		c[1][1] = a[1][1] * b[1][1];
	}
	else
	{
		float ** A11 = initMatrix(N);
		float ** A12 = initMatrix(N);
		float ** A21 = initMatrix(N);
		float ** A22 = initMatrix(N);
		float ** B11 = initMatrix(N);
		float ** B12 = initMatrix(N);
		float ** B21 = initMatrix(N);
		float ** B22 = initMatrix(N);
		float ** Ctemp;
		// Splitting the two matrices into smaller parts to be sent as parameters to the recursive call
		for(int i=1; i<=(N); ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				A11[i][j] = a[i][j];
				B11[i][j] = b[i][j];
			}
		}
		int k=1;
		for(int i=(N)+1; i<=n; ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				A21[k][j] = a[i][j];
				B21[k][j] = b[i][j];
			}
			k++;
		}
		int l = 1;
		for(int i=1; i<=(N); ++i)
		{
			for(int j=(N)+1; j<=n; ++j)
			{
				A12[i][l] = a[i][j];
				B12[i][l] = b[i][j];
				l++;
			}
			l = 1;
		}
		k = 1;
		l = 1;
		for(int i=(N)+1; i<=n; ++i)
		{
			for(int j=(N)+1; j<=n; ++j)
			{
				A22[k][l] = a[i][j];
				B22[k][l] = b[i][j];			
				l++;
			}
			l = 1;
			k++;
		}
		//All the recursive calls made here
		float ** t1 = multiply(A11, B11, (N));
		float ** t2 = multiply(A12, B21, (N));
		float ** t3 = multiply(A11, B12, (N));
		float ** t4 = multiply(A12, B22, (N));
		float ** t5 = multiply(A21, B11, (N));
		float ** t6 = multiply(A22, B21, (N));
		float ** t7 = multiply(A21, B12, (N));
		float ** t8 = multiply(A22, B22, (N));
		
		//the resultant matrices are then coupled together for further operations as specified by the algorithm
		Ctemp = addMatrix(t1, t2, (N));
		for(int i=1; i<=(N); ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				c[i][j] = Ctemp[i][j];
			}
		}
		for(int i=1; i<=N; ++i)
		{
			free(Ctemp[i]);
		}
		free(Ctemp);
		Ctemp = addMatrix(t3, t4, (N));
		l = 1;
		for(int i=1; i<=(N); ++i)
		{
			for(int j=(N)+1; j<=n; ++j)
			{
				c[i][j] = Ctemp[i][l];
				l++;
			}
			l = 1;
		}
		for(int i=1; i<=N; ++i)
		{
			free(Ctemp[i]);
		}
		free(Ctemp);
		Ctemp = addMatrix(t5, t6, (N));
		k = 1;
		for(int i=(N)+1; i<=n; ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				c[i][j] = Ctemp[k][j];
				
			}
			k++;
		}
		for(int i=1; i<=N; ++i)
		{
			free(Ctemp[i]);
		}
		free(Ctemp);
		Ctemp = addMatrix(t7, t8, (N));
		k = 1;
		l = 1;
		for(int i=(N)+1; i<=n; ++i)
		{
			for(int j=(N)+1; j<=n; ++j)
			{
				c[i][j] = Ctemp[k][l];
				l++;
			}
			l = 1;
			k++;
		}
		for(int i=1; i<=N; ++i)
		{
			free(Ctemp[i]);
		}
		free(Ctemp);
		freeThat(t1, t2, t3, t4, (N));
		freeThat(t5, t6, t7, t8, (N));
		freeThat(A11, A12, A21, A22, (N));
		freeThat(B11, B12, B21, B22, (N));
	}	
	return c; 
}

void displayMatrix(float ** m, int n)
{
	printf("\n-----------------------------------------------------------------\n");
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			printf("%f ", m[i][j]);
		}
		printf("\n");
	}
}

int main()
{
	struct timespec start, stop;	
	int n;
	scanf("%d", &n);
	float ** M1 = initMatrixV(n);
	float ** M2 = initMatrixV(n);
	//displayMatrix(M1, n);
	//displayMatrix(M2, n);
	clock_gettime(CLOCK_REALTIME,&start);
	float ** M3 = multiply(M1, M2, n);
	//displayMatrix(M3, n);
	clock_gettime(CLOCK_REALTIME,&stop);
	printf("%lf milliseconds\n",calc_time(start,stop));
	for(int i=1; i<=n; ++i)
	{
		free(M1[i]);
		free(M2[i]);
		free(M3[i]);
	}
	free(M1);
	free(M2);
	free(M3);
	return 0;
}
