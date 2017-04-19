/* Parallel Divide and Conquer Matrix multiplication */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<pthread.h>

//structure to hold pairs of matrices and their dimension
typedef struct 
{
	float ** M1;
	float ** M2;
	int n; 
} mat;

//structure to hold a single matrix and its dimension
typedef struct
{
	float ** m;
	int n;
} matS;

//function to calculate the execution time in milliseconds
double calc_time(struct timespec start, struct timespec stop) 
{
	double t;
	t = (stop.tv_sec - start.tv_sec) * 1000; 
	t += (stop.tv_nsec - start.tv_nsec) * 0.000001; 
	return t;
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

//function to allocate memory to the two-matrix structure defined above
mat * initStruct2(int N)
{
	mat * m = (mat *) malloc(sizeof(mat));
	m->M1 = initMatrix(N);
	m->M2 = initMatrix(N);
	m->n = N;
	return m;
}


//function to allocate memory to the one-matrix structure defined above
matS * initStruct(int N)
{
	matS * M = (matS *) malloc(sizeof(matS));
	M->m = initMatrix(N);
	M->n = N;
	return M;
}

//function to free allocated memory from matrix type
void freeMat(float ** m, int N)
{
	for(int i=1; i<N; ++i)
	{
		free(m[i]);
	}
	free(m);
}
 
//function to free allocated memory from mat type structures
void freeStruct2(mat * M)
{
	freeMat(M->M1, M->n);
	freeMat(M->M2, M->n);
	free(M);
}

//function to free allocated memory from matS type structures
void freeStruct(matS * M)
{
	freeMat(M->m, M->n);
	free(M);
}

//function to add two matrices
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

//Thread Routine that is called recursively
void * remul(void * M)
{
	mat * MM = (mat *) M;
	int ord = MM->n;
	int N = ord/2;
	
	matS * c = initStruct(ord); //allocating memory to the resultant matrix
	
	//Base condition 
	if(ord==1)
	{
		c->m[1][1] = MM->M1[1][1] * MM->M2[1][1] ;
	}
	else
	{
		float ** temp;
		//Creating the structures which will be used to store the matrix pairs following the algorithm
		mat * m1 = initStruct2(N);//A11, B11
	
		mat * m2 = initStruct2(N);//A12, B21
	
		mat * m3 = initStruct2(N);//A11, B12
	
		mat * m4 = initStruct2(N);//A12, B22
	
		mat * m5 = initStruct2(N);//A21, B11
	
		mat * m6 = initStruct2(N);//A22, B21
	
		mat * m7 = initStruct2(N);//A21, B12
	
		mat * m8 = initStruct2(N);//A22, B22
		for(int i=1; i<=(N); ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				m1->M1[i][j] = MM->M1[i][j];
				m1->M2[i][j] = MM->M2[i][j];
				m3->M1[i][j] = MM->M1[i][j];
				m5->M2[i][j] = MM->M2[i][j];
			}
		}
		int k=1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				m5->M1[k][j] = MM->M1[i][j];
				m2->M2[k][j] = MM->M2[i][j];
				m7->M1[k][j] = MM->M1[i][j];
				m6->M2[k][j] = MM->M2[i][j];			
			}
			k++;
		}
		int l = 1;
		for(int i=1; i<=(N); ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				m2->M1[i][l] = MM->M1[i][j];
				m3->M2[i][l] = MM->M2[i][j];
				m4->M1[i][l] = MM->M1[i][j];
				m7->M2[i][l] = MM->M2[i][j];
				l++;
			}
			l = 1;
		}
		k = 1;
		l = 1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				m6->M1[k][l] = MM->M1[i][j];
				m4->M2[k][l] = MM->M2[i][j];
				m8->M1[k][l] = MM->M1[i][j];
				m8->M2[k][l] = MM->M2[i][j];		
				l++;
			}
			l = 1;
			k++;
		}
		//Recursive calls 
		matS * t1 = remul((void *) m1);
		matS * t2 = remul((void *) m2);
		matS * t3 = remul((void *) m3);
		matS * t4 = remul((void *) m4);
		matS * t5 = remul((void *) m5);
		matS * t6 = remul((void *) m6);
		matS * t7 = remul((void *) m7);
		matS * t8 = remul((void *) m8);
		
		freeStruct2(m1);
		freeStruct2(m2);
		freeStruct2(m3);
		freeStruct2(m4);
		freeStruct2(m5);
		freeStruct2(m6);
		freeStruct2(m7);
		freeStruct2(m8);
		
		//After getting the resultant matrices from the recursive calls, 
		//pair matrices for further operations specified by the algorithm
		temp = addMatrix(t1->m, t2->m, N);
		for(int i=1; i<=(N); ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				c->m[i][j] = temp[i][j];
			}
		}	
		freeStruct(t1);
		freeStruct(t2);
		freeMat(temp, N);
		temp = addMatrix(t3->m, t4->m, N);
		l = 1;
		for(int i=1; i<=(N); ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				c->m[i][j] = temp[i][l];
				l++;
			}
			l = 1;
		}
		freeStruct(t3);
		freeStruct(t4);
		freeMat(temp, N);
		
		temp = addMatrix(t5->m, t6->m, N);
		k = 1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				c->m[i][j] = temp[k][j];
				
			}
			k++;
		}
		freeMat(temp, N);
		freeStruct(t5);
		freeStruct(t6);
		temp = addMatrix(t7->m, t8->m, N);
		k = 1;
		l = 1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				c->m[i][j] = temp[k][l];
				l++;
			}
			l = 1;
			k++;
		}
		freeMat(temp, N);
		freeStruct(t7);
		freeStruct(t8);	
	}
	return (void *) c;
}
		
//Wrapper function for the routine function,
//Here is where the threads are created and initialised
matS * multiply(mat * M)
{
	int n = M->n;
	int N = n/2;
	pthread_t tid[8];
	
	void * stat1;
	void * stat2;
	void * stat3;
	void * stat4;
	void * stat5;
	void * stat6;
	void * stat7;
	void * stat8;
	
	matS * c1;
	matS * c2;
	matS * c3;
	matS * c4;
	matS * c5;
	matS * c6;
	matS * c7;
	matS * c8;
	
	float ** temp;
	matS * c = initStruct(n);
	
	//Splitting the input matrices into their submatrices and pairing them accordingly
	mat * m1 = initStruct2(N);//A11, B11
	mat * m2 = initStruct2(N);//A12, B21
	mat * m3 = initStruct2(N);//A11, B12
	mat * m4 = initStruct2(N);//A12, B22	
	mat * m5 = initStruct2(N);//A21, B11	
	mat * m6 = initStruct2(N);//A22, B21
	mat * m7 = initStruct2(N);//A21, B12	
	mat * m8 = initStruct2(N);//A22, B22
	for(int i=1; i<=(N); ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			m1->M1[i][j] = M->M1[i][j];
			m1->M2[i][j] = M->M2[i][j];
			m3->M1[i][j] = M->M1[i][j];
			m5->M2[i][j] = M->M2[i][j];
		}
	}
	int k=1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			m5->M1[k][j] = M->M1[i][j];
			m2->M2[k][j] = M->M2[i][j];
			m7->M1[k][j] = M->M1[i][j];
			m6->M2[k][j] = M->M2[i][j];			
		}
		k++;
	}
	int l = 1;
	for(int i=1; i<=(N); ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			m2->M1[i][l] = M->M1[i][j];
			m3->M2[i][l] = M->M2[i][j];
			m4->M1[i][l] = M->M1[i][j];
			m7->M2[i][l] = M->M2[i][j];
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
			m6->M1[k][l] = M->M1[i][j];
			m4->M2[k][l] = M->M2[i][j];
			m8->M1[k][l] = M->M1[i][j];
			m8->M2[k][l] = M->M2[i][j];		
			l++;
		}
		l = 1;
		k++;
	}
	
	//Creating the threads and passing the matrix pairs to the thread routine
	pthread_create(&tid[0], NULL, remul, (void*) m1);
	pthread_create(&tid[1], NULL, remul, (void*) m2);
	pthread_create(&tid[2], NULL, remul, (void*) m3);
	pthread_create(&tid[3], NULL, remul, (void*) m4);
	pthread_create(&tid[4], NULL, remul, (void*) m5);
	pthread_create(&tid[5], NULL, remul, (void*) m6);
	pthread_create(&tid[6], NULL, remul, (void*) m7);
	pthread_create(&tid[7], NULL, remul, (void*) m8);
	
	//Collecting the threads and storing the return values in the second parameter passed
	pthread_join(tid[0], &stat1);
	
	c1 = (matS *) stat1; //A11, B11
	pthread_join(tid[1], &stat2);
	
	c2 = (matS *) stat2;//A12, B21
	pthread_join(tid[2], &stat3);
	
	c3 = (matS *) stat3;//A11, B12
	pthread_join(tid[3], &stat4);
	
	c4 = (matS *) stat4;//A12, B22
	pthread_join(tid[4], &stat5);
	
	c5 = (matS *) stat5;//A21, B11
	pthread_join(tid[5], &stat6);
	
	c6 = (matS *) stat6;//A22, B21
	pthread_join(tid[6], &stat7);
	
	c7 = (matS *) stat7;//A21, B12
	pthread_join(tid[7], &stat8);
	
	c8 = (matS *) stat8;//A22, B22
	
	freeStruct2(m1);
	freeStruct2(m2);
	freeStruct2(m3);
	freeStruct2(m4);
	freeStruct2(m5);
	freeStruct2(m6);
	freeStruct2(m7);
	freeStruct2(m8);
	
	//Resultant matrices are paired accordingly for further operations as specified by the algorithm
	temp = addMatrix(c1->m, c2->m, N);
	for(int i=1; i<=(N); ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			c->m[i][j] = temp[i][j];
		}
	}
	freeMat(temp, N);	
	freeStruct(c1);
	freeStruct(c2);
	temp = addMatrix(c3->m, c4->m, N);
	l = 1;
	for(int i=1; i<=(N); ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			c->m[i][j] = temp[i][l];
			l++;
		}
		l = 1;
	}
	freeMat(temp, N);
	freeStruct(c3);
	freeStruct(c4);
	temp = addMatrix(c5->m, c6->m, N);
	k = 1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			c->m[i][j] = temp[k][j];
				
		}
		k++;
	}
	freeMat(temp, N);
	freeStruct(c5);
	freeStruct(c6);	
	temp = addMatrix(c7->m, c8->m, N);
	k = 1;
	l = 1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			c->m[i][j] = temp[k][l];
			l++;
		}
		l = 1;
		k++;
	}
	freeStruct(c7);
	freeStruct(c8);
	freeMat(temp, N);
	return c;
}


int main()
{
	struct timespec start, stop;	
	int N;
	scanf("%d", &N);
	mat * Matrix = (mat *) malloc(sizeof(mat));
	Matrix->n = N;
	Matrix->M1 = initMatrixV(N);
	Matrix->M2 = initMatrixV(N);
	//displayMatrix(Matrix->M1, Matrix->n);
	//displayMatrix(Matrix->M2, Matrix->n);
	matS * res;
	clock_gettime(CLOCK_REALTIME,&start);
	res = multiply(Matrix);
	//displayMatrix(res->m, res->n);
	clock_gettime(CLOCK_REALTIME,&stop);
	printf("%lf milliseconds\n",calc_time(start,stop));
	freeStruct2(Matrix);
	freeStruct(res);
	return 0;
}
