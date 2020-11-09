#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double normainf(int n, double *x);

int main (void){
	int i, j, k, n;
	double tol=1e-10, suma;
	double *vect, *vectnou, *b, **mat;
	srand(time(NULL)); 
	printf("Doni la dimensio de la matriu.\n");
	scanf("%d", &n);
	mat=(double **)malloc(n*sizeof(double *));
        vect=(double *)malloc(n*sizeof(double));
        vectnou=(double *)malloc(n*sizeof(double));
        b=(double *)malloc(n*sizeof(double));
        if(mat==NULL || vect==NULL || vectnou==NULL || b==NULL){
                printf("Error de memoria.\n");
                exit(1);
        }
        for(i=0; i<n; i++){
                mat[i]=(double *)malloc(n*sizeof(double));
                if(mat[i]==NULL){
                        printf("Error de memoria.\n");
                        exit(1);
                }
        }

	/*Utilitzo el rand perque el vector sigui aleatori*/
	for(i=0; i<n; i++)
		b[i]=rand()%10;
	printf("El vector b es\n");
	for(i=0; i<n; i++)
		printf("%le   ", b[i]);
	printf("\n");
	
	/*ESCRIC MATRIU*/
	/*Omplo la diagonal*/
	for(i=0; i<n; i++)
		mat[i][i]=1+n/2.;
	
	for(k=1; k<n-n%2; k+=2) 
		for(i=k, j=0; i<n; i++, j++){
			/*Omplo la matriu per sota*/
			mat[i][j]=1;
			/*Omplo la matriu per sobre*/
			mat[j][i]=1;
		}
	/*Comprovo matriu estigui ben escrita*/
	printf("La matriu es:\n");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++)
			printf("%le  ", mat[i][j]);
		printf("\n");
	}
	
	/*APLIQUEM METODE ITERATIU DE G-S*/
	for(i=0; i<n; i++)
		vectnou[i]=1;
	do{
		for(i=0; i<n; i++)
			vect[i]=vectnou[i];
		for(i=0; i<n; i++){
			suma=0;
			for(j=0; j<=i-1; j++)	
				suma+=mat[i][j]*vectnou[j];
			for(j=i+1; j<n; j++)
				suma+=mat[i][j]*vect[j];
			vectnou[i]=(b[i]-suma)/mat[i][i];
		}
	}while(fabs(normainf(n,vectnou)-normainf(n,vect))>tol);
	
	printf("La solucio al sistema Ax=b es:\n");
	for(i=0; i<n; i++)
		printf("x[%d]=%le\n", i, vectnou[i]);
	/*ALLIBERO MEMORIA*/	
	for(i=0; i<n; i++)
		free(mat[i]);
	free(mat);
	free(vect);
	free(vectnou);
	free(b);
	return 0;
}

double normainf(int n, double *x){
        int i;
        double norma=fabs(x[0]);
        for(i=1; i<n; i++)
                if(fabs(x[i])>norma)
                        norma=fabs(x[i]);
        return norma;
}


