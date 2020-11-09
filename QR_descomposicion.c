#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void factoritzacioQR(int n, double **A, double **Q, double **R);
double elementdiag(int n, double **a);
double normavector2 (int n, double *v);
double prodvect(int n, double *v, double *u);


int main(void){
	int i, j, k, n, itermax=200, compt=0;
	double tol=1e-4, suma;
	double **A, **Q, **R;
	printf("Doni la dimensio de la matriu.\n");
	scanf("%d", &n);
	A=(double**)malloc(n*sizeof(double*));
	Q=(double**)malloc(n*sizeof(double*));
	R=(double**)malloc(n*sizeof(double*));
	if(A==NULL || Q==NULL || R==NULL){
		printf("Error de memoria.\n");
		exit(1);
	}
	for(i=0; i<n; i++){
		A[i]=(double*)malloc(n*sizeof(double));
		Q[i]=(double*)malloc(n*sizeof(double));
		R[i]=(double*)malloc(n*sizeof(double));
		if(A[i]==NULL || Q[i]==NULL || R[i]==NULL){
			printf("Error de memoria.\n");
			exit(1);
		}
	}
	printf("Doni la matriu.\n");
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			scanf("%le", &A[i][j]);
	do{
		factoritzacioQR(n,A,Q,R);
		for(i=0; i<n; i++)
			for(k=0; k<n; k++){
				suma=0;
				for(j=0; j<n; j++)
					suma+=R[i][j]*Q[j][k];
				A[i][k]=suma;
			}
		compt++;
		printf("\nIteracio %d\n", compt);
		for(i=0; i<n; i++){
			for(j=0; j<n; j++)
				printf("%8.6le ", A[i][j]);
			printf("\n");
		}
	}while(elementdiag(n,A)>tol && j<itermax);
	printf("\nEls valors propis obtinguts amb %d iteracions son:\n", compt-1);
	for(i=0; i<n; i++)
		printf("%d: %12.6le\n", i+1, A[i][i]); 
	return 0;
}

/*Funció de la factoritzacioQR*/
void factoritzacioQR (int n, double **A, double **Q, double **R){
	int k,s,i;
	double *aux1, *aux2;
	aux1=(double*)malloc(n*sizeof(double));
	aux2=(double*)malloc(n*sizeof(double));
	if(aux1==NULL || aux2==NULL){
		printf("Error de memoria.\n");
		exit(1);
	}
	for(k=0; k<n; k++){
		for(s=0; s<n; s++)
			aux1[s]=A[s][k];
		R[k][k]=normavector2(n, aux1);
		for(s=0; s<n; s++)
			A[s][k]=A[s][k]/R[k][k];
		for(s=k+1; s<n; s++){
			for(i=0; i<n; i++){
				aux1[i]=A[i][k];
				aux2[i]=A[i][s];
			}
			R[k][s]=prodvect(n,aux1,aux2);
			for(i=0; i<n; i++)
				A[i][s]=A[i][s]-R[k][s]*A[i][k];
		}
	}
	for(k=0; k<n; k++)
		for(s=0; s<n; s++)
			Q[k][s]=A[k][s];
	return;
}

/*Em serveix per comprovar que totes els elements de sota la diagonal son menors a una tolerancia*/
double elementdiag(int n, double **a){
	int i,j;
	double max=0;
	for(i=0; i<n; i++)
		for(j=0; j<i; j++)
			if(fabs(a[i][j])>max)
				max=fabs(a[i][j]);
	return max;
}

/*Funció per calcular la norma vectorial 2*/
double normavector2 (int n, double *v){
	int i;
	double norma=0;
	for(i=0; i<n; i++)
		norma+=v[i]*v[i];
	norma=sqrt(norma);
	return norma;
}

/*Funció per multiplicar vectors*/
double prodvect(int n, double *v, double *u){
	int i;
	double sol=0;
	for(i=0; i<n; i++)
		sol+=u[i]*v[i];
	return sol;
}
