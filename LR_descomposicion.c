/*Xènia Domènech*/
/*Mètode LR*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int eliminaciogauss (int n, double **a);
void prodmatRL (int n, double **a);
double elementdiag (int n, double **a);

int main(void){
	int i, j, k=1, n, itmax=200;
	double tol=1e-4;
	double **a;
	printf("Doni la dimensio.\n");
	scanf("%d", &n);
	a=(double**)calloc(n,sizeof(double*));
	if(a==NULL){
		printf("Error de memoria.\n");
		exit(1);
	}
	for(i=0; i<n; i++){
		a[i]=(double*)calloc(n, sizeof(double));
		if(a[i]==NULL){
			printf("Error de memoria.\n");
			exit(1);
		}
	}
/*	printf("Doni la matriu per files.\n");
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			scanf("%le", &a[i][j]);
*/
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			a[i][j]=1./3;
	a[0][0]=1;
	a[1][1]=3;
	a[2][2]=4;
	do{
		if(eliminaciogauss(n,a)==1){
			prodmatRL(n,a);
			printf("Pel pas %d la matriu es:\n", k);
			for(i=0; i<n; i++){
				for(j=0; j<n; j++)
					printf("%9.6le ", a[i][j]);
				printf("\n");
			}
		}else{
			printf("No s'ha pogut fer EGSP.\n");
		}
		k++;
	}while(elementdiag(n,a)>tol && k<itmax);
	printf("\n\nAmb tolerancia %6.le han calgut %d iteracions.\nEls valors propis de la matriu son:\n",tol, k-1);
	for(i=0; i<n; i++)
		printf("%d: %12.6le\n", i+1, a[i][i]);
	for(i=0; i<n; i++)
		free(a[i]);
	free(a);
	return 0;
}

/*El programa d'EGSP retorna 0 si algun pivot es proxim a zero i 1 en cas si s'ha pogut fer LU*/ 
int eliminaciogauss (int n, double **a){
	int i,j,k;
	double m, tol=1e-3;
	/*A la meva matriu A li assignaré la matriu U a la part superior i la L a la part inferior*/
	for(i=0; i<n-1; i++)
		for(k=i+1; k<n; k++){
			if(fabs(a[i][i])<tol)
				return 0;
			m=a[k][i]/a[i][i];
			for(j=i; j<n; j++)
				a[k][j]=a[k][j]-a[i][j]*m;
			a[k][i]=m;
			}
	return 1;
}	

void prodmatRL (int n, double **a){
	int i,j,k;
	double suma;
	double **aux;
	aux=(double**)calloc(n,sizeof(double*));
	if(aux==NULL){
		printf("Error de memoria.\n");
		exit(1);
	}
	for(i=0; i<n; i++){
		aux[i]=(double*)calloc(n,sizeof(double));
		if(aux[i]==NULL){
			printf("Error de memoria.\n");
			exit(1);
		}
	}
	/*Inicialment la meva matriu A té a la part superior(+diag) la matriu U i a la part inferior la matriu L*/ 
	for(i=0; i<n; i++)
		for(k=0; k<n; k++){
			suma=0;
			if(k>=i)/*Es l'equivalent a fer a[i][k]·1 */
				suma=a[i][k]; 
			for(j=k+1; j<n; j++){
				/*Nomès tinc en compte els elements per sobre la diagonal, el altres es multipliquen per 0*/
				if(j>=i)
					suma+=a[i][j]*a[j][k];
			}		
		aux[i][k]=suma;
		}
	/*Actualitzo A*/
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			a[i][j]=aux[i][j];
	for(i=0; i<n; i++)
		free(aux[i]);
	free(aux);
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
