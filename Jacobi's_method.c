/*Mètode Jacobi*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double difdiagonal(int n, double **a, double **aux);

int main(void){
	int i, j, n, l, k, signe, itermax=20, compt=0;
	double tol=1e-6, max, c, s, dif, alfa;
	double **a, **v, **aux;
	/*A serà la matriu que convergirà als vaps
	  V serà la matriu que convergirà als veps
	  AUX serà la matriu auxiliar que usaré per anar actualitzant A i V*/

	printf("Doni la dimensio de la matriu.\n");
	scanf("%d", &n);

	/*Guardo memoria*/
	a=(double **)malloc(n*sizeof(double *));
	v=(double **)malloc(n*sizeof(double *));
	aux=(double **)malloc(n*sizeof(double *));
        if(a==NULL || v==NULL || aux==NULL){
                printf("Error de memoria.\n");
                exit(1);
        }
        for(i=0; i<n; i++){
                a[i]=(double *)malloc(n*sizeof(double));
                v[i]=(double *)malloc(n*sizeof(double));
                aux[i]=(double *)malloc(n*sizeof(double));
                if(a[i]==NULL || v[i]==NULL || aux[i]==NULL){
                        printf("Error de memoria.\n");
                        exit(1);
                }
        }
	
	/*Inicialitzo la matriu que convergirà als veps com a matriu unitat*/
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			v[i][j]=0;
	for(i=0; i<n; i++)
			v[i][i]=1;
	
	/*Llegeixo només la matriu triangular i així m'asseguro que serà simètrica*/	
	printf("Doni la matriu triangular.\n");
	for(i=0; i<n; i++)
		for(j=0; j<=i; j++)
			scanf("%le", &a[i][j]);
	for(i=0; i<n; i++)
		for(j=i+1; j<n; j++)
			a[i][j]=a[j][i];
	
	/*Mètode de Jacobi clàssic*/
	do{
		/*Busco |A[r][s]|>|A[i][j]|*/
		max=0;
		for(i=0; i<n; i++)
			for(j=0; j<i; j++)
				if(fabs(a[i][j])>max){
					max=fabs(a[i][j]);
					l=i;	k=j;
				}
		
		/*Busco l'angle*/
		alfa=fabs(a[l][l]-a[k][k])/sqrt((a[l][l]-a[k][k])*(a[l][l]-a[k][k])+4*a[l][k]*a[l][k]);
		signe=fabs((a[l][l]-a[k][k])*a[l][k])/((a[l][l]-a[k][k])*a[l][k]);
		c=sqrt((1+alfa)/2.);
		s=signe*sqrt((1-alfa)/2.);
		
		/*Utilitzo AUX per poder actualitzar A*/
		/*Actualitzem les components de la diagonal per poder fer una bona comprovació de convergència*/
		for(i=0; i<n; i++)
			aux[i][i]=a[i][i];
		/*Aplico les fórmules donades a classe*/
		for(i=0; i<n; i++){
			aux[i][l]=aux[l][i]=c*a[i][l]+s*a[i][k];
			aux[i][k]=aux[k][i]=-s*a[i][l]+c*a[i][k];
		}
		aux[l][l]=c*c*a[l][l]+s*s*a[k][k]+2*s*c*a[l][k];
		aux[k][k]=s*s*a[l][l]+c*c*a[k][k]-2*s*c*a[l][k];
		aux[l][k]=aux[k][l]=-s*c*(a[l][l]-a[k][k])+(c*c-s*s)*a[l][k];

		/*Comprovo si ha convergit*/
		dif=fabs(difdiagonal(n,a,aux));

		/*Actualitzo A*/
		for(i=0; i<n; i++){
			a[i][l]=a[l][i]=aux[i][l];
			a[i][k]=a[k][i]=aux[i][k];
		}

		/*Utilitzo AUX per poder actualitzar V*/
		for(i=0; i<n; i++){
			aux[i][l]=c*v[i][l]+s*v[i][k];
			aux[i][k]=-s*v[i][l]+c*v[i][k];
		}

		/*Actualitzo V*/
		for(i=0; i<n; i++){
			v[i][l]=aux[i][l];
			v[i][k]=aux[i][k];
		}
		compt++;
	}while(dif>tol && compt<itermax);

	/*Escric els vaps i els seus respectius veps*/
	for(i=0; i<n; i++){
		printf("%d vap=%12.6le - vep=[ ", i+1, a[i][i]);
		for(j=0; j<n; j++)
			printf("%12.6le ", v[i][j]);
		printf("]\n");
	}
		
	/*Allibero memoria*/	
	for(i=0; i<n; i++){
		free(a[i]);
		free(v[i]);
		free(aux[i]);
	}
	free(a);
	free(v);
	free(aux);
	return 0;
}

/*Creo aquesta funció per calcular la diferència entre les components de la diagonal e dos matrius*/
double difdiagonal(int n, double **a, double **aux){
        int i;
        double dif=0;
        for(i=0; i<n; i++)
		if(fabs(a[i][i]-aux[i][i])>dif)
			dif=fabs(a[i][i]-aux[i][i]);
        return dif;
}

