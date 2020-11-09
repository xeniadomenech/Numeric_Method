#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 20

double fun(int n, double x[N], int i);
double der(int n, double x[N], int i, int j);
int rsl(int n, double A[N][N], double b[N], double z[N]);
void pivot(int n, double A[N][N], double b[N], int k);


int main(void){
	int i,j,k=0,n;
	double kmax=100, prec=1e-6, znorma=2*prec;
	double A[N][N], z[N], x[N], f[N];

	/*Recullo n, però m'asseguro que n<N*/
	do{
		printf("Doni la dimensio. (Recordi que ha de ser <=20)\n");
		scanf("%d", &n);
	}while(n>N);

	/*Llegeixo el vector inicial*/
	printf("Doni el vector inicial x[n]\n");
	for(i=0; i<n; i++)
		scanf("%le", &x[i]);

	while(k<=kmax && znorma>=prec){
		k++;

		/*Actualitzo f*/
		for(i=0; i<n; i++)
			f[i]=fun(n,x,i);

		/*Actualitzo A*/
		for(i=0; i<n; i++)
			for(j=0; j<n; j++)
				A[i][j]=der(n,x,i,j);
		if(rsl(n,A,f,z)==1){
			printf("No s'ha pogut dur a terme el mètode.\n");
			exit(1);
		}

		/*Faig X=X-Z*/
		for(i=0; i<n; i++)
			x[i]=x[i]-z[i];

		/*Calculo znorma*/
		znorma=0;
		for(i=0; i<n; i++)
			znorma+=(z[i]*z[i]);
		znorma=sqrt(znorma);
		
		/*Imprimeixo els resultats d'aquesta iteracio*/
		printf("ITERACIO %d\n\nEl vector z es:\n",k);
		for(i=0; i<n; i++)
			printf("z[%d]=%le\n", i, z[i]);
		printf("\nEl vector x es:\n");
		for(i=0; i<n; i++)
			printf("x[%d]=%le\n", i, x[i]);
		printf("\n||z||=%le\n\n", znorma);
	}
	return 0;
}

double fun(int n, double x[N], int i){
	double pi=4*atan(1.);
	switch(i){
		case 0 : return 3*x[0]-cos(x[1]*x[2])-1./2;
			break;
		case 1 : return x[0]*x[0]-81*(x[1]+0.1)*(x[1]+0.1)+sin(x[2])+1.06;
			break;
		case 2 : return exp(-x[0]*x[1])+20*x[2]+(10*pi-3)/3.;
			break;
		default: printf("Error\n");
			return 0;
	}
}
	

double der(int n, double x[N], int i, int j){
	double a, b, h=1e-8;
	
	/*Calculo per separat fun(n,x+hej,i) i fun(n,x-hej,i)*/
	x[j]=x[j]+h;
	a=fun(n,x,i);

	x[j]=x[j]-2*h;
	b=fun(n,x,i);

	return (a-b)/(2*h);
}

int rsl(int n, double A[N][N], double b[N], double z[N]){
	int i,j,k;
	double aux, tol=1e-10;
	for(i=0; i<n-1; i++){

		/*Fem pivotatge maximal per columnes*/
		pivot(n,A,b,i);
		if(fabs(A[i][i])<tol)
			return 1;

		/*Fem eliminació gaussina*/
		for(k=i+1; k<n; k++){
			aux=A[k][i]/A[i][i];
			b[k]=b[k]-aux*b[i];
			for(j=i; j<n; j++)
				A[k][j]=A[k][j]-A[i][j]*aux;
		}
	}

	/*Fem substitució enrere*/
	for(i=n-1; i>=0; i--){
		aux=0;
		for(j=i+1; j<n; j++)
			aux-=A[i][j]*z[j];
		aux+=b[i];
		if(fabs(A[i][i])<tol)
			return 1;
		z[i]=aux/A[i][i];
	}
	return 0;
}

void pivot(int n, double A[N][N], double b[n], int k){
	int i,r;
	double max, aux[n], baux;
	max=0;
	for(i=k; i<n; i++)
		if(max<fabs(A[i][k])){
			max=fabs(A[i][k]);
			r=i;
		}
	for(i=0; i<n; i++)
		aux[i]=A[k][i];
	for(i=0; i<n; i++)
		A[k][i]=A[r][i];
	for(i=0; i<n; i++)
		A[r][i]=aux[i];
	baux=b[k];
	b[k]=b[r];
	b[r]=baux;
}

