#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double polinp(int m, double *p, double x);
double dervE(int n, int m, int j, double *x, double *y, double *p);
int gauss (int n, double **a, double *b);

int main(void){
	int i,j,k,n,m;
	double *x,*y,*indep,**a;
	printf("Doni quants punts introduira.\n");
	scanf("%d", &m);
	printf("De quin grau vol el polinomi.\n");
	scanf("%d", &n);
	x=(double*)calloc(m+1, sizeof(double));
	y=(double*)calloc(m+1, sizeof(double));
	indep=(double*)calloc(n+1, sizeof(double));
	a=(double**)calloc(n+1, sizeof(double*));
	if(x==NULL || y==NULL || indep==NULL || a==NULL){
		printf("Error de memoria.\n");
		exit(1);
	}
	for(i=0; i<=m; i++){
		a[i]=(double*)calloc(m+1, sizeof(double));
		if(a[i]==NULL){
			printf("Error de memoria.\n");
			exit(1);
		}
	}
	printf("Doni els xk's i yk's (per parelles)\n");
	for(i=0; i<=m; i++)
		scanf("%le %le", &x[i], &y[i]);
	/*Omplim la matriu*/
	for(i=0; i<=n; i++)
		for(j=0; j<=n; j++){
			a[i][j]=0;	
			for(k=0; k<=m; k++)
				a[i][j]+=pow(x[k],i+j);
		}
	for(i=0; i<=n; i++){
		indep[i]=0;
		for(k=0; k<=m; k++)
			indep[i]+=y[k]*pow(x[k],i);
	}
	/*Resolc la matriu*/
	if(gauss(n+1,a,indep)==0){
		printf("S'ha pogut solucionar el sistema.\nEl polinomi es:\n");
		printf("p%d(x)=%10.4le ", n, indep[0]);
		for(i=1; i<=n; i++)
			printf("%+10.4e·x^%d ", indep[i], i);
	}
	printf("\n");
	/*Comprovo si la solucio verifica derivades parcials*/
	printf("E=%8.2le\n", dervE(n,m,i,x,y,indep));
	
	/*Allibero memoria*/
	for(i=0; i<=n; i++)
		free(a[i]);
	free(a);
	free(indep);
	free(x);
	free(y);
	return 0;
}

double polinp(int n, double *p, double x){
	int i;
	double pol=p[0];
	for(i=1; i<=n; i++){
		pol+=p[i]*x;
		x=x*x;
	}	
	return pol;
}

double dervE(int n, int m, int j, double *x, double *y, double *p){
	int i;
	double sol=0, pol;
	printf("Els residus yk-p2(xk) son:\n");
	for(i=0; i<=m; i++){
		pol=polinp(n,p,x[i]);
		sol+=(y[i]-pol)*(y[i]-pol);
		printf("%8.2le\t", y[i]-pol);
	}
	printf("\n");
	sol=sqrt(sol);
	return sol;
}

int gauss (int n, double **a, double *b){
        int i,j,k;
        double aux, tol=1e-10;
        for(i=0; i<n-1; i++){
		/*Comprovo que a[i][i] no és 0*/
		if(fabs(a[i][i])<tol)
			return 1;

		/*Fem eliminacio gaussiana*/
		for(k=i+1; k<n; k++){
			aux=a[k][i]/a[i][i];
			b[k]=b[k]-aux*b[i];
			for(j=i; j<n; j++)
				a[k][j]=a[k][j]-a[i][j]*aux;
		}
	}

	/*Fem substitució enrere*/
	for(i=n-1; i>=0; i--){
		aux=0;
		for(j=i+1; j<n; j++)
			aux-=a[i][j]*b[j];
		aux+=b[i];
		if(fabs(a[i][i])<tol)
			return 1;
		b[i]=aux/a[i][i];
	}
	return 0;
}

