#include <stdio.h>
#include <new>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main ( int argc, char *argv[] ){
  int K,lines,R,p;
  double delta;
  if ( 1 < argc ){R = atoi( argv[1] );}else{cout << "Enter the R number:" << endl;cin >> R;}
  if ( 2 < argc ){p = atoi( argv[2] );}else{cout << "Enter the p number:" << endl;cin >> p;}
  if ( 3 < argc ){K = atoi( argv[3] );}else{cout << "Enter the K number:" << endl;cin >> K;}
  delta = p/(2.*(p-1.));
  lines = (K+1)*R;
  double **Parametros,**Efeitos;
  double *M,*mup,*mu,*sd;
  Parametros = new double*[K];
  Efeitos = new double*[K];
  mup = new double[K];
  mu = new double[K];
  sd = new double[K];
  M = new double[lines];
  for(int k = 0; k < K; k++){
    Parametros[k] = new double[lines];
    Efeitos[k] = new double[R];
    mup[k] = mu[k] = sd[k] = 0.;
  }
  FILE *arq_data; 
  arq_data = fopen("resultados.txt","r");
  int grbg;
  for(int i = 0; i < lines; i++){
    for(int k = 0; k < K; k++)
      grbg = fscanf(arq_data, "%lf", &Parametros[k][i]); 
    grbg = fscanf(arq_data, "%lf", &M[i]); 
  }
  fclose(arq_data);
  for(int r = 0; r < R; r++)
    for(int k = 0; k < K; k++)
      for(int v = 0; v < K; v++)
	if(Parametros[v][r*(K+1) + k]!=Parametros[v][r*(K+1) + k+1]){
	  if(Parametros[v][r*(K+1) + k]>Parametros[v][r*(K+1) + k+1]){
	    Efeitos[v][r] = (M[r*(K+1) + k]-M[r*(K+1) + k+1])/delta;
	  }
	  else{
	    Efeitos[v][r] = (M[r*(K+1) + k+1]-M[r*(K+1) + k])/delta;
	  }
	}
  for(int r = 0; r < R; r++)
    for(int k = 0; k < K; k++){
      mu[k] += (1./R)*Efeitos[k][r];
      mup[k] += (1./R)*fabs(Efeitos[k][r]);
    }
  for(int r = 0; r < R; r++)
    for(int k = 0; k < K; k++)
      sd[k] += (1./(R-1.))*(Efeitos[k][r]-mu[k])*(Efeitos[k][r]-mu[k]);
  for(int k = 0; k < K; k++)
    sd[k] = sqrt(sd[k]);

  FILE *arq_out; 
  arq_out = fopen("mps.txt","w");
  for(int k = 0; k < K; k++)
    fprintf(arq_out,"%d %d %e %e %e %e\n", R, p, mup[k], mu[k], sd[k], sd[k]/mup[k]); 
  fclose(arq_out);

  for(int k = 0; k < K; k++){
    delete[] Parametros[k];
    delete[] Efeitos[k];
  }
  delete[] Parametros;
  delete[] Efeitos;
  delete[] mup;
  delete[] mu;
  delete[] sd;
  delete[] M;
  return 0;
}
