#include <iostream>
#include <new>
#include "nr3.h"
#include "ran.h"
#include "matriz.h"

using namespace std;

class TMatrix{
public:
  double*** M;
  int T,K,L;
  TMatrix(int t,int k,int l);
  void Aloca(int t,int k,int l);
  void print();
  TMatrix(){};
  ~TMatrix();
};

void combinationUtil(int arr[], int data[], int start, int end, int index, int r, double &D, int *C, double **D_ml);
 
void printCombination(int arr[], int n, int r, double &D, int *C, double **D_ml)
{
    int *data;
    data = new int[r];
    combinationUtil(arr, data, 0, n-1, 0, r, D, C, D_ml);
    delete[] data;
}

int main ( int argc, char *argv[] ){
  Ran ran(7);
  int K,R,P,Morris,M,sample,*dir,*dir_aux,*arr,*C;
  double delta,*p_val,dml,**D_ml;
  K = 6;
  Matrix B(K+1,K),D(K,K),PM(K,K),J(K+1,K,1),X(1,K),Jm(K+1,1,1);

  if ( 1 < argc ){R = atoi( argv[1] );}else{cout << "Enter the R number:" << endl;cin >> R;}
  if ( 2 < argc ){P = atoi( argv[2] );}else{cout << "Enter the p number:" << endl;cin >> P;}
  if ( 3 < argc ){Morris = atoi( argv[3] );}else{cout << "Enter 0 for Saltelli method or any other number for standard Morris:" << endl;cin >> Morris;}
  if(Morris){M = R;}else{M = 30;}

  TMatrix all_traj(M,K+1,K);
  TMatrix sel_traj(R,K+1,K);
  delta = P/(2.*(P-1.));
  p_val = new double[P];
  for(int i = 0; i < P; i++)
    p_val[i] = i/(P-1.);
  D_ml = new double*[M-1];
  for(int i = 0; i < (M-1); i++){
    D_ml[i] = new double[M-1-i];
    for(int j = 0; j < (M-1-i); j++)
      D_ml[i][j] = 0;
  }
  arr = new int[M];
  C = new int[R];
  for(int i = 0; i < M; i++)
    arr[i] = i;
  for(int i = 0; i < B.LC(0); i++){
    for(int j = i; j < B.LC(1); j++)
      B[i][j] = 0.;
    for(int j = 0; j < i; j++)
      B[i][j] = 1.;
  }
  dir = new int[K];
  dir_aux = new int[K];
  for(int m = 0; m < M; m++){
    for(int k = 0; k < K; k++){
      sample = (int) (P/2)*ran.doub();
      X[0][k] = p_val[sample];
      sample = (int) 2*ran.doub();
      D[k][k] = (sample-1)+sample;
      dir_aux[k] = k;
      for(int j = 0; j < K; j++)
	PM[k][j] = 0.;
    }
    for(int k = 0; k < K; k++){
      sample = (int) (K-k)*ran.doub();
      dir[k] = dir_aux[sample];
      dir_aux[sample] = dir_aux[K-1-k];
      PM[k][dir[k]] = 1;
    }
    Matrix Tb = (Jm*X+(delta/2.)*((2.*B-J)*D+J))*PM;
    //Tb.print();
    for(int i = 0; i < Tb.LC(0); i++)
      for(int j = 0; j < Tb.LC(1); j++)
	all_traj.M[m][i][j] = Tb[i][j];
  }
  if(Morris){
    ofstream output;
    output.open ("files.txt");
    for(int t = 0; t < all_traj.T; t++){
      for(int i = 0; i < all_traj.L; i++){
	output << scientific << all_traj.M[t][i][0];
	for(int j = 1; j < all_traj.K; j++){
	  output << scientific << " " << all_traj.M[t][i][j];
	}
	output << endl;
      }
    }
    output.close();
    int **Number;
    Number = new int*[K];
    for(int m = 0; m < K; m++)
      Number[m] = new int[P];
    for(int m = 0; m < K; m++)
      for(int n = 0; n < P; n++)
	Number[m][n] = 0;
    for(int t = 0; t < all_traj.T; t++){
      int i = 0;
      for(int j = 0; j < all_traj.K; j++){
	for(int n = 0; n < P; n++){
	  if(all_traj.M[t][i][j]==p_val[n]){
	    Number[j][n]++;
	    double val_tmp;
	    if(all_traj.M[t][i][j]+delta<=1.)
	      Number[j][n+2]++;
	    else 
	      Number[j][n-2]++;
	  }
	}
      }
    }
    for(int m = 0; m < K; m++){
      cout << Number[m][0];
      for(int n = 1; n < P; n++){
	cout << " " << Number[m][n];
      }
      cout << endl;
    }
    for(int i = 0; i < K; i++)
      delete[] Number[i];
    delete[] Number;
  }
  else{
    for(int m = 0; m < (M-1); m++){
      for(int l = (m+1) ; l < M; l++){
	for(int i = 0; i < all_traj.L; i++){
	  for(int j = 0; j < all_traj.L; j++){
	    dml = 0;
	    for(int z = 0; z < all_traj.K; z++){
	      dml+=pow(all_traj.M[m][i][z]-all_traj.M[l][j][z],2);
	    }
	    dml = sqrt(dml);
	    D_ml[m][l-(m+1)] += dml;
	  }
	}
      }
    }
    double DD = 0.;    
    printCombination(arr,M,R,DD,C,D_ml);
    
    for(int t = 0; t < sel_traj.T; t++){
      for(int i = 0; i < sel_traj.L; i++){
	for(int j = 0; j < sel_traj.K; j++){
	  sel_traj.M[t][i][j] = all_traj.M[C[t]][i][j];
	}
      }
    }
    ofstream output;
    output.open ("files.txt");
    for(int t = 0; t < sel_traj.T; t++){
      for(int i = 0; i < sel_traj.L; i++){
	output << scientific << sel_traj.M[t][i][0];
	for(int j = 1; j < sel_traj.K; j++){
	  output << scientific << " " << sel_traj.M[t][i][j];
	}
	output << endl;
      }
    }
    output.close();
    int **Number;
    Number = new int*[K];
    for(int m = 0; m < K; m++)
      Number[m] = new int[P];
    for(int m = 0; m < K; m++)
      for(int n = 0; n < P; n++)
	Number[m][n] = 0;
    for(int t = 0; t < sel_traj.T; t++){
      int i = 0;
      for(int j = 0; j < sel_traj.K; j++){
	for(int n = 0; n < P; n++){
	  if(sel_traj.M[t][i][j]==p_val[n]){
	    Number[j][n]++;
	    double val_tmp;
	    if(sel_traj.M[t][i][j]+delta<=1.)
	      Number[j][n+2]++;
	    else 
	      Number[j][n-2]++;
	  }
	}
      }
    }
    for(int m = 0; m < K; m++){
      cout << Number[m][0];
      for(int n = 1; n < P; n++){
	cout << " " << Number[m][n];
      }
      cout << endl;
    }
    for(int i = 0; i < K; i++)
      delete[] Number[i];
    delete[] Number;
    ofstream output2;
    output2.open ("traj.txt");
    for(int i = 0; i < R; i++)
      output2 << scientific << C[i] << " ";
    output2 << scientific << endl << "Distance = " << DD << endl;
    output2.close();
  }
  for(int i = 0; i < (M-1); i++)
    delete[] D_ml[i];
  delete[] D_ml;
  delete[] p_val;
  delete[] dir;
  delete[] dir_aux;
  delete[] C;
  delete[] arr;
  return 0;
}

void combinationUtil(int arr[], int data[], int start, int end, int index, int r, double &D, int *C, double **D_ml)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
      double d = 0.;
      for (int i = 0; i < (r-1); i++){
	for (int j = (i+1); j < r; j++){
	  d+=pow(D_ml[data[i]][data[j]-(i+1)],2);
	}
      }
      d = sqrt(d);
      if(d>D){
	D = d;
        for (int j=0; j<r; j++)
	  C[j] = data[j];
      }
        return;
    }
 
    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, r, D, C, D_ml);
    }
}


TMatrix::~TMatrix(){
  for(int m = 0; m < T; m++){
    for(int i = 0; i < L; i++){
      delete[] M[m][i];
    }
    delete[] M[m];
  }
  delete[] M;
}

TMatrix::TMatrix(int t,int k,int l){
  T = t;
  L = k;
  K = l;
  M = new double**[T];
  for(int m = 0; m < T; m++){
    M[m] = new double*[L];
    for(int i = 0; i < L; i++){
      M[m][i] = new double[K];
    }
  }
}

void TMatrix::Aloca(int t,int k,int l){
  T = t;
  L = k;
  K = l;
  M = new double**[T];
  for(int m = 0; m < T; m++){
    M[m] = new double*[L];
    for(int i = 0; i < L; i++){
      M[m][i] = new double[K];
    }
  }
}

void TMatrix::print(){
  //cout << fixed;
  cout << scientific;
  for(int t = 0; t < T; t++){
    for(int k = 0; k < L; k++){
      for(int l = 0; l < K; l++){
	cout << " " << M[t][k][l];
      }
      cout << endl;
    }
    cout << endl;
  }
}
