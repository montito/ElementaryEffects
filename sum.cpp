# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>

using namespace std;

int main ( int argc, char *argv[] ){
  int N;
  double X[6],Y=1,a[6]={78,12,0.5,2,97,33};
  if(1<argc){N = atoi(argv[1]);}else{cout << "Enter the output file number N" << endl;cin >> N;}
  if(2<argc){X[0] = atof(argv[2]);}else{cout << "Enter the stochastic value for X1" << endl;cin >> X[0];}
  if(3<argc){X[1] = atof(argv[3]);}else{cout << "Enter the stochastic value for X2" << endl;cin >> X[1];}
  if(4<argc){X[2] = atof(argv[4]);}else{cout << "Enter the stochastic value for X3" << endl;cin >> X[2];}
  if(5<argc){X[3] = atof(argv[5]);}else{cout << "Enter the stochastic value for X4" << endl;cin >> X[3];}
  if(6<argc){X[4] = atof(argv[6]);}else{cout << "Enter the stochastic value for X5" << endl;cin >> X[4];}
  if(7<argc){X[5] = atof(argv[7]);}else{cout << "Enter the stochastic value for X6" << endl;cin >> X[5];}
  for(int i = 0; i < 6; i++)
    Y*=(fabs(4.*X[i]-2.)+a[i])/(1.+a[i]);
  FILE *arq; 
  char name[200];
  sprintf(name,"output%d.txt",N); 
  arq = fopen(name,"w");
  fprintf(arq, "%e %e %e %e %e %e %e", X[0], X[1], X[2], X[3], X[4], X[5], Y); 
  fclose(arq);
  return 0;
}
