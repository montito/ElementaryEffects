// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

using namespace std;

class Matrix{
  friend Matrix operator*(const Matrix&, const double&);
  friend Matrix operator*(const double&, const Matrix&);
  friend Matrix operator*(const Matrix&, const Matrix&);
  friend Matrix operator+(const Matrix&, const Matrix&);
  friend Matrix operator+(const Matrix&, const double&);
  friend Matrix operator-(const Matrix&, const Matrix&);
  friend Matrix operator-(const Matrix&, const double&);
 public:
  double **M;
  int Linha,Coluna;
  Matrix(int,int,double = 0.);
  Matrix(const Matrix& r);
  void print();
  Matrix& operator=(const Matrix&);
  Matrix& operator=(const double&);
  Matrix& operator+=(const Matrix&);
  Matrix& operator+=(const double&);
  Matrix& operator-=(const Matrix&);
  Matrix& operator-=(const double&);
  Matrix& operator*=(const double&);
  Matrix& operator*=(const Matrix&);
  ~Matrix();
  int& LC(int i);
  double * operator[](int y){ return M[y]; }
};

Matrix& Matrix::operator*=(const Matrix& r){
  if (Coluna!=r.Linha) {
    cout << "ERROR: Operator (*), X cols must be equal Y rows." << endl;
    exit(1);
  }
  Matrix Z(Linha,r.Coluna);
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < r.Coluna; j++)
      for(int z = 0; z < Coluna; z++)
	Z.M[i][j] += M[i][z]*r.M[z][j];
  for(int i = 0; i < Linha; i++)
    delete[] M[i];
  delete[] M;
  M = new double*[Linha];
  for(int i = 0; i < Linha; i++)
    M[i] = new double[r.Coluna];
  Coluna = r.Coluna;
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = Z.M[i][j];
  return *this;
}

Matrix operator*(const Matrix& X, const Matrix& Y){
  if (X.Coluna!=Y.Linha) {
    cout << "ERROR: Operator (*), X cols must be equal Y rows." << endl;
    exit(1);
  }
  Matrix Z(X.Linha,Y.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < Y.Coluna; j++)
      for(int z = 0; z < X.Coluna; z++)
	Z.M[i][j] += X.M[i][z]*Y.M[z][j];
  return Z;
}

int& Matrix::LC(int i){ 
  if (i == 0) return Linha;
  else return Coluna;
}

Matrix operator*(const double& Y, const Matrix& X){
  Matrix Z(X.Linha,X.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < X.Coluna; j++)
      Z.M[i][j] = X.M[i][j]*Y;
  return Z;
}

Matrix operator*(const Matrix& X, const double& Y){
  Matrix Z(X.Linha,X.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < X.Coluna; j++)
      Z.M[i][j] = X.M[i][j]*Y;
  return Z;
}

Matrix operator-(const Matrix& X, const double& Y){
  Matrix Z(X.Linha,X.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < X.Coluna; j++)
      Z.M[i][j] = X.M[i][j] - Y;
  return Z;
}

Matrix operator-(const Matrix& X, const Matrix& Y){
  if (X.Linha!=Y.Linha || X.Coluna!=Y.Coluna) {
    cout << "ERROR: Operator (-), matrix sizes should be the same." << endl;
    exit(1);
  }
  Matrix Z(X.Linha,X.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < X.Coluna; j++)
      Z.M[i][j] = X.M[i][j] - Y.M[i][j];
  return Z;
}

Matrix operator+(const Matrix& X, const double& Y){
  Matrix Z(X.Linha,X.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < X.Coluna; j++)
      Z.M[i][j] = X.M[i][j] + Y;
  return Z;
}

Matrix operator+(const Matrix& X, const Matrix& Y){
  if (X.Linha!=Y.Linha || X.Coluna!=Y.Coluna) {
    cout << "ERROR: Operator (+), matrix sizes should be the same." << endl;
    exit(1);
  }
  Matrix Z(X.Linha,X.Coluna);
  for(int i = 0; i < X.Linha; i++)
    for(int j = 0; j < X.Coluna; j++)
      Z.M[i][j] = X.M[i][j] + Y.M[i][j];
  return Z;
}

Matrix& Matrix::operator-=(const double& r){
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = M[i][j] - r;
  return *this;
}

Matrix& Matrix::operator+=(const double& r){
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = M[i][j] + r;
  return *this;
}

Matrix& Matrix::operator*=(const double& r){
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = M[i][j] * r;
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& r){
  if (Linha!=r.Linha || Coluna!=r.Coluna) {
    cout << "ERROR: Operator (-), matrix sizes should be the same." << endl;
    exit(1);
  }
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = M[i][j] - r.M[i][j];
  return *this;
}

Matrix& Matrix::operator+=(const Matrix& r){
  if (Linha!=r.Linha || Coluna!=r.Coluna) {
    cout << "ERROR: Operator (+), matrix sizes should be the same." << endl;
    exit(1);
  }
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = M[i][j] + r.M[i][j];
  return *this;
}

Matrix& Matrix::operator=(const double& r){
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = r;  
  return *this;
}

Matrix& Matrix::operator=(const Matrix& r){
  if (Linha!=r.Linha || Coluna!=r.Coluna) {
    cout << "ERROR: Operator (=), matrix sizes should be the same." << endl;
    exit(1);
  }
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = r.M[i][j];  
  return *this;
}

Matrix::Matrix(const Matrix& r){
  Linha = r.Linha;
  Coluna = r.Coluna;
  M = new double*[Linha];
  for(int i = 0; i < Linha; i++)
    M[i] = new double[Coluna];
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = r.M[i][j];  
}

Matrix::Matrix(int l,int c,double v){
  Linha = l;
  Coluna = c;
  M = new double*[Linha];
  for(int i = 0; i < Linha; i++)
    M[i] = new double[Coluna];
  for(int i = 0; i < Linha; i++)
    for(int j = 0; j < Coluna; j++)
      M[i][j] = v;  
}

Matrix::~Matrix(){
  for(int i = 0; i < Linha; i++)
    delete[] M[i];
  delete[] M;
}

void Matrix::print(){
  cout << endl;
  cout << scientific;
  for(int l = 0; l < Linha; l++){
    cout << M[l][0];
    for(int c = 1; c < Coluna; c++)
      cout << " " << M[l][c];
    cout << endl;
  }
  cout << endl;
}
