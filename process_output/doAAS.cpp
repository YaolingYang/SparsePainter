// please compile with "make"
// g++ -I./armadillo-12.6.5/include doAAS.cpp -o doAAS -lz -fopenmp -lpthread -L./armadillo-12.6.5 -larmadillo -llapack -lblas -std=c++0x -g -O3 -Wl,-rpath=./armadillo-12.6.5
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <utility>
#include <armadillo>
#include <set>
#include "gzstream.h"
#include "gzstream.C"
using namespace std;
using namespace arma;

void readAveSNP(const string& filename, 
                vector<double>& pd, 
                vector<vector<double>>& aveSNPpainting) {
  igzstream in(filename.c_str());
  string line;
  string dummy;
  
  if (getline(in, line)) {
    // columns are space separated.
    istringstream lineStream(line);
    lineStream >> dummy;
    double number;
    while (lineStream >> number) {
      pd.push_back(number);
    }
  }
  
  cout<<pd.size()<<endl;
  
  while( getline(in, line)){
    istringstream lineStream(line);
    lineStream >> dummy;
    double number;
    vector<double> rowData;
    while (lineStream >> number) {
      rowData.push_back(number);
    }
    aveSNPpainting.push_back(rowData);
  }
  
  in.close();
}

vector<double> rowMeans(const vector<vector<double>>& data) {
  vector<double> means;
  for (const auto& row : data) {
    double sum = 0;
    for (double num : row) {
      sum += num;
    }
    means.push_back(sum / row.size());
  }
  return means;
}


void doAAS(vector<double>& pd, 
           vector<vector<double>>& aveSNPpainting,
           const string AASfile) {
  
  default_random_engine generator;
  normal_distribution<double> distribution(0.0, 1e-6);
  
  int nsnp = pd.size();
  int npop = aveSNPpainting.size();
  
  vector<double> p_values(nsnp, 1.0); // initialize p_values with 1
  vector<double> test_statistic(nsnp, 1.0); 
  
  vector<double> mu=rowMeans(aveSNPpainting);
  
  arma::mat Astar(nsnp, npop);
  for (int i = 0; i < npop; ++i) {
#pragma omp parallel for
    for (int j = 0; j < nsnp; ++j) {
      double random_error = distribution(generator);
      Astar(j, i) = aveSNPpainting[i][j] - mu[i] + random_error; // compute A*(j,k) for all j,k and store it in Astar
    }
  }
  
  // compute covariance matrix 
  arma::mat C = arma::cov(Astar);
  
  arma::mat C_inv = arma::inv(C);  // compute the inverse of the covariance matrix
  
  // Compute the SVD of the inverse covariance matrix
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd(U, s, V, C_inv);
  
  // Compute the whitening matrix W
  arma::mat W = U * arma::diagmat(arma::sqrt(s));
  
  // Compute Z = Astar * W
  arma::mat Z = Astar * W;
  
#pragma omp parallel for
  for (int j = 0; j < nsnp; ++j) {
    double t = arma::norm(Z.row(j), "fro");
    test_statistic[j]=t*t; // square of Frobenius norm
  }
  
  //output the AAS results into AASfile
  ofstream outputFile(AASfile);
  if (outputFile.is_open()) {
    outputFile.precision(15);
    outputFile << "physical_position" << " " << "test_statistic"<< "\n";
    for (int j = 0; j < nsnp; ++j) {
      outputFile << fixed<< setprecision(0) << pd[j];
      outputFile << " " << fixed<< setprecision(3) << test_statistic[j]<< "\n";
    }
    outputFile.close();
  } else {
    cerr << "Unable to open file" << AASfile;
  }
}

int main(int argc, char *argv[]){
  string aveSNPfile;
  string out;

  for (int i = 1; i < argc; i++) {
    string param = argv[i];
    if (param[0] != '-') {
      cerr << "Invalid argument format. Expected -param value or -param \n";
      return 1;
    }
    param = param.substr(1);  // Remove the -
    
    if (param == "aveSNPfile") {
      aveSNPfile = argv[++i];
    }else if (param == "out") {
      out = argv[++i];
    } else {
      cerr << "Unknown argument: " << param << ".\n";
      return 1;
    }
  }
  
  string AASfile = out+"_AAS.txt";
  
  vector<double> pd;
  vector<vector<double>> aveSNPpainting;
  
  cout<<"Begin reading aveSNPfile"<<endl;
  readAveSNP(aveSNPfile, pd, aveSNPpainting);
  cout<<"Finish reading aveSNPfile"<<endl;
  
  
  cout << "Begin calculating Ancestry Anomaly Score"<<endl;
  doAAS(pd,aveSNPpainting,AASfile);
  cout << "Finish calculating Ancestry Anomaly Score"<<endl;
  
  return 0;
}