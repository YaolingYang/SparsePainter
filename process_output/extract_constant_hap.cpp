// g++ extract_constant_hap.cpp -o extract_constant_hap -lz -std=c++0x -g -O3

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include "gzstream.h"
#include "gzstream.C"
using namespace std;

double vec_sum(const vector<double>& l){
  // calculate the sum of a vector
  double sum_l=0;
  for(int i=0;i<l.size();++i){
    sum_l=sum_l+l[i];
  }
  return(sum_l);
}

vector<int> readSNP(const string& SNPfile) {
  ifstream file(SNPfile);
  vector<int> SNPindex;
  
  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << SNPfile << endl;
    abort();
  }
  
  string line;
  int column1;
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> column1;
    SNPindex.push_back(column1);
  }
  
  file.close();
  return SNPindex;
}

int main(int argc, char *argv[]){
  int npop;
  string SNPfile;
  string probfile;
  string out;
  
  
  for (int i = 1; i < argc; i++) {
    string param = argv[i];
    if (param[0] != '-') {
      cerr << "Invalid argument format. Expected -param value or -param \n";
      return 1;
    }
    param = param.substr(1);  // Remove the -

    if (param == "npop") {
      npop=stoi(argv[++i]);
    }else if (param == "SNPfile") {
      SNPfile = argv[++i];
    }else if (param == "probfile") {
      probfile = argv[++i];
    }else if (param == "out") {
      out = argv[++i];
    } else {
      cerr << "Unknown argument: " << param << ".\n";
      return 1;
    }
  }

  vector<int> SNPidx=readSNP(SNPfile);  //start from 1
  int nsnp=SNPidx.size();

  vector<vector<double>> painting(vector<vector<double>>(nsnp, vector<double>(npop)));
  
  
  //open files
  vector<string> filename;
  for(int k=0;k<npop;++k){
    filename.push_back(out+"pop"+to_string(k)+".txt.gz");
  }
  vector<ogzstream> files(npop);
  for(int k=0;k<npop;++k){
    files[k].open(filename[k].c_str());
  }
  
  igzstream in(probfile.c_str());
  
  if (!in) {
    cerr << "Error: Unable to open file " << probfile << endl;
    abort();
  }
  
  string line;
  
  // Read and discard the first two lines
  getline(in, line);
  getline(in, line);
  
  int start, end;
  int j=0;  // SNP index
  double vl;
  int i=-1; //haplotype index
  
  // Read the data lines
  while (getline(in, line)) {
    istringstream lineStream(line);
    string firstToken;
    lineStream >> firstToken;
    
    istringstream tokenStream(firstToken);
    if (tokenStream >> start) {
      lineStream >> end;
      if(SNPidx[j]>=start && SNPidx[j]<=end){
        vector<double> values;
        for(int k=0;k<npop;++k){
          lineStream >> vl;
          values.push_back(vl);
        }
        double value_sum=vec_sum(values);
        for(int k=0;k<npop;++k){
          values[k]=values[k]/value_sum;
        }
        while(SNPidx[j]>=start && SNPidx[j]<=end){

          for(int k=0;k<npop;++k){
            painting[j][k]=values[k];
            painting[j][k]=round(painting[j][k]*1000)/1000;
          }

          j++;
          if(j>=nsnp) break;
        }
      }
    } else {
      if(i==-1){
        for(int k=0;k<npop;++k){
          files[k] << firstToken<<" ";
        }
      }else{
        for(int k=0;k<npop;++k){
          for(int j=0; j<nsnp; ++j){
            files[k] << painting[j][k];
            if(j!=nsnp) files[k] <<" ";
          }
          files[k] << "\n";
          files[k] << firstToken<<" ";
        }
      }
      j=0;
      i++;
      cout<<firstToken<<" haplotype"<<i<<endl;
    }
  }
  
  in.close();
  
  for(int k=0;k<npop;++k){
    for(int j=0; j<nsnp; ++j){
      files[k] << painting[j][k];
      if(j!=nsnp) files[k] <<" ";
    }
    files[k].close();
  }
  
  
  return 0;
}