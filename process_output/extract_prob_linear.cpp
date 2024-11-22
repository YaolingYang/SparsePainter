// g++ extract_prob_linear.cpp -o extract_linear -lz -std=c++0x -g -O3

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
  
  // Read and discard the first three lines
  getline(in, line);
  getline(in, line);
  
  int start=1;
  int end;
  int j=0;  // SNP index
  double vl;
  int i=-1; //haplotype index
  vector<double> prev_prob(npop);
  vector<double> current_prob(npop);
  
  // Read the data lines
  while (getline(in, line)) {
    istringstream lineStream(line);
    string firstToken;
    lineStream >> firstToken;
    
    istringstream tokenStream(firstToken);
    if (tokenStream >> end) {
      if(end==1){
        for(int k=0;k<npop;++k){
          lineStream >> vl;
          prev_prob[k]=vl*0.01;
          current_prob[k]=vl*0.01;
        }
      }else{
        for(int k=0;k<npop;++k){
          lineStream >> vl;
          current_prob[k]=vl*0.01;
        }
      }
      if(SNPidx[j]>=start && SNPidx[j]<=end){
        int range=end-start;
        
        while(SNPidx[j]>=start && SNPidx[j]<=end){
          vector<double> values;
          int dist=SNPidx[j]-start;
          if(range==0){
            for(int k=0;k<npop;++k){
              values.push_back(prev_prob[k]);
            }
          }else{
            for(int k=0;k<npop;++k){
              values.push_back(prev_prob[k]+(current_prob[k]-prev_prob[k])*dist/range);
            }
          }
          double value_sum=vec_sum(values);
          for(int k=0;k<npop;++k){
            values[k]=0.5*values[k]/value_sum;
          }
          
          if(i%2==0){
            for(int k=0;k<npop;++k){
              painting[j][k]=values[k];
            }
          }else{
            for(int k=0;k<npop;++k){
              painting[j][k]+=values[k];
              painting[j][k]=round(painting[j][k]*1000)/1000;
              
            }
          }
          j++;
          if(j>=nsnp) break;
        }
      }
      start=end;
      for(int k=0;k<npop;++k){
        prev_prob[k]=current_prob[k];
      }
    } else {
      // The first token is a string, add it to the indname vector
      if(i%2!=0){
        if(i>0){
          for(int k=0;k<npop;++k){
            for(int j=0; j<nsnp; ++j){
              files[k] << painting[j][k];
              if(j!=nsnp) files[k] <<" ";
            }
            files[k] << "\n";
          }
        }
        for(int k=0;k<npop;++k){
          files[k] << firstToken.substr(0, firstToken.size() - 2) <<" ";
        }
      }
      
      j=0;
      i++;
      start=1;
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