// module load languages/gcc/10.4.0
// g++ combine_chunklength.cpp -o combine.exe -lz -lpthread -llapack -lblas -std=c++0x -g -O3
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
#include "gzstream.h"
#include "gzstream.C"
using namespace std;

class hVec { // A sparse vector format
public:
  vector<int> k; // the keys that are in the vector
  unordered_map<int, double> v; // the values, stored as a map from keys to values
  int len; // nominal length of the vector; currently unused
  double x0; // default value for entries
  hVec(){
    // Create an empty vector
    len=0;
    x0=0;
  };
  hVec(int len,
       double x0){
    // Create a vector of length len filled with x0
    this->len=len;
    this->x0=x0;
  };
  hVec(vector<int> idx,
       vector<double> val,
       int len,
       double x0){
    // Create a vector of length len filled with x0 except at idx which contains val
    this->len=len;
    this->x0=x0;
    for(int i=0;i<idx.size();++i){
      k.push_back(idx[i]);
      v[idx[i]]=val[i];
    }
  };
  void setdefault(double x0){
    // Change the default value
    this->x0=x0;
  };
  void setnocheck(int p,
                  double val){
    // Set a value, should be known to be in the keys
    v[p]=val;
  };
  void set(int p,
           double val){
    //Safely set a value
    if(!in(p)) k.push_back(p);
    setnocheck(p,val);
  };
  void setall(vector<int> p,
              vector<double> val){
    for(int i=0; i<p.size();++i){
      if(!in(p[i])) k.push_back(p[i]);
      setnocheck(p[i],val[i]);
    }
  }
  bool in(int p){
    // Check if a value has a non-default entry
    if(v.find(p)==v.end()) return(false);
    return(true);
  };
  double get(int p){
    // Get a value from the vector: either its set value or the default if not present
    if(!in(p)){
      return(x0);
    }else{
      return(v[p]);
    }
  };
  
  vector<double> getall(const vector<int>& idx){
    // Get values from the vector: either its set value or the default if not present
    
    vector<double> values(idx.size());
    for(int i=0;i<idx.size();++i){
      if(!in(idx[i])) {
        values[i]=x0;
      } else {
        values[i]=v[idx[i]];
      }
    }
    return(values);
  }
  
};

class hMat {
public:
  vector<hVec> m; // sparse matrix, i.e. a vector of hVec's
  int d1; // number of rows; currently nominal
  int d2; // Number of columns; should be equal to length(m)
  hMat(int d1){
    // Empty matrix with d1 rows (can append columns)
    this->d1=d1;
    d2=0;
  };
  hMat(int d1,
       int d2,
       double x0=0.0){
    // Create a d1 by d2 matrix taking value x0
    this->d1=d1;
    this->d2=d2;
    for(int i=0;i<d2;++i){
      // each element in m is a column vector
      m.push_back(hVec(d1,x0));
    }
  };
  void appendColumn(double x0){
    // Append a default column
    ++d2;
    m.push_back(hVec(d1,x0));
  };
  void appendColumn(vector<int> idx,
                    vector<double> vals,
                    double x0){
    // Append a filled column
    ++d2;
    m.push_back(hVec(idx,vals,d1,x0));
  };
  
};

void readdatafirst(const string& filename,
                   hMat &cl,
                   double weight_use) {
  igzstream file(filename.c_str());
  
  string line;
  int ind1,ind2;
  double value;
  int q=1;
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> ind1;
    lineStream >> ind2;
    lineStream >> value;
    
    if(ind1==q){
      if(q%1000==0){
        cout<<"ind"<<q<<endl;
      }
      q++;
    }
    
    cl.m[ind1-1].set(ind2-1,value*weight_use);
  }
  
  file.close();
}



void readdata(const string& filename,
              hMat &cl,
              double weight_use) {
  igzstream file(filename.c_str());
  
  string line;
  int ind1,ind2;
  double value;
  int q=1;
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> ind1;
    lineStream >> ind2;
    lineStream >> value;
    
    if(ind1==q){
      if(q%1000==0){
        cout<<"ind"<<q<<endl;
      }
      q++;
    }
    
    cl.m[ind1-1].set(ind2-1,value*weight_use+cl.m[ind1-1].get(ind2-1));
  }
  
  file.close();
}


int main() {
  vector<double> total_SNP = {53260,52631,44263,41213,38734,45886,36061,33741,29133,32626,33134,
                              31429,22134,21337,20859,23774,22215,19250,19139,17197,9737,10911};
  vector<double> total_gd = {286.28,268.83,223.26,214.54,204.09,192.03,187.16,168.00,
                             166.35,180.91,158.22,174.68,125.70,120.20,141.35,134.05,
                             128.78,117.71,107.73,108.28,62.79,74.10};
  
  int nsnp=487409;
  
  // Compute weight
  vector<double> weight;
  for (int i = 0; i < total_SNP.size(); ++i) {
    double val = total_gd[i] / total_SNP[i];
    val = round(val * 1e6) / 1e6; // round to 6 decimal places
    weight.push_back(val);
  }
  
  hMat cl(nsnp,nsnp,0.0);
  readdatafirst("chr1_UKBall.chunklengths.s.out.gz",cl,weight[0]);
  
  for (int i = 1; i <= 21; ++i) {
    cout << "Processing chromosome " << i + 1 << endl;
    string filename = "chr" + to_string(i + 1) + "_UKBall.chunklengths.s.out.gz";
    readdata(filename,cl,weight[i]);
  }
  
  // Save the merged data to a file
  ogzstream outFile("full_chunklength_UKBall.txt.gz");
  for (int i=0; i<nsnp; ++i) {
    vector<int> twj=cl.m[i].k;
    for(int j=0;j<twj.size();++j){
      double val=cl.m[i].get(twj[j]);
      if(val>=0.000005){
        outFile << i+1 <<" "<< twj[j]+1 <<" " << fixed << setprecision(5)<<cl.m[i].get(twj[j])<<"\n";
      }
    }
  }
  outFile.close();
  
  return 0;
}