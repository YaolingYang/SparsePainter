// Compile with:
// module load languages/gcc/10.4.0
// g++ SparsePainter.cpp -o SparsePainter.exe -lz -fopenmp -lpthread -larmadillo -std=c++0x -g -O3
// on HPC
// module load libs/armadillo/12.4.0
// g++ SparsePainter.cpp -o SparsePainter.exe -lz -fopenmp -lpthread -L/mnt/storage/software/libraries/gnu/12.4.0/lib64 -larmadillo -std=c++0x -g -O3
// if the module load isn't available, please install armadillo on their official website
// using cmake . and make to install, and then
// g++ -I/mnt/storage/scratch/ip21972/1000GUKB/armadillo-9.850.1/include SparsePainter.cpp -o SparsePainter.exe -lz -fopenmp -lpthread -L/mnt/storage/scratch/ip21972/1000GUKB/armadillo-9.850.1 -larmadillo  -lopenblas -std=c++0x -g -O3
// export LD_LIBRARY_PATH=/mnt/storage/scratch/ip21972/1000GUKB/armadillo-9.850.1:$LD_LIBRARY_PATH

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


class hAnc {
public:
  //use hashmap to find the positions (rows) of each ancestry
  unordered_map<int, vector<int>> pos;
  hAnc(const vector<int>& ref) {
    // Traverse the ref vector and store the row position corresponding to each value in unordered_map
    for (int i = 0; i < ref.size(); i++) {
      int value = ref[i];
      if (pos.find(value) == pos.end()) { //if this value doesn't exist
        pos[value] = vector<int>{i};
      } else {
        pos[value].push_back(i);
      }
    }
  };
  vector<int> findrows(int value) const {
    // here the class of it is unordered_map<int, vector<int>>::const_iterator
    // we use auto to simplify
    auto it = pos.find(value);
    if (it != pos.end()) {
      return it->second;
    } else { // if the value doesn't exist, return an empty vector
      return vector<int>();
    }
  };
};


/////////////////////beginning of pbwt contents///////////////////////////

void free_PBWT_memory(vector<vector<bool>> &panel, int** &prefix, int** &divergence, int** &u, int** &v) {
  // First, delete the inner arrays of prefix, divergence, u, and v
  // Note: Since temp1, temp2, temp3, and temp4 are continuous blocks of memory 
  // you only need to delete their base pointers (the pointers originally returned by 'new').
  delete[] prefix[0];     // which is temp1
  delete[] divergence[0]; // which is temp2
  delete[] u[0];          // which is temp3
  delete[] v[0];          // which is temp4
  
  // Now, delete the outer arrays of prefix, divergence, u, and v
  delete[] prefix;
  delete[] divergence;
  delete[] u;
  delete[] v;
  
  // Clear the panel vector and minimize its memory usage
  panel.clear();
  vector<vector<bool>>().swap(panel); // This technique is used to shrink the vector's capacity to fit its size.
  
  // Nullify the pointers to ensure that they don't dangle.
  prefix = nullptr;
  divergence = nullptr;
  u = nullptr;
  v = nullptr;
}


void PBWT(vector<vector<bool>> &panel, int **prefix, int **divergence,
          int **u, int **v, int num, int N){
  for (int i = 0; i<num; ++i){
    prefix[i][0] = i;
    divergence[i][0] = 0;
  }
  for (int k = 0; k<N; ++k){
    int u2 = 0, v2 = 0, p = k+1, q = k+1;
    vector<int> a,b,d,e;
    for (int i = 0; i<num; ++i){
      u[i][k] = u2;
      v[i][k] = v2;
      if (divergence[i][k] > p) { p = divergence[i][k];}
      if (divergence[i][k] > q) { q = divergence[i][k];}
      if (!panel[prefix[i][k]][k]){
        a.push_back(prefix[i][k]);
        d.push_back(p);
        ++u2;
        p = 0;
      }
      else{
        b.push_back(prefix[i][k]);
        e.push_back(q);
        ++v2;
        q = 0;
      }
    }
    for (int i = 0; i<num; ++i){
      v[i][k] += a.size();
      if (i < a.size()){
        prefix[i][k+1] = a[i];
        divergence[i][k+1] = d[i];
      }
      else{
        prefix[i][k+1] = b[i-a.size()];
        divergence[i][k+1] = e[i-a.size()];
      }
    }
  }
}

void ReadVCF(const string inFile,
             const string qinFile,
             vector<vector<bool>> &panel, 
             const int N, 
             const int M, 
             const int qM,
             const bool haploid){
  
  cout << "Read reference data with "<<N<<" SNPs for "<<M-qM<<" haploptypes";
  if(inFile!=qinFile){
    cout<<" and target data with "<<N<<" SNPs for "<<qM<<" haploptypes" << endl;
  }else{
    cout<<endl;
  }
  
  igzstream in,qin;
  if(inFile==qinFile){
    string line = "##";
    in.open(inFile.c_str());
    if (!in) {
      cerr << "Error: unable to open file: " << inFile << endl;
      abort();
    }
    stringstream linestr;
    int x = 0;
    char y = 0;
    
    while (line[1] == '#')
      getline(in, line);
    for(int j = 0; j<N; ++j){
      getline(in, line);
      linestr.str(line);
      linestr.clear();
      for (int i = 0; i<9; ++i){
        linestr >> line;
      }
      if(!haploid){
        for (int i = 0; i<(M-qM)/2; ++i){
          linestr >> x >> y;
          panel[i*2][j] = (bool)x;
          linestr >> x;
          panel[i*2 + 1][j] = (bool)x;
        }
      }else{
        for (int i = 0; i<M-qM; ++i){
          linestr >> x;
          panel[i][j] = (bool)x;
        }
      }
    }
    in.close();
  }else{
    string line = "##", qline = "##";
    in.open(inFile.c_str());
    if (!in) {
      cerr << "Error: unable to open file: " << inFile << endl;
      abort();
    }
    qin.open(qinFile.c_str());
    if (!qin) {
      cerr << "Error: unable to open file: " << qinFile << endl;
      abort();
    }
    stringstream linestr, qlinestr;
    int x = 0;
    char y = 0;
    
    while (line[1] == '#')
      getline(in, line);
    while (qline[1] == '#')
      getline(qin, qline);
    for(int j = 0; j<N; ++j){
      getline(in, line);
      getline(qin, qline);
      linestr.str(line);
      linestr.clear();
      qlinestr.str(qline);
      qlinestr.clear();
      for (int i = 0; i<9; ++i){
        linestr >> line;
        qlinestr >> qline;
      }
      if(!haploid){
        for (int i = 0; i<(M-qM)/2; ++i){
          linestr >> x >> y;
          panel[i*2][j] = (bool)x;
          linestr >> x;
          panel[i*2 + 1][j] = (bool)x;
        }
        for (int i = (M-qM)/2; i < M/2; ++i){
          qlinestr >> x >> y;
          panel[i*2][j] = (bool)x;
          qlinestr >> x;
          panel[i*2+1][j] = (bool)x;
        }
      }else{
        for (int i = 0; i<M-qM; ++i){
          linestr >> x;
          panel[i][j] = (bool)x;
        }
        for (int i = (M-qM); i < M; ++i){
          qlinestr >> x;
          panel[i][j] = (bool)x;
        }
      }
      
    }
    in.close();
    qin.close();
  }
}

void Readphase_donor(const string inFile,
                     vector<vector<bool>> &panel, 
                     const int N, 
                     const int M, 
                     const int qM) {
  
  
  cout << "Read reference data with "<<N<<" SNPs for "<<M-qM<<" haploptypes." << endl;
  
  // read the data
  igzstream in;
  in.open(inFile.c_str());
  
  if (!in) {
    cerr << "Error: unable to open file: " << inFile << endl;
    abort();
  }
  
  string line;
  
  // Read and discard the first three lines
  for (int i = 0; i < 3; ++i) {
    getline(in, line);
  }
  
  // Read the remaining lines and store the binary data of each line in 'panel'
  // We are reading the first M-qM haplotypes
  for(int i=0; i<M-qM; ++i) {
    // i indicates which sample we are looking at
    getline(in, line);
    vector<bool> panelsnp;
    
    // convert SNP data to binary
    for (char c : line) {
      panelsnp.push_back(c == '1');
    }
    
    int Oid = i;
    // add snps to the panel
    panel.push_back(vector<bool>());
    panel[i].resize(N);
    
    for (int k = 0; k<N; ++k){ // for every SNP
      
      panel[i][k] = panelsnp[k];
      
    } // end loop over snps
  }
  
  cout<<"Finish reading reference data"<<endl;
  
  in.close();
}

vector<int> getorder(const vector<double>& vec) {
  vector<int> order(vec.size());
  iota(order.begin(), order.end(), 0);
  
  stable_sort(order.begin(), order.end(), [&vec](int i, int j) {
    return vec[i] < vec[j];
  });
  
  unordered_map<double, vector<int>> groups;
  for (int i : order) {
    groups[vec[i]].push_back(i);
  }
  
  random_device rd;
  mt19937 g(rd());
  for (auto& group : groups) {
    shuffle(group.second.begin(), group.second.end(), g);
  }
  
  vector<int> randomized_order;
  for (int i : order) {
    randomized_order.push_back(groups[vec[i]].back());
    groups[vec[i]].pop_back();
  }
  
  return randomized_order;
}

bool containsIndex(const vector<int>& fullidx, 
                   int starttemp, 
                   int endtemp) {
  bool contain=false;
  for(int i : fullidx) {
    if(i >= starttemp && i <= endtemp) {
      contain=true;
    }
  }
  return(contain);
}


tuple<vector<int>,vector<int>,vector<int>,vector<int>> longMatchpbwt(const int L_initial,
                                                                     vector<vector<bool>> &panel, 
                                                                     int **prefix, 
                                                                     int **divergence, 
                                                                     int **u, 
                                                                     int **v,
                                                                     int minmatch,
                                                                     vector<double> &gd,
                                                                     vector<int>& queryidx,
                                                                     const int N,
                                                                     const int M,
                                                                     const int qM,
                                                                     const int L_minmatch,
                                                                     const int ncores,
                                                                     const bool samefile,
                                                                     const bool phase,
                                                                     const string qinFile){
  
  igzstream in;
  
  string line;
  
  if(phase & !samefile){
    in.open(qinFile.c_str());
    if (!in) {
      cerr << "Error: unable to open file: " << qinFile << endl;
      abort();
    }
    
    // Read and discard the first three lines
    for (int i = 0; i < 3; ++i) {
      getline(in, line);
    }
  }
  
  // match of which query sample
  vector<int> queryidall={0};
  // match to which reference sample (donor)
  vector<int> donorid;
  // start position of match
  vector<int> startpos;
  // end position of match
  vector<int> endpos;
  
  struct LoopResult {
    vector<int> donorid;
    vector<int> startpos;
    vector<int> endpos;
    int queryid;
  };
  
  // define a results vector
  vector<LoopResult> allResults(queryidx.size());
  
  //read data for target haplotypes
  
  // store ncores lines of target data
  // paneltarget has ncore rows and N columns
  
  int nind=queryidx.size();
  
  int nind_left=nind;
  
  omp_set_num_threads(ncores);
  
  if(samefile){
    minmatch++;
  } 
  
  while(nind_left>0){
    int ncores_use = (ncores < nind_left) ? ncores : nind_left;
    vector<vector<bool>> panelsnp;
    if(samefile){
      panelsnp=vector<vector<bool>>(ncores_use,vector<bool>(N));
      for(int i=nind-nind_left; i<nind-nind_left+ncores_use; ++i) {
        for(int j=0;j<N;++j){
          panelsnp[i-nind+nind_left][j]=panel[i][j];
        }
      }
    }else{
      if(!phase){
        panelsnp=vector<vector<bool>>(ncores_use,vector<bool>(N));
        for(int i=nind-nind_left; i<nind-nind_left+ncores_use; ++i) {
          for(int j=0;j<N;++j){
            panelsnp[i-nind+nind_left][j]=panel[M-qM+i][j];
          }
        }
      }else{
        panelsnp=vector<vector<bool>>(ncores_use);
        for(int i=nind-nind_left; i<nind-nind_left+ncores_use; ++i) {
          getline(in, line);
          for (char c : line) {
            panelsnp[i-nind+nind_left].push_back(c == '1');
          }
        }
      }
    }
    
    cout<<"Finding matches with PBWT for target haplotypes "<<nind-nind_left<<"-"<<nind-nind_left+ncores_use-1<<endl;
    
#pragma omp parallel for
    
    for (int idx=nind-nind_left; idx<nind-nind_left+ncores_use; ++idx) {
      
      int *dZ;
      dZ = new int[M];
      for (int i = 0; i<M; i++){
        dZ[i] = 0;
      }
      
      int *t = new int[N+1];
      int *zd = new int[N+2], *bd = new int[N+2];
      
      int i = queryidx[idx];
      
      int L=L_initial;
      int prevL=L;
      
      int Oid = i;
      t[0] = 0;
      
      for (int k=0; k<N; ++k){
        if (t[k]!=M-qM)
          if (!panelsnp[Oid-nind+nind_left][k])
            t[k+1] = u[t[k]][k];
          else
            t[k+1] = v[t[k]][k];
          else
            if (!panelsnp[Oid-nind+nind_left][k])
              t[k+1] = v[0][k];
            else
              t[k+1] = M-qM;
      }
      
      
      zd[N+1] = bd[N+1] = N;
      
      for (int k = N; k>=0; --k){
        zd[k] = min(zd[k+1],k);
        bd[k] = min(bd[k+1],k);
        if (t[k]!=0)
          while(zd[k]>0 && 
                panelsnp[Oid-nind+nind_left][zd[k]-1]
                  == panel[prefix[t[k]-1][k]][zd[k]-1])
            zd[k]--;
        else
          zd[k] = k;
        if (t[k]!=M-qM)
          while(bd[k]>0 && 
                panelsnp[Oid-nind+nind_left][bd[k]-1] 
                  == panel[prefix[t[k]][k]][bd[k]-1])
            bd[k]--;
        else
          bd[k] = k;
      }
      
      // below using while loop to update L and ensure at least minmatch matches at each SNP
      
      vector<int> donoridtemp;
      
      vector<int> startpostemp;
      
      vector<int> endpostemp;
      
      vector<int> nomatchsnp;
      
      vector<int> nmatch(N,0);
      
      int times=0;
      
      bool allsnpmatches=false;
      
      vector<int> local_donorid;
      vector<int> local_startpos;
      vector<int> local_endpos;
      
      while(!allsnpmatches){
        
        vector<bool> addmatch(N,false);
        
        if(times!=0){
          for(int w=0;w<nomatchsnp.size();++w){
            for(int s=0;s<prevL;++s){
              int pos=nomatchsnp[w]+s;
              if(pos>=N) break;
              addmatch[pos]=true;
            }
          }
        }
        
        int f, g, ftemp, gtemp;
        f = g = t[0];
        
        for(int k = 0; k<N; ++k){
          if (g == M-qM){
            if (f == M-qM){
              if (!panelsnp[Oid-nind+nind_left][k]){
                ftemp = M-qM;
                f = v[0][k];
              }
              else{
                ftemp = v[0][k];
                f = M-qM;
              }
            }
            else{
              if (!panelsnp[Oid-nind+nind_left][k]){
                ftemp = v[f][k];
                f = u[f][k];
              }
              else{
                ftemp = u[f][k];
                f = v[f][k];
              }
            }
            if (!panelsnp[Oid-nind+nind_left][k]){
              gtemp = M-qM;
              g = v[0][k];
            }
            else{
              gtemp = v[0][k];
              g = M-qM;
            }
          }
          else
            if (!panelsnp[Oid-nind+nind_left][k]){
              ftemp = v[f][k];
              gtemp = v[g][k];
              f = u[f][k];
              g = u[g][k];
            }
            else{
              ftemp = u[f][k];
              gtemp = u[g][k];
              f = v[f][k];
              g = v[g][k];
            }
            
            
            while (ftemp != gtemp){
              int end=k-1;
              if(times==0){
                int start=dZ[prefix[ftemp][k+1]];
                donoridtemp.push_back(prefix[ftemp][k+1]);
                startpostemp.push_back(start);
                endpostemp.push_back(end);
                ++ftemp;
                for(int q=start;q<=end;++q){
                  nmatch[q]++;
                }
              }else{
                if(addmatch[end]){
                  int start=dZ[prefix[ftemp][k+1]];
                  //add new matches with new L
                  if(end-start+1<prevL){
                    donoridtemp.push_back(prefix[ftemp][k+1]);
                    startpostemp.push_back(start);
                    endpostemp.push_back(end);
                    ++ftemp;
                    for(int q=start;q<=end;++q){
                      nmatch[q]++;
                    }
                  }else{
                    ++ftemp;
                  }
                }else{
                  ++ftemp;
                }
              }
            }
            
            if (f==g){
              if (k+1-zd[k+1] == L){
                --f;
                //store divergence
                dZ[prefix[f][k+1]] = k+1-L;
              }
              if (k+1-bd[k+1] == L){
                //store divergence
                dZ[prefix[g][k+1]] = k+1-L;
                ++g;
              }
            }
            if (f!=g) {
              while (divergence[f][k+1] <= k+1 - L){
                --f;
                //store divergence
                dZ[prefix[f][k+1]] = k+1-L;
              }
              while (g<M-qM && divergence[g][k+1] <= k+1-L){
                //store divergence
                dZ[prefix[g][k+1]] = k+1-L;
                ++g;
              }
            }
        }
        
        while (f != g){
          int end2=N-1;
          if(times==0){
            int start2=dZ[prefix[f][N]];
            donoridtemp.push_back(prefix[f][N]);
            startpostemp.push_back(start2);
            endpostemp.push_back(end2);
            ++f;
            for(int q=start2;q<=end2;++q){
              nmatch[q]++;
            }
          }else{
            if(addmatch[end2]){
              int start2=dZ[prefix[f][N]];
              //add new matches with new L
              if(end2-start2+1<prevL){
                donoridtemp.push_back(prefix[f][N]);
                startpostemp.push_back(start2);
                endpostemp.push_back(end2);
                ++f;
                for(int q=start2;q<=end2;++q){
                  nmatch[q]++;
                }
              }else{
                ++f;
              }
            }else{
              ++f;
            }
          }
          
        }
        
        if(L<=L_minmatch){
          allsnpmatches = true;
        }else{
          if(times==0){
            // find which SNPs don't have minmatch matches
            for(int j=0;j<N;++j){
              if (nmatch[j] < minmatch) {
                nomatchsnp.push_back(j);
              }
            }
          }else{
            vector<int> nomatchsnptemp=nomatchsnp;
            nomatchsnp.clear();
            // find which SNPs still don't have minmatch matches
            for(int j=0;j<nomatchsnptemp.size();++j){
              if (nmatch[nomatchsnptemp[j]] < minmatch) {
                nomatchsnp.push_back(nomatchsnptemp[j]);
              }
            }
          }
          if(nomatchsnp.size()==0){
            // stop when all SNPs have minmatch matches
            allsnpmatches = true;
          }else{
            // update L
            prevL=L;
            L=(prevL+1)/2; 
            if(L<L_minmatch) L=L_minmatch;
            times++;
          }
        }
        
      }
      
      delete [] t;
      delete [] zd;
      delete [] bd;
      
      
      //below we remove shorter matches while ensuring at least minmatch matches at each SNP
      //information of matches is stored in donoridtemp, startpostemp and endpostemp
      //number of matches at each SNP are stored in nmatch
      //we first sort the genetic distance of each match
      
      vector<double> gdmatch(startpostemp.size());
      for(int mi=0; mi<startpostemp.size();++mi){
        gdmatch[mi]=gd[endpostemp[mi]]-gd[startpostemp[mi]];
      }
      vector<int> length_order=getorder(gdmatch);
      
      vector<int> fullidx; // record which SNP fewer than only minmatch matches
      vector<int> nmatch_output;
      for(int q=0;q<N;++q){
        fullidx.push_back(q);
        nmatch_output.push_back(0);
      }
      for(int mi=length_order.size()-1;mi>=0;--mi){
        
        int starttemp=startpostemp[length_order[mi]];
        int endtemp=endpostemp[length_order[mi]];
        
        if(containsIndex(fullidx,starttemp,endtemp)){
          local_startpos.push_back(starttemp);
          local_endpos.push_back(endtemp);
          local_donorid.push_back(donoridtemp[length_order[mi]]);
          for(int q=starttemp;q<=endtemp;++q){
            nmatch_output[q]++;
            if(nmatch_output[q]==minmatch){
              auto it = remove(fullidx.begin(), fullidx.end(), q);
              fullidx.erase(it, fullidx.end());
            }
          }
        }
        if(fullidx.size()==0) {
          break;
        }
      }
      
      LoopResult result;
      result.donorid = local_donorid;
      result.startpos = local_startpos;
      result.endpos = local_endpos;
      result.queryid = local_startpos.size();
      allResults[idx] = result;
      
      //record the position of the next start position of query haplotype
      //such that we know how many matches are there for this query haplotype
    }
    
    nind_left=nind_left-ncores_use;
  }
  
  
  for (const auto& result : allResults) {
    donorid.insert(donorid.end(), result.donorid.begin(), result.donorid.end());
    startpos.insert(startpos.end(), result.startpos.begin(), result.startpos.end());
    endpos.insert(endpos.end(), result.endpos.begin(), result.endpos.end());
    queryidall.push_back(queryidall.back() + result.queryid);
  }
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> results(queryidall,donorid,startpos,endpos);
  return(results);
}


tuple<vector<int>,vector<int>,vector<int>,vector<int>> do_pbwt(int& L_initial, 
                                                               vector<double> gd,
                                                               vector<int>& queryidx,
                                                               int ncores,
                                                               const int M,
                                                               const int N,
                                                               const int qM,
                                                               int minmatch,
                                                               int L_minmatch,
                                                               const string reffile,
                                                               const string targetfile,
                                                               const bool haploid,
                                                               const bool phase){
  
  int nrow_panel;
  bool samefile;
  if(reffile==targetfile){
    nrow_panel=M-qM;
    samefile=true;
  }else{
    nrow_panel=M;
    samefile=false;
  }
  vector<vector<bool>> panel;
  int **prefix, **divergence, **u, **v;
  
  if (!phase) {
    panel = vector<vector<bool>>(nrow_panel, vector<bool>(N)); 
    ReadVCF(reffile,targetfile,panel,N,M,qM,haploid);
  }else{
    Readphase_donor(reffile,panel,N,M,qM);
  }
  
  cout<<"Begin building PBWT for reference haplotypes"<<endl;
  
  prefix = new int*[M-qM];
  divergence = new int*[M-qM];
  u = new int*[M-qM];
  v = new int*[M-qM];
  int *temp1 = new int[(long long)(M-qM)*(N+1)];
  int *temp2 = new int[(long long)(M-qM)*(N+1)];
  int *temp3 = new int[(long long)(M-qM)*(N)];
  int *temp4 = new int[(long long)(M-qM)*(N)];
  for (long long i = 0; i<M-qM; i++){
    prefix[i] = &(temp1[i*(N+1)]);
    divergence[i] = &(temp2[i*(N+1)]);
    u[i] = &(temp3[i*(N)]);
    v[i] = &(temp4[i*(N)]);
  }
  
  PBWT(panel, prefix, divergence, u, v, M-qM, N);
  
  cout<<"Finish building PBWT for reference haplotypes"<<endl;
  
  while(L_initial>N){
    L_initial=ceil(L_initial/2);
    cout<<"Initial L cannot be greater than N, reducing L to "<<L_initial<<endl;
  }
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> matchresults=longMatchpbwt(L_initial,panel,prefix,
                                                                                    divergence,u,v,minmatch,
                                                                                    gd,queryidx,N,M,qM,
                                                                                    L_minmatch,ncores,samefile,
                                                                                    phase,targetfile);
  
  free_PBWT_memory(panel, prefix, divergence, u, v);
  
  return(matchresults);
  
}


tuple<vector<int>,vector<int>,vector<int>,vector<int>> readmatchfile(const string matchfile){
  vector<int> targetid_temp, donorid_temp, startpos_temp, endpos_temp; 
  vector<int> donorid, startpos, endpos; 
  vector<int> queryidall={0};
  
  ifstream file(matchfile);
  if (!file.is_open()) {
    cerr << "Failed to open the file: " << matchfile << endl;
    return {};
  }
  string line;
  while (getline(file, line)) {
    int t, d, s, e;
    // Assuming columns are space separated.
    istringstream lineStream(line);
    string dummy1;
    int dummy2; // Additional columns
    lineStream >> dummy1 >> t >> d >> s >> e >> dummy2;
    targetid_temp.push_back(t);
    donorid_temp.push_back(d);
    startpos_temp.push_back(s);
    endpos_temp.push_back(e);
  }
  
  file.close();
  
  vector<double> targetid_temp_double(targetid_temp.begin(), targetid_temp.end());
  vector<int> target_order = getorder(targetid_temp_double);
  
  int targetid_idx=0;
  for(int di=0;di<target_order.size();++di){
    donorid.push_back(donorid_temp[target_order[di]]);
    startpos.push_back(startpos_temp[target_order[di]]);
    endpos.push_back(endpos_temp[target_order[di]]-1);
    if(targetid_temp[target_order[di]]!=targetid_idx){
      targetid_idx++;
      queryidall.push_back(di);
    }
  }
  queryidall.push_back(target_order.size());
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> results(queryidall,donorid,startpos,endpos);
  return(results);
}


vector<vector<int>> get_matchdata(vector<int> queryidall,
                                  vector<int> donorid,
                                  vector<int> startpos,
                                  vector<int> endpos,
                                  int queryid,
                                  bool loo=false,
                                  bool haploid=false){
  int querystart=queryidall[queryid];
  int nextquerystart=queryidall[queryid+1];
  int nrow_match=nextquerystart-querystart;
  vector<vector<int>> matchinfo(nrow_match,vector<int>(3));
  // leave out the same haplotype
  if(loo){
    if(haploid){
      for(int p=querystart;p<nextquerystart;++p){
        if(donorid[p]!=queryid){
          matchinfo[p-querystart][0]=donorid[p];
          matchinfo[p-querystart][1]=startpos[p];
          matchinfo[p-querystart][2]=endpos[p];
        }
      }
    }else{
      int queryid_pair;
      if(queryid%2==0){
        queryid_pair=queryid+1;
      }else{
        queryid_pair=queryid-1;
      }
      for(int p=querystart;p<nextquerystart;++p){
        if(donorid[p]!=queryid && donorid[p]!=queryid_pair){
          matchinfo[p-querystart][0]=donorid[p];
          matchinfo[p-querystart][1]=startpos[p];
          matchinfo[p-querystart][2]=endpos[p];
        }
      }
    }
  }else{
    for(int p=querystart;p<nextquerystart;++p){
      matchinfo[p-querystart][0]=donorid[p];
      matchinfo[p-querystart][1]=startpos[p];
      matchinfo[p-querystart][2]=endpos[p];
    }
  }
  return(matchinfo);
}

/////////////////////end of pbwt contents///////////////////////////



double vec_sum(const vector<double>& l){
  // calculate the sum of a vector
  double sum_l=0;
  for(int i=0;i<l.size();++i){
    sum_l=sum_l+l[i];
  }
  return(sum_l);
}

vector<double> vec_multiply(const vector<double>& l1, 
                            const vector<double>& l2){
  // calculate the multiply of two vectors
  vector<double> multi_l(l1.size());
  for(int i=0;i<l1.size();++i){
    multi_l[i]=l1[i]*l2[i];
  }
  return(multi_l);
}

vector<int> union_vec(const vector<int>& vec1, 
                      const vector<int>& vec2) {
  // get the union of two vectors
  vector<int> result(vec1.size() + vec2.size());
  vector<int>::iterator it;
  it = set_union(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), result.begin());
  result.resize(it - result.begin());
  return(result);
}

void update_vec(vector<double>& a, 
                const vector<int>& b, 
                const vector<double>& c) {
  // change the position of b inside a to c
  for (int i = 0; i < b.size(); ++i) {
    a[b[i]] = c[i];
  }
}

double hVecSum(hVec& v,
               vector<int> idx={}){
  // calculate the sum of the values of a hVec with given indices
  if (idx.empty()) {
    idx = v.k;
  }
  return(vec_sum(v.getall(idx)));
}

hVec hVecCirc(hVec& v1, hVec& v2){
  // calculate the circ of two hVec;
  hVec ret=hVec();
  ret.setdefault(v1.x0*v2.x0);
  vector<int> v1idx=v1.k;
  vector<int> v2idx=v2.k;
  vector<int> idx=union_vec(v1idx,v2idx);
  vector<double> v1val=v1.getall(idx);
  vector<double> v2val=v2.getall(idx);
  ret.setall(idx,vec_multiply(v1val,v2val));
  return(ret);
}

void hVecScale(hVec& v,
               const double x){
  //Scale a hash vector by x
  vector<int> idx=v.k;
  vector<double> val=v.getall(idx);
  v.setdefault(x*v.x0);
  vector<double> xvec(idx.size(),x);
  v.setall(idx,vec_multiply(val,xvec));
}


vector<double> hVecdense(hVec& x,
                         const int nrow){
  //Convert a hash vector x into a dense vector with length nrow
  vector<double> vdense(nrow, x.x0);
  vector<int> k=x.k;
  vector<double> vval=x.getall(k);
  update_vec(vdense,k,vval);
  return(vdense);
}


vector<vector<double>> hMatrix2matrix(hMat& mat){
  int nrow=mat.d1;
  int ncol=mat.d2;
  vector<vector<double>> densematrix(nrow,vector<double>(ncol));
  vector<double> vecdence(nrow);
  for(int j=0;j<ncol;++j){
    vecdence=hVecdense(mat.m[j],nrow);
    for(int i=0;i<nrow;++i){
      densematrix[i][j]=vecdence[i];
    }
  }
  return(densematrix);
}

vector<int> randomsample(const vector<int>& popidx, 
                         const int number) {
  // sample number from popidx
  if(number>popidx.size()) cout<<"Number cannot be greater than the size of popidx"<<endl;
  // Initialize the random number generator
  random_device rd;
  mt19937 gen(rd());
  
  // Shuffle the elements of the vector randomly
  vector<int> shuffled_popidx = popidx;
  shuffle(shuffled_popidx.begin(), shuffled_popidx.end(), gen);
  
  // Take the first number elements of the shuffled vector to create the random sample
  vector<int> random_sample(shuffled_popidx.begin(), shuffled_popidx.begin() + number);
  
  return random_sample;
}

vector<double> cal_sameprob(const int nsnp, 
                            const double lambda, 
                            vector<double>& gd,
                            const int nref){
  //compute the sameprob
  vector<double> sameprob(nsnp);
  for(int j=0;j<nsnp-1;++j){
    sameprob[j]=exp(-lambda*(gd[j+1]-gd[j]));
    if(sameprob[j]>1-nref*0.000000000000002){
      sameprob[j]=1-nref*0.000000000000002;
    } 
  }
  return(sameprob);
}

vector<double> cal_otherprob(const int nref, 
                             const vector<double>& sameprob){
  //compute the otherprob
  int nsnp=sameprob.size();
  vector<double> otherprob(nsnp);
  for(int j=0;j<nsnp;++j){
    otherprob[j]=(1-sameprob[j])/nref;
  }
  return(otherprob);
}


tuple<hMat, vector<double>> forwardProb(const hMat& mat,
                                        const vector<double>& sameprob, 
                                        const vector<double>& otherprob){
  //compute normalised forward probability and store in hMat
  int nrow=mat.d1;
  int ncol=mat.d2;
  
  double sum_use = 0.0;
  
  for (int i = 1; i <= nrow; i++) {
    sum_use += 1.0 / i;
  }
  double mu = 0.5 / sum_use / (nrow + 1.0 / sum_use);
  
  hMat forward_prob(nrow,ncol,0);
  vector<int> twj=mat.m[0].k;
  if(twj.size()==0){
    for(int i=0;i<nrow;++i){
      twj.push_back(i);
    }
  }
  vector<double> l(twj.size(),1.0/twj.size());
  forward_prob.m[0].setall(twj,l);
  vector<double> logmultF={0};
  vector<double> fprev;
  double sameprobuse;
  double otherprobuse;
  double sumfp;
  for(int j=1;j<ncol;++j){
    twj=mat.m[j].k;
    sameprobuse=sameprob[j-1];
    otherprobuse=otherprob[j-1];
    
    if(twj.size()==0){
      twj=forward_prob.m[j-1].k;
      fprev=forward_prob.m[j-1].getall(twj);
      for(int i=0;i<twj.size();++i){
        if(log(sameprobuse*fprev[i]+otherprobuse)+log(mu)>-34.53){
          forward_prob.m[j].set(twj[i],(sameprobuse*fprev[i]+otherprobuse)*mu);
        }else{
          forward_prob.m[j].set(twj[i],exp(-34.53));
        }
      }
    }else{
      fprev=forward_prob.m[j-1].getall(twj);
      for(int i=0;i<twj.size();++i){
        forward_prob.m[j].set(twj[i],sameprobuse*fprev[i]+otherprobuse);
      }
    }
    sumfp=hVecSum(forward_prob.m[j],twj);
    logmultF.push_back(log(sumfp)+logmultF[j-1]);
    hVecScale(forward_prob.m[j],1.0/sumfp);
  }
  tuple<hMat, vector<double>> results(forward_prob, logmultF);
  return(results);
}





tuple<hMat, vector<double>> backwardProb(const hMat& mat,
                                         const vector<double>& sameprob, 
                                         const vector<double>& otherprob){
  //compute normalised backward probability and store in hMat
  int nrow=mat.d1;
  int ncol=mat.d2;
  hMat backward_prob(nrow,ncol,1.0/nrow);
  vector<int> twj;
  vector<double> Bjp1;
  double sumBjp1;
  double sameprobuse;
  double otherprobuse;
  vector<double> logmultB(ncol,log(nrow));
  
  double sum_use = 0.0;
  
  for (int i = 1; i <= nrow; i++) {
    sum_use += 1.0 / i;
  }
  double mu = 0.5 / sum_use / (nrow + 1.0 / sum_use);
  
  for(int j=ncol-2;j>=0;--j){
    sameprobuse=sameprob[j];
    otherprobuse=otherprob[j];
    twj=mat.m[j+1].k;
    if(twj.size()==0){
      twj=backward_prob.m[j+1].k;
      Bjp1=backward_prob.m[j+1].getall(twj);
      for(int i=0;i<Bjp1.size();++i){
        if(log(Bjp1[i])+log(mu)>-34.53){
          Bjp1[i]=Bjp1[i]*mu;
        }else{
          Bjp1[i]=exp(-34.53);
        }
      }
      double default_val=backward_prob.m[j+1].x0;
      sumBjp1=vec_sum(Bjp1)+default_val*mu*(nrow-Bjp1.size());
      vector<double> val(twj.size());
      for(int i=0;i<twj.size();++i){
        val[i]=sameprobuse*Bjp1[i]+otherprobuse*sumBjp1;
        backward_prob.m[j].set(twj[i],val[i]);
      }
      
      double default_backward=default_val*mu*sameprobuse+otherprobuse*sumBjp1;
      double sum_nc=vec_sum(val)+default_backward*(nrow-twj.size());
      logmultB[j]=log(sum_nc)+logmultB[j+1];
      hVecScale(backward_prob.m[j],1.0/sum_nc);
      backward_prob.m[j].setdefault(default_backward/sum_nc);
    }else{
      Bjp1=backward_prob.m[j+1].getall(twj);
      sumBjp1=vec_sum(Bjp1);
      vector<double> val(twj.size());
      for(int i=0;i<twj.size();++i){
        val[i]=sameprobuse*Bjp1[i]+otherprobuse*sumBjp1;
        backward_prob.m[j].set(twj[i],val[i]);
      }
      logmultB[j]=log(vec_sum(val)+otherprobuse*sumBjp1*(nrow-twj.size()))+logmultB[j+1];
      hVecScale(backward_prob.m[j],1.0/sumBjp1);
      backward_prob.m[j].setdefault(otherprobuse);
    }
    
  }
  tuple<hMat, vector<double>> results(backward_prob, logmultB);
  return(results);
}


hMat marginalProb(hMat& f, hMat& b){
  //compute normalised marginal probability and store in hMat
  int nrow=f.d1;
  int ncol=f.d2;
  hMat marginal_prob(nrow,ncol,1.0/nrow);
  for(int j=0;j<ncol;++j){
    marginal_prob.m[j]=hVecCirc(f.m[j],b.m[j]);
    hVecScale(marginal_prob.m[j],1.0/hVecSum(marginal_prob.m[j]));
  }
  
  return(marginal_prob);
}

void removeRowsWithValue(vector<vector<int>>& data, 
                         const vector<int>& values) {
  //remove the rows of data which first row is included in values
  auto it = remove_if(data.begin(), data.end(), [&](const vector<int>& row) {
    return find(values.begin(), values.end(), row[1]) != values.end();
  });
  data.erase(it, data.end());
}

hMat matchfiletohMat(const vector<vector<int>>& matchdata, 
                     const int& nref, 
                     const int& nsnp){
  hMat mat(nref,nsnp,0.0);
  int q=0;
  for(int i = 0; i < matchdata.size(); ++i) {
    int val = matchdata[i][0];
    for(int j = matchdata[i][1]; j <= matchdata[i][2]; ++j) {
      mat.m[j].set(val, 1.0);
    }
  }
  return(mat);
}

tuple<vector<double>,vector<double>> readmap(const string& mapfile) {
  ifstream file(mapfile);
  vector<double> pd,gd;
  
  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << mapfile << endl;
    abort();
  }
  
  string line;
  double column1,column2;
  
  // Read and discard the header line
  getline(file, line);
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> column1;
    lineStream >> column2;
    
    pd.push_back(column1);
    gd.push_back(column2);
  }
  
  file.close();
  tuple<vector<double>,vector<double>> output(pd,gd);
  return output;
}

tuple<vector<string>,vector<int>> readpopfile(const string& popfile) {
  ifstream file(popfile);
  vector<string> indnames;
  vector<int> refindex;
  
  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << popfile << endl;
    abort();
  }
  
  string line;
  string column1;
  int column2;
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> column1;
    lineStream >> column2;
    indnames.push_back(column1);
    refindex.push_back(column2);
  }
  
  file.close();
  tuple<vector<string>,vector<int>> output(indnames,refindex);
  return output;
}


vector<string> readtargetname(const string& targetname) {
  ifstream file(targetname);
  vector<string> tgnames;
  
  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << targetname << endl;
    abort();
  }
  
  string line;
  string column1;
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> column1;
    
    tgnames.push_back(column1);
  }
  
  file.close();
  return tgnames;
}



double est_lambda_EM(hMat& mat, 
                     vector<double>& gd,
                     const int ite_time){
  //estimate \lambda using EM algorithm
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob;
  vector<double> otherprob; 
  double lambda_ite=400000/nref;
  vector<double> gl(nsnp-1);
  vector<double> lambda_each((nsnp-1));
  for(int j=0;j<nsnp-1;++j){
    gl[j]=gd[j+1]-gd[j];
  }
  double totalgd=gd[nsnp-1]-gd[0];
  for(int t=0;t<ite_time;++t){
    sameprob=cal_sameprob(nsnp,lambda_ite,gd,nref);
    otherprob=cal_otherprob(nref,sameprob);
    
    tuple<hMat, vector<double>> f=forwardProb(mat,sameprob,otherprob);
    hMat forward_prob=get<0>(f);
    vector<double> logmultF=get<1>(f);
    
    tuple<hMat, vector<double>> b=backwardProb(mat,sameprob,otherprob);
    hMat backward_prob=get<0>(b);
    vector<double> logmultB=get<1>(b);
    vector<double> u((nsnp-1),0);
    for(int j=0;j<nsnp-1;++j){
      double al = exp(logmultF[j+1]+logmultB[j+1]-logmultF[nsnp-1]);
      double ar = exp(logmultF[j]+logmultB[j+1]-logmultF[nsnp-1]);
      vector<int> twj=mat.m[j+1].k;
      vector<double> valuef=forward_prob.m[j+1].getall(twj);
      vector<double> valuefprev=forward_prob.m[j].getall(twj);
      vector<double> valueb=backward_prob.m[j+1].getall(twj);
      
      for(int i=0;i<twj.size();++i){
        u[j]=u[j]+al*valuef[i]*valueb[i]-ar*valuefprev[i]*valueb[i]*sameprob[j];
      }
    }
    for(int j=0;j<nsnp-1;++j){
      lambda_each[j]=u[j]*lambda_ite*gl[j]/(1-sameprob[j]);
    }
    lambda_ite=vec_sum(lambda_each)/totalgd;
  }
  if(isnan(lambda_ite)) lambda_ite=1/totalgd;
  return(lambda_ite);
}

double est_lambda_Viterbi(const vector<int>& startpos, 
                          const vector<int>& endpos,
                          const int nsnp, 
                          const double gdall) {
  int nrec = 0;
  int j = 0;
  unordered_map<int, int> maxEndMap;
  for (size_t i = 0; i < startpos.size(); i++) {
    int currStart = startpos[i];
    int currEnd = endpos[i];
    if (maxEndMap.find(currStart) == maxEndMap.end()) {
      maxEndMap[currStart] = currEnd;
    } else {
      maxEndMap[currStart] = max(maxEndMap[currStart], currEnd);
    }
  }
  
  int maxend = -1;
  int maxkidx = 0;
  while (j < nsnp-1) {
    if(j>0){
      maxkidx = maxend+1;
    }
    int maxendnew=maxend+1;
    for (int k = j; k <= maxend+1; k++) {
      if (maxEndMap.find(k) != maxEndMap.end()) {
        int maxvalue=maxEndMap[k];
        if(maxvalue>=maxendnew){
          maxendnew=maxvalue;
          maxkidx=k;
        }
      }
    }
    maxend=maxendnew;
    if (maxend < nsnp-1) {
      nrec++;
    }
    j = maxkidx+1;
    if(maxend==nsnp-1 || maxend==nsnp-2) j=nsnp-1;
  }
  if(nrec==0){
    nrec=1;
  } 
  double lambda_est = nrec / static_cast<double>(gdall);
  return lambda_est;
}


double est_lambda_EM_average(const hAnc& refidx, 
                             const int nref, 
                             const int nsnp,
                             vector<double>& gd,
                             int L_initial,
                             int minmatch,
                             int L_minmatch,
                             int ncores,
                             const double indfrac,
                             const int ite_time, 
                             const int minsnpEM, 
                             const double EMsnpfrac,
                             bool haploid,
                             const string reffile,
                             const bool phase,
                             const bool leaveoneout,
                             tuple<vector<int>,vector<int>,vector<int>,vector<int>> pbwtall_ref=
                               make_tuple(vector<int>(), vector<int>(), vector<int>(), vector<int>()))
{
  // estimate \lambda from the reference panel
  int npop=refidx.pos.size();
  vector<double> lambda_est;
  int count=0;
  double gdall=gd[nsnp-1]-gd[0];
  
  
  vector<int> allsamples;
  vector<int> popstart={0}; //the start position of diffent population samples
  
  //stratified sampling
  for(int i=0;i<npop;++i){
    //randomly sample a percentage of indfrac reference samples
    vector<int> popidx=refidx.findrows(i);
    
    int nref_sample=static_cast<int>(ceil(popidx.size()*indfrac));  //the number of samples of this ancestry
    
    vector<int> samples=randomsample(popidx,nref_sample); // sample index of this ancestry
    
    for(int j=0;j<nref_sample;++j){
      allsamples.push_back(samples[j]);
    }
    popstart.push_back(allsamples.size());
  }
  
  if (get<0>(pbwtall_ref).empty() && 
      get<1>(pbwtall_ref).empty() && 
      get<2>(pbwtall_ref).empty() && 
      get<3>(pbwtall_ref).empty()){
    
    vector<int> queryidx;
    
    for(int i=0;i<nref;++i){
      queryidx.push_back(i);
    }
    
    
    cout<<"Begin doing PBWT and finding matches for reference haplotypes"<<endl;
    tuple<vector<int>,vector<int>,vector<int>,vector<int>> pbwtall_ref=do_pbwt(L_initial, gd,queryidx,
                                                                               ncores,nref,nsnp,0,
                                                                               minmatch,L_minmatch,
                                                                               reffile,reffile,haploid,phase);
    
    cout<<"Finish finding matches with PBWT"<<endl;
  }
  
  vector<int> queryidall=get<0>(pbwtall_ref);
  vector<int> donorid_ref=get<1>(pbwtall_ref);
  vector<int> startpos_ref=get<2>(pbwtall_ref);
  vector<int> endpos_ref=get<3>(pbwtall_ref);
  
  omp_set_num_threads(ncores); 
  for(int i=0;i<npop;++i){
    vector<int> samples;
    for(int j=popstart[i];j<popstart[i+1];++j){
      samples.push_back(allsamples[j]);
    }
    
    vector<vector<vector<int>>> match_use(samples.size());
    
    for (int k = 0; k < samples.size(); ++k) {
      match_use[k] = get_matchdata(queryidall,donorid_ref,startpos_ref,endpos_ref,samples[k], true,haploid);
    }// this reduces memory
    
#pragma omp parallel for reduction(+:count)
    for(int k=0;k<samples.size();++k){
      //leave-one-out
      vector<vector<int>> matchdata=match_use[k];
      vector<int> removeidx;
      
      if(leaveoneout){
        for(int j=0;j<npop;++j){
          if(j!=i){
            int rmidx1=randomsample(refidx.findrows(j),1)[0];
            removeidx.push_back(rmidx1);
            if(!haploid){
              int rmidx2;
              if(rmidx1%2==0){
                rmidx2=rmidx1+1;
              }else{
                rmidx2=rmidx1-1;
              }
              removeidx.push_back(rmidx2);
            }
          }
        }
        removeRowsWithValue(matchdata,removeidx);
      }
      
      hMat mat=matchfiletohMat(matchdata,nref-npop,nsnp);
      int nsnp_use=nsnp;
      if(nsnp>minsnpEM){
        if(nsnp*EMsnpfrac<minsnpEM){
          nsnp_use=minsnpEM;
        }else{
          nsnp_use=ceil(nsnp*EMsnpfrac);
        }
      }
      if(nsnp_use<nsnp){
        vector<int> allsnpidx;
        for (int i = 0; i < nsnp; i++) {
          allsnpidx.push_back(i);
        }
        vector<int> snp_use=randomsample(allsnpidx,nsnp_use);
        sort(snp_use.begin(), snp_use.end());
        hMat mat_use(nref-npop,nsnp_use,0.0);
        vector<double> gd_use(nsnp_use);
        for(int j=0; j<nsnp_use;++j){
          gd_use[j]=gd[snp_use[j]];
          vector<int> non0idx=mat.m[snp_use[j]].k;
          mat_use.m[j].setall(non0idx,vector<double>(non0idx.size(),1.0));
        }
        
        double lambda_estimated=est_lambda_EM(mat_use,gd_use,ite_time);
        count++;
        lambda_est.push_back(lambda_estimated);
      }else{
        
        double lambda_estimated=est_lambda_EM(mat,gd,ite_time);
        count++;
        lambda_est.push_back(lambda_estimated);
      }
    }
  }
  double lambda_ave=vec_sum(lambda_est)/lambda_est.size();
  return(lambda_ave);
}

hMat indpainting(const hMat& mat,
                 vector<double>& gd, 
                 const double lambda,
                 const int npop, 
                 const vector<int>& refindex,
                 tuple<hMat, vector<double>> forwardprob,
                 tuple<hMat, vector<double>> backwardprob){
  // return individual painting
  int nsnp=mat.d2;
  
  hMat marginal_prob=marginalProb(get<0>(forwardprob),get<0>(backwardprob));
  // compute individual painting in terms of different populations
  hMat marginal_prob_pop(npop,nsnp,0.0);
  
  for(int j=0; j<nsnp; ++j){
    // find the positions for non-zero marginal probability
    vector<int> margidx=marginal_prob.m[j].k;
    vector<double> popprob(npop,0);
    int popidx;
    int refnumberidx;
    for(int k=0;k<margidx.size();++k){
      // the number of the reference samples
      refnumberidx=margidx[k]; 
      // search this reference sample belongs to which population
      popidx=refindex[refnumberidx]; 
      // update the probability of this population
      popprob[popidx]=popprob[popidx]+marginal_prob.m[j].get(refnumberidx);
    }
    for(int i=0; i<npop; ++i){
      if(popprob[i]>=0.005){
        marginal_prob_pop.m[j].set(i,popprob[i]);
      }
    }
  }
  return(marginal_prob_pop);
}


vector<double> chunklength_each(vector<double>& gd, 
                                hMat& mat, 
                                const double lambda, 
                                const int npop,
                                const vector<int>& refindex,
                                tuple<hMat, vector<double>> forwardprob,
                                tuple<hMat, vector<double>> backwardprob,
                                const double gdall){
  //calculate chunk length for each haplotype
  int nsnp=mat.d2;
  vector<double> gl(nsnp-1);
  for(int j=0;j<nsnp-1;++j){
    gl[j]=gd[j+1]-gd[j];
  }
  
  hMat forward_prob=get<0>(forwardprob);
  vector<double> logmultF=get<1>(forwardprob);
  
  hMat backward_prob=get<0>(backwardprob);
  vector<double> logmultB=get<1>(backwardprob);
  
  vector<double> suml(npop,0.0);
  
  for(int j=0;j<nsnp-1;++j){
    double wl = exp(logmultF[j]+logmultB[j]-logmultF[nsnp-1]);
    double wr = exp(logmultF[j+1]+logmultB[j+1]-logmultF[nsnp-1]);
    vector<int> twj1=mat.m[j].k;
    vector<int> twj2=mat.m[j+1].k;
    vector<double> valuefprev=forward_prob.m[j].getall(twj1);
    vector<double> valuef=forward_prob.m[j+1].getall(twj2);
    vector<double> valuebprev=backward_prob.m[j].getall(twj1);
    vector<double> valueb=backward_prob.m[j+1].getall(twj2);
    
    // different chunk length for different populations
    vector<double> sumlleft(npop,0.0);
    vector<double> sumlright(npop,0.0);
    
    for(int i=0;i<twj1.size();++i){
      sumlleft[refindex[twj1[i]]]+=wl*valuefprev[i]*valuebprev[i];
    }
    for(int i=0;i<twj2.size();++i){
      sumlright[refindex[twj2[i]]]+=wr*valuef[i]*valueb[i];
    }
    for(int k=0;k<npop;++k){
      suml[k]+=(sumlleft[k]+sumlright[k])*gl[j]*0.5;
    }
  }
  
  double totalsum=0; //need to standardize when encountering 0 matches at any position
  for(int k=0;k<npop;++k){
    totalsum+=suml[k];
  }
  double standardize_factor=gdall*100/totalsum;
  for(int k=0;k<npop;++k){
    suml[k]=suml[k]*standardize_factor;
  }
  return(suml);
}


void doLDAS(hMat &LDA_result,
            const string LDASfile,
            const double window,
            const vector<double> &gd,
            const vector<double> &pd,
            const vector<int> &nsnp_left,
            const vector<int> &nsnp_right,
            const int nsnp){
  
  // calculate LDA score
  vector<double> LDAS_score(nsnp);
  vector<double> LDAS_upper(nsnp);
  vector<double> LDAS_lower(nsnp);
#pragma omp parallel for
  for(int i=0;i<nsnp;++i){
    vector<double> gdgap;
    vector<double> LDA_ave;
    vector<double> LDA_upper;
    vector<double> LDA_lower;
    for(int j=i-nsnp_left[i];j<=i+nsnp_right[i]-1;++j){
      gdgap.push_back(gd[j+1]-gd[j]);
      LDA_ave.push_back((LDA_result.m[min(i,j)].get(max(i,j))+LDA_result.m[min(i,j+1)].get(max(i,j+1)))/2);
      LDA_upper.push_back(max(LDA_result.m[min(i,j)].get(max(i,j)),LDA_result.m[min(i,j+1)].get(max(i,j+1))));
      LDA_lower.push_back(min(LDA_result.m[min(i,j)].get(max(i,j)),LDA_result.m[min(i,j+1)].get(max(i,j+1))));
    }
    double left_distance=gd[i]-gd[i-nsnp_left[i]];
    double right_distance=gd[i+nsnp_right[i]]-gd[i];
    if(i-nsnp_left[i]>0 && i+nsnp_right[i]<nsnp){
      gdgap.push_back(window-left_distance);
      gdgap.push_back(window-right_distance);
      LDA_ave.push_back((LDA_result.m[i-nsnp_left[i]].get(i)+LDA_result.m[i-nsnp_left[i]-1].get(i))/2);
      LDA_ave.push_back((LDA_result.m[i].get(i+nsnp_right[i])+LDA_result.m[i].get(i+nsnp_right[i]+1))/2);
      LDA_upper.push_back(max(LDA_result.m[i-nsnp_left[i]].get(i),LDA_result.m[i-nsnp_left[i]-1].get(i)));
      LDA_upper.push_back(max(LDA_result.m[i].get(i+nsnp_right[i]),LDA_result.m[i].get(i+nsnp_right[i]+1)));
      LDA_lower.push_back(min(LDA_result.m[i-nsnp_left[i]].get(i),LDA_result.m[i-nsnp_left[i]-1].get(i)));
      LDA_lower.push_back(min(LDA_result.m[i].get(i+nsnp_right[i]),LDA_result.m[i].get(i+nsnp_right[i]+1)));
    }
    
    if(i-nsnp_left[i]==0 && i+nsnp_right[i]<nsnp){
      // right window
      gdgap.push_back(window-right_distance);
      LDA_ave.push_back((LDA_result.m[i].get(i+nsnp_right[i])+LDA_result.m[i].get(i+nsnp_right[i]+1))/2);
      LDA_upper.push_back(max(LDA_result.m[i].get(i+nsnp_right[i]),LDA_result.m[i].get(i+nsnp_right[i]+1)));
      LDA_lower.push_back(min(LDA_result.m[i].get(i+nsnp_right[i]),LDA_result.m[i].get(i+nsnp_right[i]+1)));
      // use the right window to estimate the left window
      // window-left_distance is the distance to be estimated from the right window
      double gdright_add=0;
      int j=i+nsnp_right[i];
      
      while(gdright_add < window-left_distance){
        // how long distance from enough
        double distance_from_enough = window-left_distance-gdright_add;
        gdright_add=window+gd[i]-gd[j];
        if(gdright_add <= window-left_distance){
          // estimated distance still not enough or just enough
          if(j==i+nsnp_right[i]){
            gdgap.push_back(window-right_distance);
          }else{
            gdgap.push_back(gd[j+1]-gd[j]);
          }
        }else{
          // estimated distance is enough
          gdgap.push_back(distance_from_enough);
        }
        LDA_ave.push_back((LDA_result.m[min(i,j)].get(max(i,j))+LDA_result.m[min(i,j+1)].get(max(i,j+1)))/2);
        LDA_upper.push_back(max(LDA_result.m[min(i,j)].get(max(i,j)),LDA_result.m[min(i,j+1)].get(max(i,j+1))));
        LDA_lower.push_back(min(LDA_result.m[min(i,j)].get(max(i,j)),LDA_result.m[min(i,j+1)].get(max(i,j+1))));
        j=j-1;
        //endwhile
      }
      //endif
    }
    
    
    if(i-nsnp_left[i]>0 && i+nsnp_right[i]==nsnp){
      // left window
      gdgap.push_back(window-left_distance);
      LDA_ave.push_back((LDA_result.m[i-nsnp_left[i]].get(i)+LDA_result.m[i-nsnp_left[i]-1].get(i))/2);
      LDA_upper.push_back(max(LDA_result.m[i-nsnp_left[i]].get(i),LDA_result.m[i-nsnp_left[i]-1].get(i)));
      LDA_lower.push_back(min(LDA_result.m[i-nsnp_left[i]].get(i),LDA_result.m[i-nsnp_left[i]-1].get(i)));
      // use the left window to estimate the right window
      // window-right_distance is the distance to be estimated from the left window
      double gdleft_add=0;
      int j=i-nsnp_left[i];
      
      while(gdleft_add < window-right_distance){
        // how long distance from enough
        double distance_from_enough = window-right_distance-gdleft_add;
        gdleft_add=gd[j]-(gd[i]-window);
        if(gdleft_add <= window-right_distance){
          // estimated distance still not enough or just enough
          if(j==i-nsnp_left[i]){
            gdgap.push_back(window-left_distance);
          }else{
            gdgap.push_back(gd[j]-gd[j-1]);
          }
        }else{
          // estimated distance is enough
          gdgap.push_back(distance_from_enough);
        }
        LDA_ave.push_back((LDA_result.m[min(i,j)].get(max(i,j))+LDA_result.m[min(i,j-1)].get(max(i,j-1)))/2);
        LDA_upper.push_back(max(LDA_result.m[min(i,j)].get(max(i,j)),LDA_result.m[min(i,j-1)].get(max(i,j-1))));
        LDA_lower.push_back(min(LDA_result.m[min(i,j)].get(max(i,j)),LDA_result.m[min(i,j-1)].get(max(i,j-1))));
        j=j+1;
        //endwhile
      }
      //endif
    }
    for(int q=0;q<gdgap.size();++q){
      LDAS_score[i]+=LDA_ave[q]*gdgap[q];
      LDAS_upper[i]+=LDA_upper[q]*gdgap[q];
      LDAS_lower[i]+=LDA_lower[q]*gdgap[q];
    }
  }
  
  //output the LDAS results into LDASfile
  ofstream outputFile(LDASfile);
  if (outputFile.is_open()) {
    outputFile.precision(15);
    outputFile << "physical_position" << " " << "LDAS" << " " << "LDAS_lower" << " " << "LDAS_upper" << "\n";
    for (int i = 0; i < nsnp; ++i) {
      outputFile << fixed<< setprecision(0) << pd[i];
      outputFile << " " << fixed << setprecision(4) << LDAS_score[i] <<" "<< LDAS_lower[i] <<" " << LDAS_upper[i] << "\n";
    }
    outputFile.close();
  } else {
    cerr << "Unable to open file" << LDASfile;
  }
  
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
  
  int nsnp = pd.size();
  int npop = aveSNPpainting.size();
  
  vector<double> p_values(nsnp, 1.0); // initialize p_values with 1
  vector<double> test_statistic(nsnp, 1.0); 
  
  vector<double> mu=rowMeans(aveSNPpainting);
  
  arma::mat Astar(nsnp, npop);
  for (int i = 0; i < npop; ++i) {
#pragma omp parallel for
    for (int j = 0; j < nsnp; ++j) {
      Astar(j, i) = aveSNPpainting[i][j] - mu[i]; // compute A*(j,k) for all j,k and store it in Astar
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



void paintall(const string method,
              bool diff_lambda,
              const double fixlambda,
              const int ite_time,
              const double indfrac,
              const int minsnpEM, 
              const double EMsnpfrac,
              int L_initial,
              int nmatch,
              int L_minmatch,
              bool haploid,
              bool leaveoneout,
              const string reffile,
              const string targetfile,
              const string mapfile,
              const string popfile,
              const string targetname,
              const string matchfile,
              bool outputpainting,
              bool outputaveSNPpainting,
              bool outputaveindpainting,
              bool outputLDA,
              bool outputLDAS,
              bool outputAAS,
              const string probfile,
              const string aveSNPprobfile,
              const string aveindprobfile,
              const string chunklengthfile,
              const string LDAfile,
              const string LDASfile,
              const string AASfile,
              const string lambdafile,
              const double window,
              int ncores,
              const string run,
              bool phase){
  
  int LDAfactor=1;
  
  if(outputLDA||outputLDAS){
    LDAfactor=24/ncores+1;
  }
  
  // read the map data to get the genetic distance in Morgans
  tuple<vector<double>,vector<double>> mapinfo = readmap(mapfile);
  vector<double> gd = get<1>(mapinfo);
  vector<double> pd = get<0>(mapinfo);
  
  tuple<vector<string>,vector<int>> popinfo = readpopfile(popfile);
  
  vector<int> refindex = get<1>(popinfo);
  
  vector<string> indnames = readtargetname(targetname);
  
  int nind=indnames.size();
  
  //adjust refindex
  if(!haploid){
    vector<int> refindex_new;
    for(int i=0;i<refindex.size();++i){
      refindex_new.push_back(refindex[i]);
      refindex_new.push_back(refindex[i]);
    }
    refindex=refindex_new;
  }
  
  vector<int> allind;
  for(int i=0;i<nind;++i){
    allind.push_back(i);
  }
  int nhap_use;
  vector<int> queryidx;
  
  if(haploid){
    for(int i=0;i<nind;++i){
      queryidx.push_back(i);
    }
    nhap_use=nind;
  }else{
    for(int i=0;i<nind;++i){
      queryidx.push_back(2*i);
      queryidx.push_back(2*i+1);
    }
    nhap_use=nind*2;
  }
  
  
  //compute painting for all target individuals
  const int nsnp=gd.size();
  const int nref=refindex.size();
  hAnc refidx(refindex);
  const int npop=refidx.pos.size();
  double lambda;
  double gdall=gd[nsnp-1]-gd[0];
  
  bool loo=false;
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> pbwtall_target;
  
  if (matchfile.empty()){
    int qM;
    
    if(reffile==targetfile){
      loo = true;
      qM=0;
    }else{
      if(haploid){
        qM=nind;
      }else{
        qM=2*nind;
      }
    }
    cout<<"Begin doing PBWT and finding matches for target haplotypes"<<endl;
    pbwtall_target=do_pbwt(L_initial, gd,queryidx,ncores,nref+qM,nsnp,qM,nmatch,
                           L_minmatch,reffile,targetfile,haploid,phase);
    cout<<"Finish finding matches with PBWT"<<endl;
  }else{
    cout<<"Begin reading matches from "<<matchfile<<endl;
    pbwtall_target=readmatchfile(matchfile);
    cout<<"Finish reading matches from "<<matchfile<<endl;
  }
  
  vector<int> queryidall_target=get<0>(pbwtall_target);
  vector<int> donorid_target=get<1>(pbwtall_target);
  vector<int> startpos_target=get<2>(pbwtall_target);
  vector<int> endpos_target=get<3>(pbwtall_target);
  
  // estimate lambda
  
  if(!diff_lambda){
    if(fixlambda!=0){
      lambda=fixlambda;
      cout << "Using fixed lambda "<<lambda<<endl;
    }else{
      cout<<"Begin estimating fixed lambda"<<endl;
      if(method=="EM"){
        if(targetfile==reffile){
          lambda=est_lambda_EM_average(refidx,nref,nsnp,gd,L_initial,nmatch,L_minmatch,ncores,
                                       indfrac,ite_time,minsnpEM,EMsnpfrac,haploid,reffile,phase,leaveoneout,pbwtall_target);
        }else{
          lambda=est_lambda_EM_average(refidx,nref,nsnp,gd,L_initial,nmatch,L_minmatch,ncores,
                                       indfrac,ite_time,minsnpEM,EMsnpfrac,haploid,reffile,phase,leaveoneout);
        }
      }else{
        //estimate lambda as the average of v_nsamples target individuals
        int v_nsamples=static_cast<int>(ceil(queryidx.size()*indfrac));
        vector<int> v_samples=randomsample(queryidx,v_nsamples);
        double lambda_sum = 0.0;
        
        // get the matches before the loop
        int v_nhap_left=v_nsamples;
        while(v_nhap_left>0){
          int v_nsamples_use = (ncores < v_nhap_left) ? ncores : v_nhap_left;
          vector<vector<vector<int>>> v_targetmatch_use(v_nsamples_use);
          
          for (int ii = v_nsamples - v_nhap_left; ii < v_nsamples - v_nhap_left + v_nsamples_use; ++ii) {
            // leave one out if the donor file is the same as the target file
            v_targetmatch_use[ii - (v_nsamples - v_nhap_left)] = get_matchdata(queryidall_target,
                                                                               donorid_target,
                                                                               startpos_target,
                                                                               endpos_target,
                                                                               v_samples[ii], loo,haploid);
          }
          
#pragma omp parallel for reduction(+:lambda_sum)
          for(int ii = v_nsamples - v_nhap_left; ii < v_nsamples - v_nhap_left + v_nsamples_use; ++ii){
            vector<vector<int>> v_targetmatchdata=v_targetmatch_use[ii - (v_nsamples - v_nhap_left)];
            
            vector<int> v_startpos, v_endpos;
            for (const auto& row : v_targetmatchdata) {
              v_startpos.push_back(row[1]);
              v_endpos.push_back(row[2]);
            }
            double lambda_estimated=est_lambda_Viterbi(v_startpos,v_endpos,nsnp,gdall);
            lambda_sum += lambda_estimated;
          }
          v_nhap_left=v_nhap_left-v_nsamples_use;
        }
        
        lambda=lambda_sum/v_nsamples;
      }
      cout << "Using fixed lambda "<<lambda<<endl;
    }
    
    ofstream outputlambda(lambdafile);
    outputlambda << "The fixed lambda used for SparsePainter is "<<lambda<<".";
    outputlambda.close();
  }
  
  // begin painting
  // we only store ncores*2 samples in memory and directly output
  
  omp_set_num_threads(ncores);
  
  int nsamples_use;
  int nhap_left=nhap_use;
  
  //store data in hMat if want to compute LDA
  
  vector<int> nsnp_left(nsnp);
  vector<int> nsnp_right(nsnp);
  
  if(outputLDA || outputLDAS){
    // calculate the number of SNPs in the left and right window of each SNP
    
    for (int j = 0; j < nsnp; ++j) {
      int left_ptr = j - 1;
      int right_ptr = j + 1;
      nsnp_left[j] = 0;
      nsnp_right[j] = 0;
      // Count SNPs in the left window
      while (left_ptr >= 0 && gd[j] - gd[left_ptr] <= window) {
        nsnp_left[j]++;
        left_ptr--;
      }
      // Count SNPs in the right window
      while (right_ptr < nsnp && gd[right_ptr] - gd[j] <= window) {
        nsnp_right[j]++;
        right_ptr++;
      }
    }
  }
  
  vector<vector<double>> Dscore(nsnp - 1); // Create a vector of vectors with size nsnp - 1
  vector<vector<double>> Dprime(nsnp - 1);
  for(int i = 0; i < nsnp - 1; ++i) {
    Dscore[i].resize(nsnp_right[i], 0.0); // Resize the inner vector to nsnp_right[i] and initialize with 0.0
    Dprime[i].resize(nsnp_right[i], 0.0);
  }
  
  // the average painting for each SNP
  vector<vector<double>> aveSNPpainting(npop, vector<double>(nsnp));
  for(int j=0;j<npop;++j){
    for(int k=0;k<nsnp;++k){
      aveSNPpainting[j][k]=0;
    }
  }
  
  // the average painting for each individual
  vector<vector<double>> aveindpainting(nind, vector<double>(npop));
  for(int j=0;j<nind;++j){
    for(int k=0;k<npop;++k){
      aveindpainting[j][k]=0;
    }
  }
  
  //output the painting into probfile
  ogzstream outputFile;
  if(outputpainting && run!="chunklength"){
    outputFile.open(probfile.c_str());
    if (!outputFile) {
      cerr << "Error: unable to open file: " << probfile << endl;
      abort();
    }
    outputFile << "haplotype_name" << " ";
    //the first row is the SNP's physical position
    for (int i = 0; i < nsnp; ++i) {
      outputFile << fixed << setprecision(0) << pd[i];
      if(i != nsnp-1) outputFile << " ";
    }
    outputFile << "\n";
  }
  
  ofstream outputclFile;
  if(run!="prob"){
    outputclFile.open(chunklengthfile.c_str());
    if (!outputclFile) {
      cerr << "Error: unable to open file: " << chunklengthfile << endl;
      abort();
    }
    outputclFile << "indnames" << " ";
    for (int i = 0; i < npop; ++i) {
      outputclFile <<"pop";
      outputclFile << fixed << setprecision(0) << i;
      if(i != npop-1) outputclFile << " ";
    }
    outputclFile << "\n";
  }
  
  int looptime=0;
  
  while(nhap_left>0){
    nsamples_use = (ncores*2*LDAfactor < nhap_left) ? ncores*2*LDAfactor : nhap_left; //ensure both copies are included
    
    vector<vector<vector<double>>> painting_all(nsamples_use, 
                                                vector<vector<double>>(npop, vector<double>(nsnp)));
    
    vector<vector<double>> chunklength(nsamples_use, vector<double>(npop));
    
    // get the matches before the loop
    vector<vector<vector<int>>> targetmatch_use(nsamples_use);
    
    for (int ii = nhap_use - nhap_left; ii < nhap_use - nhap_left + nsamples_use; ++ii) {
      // leave one out if the donor file is the same as the target file
      
      targetmatch_use[ii - (nhap_use - nhap_left)] = get_matchdata(queryidall_target,
                                                                   donorid_target,
                                                                   startpos_target,
                                                                   endpos_target,ii, loo,haploid);
    }
    
    if(run=="prob"){
      cout<<"Calculating painting for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
    }else{
      if(run=="both"){
        cout<<"Calculating painting and chunk length for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
      }else{
        cout<<"Calculating chunk length for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
      }
    }
#pragma omp parallel for
    for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
      
      double lambda_use;
      
      vector<vector<int>> targetmatchdata=targetmatch_use[ii - (nhap_use - nhap_left)];
      vector<int> removeidx;
      
      if(leaveoneout){
        if(reffile==targetfile){
          int popidx=refindex[queryidx[ii]];
          for(int j=0;j<npop;++j){
            if(j!=popidx){
              int rmidx1=randomsample(refidx.findrows(j),1)[0];
              removeidx.push_back(rmidx1);
              if(!haploid){
                int rmidx2;
                if(rmidx1%2==0){
                  rmidx2=rmidx1+1;
                }else{
                  rmidx2=rmidx1-1;
                }
                removeidx.push_back(rmidx2);
              }
            }
            //removeidx contains the indices to be removed for leave-one-out
            //the same individual has already been removed
          }
          removeRowsWithValue(targetmatchdata,removeidx);
        }else{
          for(int j=0;j<npop;++j){
            int rmidx1=randomsample(refidx.findrows(j),1)[0];
            removeidx.push_back(rmidx1);
            if(!haploid){
              int rmidx2;
              if(rmidx1%2==0){
                rmidx2=rmidx1+1;
              }else{
                rmidx2=rmidx1-1;
              }
              removeidx.push_back(rmidx2);
            }
          }
          removeRowsWithValue(targetmatchdata,removeidx);
          
        }
      }
      
      hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp); 
      
      if(diff_lambda){
        vector<int> startpos, endpos;
        for (const auto& row : targetmatchdata) {
          startpos.push_back(row[1]);
          endpos.push_back(row[2]);
        }
        lambda_use=est_lambda_Viterbi(startpos,endpos,nsnp,gdall);
      }else{
        lambda_use=lambda;
      }
      
      
      vector<double> sameprob=cal_sameprob(nsnp,lambda_use,gd,nref);
      vector<double> otherprob=cal_otherprob(nref,sameprob);
      tuple<hMat, vector<double>> f=forwardProb(mat,sameprob,otherprob);
      tuple<hMat, vector<double>> b=backwardProb(mat,sameprob,otherprob);
      
      if(run!="chunklength"){
        hMat pind=indpainting(mat,gd,lambda_use,npop,refindex,f,b);
        
        vector<vector<double>> pind_dense=hMatrix2matrix(pind);
        for(int j=0;j<npop;++j){
          for(int k=0;k<nsnp;++k){
            painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
          }
        }
      }
      
      if(run!="prob"){
        vector<double> cl=chunklength_each(gd,mat,lambda_use,npop,refindex,f,b,gdall);
        for(int j=0;j<npop;++j){
          chunklength[ii-nhap_use+nhap_left][j]=cl[j];
        }
      }
    }
    
    if(run!="chunklength"){
      //compute average painting for each SNP
      if(outputaveSNPpainting||outputAAS){
        for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
          for(int j=0;j<npop;++j){
#pragma omp parallel for
            for(int k=0;k<nsnp;++k){
              aveSNPpainting[j][k]+=painting_all[ii-nhap_use+nhap_left][j][k];
            }
          }
        }
      }
      
      // compute the average painting for each individual and output.
      if(outputaveindpainting){
        if(haploid){
          for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
            for(int j=0;j<npop;++j){
              double sumpaint=0;
#pragma omp parallel for
              for(int k=0;k<nsnp;++k){
                sumpaint+=painting_all[ii-nhap_use+nhap_left][j][k];
              }
              aveindpainting[ii][j] = sumpaint/nsnp;
            }
          }
        }else{
          for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
            for (int j = 0; j < npop; ++j) {
              double sumpaint=0;
              for(int k=0;k<nsnp;++k){
                sumpaint+=painting_all[2*ii-nhap_use+nhap_left][j][k]+painting_all[2*ii-nhap_use+nhap_left+1][j][k];
              }
              aveindpainting[ii][j] = sumpaint/(2*nsnp);
            }
          }
        }
      }
      
      //compute LDA
      
      if(outputLDA || outputLDAS){
        cout<<"Calculating LDA for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        vector<int> allhaps_idx;
        for(int i=0;i<nsamples_use;++i){
          allhaps_idx.push_back(i);
        }
        vector<int> resample_idx = randomsample(allhaps_idx,nsamples_use);
#pragma omp parallel for
        for(int i=0;i<nsnp-1;++i){
          if(nsnp_right[i]!=0){
            for(int j=i+1;j<=i+nsnp_right[i];++j){
              double distance=0;
              double theo_distance=0;
              for (int nn=0; nn<nsamples_use; nn++){
                double sum_squared_diff=0;
                double sum_squared_diff_theo=0;
                for (int k=0; k<npop; k++){
                  sum_squared_diff+= pow(painting_all[nn][k][i]-painting_all[nn][k][j],2);
                  sum_squared_diff_theo+= pow(painting_all[resample_idx[nn]][k][i]-painting_all[nn][k][j],2);
                }
                distance += sqrt(sum_squared_diff/npop);
                theo_distance += sqrt(sum_squared_diff_theo/npop);
              }
              if(looptime==0){
                Dscore[i][j-i-1]=distance;
                Dprime[i][j-i-1]=theo_distance;
              }else{
                Dscore[i][j-i-1]+=distance;
                Dprime[i][j-i-1]+=theo_distance;
              }
            }
          }
        }
      }
      
      if(outputpainting){
        //output painting
        if(haploid){
          for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
            outputFile << indnames[ii] << " ";
            for (int j = 0; j < nsnp; ++j) {
              for(int k=0;k<npop;++k){
                outputFile << fixed << setprecision(2) << painting_all[ii-nhap_use+nhap_left][k][j];
                if(k!=npop-1) outputFile << ",";
              }
              if(j!=nsnp-1) outputFile << " ";
            }
            outputFile << "\n";
          }
        }else{
          for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
            outputFile << indnames[ii] << " ";
            for (int j = 0; j < nsnp; ++j) {
              for(int k=0;k<npop;++k){
                outputFile << fixed << setprecision(2) << painting_all[2*ii-nhap_use+nhap_left][k][j];
                if(k!=npop-1) outputFile << ",";
              }
              outputFile << "|";
              for(int k=0;k<npop;++k){
                outputFile << fixed << setprecision(2) << painting_all[2*ii-nhap_use+nhap_left+1][k][j];
                if(k!=npop-1) outputFile << ",";
              }
              if(j!=nsnp-1) outputFile << " ";
            }
            outputFile << "\n";
          }
        }
      }
      
      
      vector<vector<vector<double>>>().swap(painting_all);
    }
    
    if(run!="prob"){
      //output chunk length
      if(haploid){
        for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
          outputclFile << indnames[ii] << " ";
          for(int j=0;j<npop;++j){
            outputclFile << fixed << setprecision(5) << chunklength[ii-nhap_use+nhap_left][j];
            if(j!=npop-1) outputclFile << " ";
          }
          outputclFile << "\n";
        }
      }else{
        for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
          outputclFile << indnames[ii] <<"_0 ";
          for(int j=0;j<npop;++j){
            outputclFile << fixed << setprecision(3) << chunklength[2*ii-nhap_use+nhap_left][j];
            if(j!=npop-1) outputclFile << " ";
          }
          outputclFile << "\n";
          
          outputclFile << indnames[ii] <<"_1 ";
          for(int j=0;j<npop;++j){
            outputclFile << fixed << setprecision(3) << chunklength[2*ii-nhap_use+nhap_left+1][j];
            if(j!=npop-1) outputclFile << " ";
          }
          outputclFile << "\n";
        }
      }
      vector<vector<double>>().swap(chunklength);
    }
    
    nhap_left=nhap_left-nsamples_use;
    looptime++;
  }
  
  if(run!="prob"){
    outputclFile.close();
  }
  
  
  if(run!="chunklength"){
    if(outputpainting){
      outputFile.close();
    }
    
    
    // get the average painting for each SNP and output
    if(outputaveSNPpainting||outputAAS){
      for(int j=0;j<npop;++j){
#pragma omp parallel for
        for(int k=0;k<nsnp;++k){
          aveSNPpainting[j][k]=aveSNPpainting[j][k]/nhap_use;
        }
      }
      
      if(outputaveSNPpainting){
        //output the average painting for each SNP
        ofstream outputFile(aveSNPprobfile.c_str());
        
        if (outputFile) {
          outputFile << "physical_position" <<" ";
          for (int j = 0; j < npop; ++j){
            outputFile << "pop"<<j << " ";
          }
          outputFile<< "\n";
          for (int k = 0; k < nsnp; ++k) {
            outputFile << fixed << setprecision(0) << pd[k]<< " ";
            for (int j = 0; j < npop; ++j){
              outputFile << fixed << setprecision(4) <<aveSNPpainting[j][k];
              if(j != npop-1) outputFile << " ";
            }
            if(k != nsnp-1) outputFile<< "\n";
          }
          outputFile.close();
        }else {
          cerr << "Unable to open file" << aveSNPprobfile;
          abort();
        }
      }
    }
    
    //output the average painting for each individual
    if(outputaveindpainting){
      //output the average painting for each SNP
      ofstream outputFile(aveindprobfile.c_str());
      if (outputFile) {
        outputFile << "individual_name" << " ";
        //the first row is the SNP's populations
        for (int k = 0; k < npop; ++k) {
          outputFile << "pop" << k;
          if(k != npop-1) outputFile << " ";
        }
        outputFile << "\n";
        
        for (int i = 0; i < nind; ++i) {
          outputFile << indnames[i] << " ";
          for (int k = 0; k < npop; ++k) {
            outputFile << fixed << setprecision(4) <<aveindpainting[i][k];
            if(k != npop-1) outputFile << " ";
          }
          if(i != nind-1) outputFile<< "\n";
        }
        
        outputFile.close();
      } else {
        cerr << "Unable to open file" << aveindprobfile;
        abort();
      }
    }
    
    // arrange results in hMat LDA_result
    hMat LDA_result(nsnp,nsnp,0.0);
    if(outputLDA || outputLDAS){
      for(int i=0; i<nsnp; ++i){
        LDA_result.m[i].set(i,1.0);
        if(nsnp_right[i]!=0){
          for(int j=i+1;j<=i+nsnp_right[i];++j){
            LDA_result.m[i].set(j,1-Dscore[i][j-i-1]/Dprime[i][j-i-1]);
          }
        }
      }
    }
    
    
    if(outputLDA){
      //output the LDA results into LDAfile
      ogzstream outputFile(LDAfile.c_str());
      if (outputFile) {
        for (int i = 0; i < nsnp; ++i) {
          vector<int> keys = LDA_result.m[i].k;
          for (int j = 0; j < keys.size(); ++j) {
            if(i < keys[j] && LDA_result.m[i].get(keys[j])>=0.005){
              outputFile << fixed << setprecision(0) << pd[i];
              outputFile << " " << fixed << setprecision(0) << pd[keys[j]];
              outputFile << " " << fixed << setprecision(2) << LDA_result.m[i].get(keys[j]);
              outputFile << "\n";
            }
          }
        }
        
        outputFile.close();
      } else {
        cerr << "Unable to open file" << LDAfile;
        abort();
      }
    }
    
    if(outputLDAS){
      hMat Dprime(0,0,0.0);
      cout << "Begin calculating LDA score"<<endl;
      doLDAS(LDA_result,LDASfile,window,gd,pd,nsnp_left,nsnp_right,nsnp);
      cout << "Finish calculating LDA score"<<endl;
    }
    
    if(outputAAS){
      cout << "Begin calculating Ancestry Anomaly Score"<<endl;
      doAAS(pd,aveSNPpainting,AASfile);
      cout << "Finish calculating Ancestry Anomaly Score"<<endl;
    }
  }
  
}


bool ends_with(const string &value, const string &ending) {
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int main(int argc, char *argv[]){
  string run="prob";
  bool runpaint=false;
  bool chunklength=false;
  string method="Viterbi";
  bool haploid=false;
  bool leaveoneout=false;
  bool diff_lambda=false;
  double fixlambda=0;
  double indfrac=0.1;
  int minsnpEM=2000;
  double EMsnpfrac=0.1;
  int L_initial=320;
  int nmatch=10;
  int L_minmatch=20;
  int ite_time=10;
  string reffile={};
  string targetfile={};
  string mapfile={};
  string popfile={};
  string targetname={};
  string matchfile={};
  bool outputpainting=false;
  bool aveSNPpainting=false;
  bool aveindpainting=false;
  bool LDA=false;
  bool LDAS=false;
  bool AAS=false;
  bool phase=false;
  string out="SparsePainter";
  double window=4;
  int ncores=0;
  
  for (int i = 1; i < argc; i++) {
    string param = argv[i];
    if (param[0] != '-') {
      cerr << "Invalid argument format. Expected -param value or -param \n";
      return 1;
    }
    param = param.substr(1);  // Remove the -
    
    if(param=="prob" || param=="chunklength" ||
       param=="aveSNP" || param=="aveind" ||
       param=="LDA" || param=="LDAS" ||
       param=="AAS" || param=="diff_lambda" ||
       param=="haploid" || param=="loo"){
      if(i!=argc-1){
        if(argv[i+1][0]!='-'){
          cerr << "Error: No values should be given following -"<<param<<endl;
          abort();
        }
      }
    }
    
    if(param=="method" || param=="fixlambda" ||
       param=="ite_time" || param=="indfrac" ||
       param=="minsnpEM" || param=="EMsnpfrac" ||
       param=="L0" || param=="nmatch" ||
       param=="Lmin" || param=="reffile"||
       param=="targetfile" || param=="mapfile"||
       param=="popfile" || param=="targetname"||
       param=="matchfile" || param=="out"||
       param=="window" || param=="ncores"){
      if(i==argc-1){
        cerr << "Error: Values should be given following -"<<param<<endl;
        abort();
      }else if(argv[i+1][0]=='-'){
        cerr << "Error: Values should be given following -"<<param<<endl;
        abort();
      }
    }
    
    if (param == "prob") {
      runpaint=true;
      outputpainting=true;
    }else if (param == "chunklength") {
      chunklength=true;
    } else if (param == "aveSNP") {
      aveSNPpainting = true;
      runpaint=true;
    } else if (param == "aveind") {
      aveindpainting = true;
      runpaint=true;
    } else if (param == "LDA") {
      LDA = true;
      runpaint=true;
    } else if (param == "LDAS") {
      LDAS = true;
      runpaint=true;
    } else if (param == "AAS") {
      AAS = true;
      runpaint=true;
    } else if (param == "diff_lambda") {
      diff_lambda = true;
    } else if (param == "haploid") {
      haploid = true;
    } else if (param == "loo") {
      leaveoneout = true;
    }else if (param == "method") {
      method = argv[++i];
    } else if (param == "fixlambda") {
      fixlambda = stod(argv[++i]);
    } else if (param == "ite_time") {
      ite_time = stoi(argv[++i]);
    } else if (param == "indfrac") {
      indfrac = stod(argv[++i]);
    } else if (param == "minsnpEM") {
      minsnpEM = stoi(argv[++i]);
    } else if (param == "EMsnpfrac") {
      EMsnpfrac = stod(argv[++i]);
    } else if (param == "L0") {
      L_initial = stoi(argv[++i]);
    } else if (param == "nmatch") {
      nmatch = stoi(argv[++i]);
    } else if (param == "Lmin") {
      L_minmatch = stoi(argv[++i]);
    } else if (param == "reffile") {
      reffile = argv[++i];
    } else if (param == "targetfile") {
      targetfile = argv[++i];
    } else if (param == "mapfile") {
      mapfile = argv[++i];
    } else if (param == "popfile") {
      popfile = argv[++i];
    } else if (param == "namefile") {
      targetname = argv[++i];
    } else if (param == "matchfile") {
      matchfile = argv[++i];
    } else if (param == "out") {
      out = argv[++i];
    } else if (param == "window") {
      window = stod(argv[++i]);
    } else if (param == "ncores") {
      ncores = stoi(argv[++i]);
    } else {
      cerr << "Unknown parameter: " << param << ".\n";
      return 1;
    }
  }
  
  if(reffile.empty()){
    reffile="donor.vcf.gz";
    cout << "No `-reffile filename' input is found, use donor.vcf.gz as default."<<endl;
  }
  
  if(targetfile.empty()){
    targetfile="target.vcf.gz";
    cout << "No `-targetfile filename' input is found, use target.vcf.gz as default."<<endl;
  }
  
  if(mapfile.empty()){
    mapfile="map.txt";
    cout << "No `-mapfile filename' input is found, use map.txt as default."<<endl;
  }
  
  if(popfile.empty()){
    popfile="popnames.txt";
    cout << "No `-popfile filename' input is found, use popnames.txt as default."<<endl;
  }
  
  if(targetname.empty()){
    targetname="targetname.txt";
    cout << "No `-targetname filename' input is found, use targetname.txt as default."<<endl;
  }
  
  
  bool reffile_phase = ends_with(reffile, ".phase") || ends_with(reffile, ".phase.gz");
  bool targetfile_phase = ends_with(targetfile, ".phase") || ends_with(targetfile, ".phase.gz");
  bool reffile_vcf = ends_with(reffile, ".vcf") || ends_with(reffile, ".vcf.gz");
  bool targetfile_vcf = ends_with(targetfile, ".vcf") || ends_with(targetfile, ".vcf.gz");
  
  if ((reffile_phase && targetfile_phase)) {
    phase=true;
  }
  if (!((reffile_vcf && targetfile_vcf) || (reffile_phase && targetfile_phase))) {
    cerr << "The reffile and targetfile should both be vcf (including gzipped vcf) or phase (including gzipped phase) format." << "\n";
    return 1;
  }
  
  string probfile = out + "_prob.txt.gz";
  string aveSNPprobfile = out + "_aveSNPprob.txt";
  string aveindprobfile = out + "_aveindprob.txt";
  string LDAfile = out + "_LDA.txt.gz";
  string LDASfile = out + "_LDAS.txt";
  string AASfile = out + "_AAS.txt";
  string chunklengthfile = out+ "_chunklength.txt";
  string lambdafile = out+ "_fixedlambda.txt";
  
  if (!matchfile.empty()){
    targetfile=reffile;
  }
  
  if(!runpaint && !chunklength){
    cerr<<"Please specify at least one of the following command in order to run SparsePainter:"<<endl;
    cerr<<"-prob: output the local ancestry probabilities for each target sample at each SNP."<<endl;
    cerr<<"-chunklength: output the chunk length of each local ancestry for each target sample."<<endl;
    cerr<<"-aveSNP: output the average local ancestry probabilities for each SNP."<<endl;
    cerr<<"-aveind: output the average local ancestry probabilities for each target sample."<<endl;
    cerr<<"-LDA: output the LDA of each pair of SNPs."<<endl;
    cerr<<"-LDAS: output the LDAS of each SNP."<<endl;
    cerr<<"-AAS: output the AAS of each SNP."<<endl;
    return 1;
  }
  
  if(runpaint && chunklength){
    run="both";
  }
  if(!runpaint && chunklength){
    run="chunklength";
  }
  
  int ncores_temp= omp_get_num_procs();
  if(ncores==0){
    ncores = ncores_temp;
  }
  if(ncores>ncores_temp){
    ncores = ncores_temp;
    cout<<"The maximum number of cores available is "<<ncores_temp<<". SparsePainter will use "<<ncores_temp<<" cores for parallel programming only."<<endl;
  }
  
  paintall(method, diff_lambda, fixlambda, ite_time, indfrac, minsnpEM, EMsnpfrac, L_initial, nmatch, 
           L_minmatch, haploid, leaveoneout, reffile, targetfile, mapfile, popfile, targetname, matchfile, 
           outputpainting,aveSNPpainting,aveindpainting,LDA, LDAS, 
           AAS,probfile, aveSNPprobfile,aveindprobfile, chunklengthfile,
           LDAfile, LDASfile, AASfile, lambdafile, window, ncores,run,phase);
  
  return 0;
} 
