// please compile with "make"
// g++ -I./armadillo-12.6.5/include SparsePainter.cpp -o SparsePainter -lz -fopenmp -lpthread -L./armadillo-12.6.5 -larmadillo -llapack -lblas -std=c++0x -g -O3 -Wl,-rpath=./armadillo-12.6.5
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
          if(x==0 || x==1){
            panel[i*2][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << inFile << endl;
            abort();
          }
          linestr >> x;
          if(x==0 || x==1){
            panel[i*2 + 1][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << inFile << endl;
            abort();
          }
        }
      }else{
        for (int i = 0; i<M-qM; ++i){
          linestr >> x;
          if(x==0 || x==1){
            panel[i][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << inFile << endl;
            abort();
          }
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
          if(x==0 || x==1){
            panel[i*2][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << inFile << endl;
            abort();
          }
          linestr >> x;
          if(x==0 || x==1){
            panel[i*2 + 1][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << inFile << endl;
            abort();
          }
        }
        for (int i = (M-qM)/2; i < M/2; ++i){
          qlinestr >> x >> y;
          if(x==0 || x==1){
            panel[i*2][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << qinFile << endl;
            abort();
          }
          qlinestr >> x;
          if(x==0 || x==1){
            panel[i*2 + 1][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << qinFile << endl;
            abort();
          }
        }
      }else{
        for (int i = 0; i<M-qM; ++i){
          linestr >> x;
          if(x==0 || x==1){
            panel[i][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << inFile << endl;
            abort();
          }
        }
        for (int i = (M-qM); i < M; ++i){
          qlinestr >> x;
          if(x==0 || x==1){
            panel[i][j] = (bool)x;
          }else{
            cerr << "Error: Genotypes are not represented as 0 or 1 in " << qinFile << endl;
            abort();
          }
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
      break;
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
                for(int q=start;q<=end;++q){
                  nmatch[q]++;
                }
                ++ftemp;
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
            for(int q=start2;q<=end2;++q){
              nmatch[q]++;
            }
            ++f;
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
      vector<int> nmatch_output(N, 0);
      for(int q=0;q<N;++q){
        fullidx.push_back(q);
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

      //record the position of the next start position of query haplotype
      //such that we know how many matches are there for this query haplotype

      LoopResult result;
      result.donorid = local_donorid;
      result.startpos = local_startpos;
      result.endpos = local_endpos;
      result.queryid = local_startpos.size();
      allResults[idx] = result;

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

/////////////////////end of pbwt contents///////////////////////////


pair<vector<int>, vector<vector<double>>> kMeans(const mat& data, int ncluster, int max_ite) {
  int n_samples = data.n_rows;
  int n_features = data.n_cols;

  mat centroids(ncluster, n_features);
  for (int j = 0; j < ncluster; ++j) {
    centroids.row(j) = data.row(randi<uword>(distr_param(0, n_samples - 1)));
  }

  vec labels = zeros<vec>(n_samples);
  bool converged = false;
  int iterations = 0;

  while (!converged && iterations < max_ite) {
    converged = true;
    ++iterations;

    for (int i = 0; i < n_samples; ++i) {
      rowvec sample = data.row(i);
      vec distances = sum(square(centroids.each_row() - sample), 1);
      uword min_index;
      double min_distance = distances.min(min_index);

      if (min_index != labels[i]) {
        labels[i] = min_index;
        converged = false;
      }
    }

    vector<int> empty_clusters;
    for (int j = 0; j < ncluster; ++j) {
      uvec cluster_indices = find(labels == j);
      if (!cluster_indices.is_empty()) {
        centroids.row(j) = mean(data.rows(cluster_indices), 0);
      } else {
        empty_clusters.push_back(j);
      }
    }

    for (int j = empty_clusters.size() - 1; j >= 0; --j) {
      centroids.shed_row(empty_clusters[j]);
    }
    ncluster -= empty_clusters.size();
  }

  map<int, int> new_labels;
  int new_label = 0;
  for (int j = 0; j < centroids.n_rows; ++j) {
    new_labels[j] = new_label++;
  }

  vector<int> stdVec(n_samples);
  for (size_t i = 0; i < labels.n_elem; ++i) {
    stdVec[i] = new_labels[labels[i]];
  }

  vector<vector<double>> centroidsVec;
  centroidsVec.reserve(centroids.n_rows);
  for (size_t i = 0; i < centroids.n_rows; ++i) {
    centroidsVec.push_back(conv_to<vector<double>>::from(centroids.row(i)));
  }

  return make_pair(stdVec, centroidsVec);
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


vector<vector<int>> get_matchdata(vector<int>& queryidall,
                                  vector<int>& donorid,
                                  vector<int>& startpos,
                                  vector<int>& endpos,
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
  if(number>popidx.size()) cout<<"Number cannot be greater than the size of popidx. Please check the populations' indices are continuous integers start from 0, as provided by the popfile."<<endl;
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
                            const vector<double> &gd,
                            const int nref){
  //compute the sameprob
  vector<double> sameprob(nsnp);
  for(int j=0;j<nsnp-1;++j){
    sameprob[j]=exp(-lambda*(gd[j+1]-gd[j]));
    if(sameprob[j]>1-nref*2e-10){ //control the value within the limit of C++
      sameprob[j]=1-nref*2e-10;
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


pair<hMat, vector<double>> forwardProb(const hMat& mat,
                                       const vector<double>& sameprob,
                                       const vector<double>& otherprob){
  //compute normalised forward probability and store in hMat
  int nrow=mat.d1;
  int ncol=mat.d2;

  hMat forward_prob(nrow,ncol,0);
  vector<int> twj=mat.m[0].k;
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

    fprev=forward_prob.m[j-1].getall(twj);
    for(int i=0;i<twj.size();++i){
      forward_prob.m[j].set(twj[i],sameprobuse*fprev[i]+otherprobuse);
    }

    sumfp=hVecSum(forward_prob.m[j],twj);
    logmultF.push_back(log(sumfp)+logmultF[j-1]);
    hVecScale(forward_prob.m[j],1.0/sumfp);
  }
  return(make_pair(forward_prob, logmultF));
}

pair<hMat, vector<double>> backwardProb(const hMat& mat,
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

  for(int j=ncol-2;j>=0;--j){
    sameprobuse=sameprob[j];
    otherprobuse=otherprob[j];
    twj=mat.m[j+1].k;
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
  return(make_pair(backward_prob, logmultB));
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

pair<hMat,vector<int>> matchfiletohMat(const vector<vector<int>>& matchdata,
                                       const int& nref,
                                       const int& nsnp,
                                       const vector<double>& gd){
  hMat mat(nref,nsnp,0.0);
  vector<int> nmatch;
  for(int i = 0; i < matchdata.size(); ++i) {
    int val = matchdata[i][0];
    for(int j = matchdata[i][1]; j <= matchdata[i][2]; ++j) {
      mat.m[j].set(val, 1.0);
    }
  }

  //impute matches
  vector<int> nomatch;
  vector<bool> withmatch;
  for(int i=0;i<nsnp;++i){
    int matchsize=mat.m[i].k.size();
    nmatch.push_back(matchsize);
    if(matchsize==0){
      nomatch.push_back(i);
      withmatch.push_back(false);
    }else{
      withmatch.push_back(true);
    }
  }

  for(int i=0;i<nomatch.size();++i){
    double gd_this=gd[nomatch[i]];
    int left = nomatch[i] - 1;
    int right = nomatch[i] + 1;
    while (left >= 0 && !withmatch[left]) {
      left--;
    }

    while (right <=nsnp-1 && !withmatch[right]) {
      right++;
    }
    if(left<0){
      vector<int> twj=mat.m[right].k;
      mat.m[nomatch[i]].setall(twj,vector<double>(1.0,twj.size()));
    }else if(right>nsnp-1){
      vector<int> twj=mat.m[left].k;
      mat.m[nomatch[i]].setall(twj,vector<double>(1.0,twj.size()));
    }else{
      double gd_left=gd_this-gd[left];
      double gd_right=gd[right]-gd_this;
      if(gd_left>gd_right){
        vector<int> twj=mat.m[right].k;
        mat.m[nomatch[i]].setall(twj,vector<double>(1.0,twj.size()));
      }else{
        vector<int> twj=mat.m[left].k;
        mat.m[nomatch[i]].setall(twj,vector<double>(1.0,twj.size()));
      }
    }
  }

  return(make_pair(mat,nmatch));
}

pair<vector<double>,vector<double>> readmap(const string& mapfile) {
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

    if(!pd.empty() && column1<pd.back()){
      cerr<<"Error: The physical distance of SNPfile is of the wrong order!"<<endl;
      abort();
    }

    if(!gd.empty() && column2<gd.back()){
      cerr<<"Error: The genetic distance of SNPfile should be increasing!"<<endl;
      abort();
    }

    pd.push_back(column1);
    gd.push_back(column2/100);
  }

  file.close();
  return make_pair(pd,gd);
}


pair<vector<double>,vector<int>> readSNP(const string& SNPfile, vector<double>& pd) {
  ifstream file(SNPfile);
  vector<double> oppd;
  vector<int> oppd_idx;

  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << SNPfile << endl;
    abort();
  }

  string line;
  double column1;

  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> column1;
    auto it = find(pd.begin(), pd.end(), column1);
    if (it != pd.end()) {
      int index = distance(pd.begin(), it);
      oppd.push_back(column1);
      oppd_idx.push_back(index);
    } else {
      cout << "Unable to find SNP "<<column1<<", skip this SNP." << endl;
    }
  }

  file.close();
  return make_pair(oppd,oppd_idx);
}

pair<vector<int>,vector<string>> readpopfile(const string popfile) {
  ifstream file(popfile);
  vector<string> refindex_raw;

  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << popfile << endl;
    abort();
  }

  string line;
  string column1;
  string column2;

  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    lineStream >> column1;
    lineStream >> column2;
    refindex_raw.push_back(column2);
  }

  file.close();

  map<string, int> ref_map;
  vector<int> refindex;
  vector<string> refidmatch;
  int current_index = 0;
  for (const string& ref : refindex_raw) {
    if (ref_map.find(ref) == ref_map.end()) {
      ref_map[ref] = current_index;
      refidmatch.push_back(ref);
      current_index++;
    }
    refindex.push_back(ref_map[ref]);
  }

  return make_pair(refindex,refidmatch);
}


vector<string> readtargetname(const string& namefile) {
  ifstream file(namefile);
  vector<string> tgnames;

  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << namefile << endl;
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
                     const vector<double>& gd,
                     const int EM_ite){
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
  for(int t=0;t<EM_ite;++t){
    sameprob=cal_sameprob(nsnp,lambda_ite,gd,nref);
    otherprob=cal_otherprob(nref,sameprob);

    pair<hMat, vector<double>> f=forwardProb(mat,sameprob,otherprob);
    hMat forward_prob=f.first;
    vector<double> logmultF=f.second;

    pair<hMat, vector<double>> b=backwardProb(mat,sameprob,otherprob);
    hMat backward_prob=b.first;
    vector<double> logmultB=b.second;
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
                             const vector<double>& gd,
                             int L_initial,
                             const int nmatch,
                             const int L_minmatch,
                             const int ncores,
                             const double indfrac,
                             const int EM_ite,
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
  int qM=0;

  int lambda_max=static_cast<int>(nsnp/gdall);


  vector<int> allsamples;
  vector<int> popstart={0}; //the start position of different population samples

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

      pair<hMat,vector<int>> matall=matchfiletohMat(matchdata,nref-npop,nsnp,gd);
      hMat mat=matall.first;
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

        double lambda_estimated=est_lambda_EM(mat_use,gd_use,EM_ite);
        count++;
        if(lambda_estimated<lambda_max){
          lambda_est.push_back(lambda_estimated);
        }
      }else{

        double lambda_estimated=est_lambda_EM(mat,gd,EM_ite);
        count++;
        if(lambda_estimated<lambda_max){
          lambda_est.push_back(lambda_estimated);
        }
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
                 const int dp,
                 pair<hMat, vector<double>>& forwardprob,
                 pair<hMat, vector<double>>& backwardprob)
{

  int precision=pow(10,dp);

  // return individual painting
  int nsnp=mat.d2;

  hMat marginal_prob=marginalProb(forwardprob.first,backwardprob.first);
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

      popprob[popidx] += marginal_prob.m[j].get(refnumberidx);
    }

    double probsum=vec_sum(popprob);
    for(int i=0; i<npop; ++i){
      marginal_prob_pop.m[j].set(i,round(popprob[i]/probsum* precision)/precision);
    }


  }
  return(marginal_prob_pop);
}


vector<double> chunklength_each(vector<double>& gd,
                                hMat& mat,
                                const int npop,
                                const vector<int>& refindex,
                                pair<hMat, vector<double>>& forwardprob,
                                pair<hMat, vector<double>>& backwardprob,
                                const double gdall){
  //calculate chunk length for each haplotype
  int nsnp=mat.d2;
  vector<double> gl(nsnp-1);
  for(int j=0;j<nsnp-1;++j){
    gl[j]=gd[j+1]-gd[j];
  }

  hMat forward_prob=forwardprob.first;
  vector<double> logmultF=forwardprob.second;

  hMat backward_prob=backwardprob.first;
  vector<double> logmultB=backwardprob.second;

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

vector<double> chunkcount_each(hMat& mat,
                               const int npop,
                               const vector<int>& refindex,
                               pair<hMat, vector<double>>& forwardprob,
                               pair<hMat, vector<double>>& backwardprob,
                               const vector<double>& sameprob){
  //calculate chunk length for each haplotype
  int nsnp=mat.d2;

  hMat forward_prob=forwardprob.first;
  vector<double> logmultF=forwardprob.second;

  hMat backward_prob=backwardprob.first;
  vector<double> logmultB=backwardprob.second;

  vector<double> sumc(npop,0.0);

  double sF = exp(logmultF[1]+logmultB[1]-logmultF[nsnp-1]);
  vector<int> twj1=mat.m[1].k;
  vector<double> valuef1=forward_prob.m[1].getall(twj1);
  vector<double> valueb1=backward_prob.m[1].getall(twj1);

  for(int i=0;i<twj1.size();++i){
    sumc[refindex[twj1[i]]]+=sF*valuef1[i]*valueb1[i];
  }

  for(int j=0;j<nsnp-1;++j){

    double al = exp(logmultF[j+1]+logmultB[j+1]-logmultF[nsnp-1]);
    double ar = exp(logmultF[j]+logmultB[j+1]-logmultF[nsnp-1]);
    vector<int> twj=mat.m[j+1].k;
    vector<double> valuefprev=forward_prob.m[j].getall(twj);
    vector<double> valuef=forward_prob.m[j+1].getall(twj);
    vector<double> valueb=backward_prob.m[j+1].getall(twj);

    // different chunk length for different populations
    vector<double> sumcright(npop,0.0);

    for(int i=0;i<twj.size();++i){
      sumcright[refindex[twj[i]]]+=al*valuef[i]*valueb[i]-ar*valuefprev[i]*valueb[i]*sameprob[j];
    }

    for(int k=0;k<npop;++k){
      sumc[k]+=sumcright[k];
    }
  }
  return(sumc);
}

vector<vector<int>> sample_each(pair<hMat, vector<double>>& forwardprob,
                                const vector<double>& sameprob,
                                const vector<double>& otherprob,
                                const int nsample){

  hMat forward_prob=forwardprob.first;

  int nsnp=forward_prob.d2;

  vector<vector<int>> sample_state(nsample,vector<int>(nsnp));

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> dis(0.0, 1.0);
  double random_value;
  double temp_sum=0;

  vector<int> twj;
  vector<double> valuef;
  double f_use;
  double total_prob;

  for(int round=0;round<nsample;round++){
    //sample the reference haplotype for the last SNP
    temp_sum=0;
    random_value = dis(gen);
    twj=forward_prob.m[nsnp-1].k;
    valuef=forward_prob.m[nsnp-1].getall(twj);
    for(int i=0;i<twj.size();i++){
      temp_sum+=valuef[i];
      if(temp_sum>=random_value){
        sample_state[round][nsnp-1]=twj[i];
        break;
      }
    }
    for(int j=nsnp-2;j>=0;j--){
      temp_sum=0;
      f_use=forward_prob.m[j].get(sample_state[round][j+1]);
      total_prob= f_use*sameprob[j]+otherprob[j];
      random_value = dis(gen)*total_prob;
      temp_sum=f_use*sameprob[j];
      // do sampling to determine if switch happens
      if(temp_sum>=random_value){
        // no switch, the previous sampled haplotype carries on
        sample_state[round][j]=sample_state[round][j+1];
      }else{
        // switch happens (still possible to switch to the previous haplotype if there is a match)
        twj=forward_prob.m[j].k;
        valuef=forward_prob.m[j].getall(twj);
        for(int i=0;i<twj.size();i++){
          temp_sum+=valuef[i]*otherprob[j];
          if(temp_sum>=random_value){
            sample_state[round][j]=twj[i];
            break;
          }
        }
      }
    }
  }

  return sample_state;
}


void doLDAS(hMat &LDA_result,
            const string LDASfile,
            const double window,
            const vector<double>& gd,
            const vector<double>& pd,
            const vector<int>& nsnp_left,
            const vector<int>& nsnp_right,
            const int nsnp){

  // calculate LDA score
  vector<double> LDAS_score(nsnp,0);
  vector<double> LDAS_upper(nsnp,0);
  vector<double> LDAS_lower(nsnp,0);
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
    if(i-nsnp_left[i]>0 && i+nsnp_right[i]<nsnp-1){
      gdgap.push_back(window-left_distance);
      gdgap.push_back(window-right_distance);
      LDA_ave.push_back((LDA_result.m[i-nsnp_left[i]].get(i)+LDA_result.m[i-nsnp_left[i]-1].get(i))/2);
      LDA_ave.push_back((LDA_result.m[i].get(i+nsnp_right[i])+LDA_result.m[i].get(i+nsnp_right[i]+1))/2);
      LDA_upper.push_back(max(LDA_result.m[i-nsnp_left[i]].get(i),LDA_result.m[i-nsnp_left[i]-1].get(i)));
      LDA_upper.push_back(max(LDA_result.m[i].get(i+nsnp_right[i]),LDA_result.m[i].get(i+nsnp_right[i]+1)));
      LDA_lower.push_back(min(LDA_result.m[i-nsnp_left[i]].get(i),LDA_result.m[i-nsnp_left[i]-1].get(i)));
      LDA_lower.push_back(min(LDA_result.m[i].get(i+nsnp_right[i]),LDA_result.m[i].get(i+nsnp_right[i]+1)));
    }

    if(i-nsnp_left[i]==0 && i+nsnp_right[i]<nsnp-1){
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


    if(i-nsnp_left[i]>0 && i+nsnp_right[i]==nsnp-1){
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
      LDAS_score[i]+=LDA_ave[q]*gdgap[q]*100;
      LDAS_upper[i]+=LDA_upper[q]*gdgap[q]*100;
      LDAS_lower[i]+=LDA_lower[q]*gdgap[q]*100;
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



void paintall(const string method,
              bool diff_lambda,
              const double fixlambda,
              const int EM_ite,
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
              const string namefile,
              const string matchfile,
              const string SNPfile,
              const string nmatchfile,
              bool outputpainting,
              bool clength,
              bool ccount,
              bool csample,
              bool outputaveSNPpainting,
              bool outputaveindpainting,
              bool outputLDA,
              bool outputLDAS,
              bool outputAAS,
              bool outputnmatch,
              bool outputallSNP,
              bool rmrelative,
              const string probfile,
              const string aveSNPprobfile,
              const string aveindprobfile,
              const string chunklengthfile,
              const string chunkcountfile,
              const string samplefile,
              const string LDAfile,
              const string LDASfile,
              const string AASfile,
              const string lambdafile,
              const string probstore,
              const double window,
              const int dp,
              const int nsample,
              const double rmsethre,
              const double relafrac,
              const int ncluster,
              const int max_ite,
              int ncores,
              const string run,
              bool phase){

  int LDAfactor=1;
  int precision=pow(10,dp);

  if(outputLDA||outputLDAS){
    LDAfactor=24/ncores+1;
    if(LDAfactor<10) LDAfactor=10;
  }

  // read the map data to get the genetic distance in Morgans
  pair<vector<double>,vector<double>> mapinfo = readmap(mapfile);
  vector<double> gd = mapinfo.second;
  vector<double> pd = mapinfo.first;

  // read the SNPs to be output if given
  int nsnp_op;
  vector<double> pd_op;
  vector<int> SNPidx_op;
  if(!outputallSNP){
    pair<vector<double>,vector<int>> SNPinfo=readSNP(SNPfile,pd);
    pd_op = SNPinfo.first;
    SNPidx_op = SNPinfo.second;
    nsnp_op=SNPidx_op.size();
  }

  pair<vector<int>,vector<string>> refindex_read = readpopfile(popfile);

  vector<int> refindex=refindex_read.first;
  vector<string> refmatch=refindex_read.second;

  int nref_ind=refindex.size();

  vector<string> indnames = readtargetname(namefile);

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
                                       indfrac,EM_ite,minsnpEM,EMsnpfrac,haploid,reffile,phase,leaveoneout,pbwtall_target);
        }else{
          cout<<"Begin doing PBWT and finding matches for reference haplotypes"<<endl;
          vector<int> queryidx_ref;

          for(int i=0;i<nref;++i){
            queryidx_ref.push_back(i);
          }

          tuple<vector<int>,vector<int>,vector<int>,vector<int>> pbwtall_ref=do_pbwt(L_initial, gd,queryidx_ref,ncores,nref,nsnp,0,nmatch,
                                                                                     L_minmatch,reffile,reffile,haploid,phase);
          cout<<"Finish finding matches with PBWT"<<endl;
          lambda=est_lambda_EM_average(refidx,nref,nsnp,gd,L_initial,nmatch,L_minmatch,ncores,
                                       indfrac,EM_ite,minsnpEM,EMsnpfrac,haploid,reffile,phase,leaveoneout,pbwtall_ref);
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
  if(outputpainting && run!="chunk"){
    outputFile.open(probfile.c_str());
    if (!outputFile) {
      cerr << "Error: unable to open file: " << probfile << endl;
      abort();
    }
    if(outputallSNP){
      if(probstore=="constant"){
        outputFile << "#Storage Mode: Constant" <<"\n";
        outputFile << "SNPidx_start"<< " "<<"SNPidx_end"<< " ";
        for (int j = 0; j < npop; ++j){
          outputFile << refmatch[j] << " ";
        }
      }else if(probstore=="raw"){
        outputFile << "#Storage Mode: Raw" <<"\n";
        outputFile << "ind_name" << " ";
        //the first row is the SNP's physical position
        for (int i = 0; i < nsnp; ++i) {
          outputFile << static_cast<int>(pd[i]);
          if(i != nsnp-1) outputFile << " ";
        }
      }else if (probstore=="linear"){
        outputFile << "#Storage Mode: Linear" <<"\n";
        outputFile << "SNPidx"<< " ";
        for (int j = 0; j < npop; ++j){
          outputFile << refmatch[j] << " ";
        }
      }else if (probstore=="cluster"){
        outputFile << "#Storage Mode: Cluster" <<"\n";
        for (int j = 0; j < npop; ++j){
          outputFile << refmatch[j] << " ";
        }
      }
    }else{
      outputFile << "#Storage Mode: Raw" <<"\n";
      outputFile << "ind_name" << " ";
      //the first row is the SNP's physical position
      for (int i = 0; i < nsnp_op; ++i) {
        outputFile << static_cast<int>(pd_op[i]);
        if(i != nsnp_op-1) outputFile << " ";
      }
    }
    outputFile <<"\n";
    outputFile.unsetf(std::ios_base::fixed);
  }

  ogzstream outputclFile;
  ogzstream outputccFile;
  ogzstream outputcsFile;
  if(run!="prob"){
    if(clength){
      outputclFile.open(chunklengthfile.c_str());
      if (!outputclFile) {
        cerr << "Error: unable to open file: " << chunklengthfile << endl;
        abort();
      }
      outputclFile << "indnames" << " ";
    }
    if(ccount){
      outputccFile.open(chunkcountfile.c_str());
      if (!outputccFile) {
        cerr << "Error: unable to open file: " << chunkcountfile << endl;
        abort();
      }
      outputccFile << "indnames" << " ";
    }
    if(csample){
      outputcsFile.open(samplefile.c_str());
      if (!outputcsFile) {
        cerr << "Error: unable to open file: " << samplefile << endl;
        abort();
      }
    }

    for (int i = 0; i < npop; ++i) {
      if(clength){
        outputclFile << refmatch[i];
        if(i != npop-1){
          outputclFile << " ";
        }
      }
      if(ccount){
        outputccFile << refmatch[i];
        if(i != npop-1){
          outputccFile << " ";
        }
      }
    }
    if(clength) outputclFile << "\n";
    if(ccount) outputccFile << "\n";
  }

  //output nmatch file
  ogzstream outputnmatchFile;
  if(outputnmatch){
    outputnmatchFile.open(nmatchfile.c_str());
    outputnmatchFile << "haplotype_name" << " ";
    for (int k = 0; k < nsnp; ++k) {
      outputnmatchFile << fixed << setprecision(0) << pd[k];
      if(k != nsnp-1) outputnmatchFile << " ";
    }
    outputnmatchFile << "\n";
  }

  int looptime=0;

  while(nhap_left>0){
    nsamples_use = (ncores*2*LDAfactor < nhap_left) ? ncores*2*LDAfactor : nhap_left; //ensure both copies are included

    vector<vector<vector<double>>> painting_all(nsamples_use,
                                                vector<vector<double>>(npop, vector<double>(nsnp)));

    vector<vector<double>> chunklength(nsamples_use, vector<double>(npop));
    vector<vector<double>> chunkcount(nsamples_use, vector<double>(npop));
    vector<vector<vector<int>>> samplestate(nsamples_use, vector<vector<int>>(nsample, vector<int>(nsnp)));


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
        if(clength&&ccount&&csample){
          cout<<"Calculating painting, chunk length, chunk counts, and sampling reference haplotypes for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(clength && !ccount && !csample){
          cout<<"Calculating painting, chunk length and chunk counts for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(clength && !ccount && !csample){
          cout<<"Calculating painting and chunk length for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(!clength && ccount && !csample){
          cout<<"Calculating painting and chunk counts for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(!clength && !ccount && csample){
          cout<<"Calculating painting and sampling reference haplotypes for target haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(clength&&!ccount&&csample){
          cout<<"Calculating painting, chunk length, and sampling reference haplotypes for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(!clength&&ccount&&csample){
          cout<<"Calculating painting, chunk counts, and sampling reference haplotypes for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }
      }else{
        if(clength&&ccount&&csample){
          cout<<"Calculating chunk length, chunk counts, and sampling reference haplotypes for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(clength && !ccount && !csample){
          cout<<"Calculating chunk length and chunk counts for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(clength && !ccount && !csample){
          cout<<"Calculating chunk length for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(!clength && ccount && !csample){
          cout<<"Calculating chunk counts for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(!clength && !ccount && csample){
          cout<<"Sampling reference haplotypes for target haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(clength&&!ccount&&csample){
          cout<<"Calculating chunk length and sampling reference haplotypes for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }else if(!clength&&ccount&&csample){
          cout<<"Calculating chunk counts and sampling reference haplotypes for haplotypes "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
        }
      }
    }

    vector<vector<int>> nmatch_use(nsamples_use);

#pragma omp parallel for
    for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){

      double lambda_use;

      vector<vector<int>> targetmatchdata=targetmatch_use[ii - (nhap_use - nhap_left)];
      vector<int> removeidx;

      vector<int> maxind(npop,0);
      vector<double> maxPHAT(npop,0.0);
      if(leaveoneout && rmrelative){
        vector<int> nsameSNP(nref_ind, 0);
        for(int j = 0; j < targetmatchdata.size(); ++j) {
          int val = targetmatchdata[j][0];
          if(haploid){
            nsameSNP[val]+=targetmatchdata[j][2]-targetmatchdata[j][1]+1;
          }else{
            nsameSNP[val/2]+=targetmatchdata[j][2]-targetmatchdata[j][1]+1;
          }
        }
        int refpopidx;
        for(int j=0;j<nref_ind;++j){
          if(haploid){
            refpopidx=refindex[j];
          }else{
            refpopidx=refindex[2*j];
          }
          if(nsameSNP[j]>maxPHAT[refpopidx]){
            maxind[refpopidx]=j;
            maxPHAT[refpopidx]=nsameSNP[j];
          }
        }
      }


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
            if(rmrelative && maxPHAT[j]/nsnp/2>=relafrac){
              removeidx.push_back(maxind[j]*2);
              removeidx.push_back(maxind[j]*2+1);
            }else{
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
      }

      pair<hMat,vector<int>> matall=matchfiletohMat(targetmatchdata,nref,nsnp,gd);
      hMat mat=matall.first;
      nmatch_use[ii - (nhap_use - nhap_left)]=matall.second;

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
      pair<hMat, vector<double>> f=forwardProb(mat,sameprob,otherprob);
      pair<hMat, vector<double>> b=backwardProb(mat,sameprob,otherprob);

      if(run!="chunk"){
        hMat pind=indpainting(mat,gd,lambda_use,npop,refindex,dp,f,b);

        vector<vector<double>> pind_dense=hMatrix2matrix(pind);
        for(int j=0;j<npop;++j){
          for(int k=0;k<nsnp;++k){
            painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
          }
        }
      }

      if(run!="prob"){
        if(clength){
          vector<double> cl=chunklength_each(gd,mat,npop,refindex,f,b,gdall);
          for(int j=0;j<npop;++j){
            chunklength[ii-nhap_use+nhap_left][j]=cl[j];
          }
        }
        if(ccount){
          vector<double> cc=chunkcount_each(mat,npop,refindex,f,b,sameprob);
          for(int j=0;j<npop;++j){
            chunkcount[ii-nhap_use+nhap_left][j]=cc[j];
          }
        }
       if(csample){
          vector<vector<int>> cs=sample_each(f,sameprob,otherprob,nsample);
          for(int j=0;j<nsample;++j){
            for(int z=0;z<nsnp;++z){
              samplestate[ii-nhap_use+nhap_left][j][z]=cs[j][z];
            }
          }
        }

      }
    }

    if(run!="chunk"){
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
        if(outputallSNP){
          if(probstore=="constant"){
            if(haploid){
              for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
                outputFile << indnames[ii] << " "<<"\n";
                int snpidx=1;
                bool same=true;
                for (int j = 1; j < nsnp; ++j) {
                  for(int k=0;k<npop;++k){
                    if(painting_all[ii-nhap_use+nhap_left][k][j]!=painting_all[ii-nhap_use+nhap_left][k][j-1]){
                      same=false;
                      break;
                    }
                  }
                  if(!same || j==nsnp-1){
                    if(j==nsnp-1){
                      outputFile << snpidx <<" "<<j+1<<" ";
                    }else{
                      outputFile << snpidx <<" "<<j<<" ";
                    }

                    for(int k=0;k<npop;++k){
                      if(j==nsnp-1){
                        outputFile << painting_all[ii-nhap_use+nhap_left][k][j];
                      }else{
                        outputFile << painting_all[ii-nhap_use+nhap_left][k][j-1];
                      }
                      if(k!=npop-1) outputFile << " ";
                    }
                    same=true;
                    snpidx=j+1;
                    outputFile <<"\n";
                  }
                  if(j==nsnp-1) snpidx=1;
                }
              }
            }else{
              for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
                outputFile << indnames[ii] << "_0 "<<"\n";
                int snpidx=1;
                bool same=true;

                for (int j = 1; j < nsnp; ++j) {
                  for(int k=0;k<npop;++k){
                    if(painting_all[2*ii-nhap_use+nhap_left][k][j]!=painting_all[2*ii-nhap_use+nhap_left][k][j-1]){
                      same=false;
                      break;
                    }
                  }
                  if(!same || j==nsnp-1){
                    if(j==nsnp-1){
                      outputFile << snpidx <<" "<<j+1<<" ";
                    }else{
                      outputFile << snpidx <<" "<<j<<" ";
                    }

                    for(int k=0;k<npop;++k){
                      if(j==nsnp-1){
                        outputFile << painting_all[2*ii-nhap_use+nhap_left][k][j];
                      }else{
                        outputFile << painting_all[2*ii-nhap_use+nhap_left][k][j-1];
                      }
                      if(k!=npop-1) outputFile << " ";
                    }
                    same=true;
                    snpidx=j+1;
                    outputFile <<"\n";
                  }
                  if(j==nsnp-1) snpidx=1;
                }

                outputFile << indnames[ii] << "_1 "<<"\n";

                for (int j = 1; j < nsnp; ++j) {
                  for(int k=0;k<npop;++k){
                    if(painting_all[2*ii-nhap_use+nhap_left+1][k][j]!=painting_all[2*ii-nhap_use+nhap_left+1][k][j-1]){
                      same=false;
                      break;
                    }
                  }
                  if(!same || j==nsnp-1){
                    if(j==nsnp-1){
                      outputFile << snpidx <<" "<<j+1<<" ";
                    }else{
                      outputFile << snpidx <<" "<<j<<" ";
                    }

                    for(int k=0;k<npop;++k){
                      if(j==nsnp-1){
                        outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][j];
                      }else{
                        outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][j-1];
                      }
                      if(k!=npop-1) outputFile << " ";
                    }
                    same=true;
                    snpidx=j+1;
                    outputFile <<"\n";
                  }
                  if(j==nsnp-1) snpidx=1;
                }
              }
            }
          }else if(probstore=="raw"){
            if(haploid){
              for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
                outputFile << indnames[ii] << " ";
                for (int j = 0; j < nsnp; ++j) {
                  for(int k=0;k<npop;++k){
                    outputFile << painting_all[ii-nhap_use+nhap_left][k][j];
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
                    outputFile << painting_all[2*ii-nhap_use+nhap_left][k][j];
                    if(k!=npop-1) outputFile << ",";
                  }
                  outputFile << "|";
                  for(int k=0;k<npop;++k){
                    outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][j];
                    if(k!=npop-1) outputFile << ",";
                  }
                  if(j!=nsnp-1) outputFile << " ";
                }
                outputFile << "\n";
              }
            }
          }else if(probstore=="linear"){
            // piecewise linear output which is the most sparse version
            if(haploid){
              for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
                outputFile << indnames[ii]<<"\n";
                //output the first SNP's painting
                outputFile << 1 <<" ";
                for(int k=0;k<npop;++k){
                  outputFile << painting_all[ii-nhap_use+nhap_left][k][0];
                  if(k!=npop-1) outputFile << " ";
                }
                outputFile <<"\n";

                int startidx=0;
                for (int j = 1; j < nsnp; ++j) {
                  double rmse=0;
#pragma omp parallel for reduction(+:rmse)
                  for(int k=0;k<npop;++k){
                    double slope=(painting_all[ii-nhap_use+nhap_left][k][j]-painting_all[ii-nhap_use+nhap_left][k][startidx])/(j-startidx);
                    if(j-startidx>=2){
                      for(int jj=startidx+1;jj<j;++jj){
                        double se=painting_all[ii-nhap_use+nhap_left][k][startidx]+slope*(jj-startidx)-painting_all[ii-nhap_use+nhap_left][k][jj];
                        rmse+=se*se;
                      }
                    }
                  }
                  if(j-startidx>=2){
                    rmse=sqrt(rmse/npop/(j-startidx-1));
                  }

                  if(rmse>rmsethre){
                    startidx=j;
                    // output the painting of the previous SNP
                    outputFile << j <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[ii-nhap_use+nhap_left][k][j-1];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                    // output the painting of the this SNP
                    outputFile << j+1 <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[ii-nhap_use+nhap_left][k][j];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                  }else if(j==nsnp-1){
                    outputFile << j+1 <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[ii-nhap_use+nhap_left][k][j];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                  }
                }
              }
            }else{
              for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
                outputFile << indnames[ii] << "_0 "<<"\n";
                //output the first SNP's painting
                outputFile << 1 <<" ";
                for(int k=0;k<npop;++k){
                  outputFile << painting_all[2*ii-nhap_use+nhap_left][k][0];
                  if(k!=npop-1) outputFile << " ";
                }
                outputFile <<"\n";

                int startidx=0;
                for (int j = 1; j < nsnp; ++j) {
                  double rmse=0;
                  for(int k=0;k<npop;++k){
                    double slope=(painting_all[2*ii-nhap_use+nhap_left][k][j]-painting_all[2*ii-nhap_use+nhap_left][k][startidx])/(j-startidx);
                    if(j-startidx>=2){
                      for(int jj=startidx+1;jj<j;++jj){
                        double se=painting_all[2*ii-nhap_use+nhap_left][k][startidx]+slope*(jj-startidx)-painting_all[2*ii-nhap_use+nhap_left][k][jj];
                        rmse+=se*se;
                      }
                    }
                  }
                  if(j-startidx>=2){
                    rmse=sqrt(rmse/npop/(j-startidx-1));
                  }

                  if(rmse>rmsethre){
                    startidx=j;
                    // output the painting of the previous SNP
                    outputFile << j <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[2*ii-nhap_use+nhap_left][k][j-1];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                    // output the painting of the this SNP
                    outputFile << j+1 <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[2*ii-nhap_use+nhap_left][k][j];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                  }else if(j==nsnp-1){
                    outputFile << j+1 <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[2*ii-nhap_use+nhap_left][k][j];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                  }
                }


                outputFile << indnames[ii] << "_1 "<<"\n";
                //output the first SNP's painting
                outputFile << 1 <<" ";
                for(int k=0;k<npop;++k){
                  outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][0];
                  if(k!=npop-1) outputFile << " ";
                }
                outputFile <<"\n";

                startidx=0;
                for (int j = 1; j < nsnp; ++j) {
                  double rmse=0;
                  for(int k=0;k<npop;++k){
                    double slope=(painting_all[2*ii-nhap_use+nhap_left+1][k][j]-painting_all[2*ii-nhap_use+nhap_left+1][k][startidx])/(j-startidx);
                    if(j-startidx>=2){
                      for(int jj=startidx+1;jj<j;++jj){
                        double se=painting_all[2*ii-nhap_use+nhap_left+1][k][startidx]+slope*(jj-startidx)-painting_all[2*ii-nhap_use+nhap_left+1][k][jj];
                        rmse+=se*se;
                      }
                    }
                  }
                  if(j-startidx>=2){
                    rmse=sqrt(rmse/npop/(j-startidx-1));
                  }


                  if(rmse>rmsethre){
                    startidx=j;
                    // output the painting of the previous SNP
                    outputFile << j <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][j-1];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                    // output the painting of the this SNP
                    outputFile << j+1 <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][j];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                  }else if(j==nsnp-1){
                    outputFile << j+1 <<" ";
                    for(int k=0;k<npop;++k){
                      outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][j];
                      if(k!=npop-1) outputFile << " ";
                    }
                    outputFile <<"\n";
                  }
                }

              }
            }
          }else if(probstore=="cluster"){
            // cluster the SNPs with K-means clustering for each haplotype
            if(haploid){
              vector<vector<vector<double>>> aveclus(nsamples_use);
              vector<vector<int>> stdVec(nsamples_use);
              vector<int> nclus(nsamples_use);
#pragma omp parallel for
              for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
                arma::mat data(nsnp, npop);
                for (int k = 0; k < npop; ++k) {
                  for (int j = 0; j < nsnp; ++j) {
                    data(j, k) = painting_all[ii-nhap_use+nhap_left][k][j];
                  }
                }

                auto result = kMeans(data, ncluster,max_ite);
                stdVec[ii-nhap_use+nhap_left] = result.first;
                aveclus[ii-nhap_use+nhap_left] = result.second;
                nclus[ii-nhap_use+nhap_left] = aveclus[ii-nhap_use+nhap_left].size();
              }

              for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
                outputFile <<"ID "<< indnames[ii] << "\n";

                for(int cl=0; cl<nclus[ii-nhap_use+nhap_left]; ++cl){
                  for(int k=0;k<npop;++k){
                    outputFile << round(aveclus[ii-nhap_use+nhap_left][cl][k]*precision)/precision;
                    if(k!=npop-1) outputFile << " ";
                  }
                  outputFile <<"\n";
                }
                outputFile <<"SNP_cluster ";
                for(int j=0;j<nsnp;++j){
                  outputFile <<stdVec[ii-nhap_use+nhap_left][j];
                  if(j!=nsnp-1) outputFile << ",";
                }
                outputFile <<"\n";
              }
            }else{
              vector<vector<vector<double>>> aveclus(nsamples_use);
              vector<vector<int>> stdVec(nsamples_use);
              vector<int> nclus(nsamples_use);
#pragma omp parallel for
              for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
                arma::mat data(nsnp, npop);
                for (int k = 0; k < npop; ++k) {
                  for (int j = 0; j < nsnp; ++j) {
                    data(j, k) = painting_all[2*ii-nhap_use+nhap_left][k][j];
                  }
                }

                auto result = kMeans(data, ncluster,max_ite);
                stdVec[2*ii-nhap_use+nhap_left] = result.first;
                aveclus[2*ii-nhap_use+nhap_left] = result.second;
                nclus[2*ii-nhap_use+nhap_left] = aveclus[2*ii-nhap_use+nhap_left].size();

                for (int k = 0; k < npop; ++k) {
                  for (int j = 0; j < nsnp; ++j) {
                    data(j, k) = painting_all[2*ii-nhap_use+nhap_left+1][k][j];
                  }
                }

                result = kMeans(data, ncluster,max_ite);
                stdVec[2*ii-nhap_use+nhap_left+1] = result.first;
                aveclus[2*ii-nhap_use+nhap_left+1] = result.second;
                nclus[2*ii-nhap_use+nhap_left+1] = aveclus[2*ii-nhap_use+nhap_left+1].size();
              }
              for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
                outputFile <<"ID "<< indnames[ii] << "_0" <<"\n";

                for(int cl=0; cl<nclus[2*ii-nhap_use+nhap_left]; ++cl){
                  for(int k=0;k<npop;++k){
                    outputFile << round(aveclus[2*ii-nhap_use+nhap_left][cl][k]*precision)/precision;
                    if(k!=npop-1) outputFile << " ";
                  }
                  outputFile <<"\n";
                }
                outputFile <<"SNP_cluster ";
                for(int j=0;j<nsnp;++j){
                  outputFile <<stdVec[2*ii-nhap_use+nhap_left][j];
                  if(j!=nsnp-1) outputFile << ",";
                }
                outputFile <<"\n";

                outputFile <<"ID "<< indnames[ii] << "_1" <<"\n";

                for(int cl=0; cl<nclus[2*ii-nhap_use+nhap_left+1]; ++cl){
                  for(int k=0;k<npop;++k){
                    outputFile << round(aveclus[2*ii-nhap_use+nhap_left+1][cl][k]*precision)/precision;
                    if(k!=npop-1) outputFile << " ";
                  }
                  outputFile <<"\n";
                }
                outputFile <<"SNP_cluster ";
                for(int j=0;j<nsnp;++j){
                  outputFile <<stdVec[2*ii-nhap_use+nhap_left+1][j];
                  if(j!=nsnp-1) outputFile << ",";
                }
                outputFile <<"\n";
              }
            }
          }
        }else{
          // output specific SNPs' local ancestry probabilities
          if(haploid){
            for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
              outputFile << indnames[ii] << " ";
              for (int j = 0; j < nsnp_op; ++j) {
                for(int k=0;k<npop;++k){
                  outputFile << painting_all[ii-nhap_use+nhap_left][k][SNPidx_op[j]];
                  if(k!=npop-1) outputFile << ",";
                }
                if(j!=nsnp_op-1) outputFile << " ";
              }
              outputFile << "\n";
            }
          }else{
            for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
              outputFile << indnames[ii] << " ";
              for (int j = 0; j < nsnp_op; ++j) {
                for(int k=0;k<npop;++k){
                  outputFile << painting_all[2*ii-nhap_use+nhap_left][k][SNPidx_op[j]];
                  if(k!=npop-1) outputFile << ",";
                }
                outputFile << "|";
                for(int k=0;k<npop;++k){
                  outputFile << painting_all[2*ii-nhap_use+nhap_left+1][k][SNPidx_op[j]];
                  if(k!=npop-1) outputFile << ",";
                }
                if(j!=nsnp_op-1) outputFile << " ";
              }
              outputFile << "\n";
            }
          }
        }

      }


      vector<vector<vector<double>>>().swap(painting_all);
    }

    if(run!="prob"){
      //output chunk length

      if(haploid){
        //haploid output
        for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
          if(clength){
            outputclFile << indnames[ii] << " ";
            for(int j=0;j<npop;++j){
              outputclFile << fixed << setprecision(5) << chunklength[ii-nhap_use+nhap_left][j];
              if(j!=npop-1){
                outputclFile << " ";
              }
            }
            outputclFile << "\n";
          }
          if(ccount){
            outputccFile << indnames[ii] << " ";
            for(int j=0;j<npop;++j){
              outputccFile << fixed << setprecision(5) << chunkcount[ii-nhap_use+nhap_left][j];
              if(j!=npop-1){
                outputccFile << " ";
              }
            }
            outputccFile << "\n";
          }
          if(csample){
            outputcsFile << indnames[ii] << "\n";
            for(int j=0;j<nsample;++j){
              for(int z=0;z<nsnp;++z){
                outputcsFile << fixed << samplestate[ii-nhap_use+nhap_left][j][z];
                if(z!=nsnp-1){
                  outputcsFile << " ";
                }
              }
              outputcsFile << "\n";
            }
          }
        }
      }else{
        //diploid output
        for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
          if(clength){
            outputclFile << indnames[ii] <<"_0 ";
            for(int j=0;j<npop;++j){
              outputclFile << fixed << setprecision(3) << chunklength[2*ii-nhap_use+nhap_left][j];
              if(j!=npop-1){
                outputclFile << " ";
              }
            }
            outputclFile << "\n";
            outputclFile << indnames[ii] <<"_1 ";
            for(int j=0;j<npop;++j){
              outputclFile << fixed << setprecision(3) << chunklength[2*ii-nhap_use+nhap_left+1][j];
              if(j!=npop-1){
                outputclFile << " ";
              }
            }
            outputclFile << "\n";
          }
          if(ccount){
            outputccFile << indnames[ii] <<"_0 ";
            for(int j=0;j<npop;++j){
              outputccFile << fixed << setprecision(3) << chunkcount[2*ii-nhap_use+nhap_left][j];
              if(j!=npop-1){
                outputccFile << " ";
              }
            }
            outputccFile << "\n";
            outputccFile << indnames[ii] <<"_1 ";
            for(int j=0;j<npop;++j){
              outputccFile << fixed << setprecision(3) << chunkcount[2*ii-nhap_use+nhap_left+1][j];
              if(j!=npop-1){
                outputccFile << " ";
              }
            }
            outputccFile << "\n";
          }
          if(csample){
            outputcsFile << indnames[ii] <<"_0"<< "\n";
            for(int j=0;j<nsample;++j){
              for(int z=0;z<nsnp;++z){
                outputcsFile << samplestate[2*ii-nhap_use+nhap_left][j][z];
                if(z!=nsnp-1){
                  outputcsFile << " ";
                }
              }
              outputcsFile << "\n";
            }
            outputcsFile << indnames[ii] <<"_1"<< "\n";
            for(int j=0;j<nsample;++j){
              for(int z=0;z<nsnp;++z){
                outputcsFile << samplestate[2*ii-nhap_use+nhap_left+1][j][z];
                if(z!=nsnp-1){
                  outputcsFile << " ";
                }
              }
              outputcsFile << "\n";
            }
          }
        }
      }
      if(clength) vector<vector<double>>().swap(chunklength);
      if(ccount) vector<vector<double>>().swap(chunkcount);
      if(csample) vector<vector<vector<int>>>().swap(samplestate);
    }

    //output nmatch file
    if(outputnmatch){
      if(haploid){
        for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
          outputnmatchFile << indnames[ii] << " ";
          for(int j=0;j<nsnp;++j){
            outputnmatchFile << fixed << setprecision(0) << nmatch_use[ii-nhap_use+nhap_left][j];
            if(j!=nsnp-1) outputnmatchFile << " ";
          }
          outputnmatchFile << "\n";
        }
      }else{
        for(int ii=(nhap_use-nhap_left)/2; ii<(nhap_use-nhap_left+nsamples_use)/2; ++ii){
          outputnmatchFile << indnames[ii] <<"_0 ";
          for(int j=0;j<nsnp;++j){
            outputnmatchFile << fixed << setprecision(0) << nmatch_use[2*ii-nhap_use+nhap_left][j];
            if(j!=nsnp-1) outputnmatchFile << " ";
          }
          outputnmatchFile << "\n";

          outputnmatchFile << indnames[ii] <<"_1 ";
          for(int j=0;j<nsnp;++j){
            outputnmatchFile << fixed << setprecision(0) << nmatch_use[2*ii-nhap_use+nhap_left+1][j];
            if(j!=nsnp-1) outputnmatchFile << " ";
          }
          outputnmatchFile << "\n";
        }
      }
    }

    nhap_left=nhap_left-nsamples_use;
    looptime++;
  }
  if(outputnmatch){
    outputnmatchFile.close();
  }

  if(run!="prob"){
    if(clength) outputclFile.close();
    if(ccount) outputccFile.close();
    if(csample) outputcsFile.close();
  }


  if(run!="chunk"){
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
            outputFile << refmatch[j] << " ";
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
          outputFile << refmatch[k];
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
            double val=1-Dscore[i][j-i-1]/Dprime[i][j-i-1];
            if(val>0){
              LDA_result.m[i].set(j,val);
            }
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

  // Check if no commands are provided
  if (argc == 1 || string(argv[1]) == "-help" || string(argv[1]) == "-h") {
    cout << "Program: SparsePainter" << endl<< endl;
    cout << "Version: 1.3.2" << endl<< endl;
    cout << "SparsePainter reference: Yang, Y., Durbin, R., Iversen, A.K.N & Lawson, D.J. Sparse haplotype-based fine-scale local ancestry inference at scale reveals recent selection on immune responses. Nature Communications 16, 2742 (2025)." << endl<< endl;
    cout << "Contact: Yaoling Yang [yaoling.yang@bristol.ac.uk] or Daniel Lawson [dan.lawson@bristol.ac.uk]" << endl<< endl;
    cout << "Usage: ./SparsePainter [-command1 -command2 ...... -command3 parameter3 -command4 parameter4 ......]" << endl<< endl;

    cout << "Type ./SparsePainter, ./SparsePainter -h or ./SparsePainter -help to see this help file." << endl<< endl;

    cout << "Required Commands" << endl<< endl;
    cout << "SparsePainter has the following 6 required commands together with additional commands that specify the desired output." << endl<< endl;

    cout << "  -reffile [file]: Reference vcf (including gzipped vcf), or phase (including gzipped phase) file that contains the (phased non-missing) genotype data for the reference samples." << endl<< endl;

    cout << "  -targetfile [file]: Reference vcf (including gzipped vcf), or phase (including gzipped phase) file that contains the (phased non-missing) genotype data for target samples. To paint reference samples against themselves, please set [targetfile] to be the same as [reffile]. The file type of [targetfile] and [reffile] should be the same." << endl<< endl;

    cout << "  -mapfile [file]: Genetic map file that contains two columns with headers. The first column is the SNP position (in base) and the second column is the genetic distance of each SNP (in centiMorgan). The SNPs must be the same and of the same order as those in [reffile] and [targetfile]." << endl<< endl;

    cout << "  -popfile [file]: Population file of reference individuals that contains two columns without headers. The first column is the names of all the reference samples (must be in the same order as [reffile]). The second column is the population labels of the reference samples, which can be either strings or numbers." << endl<< endl;

    cout << "  -namefile [file]: Name file that contains the names of samples to be painted, following the same order as they appear in [targetfile].." << endl<< endl;

    cout << "  -out [string]: Prefix of the output file names (default=SparsePainter)." << endl<< endl;

    cout << "At least one of the below commands should also be given in order to run SparsePainter" << endl<< endl;

    cout << "  -prob: Output the local ancestry probabilities for each target sample at each SNP. The output is a gzipped text file (.txt.gz) with format specified in [probstore]." << endl<< endl;

    cout << "  -chunklength: Output the expected length (in centiMorgan) of copied chunks of each local ancestry for each target sample. The output is a gzipped text file (.txt.gz)." << endl<< endl;

    cout << "  -chunkcount: Output the expected number of copied chunks of each local ancestry for each target sample. The output is a gzipped text file (.txt.gz)." << endl<< endl;

    cout << "  -sample: Output the sampled reference haplotypes' indices for each target sample at each SNP. The output is a gzipped text file (.txt.gz), which is the same format as the .samples.out file of ChromoPainter, and is the required input file to run GLOBETROTTER and fastGLOBETROTTER."<< endl<< endl;

    cout << "  -aveSNP: Output the average local ancestry probabilities for each SNP. The output is a text file (.txt)." << endl<< endl;

    cout << "  -aveind: Output the average local ancestry probabilities for each target individual. The output is a text file (.txt)." << endl<< endl;

    cout << "  -LDA: Output the Linakage Disequilibrium of Ancestry (LDA) of each pair of SNPs. The output is a gzipped text file (.txt.gz). It might be slow: the computational time is proportional to the number of local ancestries and the density of SNPs in the chromosome." << endl<< endl;

    cout << "  -LDAS: Output the Linakage Disequilibrium of Ancestry Score (LDAS) of each SNP. The output is a text file (.txt), including the LDAS and its lower and upper bound, which can be used for quality control. It might be slow: the computational time is proportional to the number of local ancestries and the density of SNPs in the genome." << endl<< endl;

    cout << "  -AAS: Output the test statistic of Ancestry Anomaly Score (AAS) of each SNP. The output is a text file (.txt). The AAS test statistic follows chi-squared distribution with K degrees of freedom under the null, where K is the number of reference populations." << endl<< endl;

    cout << "Optional Commands" << endl<< endl;
    cout << "(a) Commands without parameters" << endl<< endl;

    cout << "  -haploid: The individuals are haploid." << endl<< endl;

    cout << "  -diff_lambda: Use different recombination scaling constants for each target sample. If this parameter is not given, the fixed lambda will be output in a text file (.txt) for future reference." << endl<< endl;

    cout << "  -loo: Paint with leave-one-out strategy: one individual is left out of each population (self from own population). If [-loo] is not specified under reference-vs-reference painting ([reffile] = [targetfile]), each individual will be automatically left out of painting. For accuracy, please do not use this command if any of the reference populations has very few (e.g. <=5) samples." << endl<< endl;

    cout << "  -rmrelative: Leave out the reference sample that is the most related to the target sample under leave-one-out mode [-loo], if they share at least [relafrac] proportion of SNPs of a continuous segment. Please do not use this command for reference-vs-reference painting." << endl<< endl;

    cout << "  -outmatch: Output the number of matches at each SNP for each target haplotype. The output file format is a gzipped text file (.txt.gz)." << endl<< endl;

    cout << "(b) Commands with parameters" << endl<< endl;

    cout << "  -ncores [integer>=0]: The number of CPU cores used for the analysis (default=0). The default ncores uses all the available CPU cores of your device." << endl<< endl;

    cout << "  -fixlambda [number>=0]: The value of the fixed recombination scaling constant (default=0). SparsePainter will estimate lambda as the average recombination scaling constant of [indfrac] target samples under the default [fixlambda] and [diff_lambda]." << endl<< endl;

    cout << "  -nmatch [integer>=1]: The number of haplotype matches of at least [Lmin] SNPs that SparsePainter searches for (default=10). Positions with more than [nmatch] matches of at least [Lmin] SNPs will retain at least the longest [nmatch] matches. A larger [nmatch] slightly improves accuracy but significantly increases the computational time." << endl<< endl;

    cout << "  -L0 [integer>0]: The initial length of matches (the number of SNPs) that SparsePainter searches for (default=320). [L0] must be bigger than [Lmin] and preferrably be a power of 2 of [Lmin] for computational efficiency." << endl<< endl;

    cout << "  -Lmin [integer>0]: The minimal length of matches that SparsePainter searches for (default=20). Positions with fewer than [nmatch] matches of at least [Lmin] SNPs will retain all the matches of at least [Lmin]. A larger [Lmin] increases both the accuracy and the computational time." << endl<< endl;

    cout << "  -method [Viterbi/EM]: The algorithm used for estimating the recombination scaling constant (default=Viterbi)." << endl<< endl;

    cout << "  -probstore [raw/constant/linear/ASCII/cluster]: Output the local ancestry probabilities in [raw], [constant], [linear], [ASCII] or [cluster] form (default=constant). For each haplotype, in [raw] form, we output the probabilities of each SNP with the SNP name being their physical positions in base; in [constant] form, we output the range of SNP index, and the painting probabilities that those SNPs share; in [linear] form, we output the range of SNP index, and the painting probabilities of the start SNP and the end SNP, while the intermediate SNPs are estimated by the simple linear regression with root mean squared error smaller than [rmsethre]; in [ASCII] form, we store the same results as [linear] form but converted probabilities with ASCII characters; in [cluster] form, we perform K-means clustering on the painting of each haplotype with [ncluster] clusters and maximum [max_ite] iterations, and output the average probabilities of each cluster and the cluster of each SNP. Storing in [constant] considerably reduces the file size while has the same accuracy compared with storing in [raw]; storing in [linear] has an even smaller file size but becomes slightly slower and loses some accuracy; storing in [cluster] has the smallest file size but with slowest speed." << endl<< endl;

    cout << "  -dp [integer>0]: The decimal places of the output of local ancestry probabilities (default=2). This also controls the size of the output file for local ancestry probabilities." << endl<< endl;

    cout << "  -nsample [integer>0]: The number of different sampled reference haplotypes for each target haplotype at each SNP (default=10) implemented by command [-sample]." << endl<< endl;

    cout << "  -rmsethre [number∈(0,1)]: The upper bound that the root mean squared error of the estimated local ancestry probabilities (default=0.01) when storing them in linear form by argument, i.e. [-probstore linear]." <<endl<< endl;

    cout << "  -relafrac [number∈(0,1)]: The proportion of total number of SNPs shared between a reference and target haplotype sample (default=0.2). The reference sample will be removed under the leave-one-out [-loo] and remove relative [-rmrelative] modes. " << endl<< endl;

    cout << "  -ncluster [integer>0]: The number of clusters (default=100) for K-means clustering under [-probstore cluster] mode." << endl<< endl;

    cout << "  -kmeans_ite [integer>0]: The number of maximum iterations (default=30) for K-means clustering under [-probstore cluster] mode." << endl<< endl;

    cout << "  -SNPfile [file]: File contains the specific physical position (in base) of the SNPs whose local ancestry probabilities are output in the raw form. If this file is not specified (default), then all the SNPs' local ancestry probabilities will be output in the form specified by [probstore]." <<endl<< endl;

    cout << "  -indfrac [number∈(0,1]]: The proportion of individuals used to estimate the recombination scaling constant (default=0.1)." << endl<< endl;

    cout << "  -minsnpEM [integer>0]: The minimum number of SNPs used for EM algorithm if [-method EM] is specified (default=2000)." << endl<< endl;

    cout << "  -EMsnpfrac [number∈(0,1]]: The proportion of SNPs used for EM algorithm if [-method EM] is specified (default=0.1). Note that if nsnp * [-EMsnpfrac] < [-minsnpEM], [-minsnpEM] SNPs will be used for EM algorithm." << endl<< endl;

    cout << "  -EM_ite [integer>0]: The iteration times for EM algorithm if [-method EM] is specified (default=10)." << endl<< endl;

    cout << "  -window [number>0]: The window for calculating LDA score (LDAS) in centiMorgan (default=4)." << endl<< endl;

    cout << "  -matchfile [file]: The file name of the set-maximal match file which is the output of pbwt -maxWithin (https://github.com/richarddurbin/pbwt). This can only be used for painting reference samples against themselves. When [-matchfile] is given, there is no need to provide [-reffile] and [-targetfile], because all the match information required for painting is contained in [-matchfile]. Using set-maximal matches is not recommended because set-maximal matches are extremely sparse and will significantly reduce the accuracy, despite saving compute time." << endl<< endl;

    return 0;
  }

  string run="prob";
  bool runpaint=false;
  bool chunk=false;
  bool clength=false;
  bool ccount=false;
  bool csample=false;
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
  int EM_ite=10;
  string reffile={};
  string targetfile={};
  string mapfile={};
  string popfile={};
  string namefile={};
  string matchfile={};
  string SNPfile={};
  bool outputpainting=false;
  bool aveSNPpainting=false;
  bool aveindpainting=false;
  bool LDA=false;
  bool LDAS=false;
  bool AAS=false;
  bool phase=false;
  bool outputallSNP=true;
  bool outputnmatch=false;
  bool rmrelative=false;
  string out="SparsePainter";
  double window=4;
  string probstore="constant";
  int dp=2;
  int nsample=10;
  double rmsethre=0.01;
  double relafrac=0.2;
  int ncluster=100;
  int max_ite=30;
  int ncores=0;

  for (int i = 1; i < argc; i++) {
    string param = argv[i];
    if (param[0] != '-') {
      cerr << "Error: Invalid argument format. Expected -param value or -param. \n";
      cerr<<"Type -h or -help to see the help file."<<endl;
      return 1;
    }
    param = param.substr(1);  // Remove the -

    if(param=="prob" || param=="chunklength" || param=="chunkcount" ||
       param=="sample" || param=="aveSNP" || param=="aveind" ||
       param=="LDA" || param=="LDAS" || param=="outmatch" ||
       param=="AAS" || param=="diff_lambda" || param=="rmrelative" ||
       param=="haploid" || param=="loo"){
      if(i!=argc-1){
        if(argv[i+1][0]!='-'){
          cerr << "Error: No values should be given following -"<<param<<"."<<endl;
          cerr<<"Type -h or -help to see the help file."<<endl;
          return 0;
        }
      }
    }

    if(param=="method" || param=="fixlambda" ||
       param=="EM_ite" || param=="indfrac" || param=="max_ite" ||
       param=="minsnpEM" || param=="EMsnpfrac" || param=="ncluster" ||
       param=="L0" || param=="nmatch" || param=="relafrac" ||
       param=="Lmin" || param=="reffile"||
       param=="targetfile" || param=="mapfile"|| param=="rmsethre"||
       param=="popfile" || param=="namefile"|| param=="SNPfile"||
       param=="matchfile" || param=="out" || param=="probstore" ||
       param=="window" || param=="ncores" || param=="dp" || param=="nsample"){
      if(i==argc-1){
        cerr << "Error: Parameters should be given following -"<<param<<"."<<endl;
        cerr<<"Type -h or -help to see the help file."<<endl;
        return 0;
      }else if(argv[i+1][0]=='-'){
        cerr << "Error: Parameters should be given following -"<<param<<"."<<endl;
        cerr<<"Type -h or -help to see the help file."<<endl;
        return 0;
      }
    }

    if (param == "prob") {
      runpaint=true;
      outputpainting=true;
    } else if (param == "chunklength") {
      clength=true;
      chunk=true;
    } else if (param == "chunkcount") {
      ccount=true;
      chunk=true;
    } else if (param == "sample") {
      csample=true;
      chunk=true;
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
    } else if (param == "outmatch") {
      outputnmatch = true;
    } else if (param == "loo") {
      leaveoneout = true;
    } else if (param == "rmrelative") {
      rmrelative = true;
    }else if (param == "method") {
      method = argv[++i];
    } else if (param == "fixlambda") {
      fixlambda = stod(argv[++i]);
    } else if (param == "EM_ite") {
      EM_ite = stoi(argv[++i]);
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
      namefile = argv[++i];
    } else if (param == "matchfile") {
      matchfile = argv[++i];
    } else if (param == "SNPfile") {
      SNPfile = argv[++i];
    } else if (param == "out") {
      out = argv[++i];
    } else if (param == "probstore") {
      probstore = argv[++i];
    } else if (param == "window") {
      window = stod(argv[++i]);
    } else if (param == "dp") {
      dp = stoi(argv[++i]);
    } else if (param == "nsample") {
      nsample = stoi(argv[++i]);
    } else if (param == "rmsethre") {
      rmsethre = stod(argv[++i]);
    } else if (param == "relafrac") {
      relafrac = stod(argv[++i]);
    } else if (param == "ncluster") {
      ncluster = stoi(argv[++i]);
    } else if (param == "max_ite") {
      max_ite = stoi(argv[++i]);
    } else if (param == "ncores") {
      ncores = stoi(argv[++i]);
    } else {
      cerr << "Error: Unknown argument: " << param << ".\n";
      cerr<<"Type -h or -help to see the help file."<<endl;
      return 1;
    }
  }

  if(reffile.empty()){
    cerr << "Warning: No `-reffile filename' input is found, please check your command." <<endl;
    abort();
  }

  if(targetfile.empty()){
    cerr << "Warning: No `-targetfile filename' input is found, please check your command." <<endl;
  }

  if(mapfile.empty()){
    cerr << "Warning: No `-mapfile filename' input is found, please check your command." <<endl;
  }

  if(popfile.empty()){
    cerr << "Warning: No `-popfile filename' input is found, please check your command." <<endl;
  }

  if(namefile.empty()){
    cerr << "Warning: No `-namefile filename' input is found, please check your command." <<endl;
  }

  if(rmrelative && !leaveoneout){
    cerr<<"Error: Argument rmrelative is invalid because it only works when -loo is specified." <<endl;
    abort();
  }


  bool reffile_phase = ends_with(reffile, ".phase") || ends_with(reffile, ".phase.gz");
  bool targetfile_phase = ends_with(targetfile, ".phase") || ends_with(targetfile, ".phase.gz");
  bool reffile_vcf = ends_with(reffile, ".vcf") || ends_with(reffile, ".vcf.gz");
  bool targetfile_vcf = ends_with(targetfile, ".vcf") || ends_with(targetfile, ".vcf.gz");

  if ((reffile_phase && targetfile_phase)) {
    phase=true;
  }
  if (!((reffile_vcf && targetfile_vcf) || (reffile_phase && targetfile_phase))) {
    cerr << "The reffile and targetfile should both be vcf (including gzipped vcf) or phase (including gzipped phase) format." << endl;
    cerr<<"Type -h or -help to see the help file."<<endl;
    return 1;
  }

  string probfile;

  probfile=out + "_prob.txt.gz";

  string aveSNPprobfile = out + "_aveSNPprob.txt";
  string aveindprobfile = out + "_aveindprob.txt";
  string LDAfile = out + "_LDA.txt.gz";
  string LDASfile = out + "_LDAS.txt";
  string AASfile = out + "_AAS.txt";
  string chunklengthfile = out+ "_chunklength.txt.gz";
  string chunkcountfile = out+ "_chunkcount.txt.gz";
  string samplefile = out+ "_samples.txt.gz";
  string lambdafile = out+ "_fixedlambda.txt";
  string nmatchfile = out + "_nmatches.txt.gz";

  if (!matchfile.empty()){
    targetfile=reffile;
  }

  if (!SNPfile.empty()){
    outputallSNP=false;
  }

  if(!runpaint && !clength && !ccount && !csample){
    cerr<<"Please give at least one of the following command in order to run SparsePainter:"<<endl;
    cerr<<"-prob: output the local ancestry probabilities for each target sample at each SNP."<<endl;
    cerr<<"-chunklength: output the expected chunk length of each local ancestry for each target sample."<<endl;
    cerr<<"-chunkcount: output the expected chunk connts of each local ancestry for each target sample."<<endl;
    cerr<<"-sample: output the sampled reference haplotypes for each target sample at each SNP."<<endl;
    cerr<<"-aveSNP: output the average local ancestry probabilities for each SNP."<<endl;
    cerr<<"-aveind: output the average local ancestry probabilities for each target sample."<<endl;
    cerr<<"-LDA: output the LDA of each pair of SNPs."<<endl;
    cerr<<"-LDAS: output the LDAS of each SNP."<<endl;
    cerr<<"-AAS: output the AAS of each SNP."<<endl;
    cerr<<"Type -h or -help to see the help file."<<endl;
    return 1;
  }

  if(runpaint && chunk){
    run="both";
  }
  if(!runpaint && chunk){
    run="chunk";
  }

  int ncores_temp= omp_get_num_procs();
  if(ncores==0){
    ncores = ncores_temp;
  }
  if(ncores>ncores_temp){
    ncores = ncores_temp;
    cout<<"Warning: The maximum number of cores available is "<<ncores_temp<<". SparsePainter will use "<<ncores_temp<<" cores for parallel programming only."<<endl;
  }

  if(probstore!="raw" && probstore!="constant" && probstore!="linear" && probstore!="cluster"){
    cerr<<"Warning: Invalid parameter given to probstore, using probstore=constant as default."<<endl;
    probstore="constant";
  }

  if(method!="Viterbi" && method!="EM"){
    cerr<<"Warning: Invalid parameter given to method, using method=Viterbi as default."<<endl;
    probstore="Viterbi";
  }

  paintall(method, diff_lambda, fixlambda, EM_ite, indfrac, minsnpEM, EMsnpfrac, L_initial, nmatch, L_minmatch, haploid,
           leaveoneout, reffile, targetfile, mapfile, popfile, namefile, matchfile, SNPfile, nmatchfile,
           outputpainting, clength, ccount, csample, aveSNPpainting, aveindpainting, LDA, LDAS, AAS, outputnmatch,outputallSNP,
           rmrelative, probfile, aveSNPprobfile, aveindprobfile, chunklengthfile, chunkcountfile, samplefile, LDAfile, LDASfile, AASfile,
           lambdafile, probstore, window/100, dp, nsample, rmsethre, relafrac, ncluster,max_ite, ncores, run, phase);
  cout<<"SparsePainter completed successfully!"<<endl;
  return 0;
}
