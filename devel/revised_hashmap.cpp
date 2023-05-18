// [[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]
// Compile with:
// module load languages/gcc/10.4.0
// g++ revised_hashmap.cpp -o test.exe -lz -fopenmp -lpthread
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
//#include <Rcpp.h>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cmath>

#include <iomanip>
#include <sstream>
#include <utility>

#include "gzstream.h"
#include "gzstream.C"

//using namespace Rcpp;
using namespace std;

//namespace hMatRcpp {

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




/////////////////////beginning of dpbwt contents///////////////////////////



struct dpbwtnode{
  dpbwtnode *below, *above, *u, *v;
  int divergence, id, originalid;
};

struct dpbwt{
  dpbwtnode bottomfirstcol;
  vector<dpbwtnode*> firstcol;
  vector<vector<bool>> panel;
  int size, count;
};

  
  void free_dpbwt_memory(dpbwt& x) {
    for (dpbwtnode* node : x.firstcol) {
      delete[] node;
    }
    
    for (auto& inner_vector : x.panel) {
      inner_vector.clear();
      inner_vector.shrink_to_fit();
    }
    
    x.panel.clear();
    x.panel.shrink_to_fit();
  }

void Readphase_donor(const string inFile,
                     dpbwt & x,
                     const int N, 
                     const int M, 
                     const int qM, 
                     bool haploid=false) {
  
  cout << "Read data and do dpbwt for reference panel" << endl;
  
  x.bottomfirstcol.divergence = 0;
  x.bottomfirstcol.id = -1;
  x.bottomfirstcol.below = x.bottomfirstcol.above = nullptr;
  dpbwtnode *prev = &x.bottomfirstcol, *current;
  for (int j = 0; j<N; ++j){
    current = new dpbwtnode;
    current->divergence = j+1;
    current->id = -1;
    current->below = current->above = nullptr;
    prev->u = prev->v = current;
    prev = current;
  }
  
  igzstream in(inFile.c_str());
  
  if (!in) {
    cerr << "Error opening file: " << inFile << endl;
  }
  
  string line;
  
  // Read and discard the first three lines
  for (int i = 0; i < 3; ++i) {
    getline(in, line);
  }
  
  // Read the remaining lines and store the binary data of each line in 'panel'
  for(int i=0; i<M-qM; ++i) {
    // i indicates which sample we are looking at
    //cout <<"dpbwt for donor" <<i<<endl;
    getline(in, line);
    vector<bool> panelsnp;
    
    for (char c : line) {
      panelsnp.push_back(c == '1');
    }
  
    int Oid = i;
    
    x.panel.push_back(vector<bool>());
    x.panel[i].resize(N);
    
    dpbwtnode *t = &(x.bottomfirstcol),
      *botk = &(x.bottomfirstcol),
      *z = new dpbwtnode[N+1];
      
      x.firstcol.push_back(&z[0]);
      
      
      z[0].id = i;
      z[0].originalid = Oid; //!!!! different from indel benchmark
      z[0].divergence = 0;
      z[0].below = t;
      z[0].above = t->above;
      t->above = &z[0];
      if (z[0].above != nullptr)
        z[0].above->below = &z[0];
      // open vcf as file handle
      for (int k = 0; k<N; ++k){
        // read snpk from vcf file handle
        dpbwtnode* temp = z[k].above;
        while (temp != nullptr && x.panel[temp->id][k] != panelsnp[k]){
          if (!panelsnp[k])
            temp->u = &z[k+1];
          else 
            temp->v = &z[k+1];
          temp = temp->above;
        }
        if (temp == nullptr && panelsnp[k]){
          temp = botk->above;
          botk->u = &z[k+1];
          while (temp != &z[k] && x.panel[temp->id][k]){
            temp->u = &z[k+1];
            temp = temp->above;
          }
        }
        if (!panelsnp[k]){
          z[k].u = &z[k+1];
          z[k].v = t->v;
        }
        else{
          z[k].u = t->u;
          z[k].v = &z[k+1];
        }
        t = (panelsnp[k])? t->v : t->u;
        z[k+1].id = i;
        z[k+1].originalid = Oid;
        z[k+1].below = t;
        z[k+1].above = t->above;
        t->above = &z[k+1];
        if (z[k+1].above != nullptr)
          z[k+1].above->below = &z[k+1];
        botk = botk->v;
        
        x.panel[i][k] = panelsnp[k];
        
      }
      
      int zdtemp = N,
        bdtemp = N;
      
      for (int k = N; k>= 0; --k){
        zdtemp = min(zdtemp, k);
        bdtemp = min(bdtemp, k);
        if (z[k].above!=nullptr)
          while (zdtemp>0 && panelsnp[zdtemp-1] == x.panel[z[k].above->id][zdtemp-1])
            zdtemp--;
        else 
          zdtemp = k;
        if (z[k].below->id!=-1)
          while (bdtemp>0 && panelsnp[bdtemp-1] == x.panel[z[k].below->id][bdtemp-1])
            bdtemp--;
        else 
          bdtemp = k;
        z[k].divergence = zdtemp;
        z[k].below->divergence = bdtemp;
      }
      
      
  }
  
  cout<<"finish dpbwt for reference panel"<<endl;
  
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


tuple<vector<int>,vector<int>,vector<int>,vector<int>> longMatchdpbwt(const int L_initial,
                                                                      dpbwt & x,
                                                                      const int minmatch,
                                                                      vector<double> &gd,
                                                                      vector<int>& queryidx,
                                                                      const int N,
                                                                      const int M,
                                                                      const int qM,
                                                                      const int L_minmatch,
                                                                      const int ncores,
                                                                      const string qinFile){

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
  
  igzstream in(qinFile.c_str());
  
  if (!in) {
    cerr << "Error opening file: " << qinFile << endl;
    tuple<vector<int>,vector<int>,vector<int>,vector<int>> a;
    return a;
  }
  
  string line;
  
  // Read and discard the first three lines
  for (int i = 0; i < 3; ++i) {
    getline(in, line);
  }
  
  // Read ncores lines of target data
  // paneltarget has ncore rows and N columns
  
  int nind=queryidx.size();
  
  int nind_left=nind;
  
  omp_set_num_threads(ncores);
  
  while(nind_left>0){
    int ncores_use = (ncores < nind_left) ? ncores : nind_left;
    vector<vector<bool>> panelsnp(ncores_use);
    for(int i=nind-nind_left; i<nind-nind_left+ncores_use; ++i) {
      // i indicates which sample we are looking at
      getline(in, line);
      for (char c : line) {
        panelsnp[i-nind+nind_left].push_back(c == '1');
      }
    }
    
    #pragma omp parallel for
    
    cout<<"Building dpbwt for target samples "<<nind-nind_left<<"-"<<nind-nind_left+ncores_use-1<<endl;
    
    for (int idx=nind-nind_left; idx<nind-nind_left+ncores_use; ++idx) {
      
      int *dZ;
      dZ = new int[M];
      for (int i = 0; i<M; i++){
        dZ[i] = 0;
      }
      
      int i = queryidx[idx];
      
      //cout<<i<<endl;
      int L=L_initial;
      int prevL=L;
      
      int Oid = i;
      
      dpbwtnode *t = &(x.bottomfirstcol),
        *botk = &(x.bottomfirstcol),
        *z = new dpbwtnode[N+1];
        
        z[0].id = i;
        z[0].originalid = Oid; //!!!! different from indel benchmark
        z[0].divergence = 0;
        z[0].below = t;
        z[0].above = t->above;
        
        for (int k = 0; k<N; ++k){
          
          t = (panelsnp[Oid-nind+nind_left][k])? t->v : t->u;
          z[k+1].id = i;
          z[k+1].originalid = Oid;
          z[k+1].below = t;
          z[k+1].above = t->above;
        }
        
        int zdtemp = N,
          bdtemp = N;
        int zd[N+1], bd[N+1];
        for (int k = N; k>= 0; --k){
          zdtemp = min(zdtemp, k);
          bdtemp = min(bdtemp, k);
          if (z[k].above!=nullptr)
            while(zdtemp>0 && panelsnp[Oid-nind+nind_left][zdtemp-1] == x.panel[z[k].above->id][zdtemp-1])
              zdtemp--;
          else
            zdtemp = k;
          if (z[k].below->id!=-1)
            while (bdtemp>0 && panelsnp[Oid-nind+nind_left][bdtemp-1] == x.panel[z[k].below->id][bdtemp-1])
              bdtemp--;
          else 
            bdtemp = k;
          
          zd[k] = zdtemp;
          bd[k] = bdtemp;
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
          
          
          dpbwtnode *f, *g, *ftemp, *gtemp;
          f = g = z[0].below;
          for(int k = 0; k<N; ++k){
            ftemp = (panelsnp[Oid-nind+nind_left][k])? f->u : f->v;
            gtemp = (panelsnp[Oid-nind+nind_left][k])? g->u : g->v;
            f = (panelsnp[Oid-nind+nind_left][k])? f->v : f->u;
            g = (panelsnp[Oid-nind+nind_left][k])? g->v : g->u;
            while (ftemp != gtemp){
              //matchOut << ftemp->originalid << " = q" << i-M+qM << " at [" << dZ[ftemp->id] << ", " << k << ")\n";
              int end=k-1;
              if(times==0){
                int start=dZ[ftemp->id];
                donoridtemp.push_back(ftemp->originalid);
                startpostemp.push_back(start);
                endpostemp.push_back(end);
                ftemp = ftemp->below;
                for(int q=start;q<=end;++q){
                  nmatch[q]++;
                }
              }else{
                if(addmatch[end]){
                  int start=dZ[ftemp->id];
                  //add new matches with new L
                  if(end-start+1<prevL){
                    donoridtemp.push_back(ftemp->originalid);
                    startpostemp.push_back(start);
                    endpostemp.push_back(end);
                    ftemp = ftemp->below;
                    for(int q=start;q<=end;++q){
                      nmatch[q]++;
                    }
                  }else{
                    ftemp = ftemp->below;
                  }
                }else{
                  ftemp = ftemp->below;
                }
              }
            }
            
            if (f==g){
              if (k+1-zd[k+1] == L){
                f = f->above;
                //store divergence
                dZ[f->id] = k+1-L;
              }
              if (k+1-bd[k+1] == L){
                //store divergence
                dZ[g->id] = k+1-L;
                g = g->below;
              }
            }
            if (f!=g) {
              while (f->divergence <= k+1 - L){
                f = f->above;
                //store divergence
                dZ[f->id] = k+1-L;
              }
              while (g->divergence <= k+1 - L){
                //store divergence
                dZ[g->id] = k+1-L;
                g = g->below;
              }
            }
          }
          
          while (f != g){
            //output match
            //matchOut << f->originalid << " = q" << i-M+qM << " at [" << dZ[f->id] << ", " << N << ")\n";
            int end2=N-1;
            if(times==0){
              int start2=dZ[f->id];
              donoridtemp.push_back(f->originalid);
              startpostemp.push_back(start2);
              endpostemp.push_back(end2);
              f = f->below;
              for(int q=start2;q<=end2;++q){
                nmatch[q]++;
              }
            }else{
              if(addmatch[end2]){
                int start2=dZ[f->id];
                //add new matches with new L
                if(end2-start2+1<prevL){
                  donoridtemp.push_back(f->originalid);
                  startpostemp.push_back(start2);
                  endpostemp.push_back(end2);
                  f = f->below;
                  for(int q=start2;q<=end2;++q){
                    nmatch[q]++;
                  }
                }else{
                  f = f->below;
                }
              }else{
                f = f->below;
              }
            }
            
          }
          
          if(L<=L_minmatch){
            // we don't need to check whether we have all positions with minmatch matches
            // if L is already the minimum value allowed
            //if(times==0){
              // find which SNPs don't have minmatch matches
              //for(int j=0;j<N;++j){
               //if (nmatch[j] ==0) {
                  //nomatchsnp.push_back(j);
                //}
              //}
            //}else{
              //cout<<"here"<<endl;
              //vector<int> nomatchsnptemp=nomatchsnp;
              //nomatchsnp.clear();
              // find which SNPs still don't have any matches
             //for(int j=0;j<nomatchsnptemp.size();++j){
              //  if (nmatch[nomatchsnptemp[j]] ==0) {
              //    nomatchsnp.push_back(nomatchsnptemp[j]);
              //  }
              //}
            //}
            //cout<<L<<" and "<<nomatchsnp.size()<<endl;
            //if(nomatchsnp.size()==0){
              // stop when all SNPs have at least one match
              allsnpmatches = true;
            //}else{
              // update L
             // prevL=L;
              //L=(prevL+1)/2; // this is equal to ceil(L/2) when L is double type
              //if(L==1) allsnpmatches = true;
              //times++;
            //}
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
              L=(prevL+1)/2; // this is equal to ceil(L/2) when L is double type
              if(L<L_minmatch) L=L_minmatch;
              times++;
            }
          }
          
        }
        
        delete [] z;
        delete[] dZ;
        
        
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
                auto it = std::remove(fullidx.begin(), fullidx.end(), q);
                fullidx.erase(it, fullidx.end());
              }
            }
          }
        }
        
        
        
       // vector<int> fullidx; // record which SNP has only minmatch matches
       // for(int q=0;q<N;++q){
       //   if(nmatch[q]<=minmatch) fullidx.push_back(q);
       // }
       // for(int mi=0;mi<length_order.size();++mi){
          
        //  int starttemp=startpostemp[length_order[mi]];
       //   int endtemp=endpostemp[length_order[mi]];
          
       //   if(!containsIndex(fullidx,starttemp,endtemp)){
       //     for(int q=starttemp;q<=endtemp;++q){
       //       nmatch[q]--;
        //      if(nmatch[q]==minmatch) fullidx.push_back(q);
       //     }
       //   }else{
       //     local_startpos.push_back(starttemp);
       //     local_endpos.push_back(endtemp);
       //     local_donorid.push_back(donoridtemp[length_order[mi]]);
       //   }
       // }
        
        LoopResult result;
        result.donorid = local_donorid;
        result.startpos = local_startpos;
        result.endpos = local_endpos;
        result.queryid = local_startpos.size();
        allResults[idx] = result;
        
        //record the position of the next start position of query haplotype
        //such that we know how many matches are there for this query haplotype
        //queryidall.push_back(startpos.size());
    }
    
    nind_left=nind_left-ncores_use;
  }
  
  in.close();
  
  
  for (const auto& result : allResults) {
    donorid.insert(donorid.end(), result.donorid.begin(), result.donorid.end());
    startpos.insert(startpos.end(), result.startpos.begin(), result.startpos.end());
    endpos.insert(endpos.end(), result.endpos.begin(), result.endpos.end());
    queryidall.push_back(queryidall.back() + result.queryid);
  }
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> results(queryidall,donorid,startpos,endpos);
  return(results);
}


tuple<vector<int>,vector<int>,vector<int>,vector<int>> do_dpbwt(int& L_initial, 
                                                                vector<double> gd,
                                                                vector<int>& queryidx,
                                                                int ncores,
                                                                const int M=0,
                                                                const int N=0,
                                                                const int qM=0,
                                                                int minmatch=100,
                                                                int L_minmatch=20,
                                                                bool haploid=false,
                                                                const string donorfile="donor.phase.gz",
                                                                const string targetfile="target.phase.gz"){
  
  dpbwt x;
  
  Readphase_donor(donorfile,x,N,M,qM,haploid);
  
  while(L_initial>N){
    L_initial=ceil(L_initial/2);
    cout<<"Initial L cannot be greater than N, reducing L to "<<L_initial<<endl;
  }
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> matchresults=longMatchdpbwt(L_initial,x,minmatch,
                                                                                     gd,queryidx,N,M,qM,
                                                                                     L_minmatch,ncores,targetfile);
  
  free_dpbwt_memory(x);
    
  return(matchresults);
  
}


vector<vector<int>> get_matchdata(vector<int> queryidall,
                                  vector<int> donorid,
                                  vector<int> startpos,
                                  vector<int> endpos,
                                  int queryid,
                                  bool loo=false){
  int querystart=queryidall[queryid];
  int nextquerystart=queryidall[queryid+1];
  int nrow_match=nextquerystart-querystart;
  vector<vector<int>> matchinfo(nrow_match,vector<int>(3));
  // leave out the same haplotype
  if(loo){
    for(int p=querystart;p<nextquerystart;++p){
      if(donorid[p]!=queryid){
        matchinfo[p-querystart][0]=donorid[p];
        matchinfo[p-querystart][1]=startpos[p];
        matchinfo[p-querystart][2]=endpos[p];
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

/////////////////////end of dpbwt contents///////////////////////////



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

vector<int> intesect_vec(const vector<int>& vec1, 
                         const vector<int>& vec2) {
  vector<int> result;
  // Sort the vectors to use set_intersection algorithm
  vector<int> s1(vec1);
  vector<int> s2(vec2);
  sort(s1.begin(), s1.end());
  sort(s2.begin(), s2.end());
  // Use set_intersection algorithm to find intersect elements
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), back_inserter(result));
  return result;
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


hMat matchmat2hMatrix(const vector<vector<bool>>& matchmat,
                      const double default_val=0.0){
  //Convert a logical matrix into a hash matrix
  int nrow=matchmat.size();
  int ncol=matchmat[0].size();
  hMat mat(nrow,ncol,default_val);
  
  for(int j=0;j<ncol;++j){
    for(int i=0;i<nrow;++i){
      if(matchmat[i][j]){
        mat.m[j].set(i,1.0);
      }
    }
  }
  return(mat);
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
                            const double rho, 
                            vector<double>& gd){
  //compute the sameprob
  vector<double> sameprob(nsnp);
  for(int j=0;j<nsnp-1;++j){
    sameprob[j]=exp(-rho*(gd[j+1]-gd[j]));
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
        forward_prob.m[j].set(twj[i],(sameprobuse*fprev[i]+otherprobuse)*mu);
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
        Bjp1[i]=Bjp1[i]*mu;
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
  //compute normalised backward probability and store in hMat
  int nrow=f.d1;
  int ncol=f.d2;
  hMat marginal_prob(nrow,ncol,1.0/nrow);
  for(int j=0;j<ncol;++j){
    marginal_prob.m[j]=hVecCirc(f.m[j],b.m[j]);
    hVecScale(marginal_prob.m[j],1.0/hVecSum(marginal_prob.m[j]));
  }
  
  return(marginal_prob);
}

hMat forwardBackward(const hMat& mat,
                     const vector<double>& sameprob,
                     const vector<double>& otherprob){
  //compute marginal probability
  tuple<hMat, vector<double>> f=forwardProb(mat,sameprob,otherprob);
  hMat forward_prob=get<0>(f);
  tuple<hMat, vector<double>> b=backwardProb(mat,sameprob,otherprob);
  hMat backward_prob=get<0>(b);
  hMat marginal_prob=marginalProb(forward_prob,backward_prob);
  return(marginal_prob);
}




void removeRowsWithValue(vector<vector<int>>& data, 
                         const vector<int>& values) {
  //remove the rows of data which first row is included in values
  auto it = remove_if(data.begin(), data.end(), [&](const vector<int>& row) {
    return find(values.begin(), values.end(), row.front()) != values.end();
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
    tuple<vector<double>,vector<double>> output(pd,gd);
    return output;
  }
  
  string line;
  double column1,column2;
  
  // Read and discard the header line
  getline(file, line);
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    
    // Read the first column and store it in 'column1'
    lineStream >> column1;
    
    // Read the second column and store it in 'column2'
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
    tuple<vector<string>,vector<int>> output(indnames,refindex);
    return output;
  }
  
  string line;
  string column1;
  int column2;
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    
    // Read the first column and store it in 'column1'
    lineStream >> column1;
    
    // Read the second column and store it in 'column2'
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
    return tgnames;
  }
  
  string line;
  string column1;
  
  // Read the data lines
  while (getline(file, line)) {
    istringstream lineStream(line);
    
    // Read the first column and store it in 'column1'
    lineStream >> column1;
    
    tgnames.push_back(column1);
  }
  
  file.close();
  return tgnames;
}



double est_rho_EM(hMat& mat, 
                  vector<double>& gd,
                  const int ite_time=10){
  //estimate \rho using EM algorithm
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob;
  vector<double> otherprob; 
  double rho_ite=400000/nref;
  vector<double> gl(nsnp-1);
  vector<double> rho_each((nsnp-1));
  for(int j=0;j<nsnp-1;++j){
    gl[j]=gd[j+1]-gd[j];
  }
  double totalgd=gd[nsnp-1]-gd[0];
  for(int t=0;t<ite_time;++t){
    sameprob=cal_sameprob(nsnp,rho_ite,gd);
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
      rho_each[j]=u[j]*rho_ite*gl[j]/(1-sameprob[j]);
    }
    rho_ite=vec_sum(rho_each)/totalgd;
  }
  if(isnan(rho_ite)) rho_ite=1/totalgd;
  return(rho_ite);
}

// [[Rcpp::export]]
double est_rho_Viterbi(const vector<int>& startpos, 
                       const vector<int>& endpos,
                       const int nsnp, 
                       const double gdtotal) {
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
      //maxend = maxEndMap[j-1];
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
    cout<<"There is a perfect match but our model assumes at least one match"<<endl;
    nrec=1;
  } 
  double rho_est = nrec / static_cast<double>(gdtotal);
  return rho_est;
}



double est_rho_average(const hAnc& refidx, 
                       const int nref, 
                       const int nsnp,
                       vector<double>& gd,
                       int L_initial,
                       int minmatch,
                       int L_minmatch,
                       int ncores,
                       const double indfrac=0.1,
                       const int ite_time=10, 
                       const string method="Viterbi",
                       const int minsnpEM=10000, 
                       const double EMsnpfrac=0.1,
                       bool haploid=false,
                       const string donorfile="donor.vcf",
                       tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_ref=
                         make_tuple(vector<int>(), vector<int>(), vector<int>(), vector<int>()))
  {
  // estimate \rho from the reference panel
  int npop=refidx.pos.size();
  vector<double> rho_est;
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
  
  if (get<0>(dpbwtall_ref).empty() && 
      get<1>(dpbwtall_ref).empty() && 
      get<2>(dpbwtall_ref).empty() && 
      get<3>(dpbwtall_ref).empty()){
    
    vector<int> queryidx;
    
    for(int i=0;i<nref;++i){
      queryidx.push_back(i);
    }
    
    
    cout<<"Do dPBWT for donor haplotypes"<<endl;
    tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_ref=do_dpbwt(L_initial, gd,queryidx,
                                                                                 ncores,nref,nsnp,0,
                                                                                 minmatch,L_minmatch,
                                                                                 haploid,donorfile,donorfile);
    
    cout<<"dPBWT works successfully"<<endl;
  }
  
  vector<int> queryidall=get<0>(dpbwtall_ref);
  vector<int> donorid_ref=get<1>(dpbwtall_ref);
  vector<int> startpos_ref=get<2>(dpbwtall_ref);
  vector<int> endpos_ref=get<3>(dpbwtall_ref);
  
  
  for(int i=0;i<npop;++i){
    vector<int> samples;
    for(int j=popstart[i];j<popstart[i+1];++j){
      samples.push_back(allsamples[j]);
    }
    
    for(int k=0;k<samples.size();++k){
      //leave-one-out
      vector<vector<int>> matchdata=get_matchdata(queryidall,donorid_ref,startpos_ref,endpos_ref,samples[k]);
      vector<int> removeidx;
      for(int j=0;j<npop;++j){
        if(j==i){
          removeidx.push_back(samples[k]);
        }else{
          removeidx.push_back(randomsample(refidx.findrows(j),1)[0]);
        }
        //removeidx contains the indices to be removed for leave-one-out
      }
      removeRowsWithValue(matchdata,removeidx);
      
      
      if(method=="Viterbi"){
        vector<int> startpos, endpos;
        for (const auto& row : matchdata) {
          startpos.push_back(row[1]);
          endpos.push_back(row[2]);
        }
        double rho_estimated=est_rho_Viterbi(startpos,endpos,nsnp,gdall);
        count=count+1;
        //cout<<"Estimated rho for sample "<<count<<" is "<<rho_estimated<<endl;
        rho_est.push_back(rho_estimated);
      }else{
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
          
          double rho_estimated=est_rho_EM(mat_use,gd_use,ite_time);
          count=count+1;
          //cout<<"Estimated rho for sample "<<count<<" is "<<rho_estimated<<endl;
          rho_est.push_back(rho_estimated);
        }else{
          
          double rho_estimated=est_rho_EM(mat,gd,ite_time);
          count=count+1;
          //cout<<"Estimated rho for sample "<<count<<" is "<<rho_estimated<<endl;
          rho_est.push_back(rho_estimated);
        }
      }
    }
  }
  double rho_ave=vec_sum(rho_est)/rho_est.size();
  //cout<<"Average estimated rho is "<<rho_ave<<endl;
  return(rho_ave);
}


hMat indpainting(const hMat& mat,
                 vector<double>& gd, 
                 const double rho,
                 const int npop, 
                 const vector<int>& refindex){
  // return individual painting
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob=cal_sameprob(nsnp,rho,gd);
  vector<double> otherprob=cal_otherprob(nref,sameprob);
  hMat marginal_prob=forwardBackward(mat,sameprob,otherprob);
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
      if(popprob[i]>=0.01){
        marginal_prob_pop.m[j].set(i,popprob[i]);
      }
    }
  }
  return(marginal_prob_pop);
}


vector<double> chunklength_each(vector<double>& gd, 
                                hMat& mat, 
                                const double rho, 
                                const int npop,
                                const vector<int>& refindex){
  //calculate chunk length for each reference sample
  //chunk length is the elements of coancestry matrix
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob;
  vector<double> otherprob; 
  vector<double> gl(nsnp-1);
  for(int j=0;j<nsnp-1;++j){
    gl[j]=gd[j+1]-gd[j];
  }
  sameprob=cal_sameprob(nsnp,rho,gd);
  otherprob=cal_otherprob(nref,sameprob);
  
  tuple<hMat, vector<double>> f=forwardProb(mat,sameprob,otherprob);
  hMat forward_prob=get<0>(f);
  vector<double> logmultF=get<1>(f);
  
  tuple<hMat, vector<double>> b=backwardProb(mat,sameprob,otherprob);
  hMat backward_prob=get<0>(b);
  vector<double> logmultB=get<1>(b);
  
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
  return(suml);
}


// [[Rcpp::export]]
vector<vector<double>> chunklengthall(const string method="Viterbi", 
                                      const int ite_time=10,
                                      const double indfrac=0.1,
                                      const int minsnpEM=10000, 
                                      const double EMsnpfrac=0.1,
                                      int L_initial=200,
                                      double minmatchfrac=0.002,
                                      int L_minmatch=25,
                                      bool haploid=false,
                                      const string donorfile="donor.phase.gz",
                                      const string mapfile="map.txt",
                                      const string popfile="popnames.txt",
                                      const string chunklengthfile="chunklength.txt",
                                      int ncores=0){
  
  //detect cores
  if(ncores==0){
    ncores = omp_get_num_procs();
  }
  
  // read the map data to get the genetic distance in Morgans
  tuple<vector<double>,vector<double>> mapinfo = readmap(mapfile);
  vector<double> gd = get<1>(mapinfo);
  
  tuple<vector<string>,vector<int>> popinfo = readpopfile(popfile);
  vector<string> indnames = get<0>(popinfo);
  vector<int> refindex = get<1>(popinfo);
  
  if(!haploid){
    vector<int> refindex_new;
    for(int i=0;i<refindex.size();++i){
      refindex_new.push_back(refindex[i]);
      refindex_new.push_back(refindex[i]);
    }
    refindex=refindex_new;
  }
  
  //compute coancestry matrix for all reference individuals
  const int nsnp=gd.size();
  const int nref=refindex.size();
  hAnc refidx(refindex);
  const int npop=refidx.pos.size();
  double rho;
  int minmatch=static_cast<int>(ceil(nref*minmatchfrac));
  
  
  vector<int> queryidx;
  

  for(int i=0;i<nref;++i){
    queryidx.push_back(i);
  }
  
  cout<<"Do dPBWT on the reference"<<endl;
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_ref=do_dpbwt(L_initial, gd,queryidx,
                                                                               ncores,nref,nsnp,0,
                                                                               minmatch,L_minmatch,
                                                                               haploid,donorfile,donorfile);
  vector<int> queryidall_ref=get<0>(dpbwtall_ref);
  vector<int> donorid_ref=get<1>(dpbwtall_ref);
  vector<int> startpos_ref=get<2>(dpbwtall_ref);
  vector<int> endpos_ref=get<3>(dpbwtall_ref);
  cout<<"dPBWT works successfully"<<endl;
  
  cout<<"Begin estimating fixed rho"<<endl;
  
  rho=est_rho_average(refidx,nref,nsnp,gd,L_initial,minmatch,L_minmatch,ncores,indfrac,
                      ite_time,method,minsnpEM,EMsnpfrac,haploid,donorfile,dpbwtall_ref);
  
  int nrefpaint=queryidx.size();
  
  vector<vector<double>> chunklength(nrefpaint, vector<double>(npop));
  
  omp_set_num_threads(ncores);
  
  int nhap_left=nrefpaint;
  
  while(nhap_left>0){
    int nsamples_use = (ncores < nhap_left) ? ncores : nhap_left;
    // get the matches before the loop
    vector<vector<vector<int>>> matchdata_use(nsamples_use);
    
    for (int ii = nrefpaint - nhap_left; ii < nrefpaint - nhap_left + nsamples_use; ++ii) {
      // leave one out if the donor file is the same as the target file
      vector<vector<int>> match_data = get_matchdata(queryidall_ref,
                                                     donorid_ref,
                                                     startpos_ref,
                                                     endpos_ref,
                                                     ii);
      
      matchdata_use[ii - (nrefpaint - nhap_left)] = match_data;
    }
    
    #pragma omp parallel for
    cout<<"Calculating chunk length for donor samples "<<nrefpaint-nhap_left<<"-"<<nrefpaint-nhap_left+nsamples_use-1<<endl;
    for(int i=nrefpaint-nhap_left;i<nrefpaint-nhap_left+nsamples_use;++i){
      //leave-one-out
      vector<vector<int>> matchdata=matchdata_use[i-nrefpaint+nhap_left];
      vector<int> removeidx;
      int popidx=refindex[queryidx[i]];
      for(int j=0;j<npop;++j){
        if(j==popidx){
          removeidx.push_back(queryidx[i]);
        }else{
          removeidx.push_back(randomsample(refidx.findrows(j),1)[0]);
        }
        //removeidx contains the indices to be removed for leave-one-out
      }
      removeRowsWithValue(matchdata,removeidx);
      hMat mat=matchfiletohMat(matchdata,nref-npop,nsnp);
      vector<double> cl=chunklength_each(gd,mat,rho,npop,refindex);
      for(int j=0;j<npop;++j){
        chunklength[i][j]=cl[j];
      }
    }
    nhap_left=nhap_left-nsamples_use;
    
  }
  
  
  

  
  
  ofstream outputFile(chunklengthfile.c_str());
  if (outputFile) {
    outputFile << "indnames" << " ";
    //the first row is the SNP's physical position
    for (int i = 0; i < npop; ++i) {
      outputFile <<"pop";
      outputFile << fixed << setprecision(0) << i;
      if(i != npop-1) outputFile << " ";
    }
    outputFile << "\n";
    
    if(haploid){
      for(int ii=0;ii<nrefpaint;++ii){
        outputFile << indnames[ii] << " ";
        for(int j=0;j<npop;++j){
          outputFile << fixed << setprecision(5) << chunklength[ii][j];
          if(j!=npop-1) outputFile << " ";
        }
        outputFile << "\n";
      }
    }else{
      for(int ii=0;ii<nrefpaint/2;++ii){
        outputFile << indnames[ii] <<"_0 ";
        for(int j=0;j<npop;++j){
          outputFile << fixed << setprecision(5) << chunklength[2*ii][j];
          if(j!=npop-1) outputFile << " ";
        }
        outputFile << "\n";
        
        outputFile << indnames[ii] <<"_1 ";
        for(int j=0;j<npop;++j){
          outputFile << fixed << setprecision(5) << chunklength[2*ii+1][j];
          if(j!=npop-1) outputFile << " ";
        }
        outputFile << "\n";
      }
    }
    
    outputFile.close();
  } else {
    cerr << "Unable to open file" << chunklengthfile;
  }
  
  
  
  return(chunklength);
}


void LDAS(hMat &LDA_result,
          const string LDASfile,
          const double window,
          const vector<double> gd,
          const vector<double> pd,
          const vector<int> nsnp_left,
          const vector<int> nsnp_right,
          const int nsnp){
    
    // calculate LDA score
    vector<double> LDAS(nsnp);
    vector<double> LDAS_upper(nsnp);
    vector<double> LDAS_lower(nsnp);
#pragma omp parallel for
    for(int i=0;i<nsnp;++i){
      //cout<<"Computing LDAS for SNP "<<i<<endl;
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
        LDAS[i]+=LDA_ave[q]*gdgap[q];
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
        outputFile << " " << fixed << setprecision(4) << LDAS[i] <<" "<< LDAS_lower[i] <<" " << LDAS_upper[i] << "\n";
      }
      outputFile.close();
    } else {
      cerr << "Unable to open file" << LDASfile;
    }
    
  }



void paintingalldense(const string method="Viterbi",
                      bool fixrho=false,
                      const int ite_time=10,
                      const double indfrac=0.1,
                      const int minsnpEM=10000, 
                      const double EMsnpfrac=0.1,
                      int L_initial=200,
                      double minmatchfrac=0.002,
                      int L_minmatch=25,
                      bool haploid=false,
                      const string donorfile="donor.phase.gz",
                      const string targetfile="target.phase.gz",
                      const string mapfile="map.txt",
                      const string popfile="popnames.txt",
                      const string targetname="targetname.txt",
                      bool outputpainting=true,
                      bool outputLDA=true,
                      bool outputLDAS=true,
                      const string paintingfile="painting.txt.gz",
                      const string LDAfile="LDA.txt.gz",
                      const string LDASfile="LDAS.txt",
                      const double window=0.05,
                      const int LDAfactor=1,
                      int ncores=0){
  //detect cores
  if(ncores==0){
    ncores = omp_get_num_procs();
  }
  
  // read the map data to get the genetic distance in Morgans
  tuple<vector<double>,vector<double>> mapinfo = readmap(mapfile);
  vector<double> gd = get<1>(mapinfo);
  vector<double> pd = get<0>(mapinfo);
  
  tuple<vector<string>,vector<int>> popinfo = readpopfile(popfile);
  //vector<string> indnames = get<0>(popinfo);
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
  double rho_use;
  double gdall=gd[nsnp-1]-gd[0];
  int minmatch=static_cast<int>(ceil(nref*minmatchfrac));
  
  if(fixrho){
    
    cout<<"Begin estimating fixed rho"<<endl;
    rho_use=est_rho_average(refidx,nref,nsnp,gd,L_initial,minmatch,L_minmatch,ncores,
                            indfrac,ite_time,method,minsnpEM,EMsnpfrac,haploid,donorfile);
  }
  
  
  cout<<"Do dPBWT for target haplotypes"<<endl;
  
  bool loo=false;
  
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_target;
  
  int qM;
  
  if(donorfile==targetfile){
    loo = true;
    qM=0;
  }else{
    if(haploid){
      qM=nind;
    }else{
      qM=2*nind;
    }
  }
  
  dpbwtall_target=do_dpbwt(L_initial, gd,queryidx,ncores,nref+qM,nsnp,qM,minmatch,
                           L_minmatch,haploid,donorfile,targetfile);
  
  vector<int> queryidall_target=get<0>(dpbwtall_target);
  vector<int> donorid_target=get<1>(dpbwtall_target);
  vector<int> startpos_target=get<2>(dpbwtall_target);
  vector<int> endpos_target=get<3>(dpbwtall_target);
  cout<<"dPBWT works successfully"<<endl;
  
  // begin painting
  // we only store ncores*2 samples in memory and directly output
  
  omp_set_num_threads(ncores);
  
  int nsamples_use;
  int nhap_left=nhap_use;

  
  //store data in hMat if want to compute LDA
  //vector<hMat> painting_all_hmat(nhap_use, hMat(npop, nsnp));
  
  
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
  
  // calculate LDA and store in LDA_result (hMat)
  // should be defined outside of the if command
  hMat LDA_result(nsnp,nsnp,0.0);
  hMat Dprime(nsnp,nsnp,0.0);
  
  if(outputpainting){
    //output the painting into paintingfile
    //ofstream outputFile(paintingfile);
    ogzstream outputFile(paintingfile.c_str());
    outputFile << "indnames" << " ";
    //the first row is the SNP's physical position
    for (int i = 0; i < nsnp; ++i) {
      outputFile << fixed << setprecision(0) << pd[i];
      if(i != nsnp-1) outputFile << " ";
    }
    outputFile << "\n";
    
    int looptime=0;
    
    while(nhap_left>0){
      nsamples_use = (ncores*2*LDAfactor < nhap_left) ? ncores*2*LDAfactor : nhap_left; //ensure both copies are included
      
      vector<vector<vector<double>>> painting_all(nsamples_use, 
                                                  vector<vector<double>>(npop, vector<double>(nsnp)));
      
      // get the matches before the loop
      vector<vector<vector<int>>> targetmatch_use(nsamples_use);
      
      for (int ii = nhap_use - nhap_left; ii < nhap_use - nhap_left + nsamples_use; ++ii) {
        // leave one out if the donor file is the same as the target file
        vector<vector<int>> match_data = get_matchdata(queryidall_target,
                                                       donorid_target,
                                                       startpos_target,
                                                       endpos_target,
                                                       ii, loo);
        
        targetmatch_use[ii - (nhap_use - nhap_left)] = match_data;
      }
      
      
#pragma omp parallel for
      cout<<"Calculating painting for samples "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
      for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
        
        
        vector<vector<int>> targetmatchdata=targetmatch_use[ii - (nhap_use - nhap_left)];
        
        hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
        if(fixrho){
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
          
          //if want to compute LDA
          //painting_all_hmat[ii]=pind;
          
          
          vector<vector<double>> pind_dense=hMatrix2matrix(pind);
          for(int j=0;j<npop;++j){
            for(int k=0;k<nsnp;++k){
              painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
            }
          }
        }else{
          vector<int> startpos, endpos;
          for (const auto& row : targetmatchdata) {
            startpos.push_back(row[1]);
            endpos.push_back(row[2]);
          }
          rho_use=est_rho_Viterbi(startpos,endpos,nsnp,gdall);
          //cout<<"Estimated rho is "<<rho_use<<endl;
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
          
          //if want to compute LDA
          //painting_all_hmat[ii]=pind;
          
          
          vector<vector<double>> pind_dense=hMatrix2matrix(pind);
          for(int j=0;j<npop;++j){
            for(int k=0;k<nsnp;++k){
              painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
            }
          }
        }
      }
      
      //compute LDA
      cout<<"Calculating LDA for samples "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
      if(outputLDA || outputLDAS){
        vector<int> allhaps_idx;
        for(int i=0;i<nsamples_use;++i){
          allhaps_idx.push_back(i);
        }
        vector<int> resample_idx = randomsample(allhaps_idx,nsamples_use);
#pragma omp parallel for
        for(int i=0;i<nsnp;++i){
          LDA_result.m[i].set(i,1.0);
          Dprime.m[i].set(i,1.0);
          //cout<<"Computing LDA for SNP "<<i<<endl;
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
                LDA_result.m[i].set(j,distance);
                Dprime.m[i].set(j,theo_distance);
              }else{
                LDA_result.m[i].set(j,distance+LDA_result.m[i].get(j));
                Dprime.m[i].set(j,theo_distance+Dprime.m[i].get(j));
              }
            }
          }
        }
      }
      
      
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
      vector<vector<vector<double>>>().swap(painting_all);
      nhap_left=nhap_left-nsamples_use;
      looptime++;
    }
    
    outputFile.close();
  }else{
    //don't output the painting into paintingfile
    
    int looptime=0;
    
    while(nhap_left>0){
      nsamples_use = (ncores*2*LDAfactor < nhap_left) ? ncores*2*LDAfactor : nhap_left; //ensure both copies are included
      
      vector<vector<vector<double>>> painting_all(nsamples_use, 
                                                  vector<vector<double>>(npop, vector<double>(nsnp)));
      
      // get the matches before the loop
      vector<vector<vector<int>>> targetmatch_use(nsamples_use);
      
      for (int ii = nhap_use - nhap_left; ii < nhap_use - nhap_left + nsamples_use; ++ii) {
        // leave one out if the donor file is the same as the target file
        vector<vector<int>> match_data = get_matchdata(queryidall_target,
                                                       donorid_target,
                                                       startpos_target,
                                                       endpos_target,
                                                       ii, loo);
        
        targetmatch_use[ii - (nhap_use - nhap_left)] = match_data;
      }
      
      
#pragma omp parallel for
      cout<<"Calculating painting for samples "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
      for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
        
        vector<vector<int>> targetmatchdata=targetmatch_use[ii - (nhap_use - nhap_left)];
        
        hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
        if(fixrho){
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
          
          //if want to compute LDA
          //painting_all_hmat[ii]=pind;
          
          
          vector<vector<double>> pind_dense=hMatrix2matrix(pind);
          for(int j=0;j<npop;++j){
            for(int k=0;k<nsnp;++k){
              painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
            }
          }
        }else{
          vector<int> startpos, endpos;
          for (const auto& row : targetmatchdata) {
            startpos.push_back(row[1]);
            endpos.push_back(row[2]);
          }
          rho_use=est_rho_Viterbi(startpos,endpos,nsnp,gdall);
          //cout<<"Estimated rho is "<<rho_use<<endl;
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
          
          //if want to compute LDA
          //painting_all_hmat[ii]=pind;
          
          
          vector<vector<double>> pind_dense=hMatrix2matrix(pind);
          for(int j=0;j<npop;++j){
            for(int k=0;k<nsnp;++k){
              painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
            }
          }
        }
      }
      
      //compute LDA
      cout<<"Calculating LDA for samples "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
      if(outputLDA || outputLDAS){
        vector<int> allhaps_idx;
        for(int i=0;i<nsamples_use;++i){
          allhaps_idx.push_back(i);
        }
        vector<int> resample_idx = randomsample(allhaps_idx,nsamples_use);
#pragma omp parallel for
        for(int i=0;i<nsnp;++i){
          LDA_result.m[i].set(i,1.0);
          Dprime.m[i].set(i,1.0);
          //cout<<"Computing LDA for SNP "<<i<<endl;
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
                LDA_result.m[i].set(j,distance);
                Dprime.m[i].set(j,theo_distance);
              }else{
                LDA_result.m[i].set(j,distance+LDA_result.m[i].get(j));
                Dprime.m[i].set(j,theo_distance+Dprime.m[i].get(j));
              }
            }
          }
        }
      }

      vector<vector<vector<double>>>().swap(painting_all);
      nhap_left=nhap_left-nsamples_use;
      looptime++;
    }
    
  }
  
  if(outputLDA){
    //output the LDA results into LDAfile
    ogzstream outputFile(LDAfile.c_str());
    if (outputFile) {
      for (int i = 0; i < nsnp; ++i) {
        vector<int> keys = LDA_result.m[i].k;
        for (int j = 0; j < keys.size(); ++j) {
          LDA_result.m[i].set(keys[j],1-LDA_result.m[i].get(keys[j])/Dprime.m[i].get(keys[j]));
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
    }
  }
  
  if(outputLDAS){
    hMat Dprime(0,0,0.0);
    cout << "Begin calculating LDA score"<<endl;
    LDAS(LDA_result,LDASfile,window,gd,pd,nsnp_left,nsnp_right,nsnp);
    cout << "Finish calculating LDA score"<<endl;
  }
  
}
  
  
int main(int argc, char *argv[]){
    std::string run="paint";
    std::string method="Viterbi";
    bool fixrho=false;
    int ite_time=10;
    double indfrac=0.1;
    int minsnpEM=10000;
    double EMsnpfrac=0.1;
    int L_initial=200;
    double minmatchfrac=0.002;
    int L_minmatch=25;
    bool haploid=false;
    std::string donorfile="donor.phase.gz";
    std::string targetfile="target.phase.gz";
    std::string mapfile="map.txt";
    std::string popfile="popnames.txt";
    std::string targetname="targetname.txt";
    bool outputpainting=true;
    bool outputLDA=true;
    bool outputLDAS=true;
    std::string paintingfile="painting.txt.gz";
    std::string LDAfile="LDA.txt.gz";
    std::string LDASfile="LDAS.txt";
    double window=0.05;
    int LDAfactor=1;
    int ncores=0;
    std::string chunklengthfile="chunklength.txt";
    
    for (int i = 1; i < argc; i+=2) {
        std::string param = argv[i];
        if (param[0] != '-') {
            std::cerr << "Invalid argument format. Expected -param value\n";
            return 1;
        }
        param = param.substr(1);  // Remove the -
        
        if (param == "run") {
          run = argv[i+1];
        }else if (param == "method") {
            method = argv[i+1];
        } else if (param == "fixrho") {
            fixrho = std::stoi(argv[i+1]);
        } else if (param == "ite_time") {
            ite_time = std::stoi(argv[i+1]);
        } else if (param == "indfrac") {
            indfrac = std::stod(argv[i+1]);
        } else if (param == "minsnpEM") {
            minsnpEM = std::stoi(argv[i+1]);
        } else if (param == "EMsnpfrac") {
            EMsnpfrac = std::stod(argv[i+1]);
        } else if (param == "L_initial") {
            L_initial = std::stoi(argv[i+1]);
        } else if (param == "minmatchfrac") {
            minmatchfrac = std::stod(argv[i+1]);
        } else if (param == "L_minmatch") {
            L_minmatch = std::stoi(argv[i+1]);
        } else if (param == "haploid") {
            haploid = std::stoi(argv[i+1]);
        } else if (param == "donorfile") {
            donorfile = argv[i+1];
        } else if (param == "targetfile") {
            targetfile = argv[i+1];
        } else if (param == "mapfile") {
            mapfile = argv[i+1];
        } else if (param == "popfile") {
            popfile = argv[i+1];
        } else if (param == "targetname") {
            targetname = argv[i+1];
        } else if (param == "outputpainting") {
            outputpainting = std::stoi(argv[i+1]);
        } else if (param == "outputLDA") {
            outputLDA = std::stoi(argv[i+1]);
        } else if (param == "outputLDAS") {
            outputLDAS = std::stoi(argv[i+1]);
        } else if (param == "paintingfile") {
            paintingfile = argv[i+1];
        } else if (param == "LDAfile") {
          LDAfile = argv[i+1];
        } else if (param == "LDASfile") {
          LDASfile = argv[i+1];
        } else if (param == "window") {
          window = std::stod(argv[i+1]);
        } else if (param == "LDAfactor") {
          LDAfactor = std::stoi(argv[i+1]);
        } else if (param == "ncores") {
          ncores = std::stoi(argv[i+1]);
        } else if (param == "chunklengthfile") {
          chunklengthfile = argv[i+1];
        } else {
          std::cerr << "Unknown parameter: " << param << "\n";
          return 1;
        }
    }
    if(run=="paint"){
      paintingalldense(method, fixrho, ite_time, indfrac, minsnpEM, EMsnpfrac, L_initial, minmatchfrac, 
                       L_minmatch, haploid, donorfile, targetfile, mapfile, popfile, targetname, outputpainting,
                       outputLDA, outputLDAS, paintingfile, LDAfile, LDASfile, window, LDAfactor, ncores);
    }else if (run=="chunklength"){
      chunklengthall(method,ite_time,indfrac,minsnpEM,EMsnpfrac,
                     L_initial,minmatchfrac,L_minmatch,haploid,
                     donorfile,mapfile,popfile,chunklengthfile,ncores);
    }else if (run=="both"){
      paintingalldense(method, fixrho, ite_time, indfrac, minsnpEM, EMsnpfrac, L_initial, minmatchfrac, 
                       L_minmatch, haploid, donorfile, targetfile, mapfile, popfile, targetname, outputpainting,
                       outputLDA, outputLDAS, paintingfile, LDAfile, LDASfile, window, LDAfactor, ncores);
      chunklengthall(method,ite_time,indfrac,minsnpEM,EMsnpfrac,
                     L_initial,minmatchfrac,L_minmatch,haploid,
                     donorfile,mapfile,popfile,chunklengthfile,ncores);
    }else{
      std::cerr << "Unknown argument given to run" << "\n";
    }
    return 0;
} 
