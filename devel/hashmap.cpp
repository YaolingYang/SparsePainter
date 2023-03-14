// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <Rcpp.h>
#include <unordered_map>
#include <algorithm>
#include <random>

#include<iomanip>
#include<sstream>
#include<utility>

using namespace Rcpp;
using namespace std;

int M = 0;
int qM = 0;
int N = 0;
int *dZ;

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

void ReadVCF(string inFile, string qinFile, bool ** & panel){
  using namespace std;
  {//count M and N
  string line = "1#";
  ifstream in(inFile);
  stringstream linestr;
  while (line[1] == '#')
    getline(in, line);
  
  
  linestr.str(line);
  linestr.clear();
  for (int i = 0; i<9;++i)
    linestr>>line;
  N = M = 0;
  while (!linestr.eof()){
    ++M;
    ++M;
    linestr >> line;
  }
  
  while (getline(in, line))
    ++N;
  in.close();
  }
  {//count qM, finish M
    string line = "1#";
    ifstream qin(qinFile);
    while(line[1] == '#')
      getline(qin, line);
    stringstream linestr;
    
    linestr.str(line);
    linestr.clear();
    for (int i = 0; i<9; ++i)
      linestr >> line;
    qM = 0;
    while (!linestr.eof()){
      ++qM;
      ++qM;
      linestr >> line;
    }
    M +=qM;
    int qN = 0;
    while (getline(qin, line)){
      ++qN;
    }
    if (qN != N){
      cout << "Query file and input file have different numbers of sites. Query has " << qN << ". Panel has " << N << endl;
      throw "Query and input file have different numbers of sites";
    }
    qin.close();
  }
  bool *temp = new bool[(long long)M * N];
  panel = new bool*[M];
  for (int i = 0; i<M; i++){
    panel[i] = &(temp[i*N]);
  }
  string line = "##", qline = "##";
  ifstream in (inFile);
  ifstream qin(qinFile);
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
  }
  in.close();
  qin.close();
}


vector<int> getorder(const vector<double>& vec) {
  
  vector<int> order(vec.size());
  iota(order.begin(), order.end(), 0);
  
  sort(order.begin(), order.end(), [&vec](int i, int j) {
    if (vec[i] == vec[j]) {
      return i < j;
    }
    return vec[i] < vec[j];
  });
  
  return order;
}

bool containsIndex(const vector<int>& fullidx, int starttemp, int endtemp) {
  bool contain=false;
  for(int i : fullidx) {
    if(i >= starttemp && i <= endtemp) {
      contain=true;
    }
  }
  return(contain);
}


tuple<vector<int>,vector<int>,vector<int>,vector<int>> longMatchdpbwt(int& L,
                                                                      bool **panel, 
                                                                      dpbwt & x,
                                                                      const int minmatch,
                                                                      vector<double> &gd){
  const int L0=L;
  
  dZ = new int[M];
  for (int i = 0; i<M; i++){
    dZ[i] = 0;
  }
  
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
  
  for (int i = 0; i < M-qM; i++){
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
      for (int k = 0; k<N; ++k){
        dpbwtnode* temp = z[k].above;
        while (temp != nullptr && x.panel[temp->id][k] != panel[Oid][k]){
          if (!panel[Oid][k])
            temp->u = &z[k+1];
          else 
            temp->v = &z[k+1];
          temp = temp->above;
        }
        if (temp == nullptr && panel[Oid][k]){
          temp = botk->above;
          botk->u = &z[k+1];
          while (temp != &z[k] && x.panel[temp->id][k]){
            temp->u = &z[k+1];
            temp = temp->above;
          }
        }
        if (!panel[Oid][k]){
          z[k].u = &z[k+1];
          z[k].v = t->v;
        }
        else{
          z[k].u = t->u;
          z[k].v = &z[k+1];
        }
        t = (panel[Oid][k])? t->v : t->u;
        z[k+1].id = i;
        z[k+1].originalid = Oid;
        z[k+1].below = t;
        z[k+1].above = t->above;
        t->above = &z[k+1];
        if (z[k+1].above != nullptr)
          z[k+1].above->below = &z[k+1];
        botk = botk->v;
        x.panel[i][k] = panel[Oid][k];
      }
      
      int zdtemp = N,
        bdtemp = N;
      
      for (int k = N; k>= 0; --k){
        zdtemp = min(zdtemp, k);
        bdtemp = min(bdtemp, k);
        if (z[k].above!=nullptr)
          while (zdtemp>0 && panel[Oid][zdtemp-1] == x.panel[z[k].above->id][zdtemp-1])
            zdtemp--;
        else 
          zdtemp = k;
        if (z[k].below->id!=-1)
          while (bdtemp>0 && panel[Oid][bdtemp-1] == x.panel[z[k].below->id][bdtemp-1])
            bdtemp--;
        else 
          bdtemp = k;
        z[k].divergence = zdtemp;
        z[k].below->divergence = bdtemp;
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
  
  
  for (int i = M-qM; i<M; i++){
    cout<<i<<endl;
    L=L0;
    int prevL=L0;
    
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
        
        t = (panel[Oid][k])? t->v : t->u;
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
          while(zdtemp>0 && panel[Oid][zdtemp-1] == x.panel[z[k].above->id][zdtemp-1])
            zdtemp--;
        else
          zdtemp = k;
        if (z[k].below->id!=-1)
          while (bdtemp>0 && panel[Oid][bdtemp-1] == x.panel[z[k].below->id][bdtemp-1])
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
      
      while(!allsnpmatches){
        

        
        vector<bool> addmatch(N,false);
        
        //int prevL=L*2;
        
        if(times!=0){
          for(int w=0;w<nomatchsnp.size();++w){
            for(int s=0;s<prevL;++s){
              int pos=nomatchsnp[w]+s;
              if(pos>=N) break;
              addmatch[pos]=true;
            }
          }
        }


        
        //donoridtemp.clear();
        //startpostemp.clear();
        //endpostemp.clear();
        
        dpbwtnode *f, *g, *ftemp, *gtemp;
        f = g = z[0].below;
        for(int k = 0; k<N; ++k){
          ftemp = (panel[Oid][k])? f->u : f->v;
          gtemp = (panel[Oid][k])? g->u : g->v;
          f = (panel[Oid][k])? f->v : f->u;
          g = (panel[Oid][k])? g->v : g->u;
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
        
        if(times==0){
          for(int j=0;j<N;++j){
            if (nmatch[j] < minmatch) {
              nomatchsnp.push_back(j);
            }
          }
        }else{
          vector<int> nomatchsnptemp=nomatchsnp;
          nomatchsnp.clear();
          for(int j=0;j<nomatchsnptemp.size();++j){
            if (nmatch[nomatchsnptemp[j]] < minmatch) {
              nomatchsnp.push_back(nomatchsnptemp[j]);
            }
          }
        }
        if(nomatchsnp.size()==0 || L==1){
          allsnpmatches = true;
        }else{
          prevL=L;
          L=(prevL+1)/2; // this is equal to ceil(L/2) is L is double type
        }
        times++;
      }
      
      delete [] z;
      
      
      //below we remove shorter matches while ensuring at least minmatch matches at each SNP
      //information of matches is stored in donoridtemp, startpostemp and endpostemp
      //number of matches at each SNP are stored in nmatch
      //we first sort the genetic distance of each match
      
      vector<double> gdmatch(startpostemp.size());
      for(int mi=0; mi<startpostemp.size();++mi){
        gdmatch[mi]=gd[endpostemp[mi]]-gd[startpostemp[mi]];
      }
      vector<int> length_order=getorder(gdmatch);
      
      vector<int> fullidx; // record which SNP has only minmatch matches
      for(int q=0;q<N;++q){
        if(nmatch[q]<=minmatch) fullidx.push_back(q);
      }
      for(int mi=0;mi<length_order.size();++mi){
        
        int starttemp=startpostemp[length_order[mi]];
        int endtemp=endpostemp[length_order[mi]];
        
        if(!containsIndex(fullidx,starttemp,endtemp)){
          for(int q=starttemp;q<=endtemp;++q){
            nmatch[q]--;
            if(nmatch[q]==minmatch) fullidx.push_back(q);
          }
        }else{
          startpos.push_back(starttemp);
          endpos.push_back(endtemp);
          donorid.push_back(donoridtemp[length_order[mi]]);
        }
      }
      //record the position of the next start position of query haplotype
      //such that we know how many matches are there for this query haplotype
      queryidall.push_back(startpos.size());
  }
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> results(queryidall,donorid,startpos,endpos);
  return(results);
}



tuple<vector<int>,vector<int>,vector<int>,vector<int>> do_dpbwt(int& L, 
                                                                vector<double> gd,
                                                                string query="target",
                                                                int minmatch=100){
  string qin = "p_" + query + ".vcf";
  string in = "p_donor.vcf";
  
  bool **panel;
  dpbwt x;
  ReadVCF(in, qin, panel);
  
  while(L>N){
    L=ceil(L/2);
    cout<<"L cannot be greater than N, reducing L to "<<L<<endl;
  }
  
  return(longMatchdpbwt(L,panel,x,minmatch,gd));
  
}

vector<vector<int>> get_matchdata(vector<int> queryidall,
                                  vector<int> donorid,
                                  vector<int> startpos,
                                  vector<int> endpos,
                                  int queryid){
  int querystart=queryidall[queryid];
  int nextquerystart=queryidall[queryid+1];
  int nrow_match=nextquerystart-querystart;
  vector<vector<int>> matchinfo(nrow_match,vector<int>(3));
  for(int p=querystart;p<nextquerystart;++p){
    matchinfo[p-querystart][0]=donorid[p];
    matchinfo[p-querystart][1]=startpos[p];
    matchinfo[p-querystart][2]=endpos[p];
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
    fprev=forward_prob.m[j-1].getall(twj);
    sameprobuse=sameprob[j-1];
    otherprobuse=otherprob[j-1];
    for(int i=0;i<twj.size();++i){
      forward_prob.m[j].set(twj[i],sameprobuse*fprev[i]+otherprobuse);
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
  for(int j=ncol-2;j>=0;--j){
    twj=mat.m[j+1].k;
    Bjp1=backward_prob.m[j+1].getall(twj);
    sumBjp1=vec_sum(Bjp1);
    sameprobuse=sameprob[j];
    otherprobuse=otherprob[j];
    vector<double> val(twj.size());
    for(int i=0;i<twj.size();++i){
      val[i]=sameprobuse*Bjp1[i]+otherprobuse*sumBjp1;
      backward_prob.m[j].set(twj[i],val[i]);
    }
    logmultB[j]=log(vec_sum(val)+otherprob[j]*sumBjp1*(nrow-twj.size()))+logmultB[j+1];
    hVecScale(backward_prob.m[j],1.0/sumBjp1);
    backward_prob.m[j].setdefault(otherprobuse);
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
                       const vector<int>& queryidall, 
                       const vector<int>& donorid_ref,
                       const vector<int>& startpos_ref,
                       const vector<int>& endpos_ref,
                       const double indfrac=0.1,
                       const int ite_time=10, 
                       const string method="Viterbi",
                       const int minsnpEM=10000, 
                       const double EMsnpfrac=0.1){
  // estimate \rho from the reference panel
  int npop=refidx.pos.size();
  vector<double> rho_est;
  int count=0;
  double gdall=gd[nsnp-1]-gd[0];
  
  for(int i=0;i<npop;++i){
    //randomly sample a percentage of indfrac reference samples
    vector<int> popidx=refidx.findrows(i);
    
    int nref_sample=static_cast<int>(ceil(popidx.size()*indfrac));
    
    vector<int> samples=randomsample(popidx,nref_sample);
    
    for(int k=0;k<samples.size();k++){
      //leave-one-out
      //vector<vector<int>> matchdata=read_data("donor",samples[k]);
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
        cout<<"Estimated rho for sample "<<count<<" is "<<rho_estimated<<endl;
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
          cout<<"Estimated rho for sample "<<count<<" is "<<rho_estimated<<endl;
          rho_est.push_back(rho_estimated);
        }else{
          
          double rho_estimated=est_rho_EM(mat,gd,ite_time);
          count=count+1;
          cout<<"Estimated rho for sample "<<count<<" is "<<rho_estimated<<endl;
          rho_est.push_back(rho_estimated);
        }
      }
    }
  }
  double rho_ave=vec_sum(rho_est)/rho_est.size();
  cout<<"Average estimated rho is "<<rho_ave<<endl;
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
      marginal_prob_pop.m[j].set(i,popprob[i]);
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
vector<vector<double>> chunklengthall(vector<double>& gd,
                                      const vector<int>& refindex, 
                                      const string method="Viterbi", 
                                      const int ite_time=10,
                                      const double indfrac=0.1,
                                      const int minsnpEM=10000, 
                                      const double EMsnpfrac=0.1,
                                      int L=500,
                                      double minmatchfrac=0.001){
  //compute coancestry matrix for all reference individuals
  const int nsnp=gd.size();
  const int nref=refindex.size();
  hAnc refidx(refindex);
  const int npop=refidx.pos.size();
  vector<vector<double>> chunklength(nref, vector<double>(npop));
  double rho;
  int minmatch=static_cast<int>(ceil(nref*minmatchfrac));
  cout<<"Do dPBWT for donor haplotypes"<<endl;
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_ref=do_dpbwt(L, gd,"donor",minmatch);
  vector<int> queryidall_ref=get<0>(dpbwtall_ref);
  vector<int> donorid_ref=get<1>(dpbwtall_ref);
  vector<int> startpos_ref=get<2>(dpbwtall_ref);
  vector<int> endpos_ref=get<3>(dpbwtall_ref);
  cout<<"dPBWT works successfully"<<endl;

  cout<<"Begin estimating fixed rho"<<endl;
  rho=est_rho_average(refidx,nref,nsnp,gd,queryidall_ref,
                      donorid_ref,startpos_ref,endpos_ref,indfrac,
                      ite_time,method,minsnpEM,EMsnpfrac);
  
  for(int i=0;i<nref;++i){
    cout<<"Calculating chunk length for donor sample "<<i+1<<endl;
    //leave-one-out
    vector<vector<int>> matchdata=get_matchdata(queryidall_ref,
                                                donorid_ref,
                                                startpos_ref,
                                                endpos_ref,
                                                i);
    vector<int> removeidx;
    int popidx=refindex[i];
    for(int j=0;j<npop;++j){
      if(j==popidx){
        removeidx.push_back(i);
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
  return(chunklength);
}



// [[Rcpp::export]]
vector<vector<vector<double>>> paintingalldense(vector<double>& gd,
                                                const vector<int>& refindex,
                                                const int nind,
                                                const string method="Viterbi",
                                                bool fixrho=true,
                                                const int ite_time=10,
                                                const double indfrac=0.1,
                                                const int minsnpEM=10000, 
                                                const double EMsnpfrac=0.1,
                                                int L=500,
                                                double minmatchfrac=0.001){
  //compute painting for all target individuals
  const int nsnp=gd.size();
  const int nref=refindex.size();
  hAnc refidx(refindex);
  const int npop=refidx.pos.size();
  vector<vector<vector<double>>> painting_all(nind, vector<vector<double>>(npop, vector<double>(nsnp)));
  double rho_use;
  double gdall=gd[nsnp-1]-gd[0];
  int minmatch=static_cast<int>(ceil(nref*minmatchfrac));
  if(fixrho){
    cout<<"do dPBWT for donor haplotypes"<<endl;
    tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_ref=do_dpbwt(L, gd,"donor",minmatch);
    vector<int> queryidall_ref=get<0>(dpbwtall_ref);
    vector<int> donorid_ref=get<1>(dpbwtall_ref);
    vector<int> startpos_ref=get<2>(dpbwtall_ref);
    vector<int> endpos_ref=get<3>(dpbwtall_ref);
    cout<<"dPBWT works successfully"<<endl;
    
    cout<<"Begin estimating fixed rho"<<endl;
    rho_use=est_rho_average(refidx,nref,nsnp,gd,queryidall_ref,donorid_ref,
                            startpos_ref,endpos_ref,indfrac,ite_time,
                            method,minsnpEM,EMsnpfrac);
    // clear the space occupied
    donorid_ref.clear();
    donorid_ref.shrink_to_fit();
    startpos_ref.clear();
    startpos_ref.shrink_to_fit();
    endpos_ref.clear();
    endpos_ref.shrink_to_fit();
  }
  
  
  cout<<"Do dPBWT for target haplotypes"<<endl;
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_target=do_dpbwt(L, gd,"target",minmatch);
  vector<int> queryidall_target=get<0>(dpbwtall_target);
  vector<int> donorid_target=get<1>(dpbwtall_target);
  vector<int> startpos_target=get<2>(dpbwtall_target);
  vector<int> endpos_target=get<3>(dpbwtall_target);
  cout<<"dPBWT works successfully"<<endl;
  for(int ii=0;ii<nind;++ii){
    cout<<"Calculating painting for individual "<<ii+1<<endl;
    vector<vector<int>> targetmatchdata=get_matchdata(queryidall_target,
                                                      donorid_target,
                                                      startpos_target,
                                                      endpos_target,
                                                      ii);
    hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
    if(fixrho){
      hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
      vector<vector<double>> pind_dense=hMatrix2matrix(pind);
      for(int j=0;j<npop;++j){
        for(int k=0;k<nsnp;++k){
          painting_all[ii][j][k]=pind_dense[j][k];
        }
      }
    }else{
      vector<int> startpos, endpos;
      for (const auto& row : targetmatchdata) {
        startpos.push_back(row[1]);
        endpos.push_back(row[2]);
      }
      rho_use=est_rho_Viterbi(startpos,endpos,nsnp,gdall);
      cout<<"Estimated rho is "<<rho_use<<endl;
      hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
      vector<vector<double>> pind_dense=hMatrix2matrix(pind);
      for(int j=0;j<npop;++j){
        for(int k=0;k<nsnp;++k){
          painting_all[ii][j][k]=pind_dense[j][k];
        }
      }
    }
  }
  return(painting_all);
}




// [[Rcpp::export]]
vector<double> checksparsitytarget(const int nind,
                                   const int nref,
                                   const int nsnp, 
                                   vector<double> gd,
                                   int L=500,
                                   double minmatchfrac=0.001){
  int minmatch=static_cast<int>(ceil(nref*minmatchfrac));
  vector<double> sparsity(nind,0);
  cout<<"Do dPBWT for target haplotypes"<<endl;
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_target=do_dpbwt(L, gd,"target",minmatch);
  vector<int> queryidall_target=get<0>(dpbwtall_target);
  vector<int> donorid_target=get<1>(dpbwtall_target);
  vector<int> startpos_target=get<2>(dpbwtall_target);
  vector<int> endpos_target=get<3>(dpbwtall_target);
  cout<<"dPBWT works successfully"<<endl;
  for(int ii=0;ii<nind;++ii){
    int count=0;
    vector<vector<int>> matchdata=get_matchdata(queryidall_target,
                                                donorid_target,
                                                startpos_target,
                                                endpos_target,
                                                ii);
    for(int j=0;j<matchdata.size();++j){
      count=count+matchdata[j][2]-matchdata[j][1];
    }
    sparsity[ii]=static_cast<double>(count)/nref/nsnp;
  }
  return(sparsity);
}

// [[Rcpp::export]]
vector<double> checksparsitydonor(const int nref,
                                  const int nsnp,
                                  vector<double> gd,
                                  int L=500,
                                  double minmatchfrac=0.001){
  int minmatch=static_cast<int>(ceil(nref*minmatchfrac));
  vector<double> sparsity(nref,0);
  cout<<"Do dPBWT for donor haplotypes"<<endl;
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> dpbwtall_ref=do_dpbwt(L, gd,"donor",minmatch);
  vector<int> queryidall_ref=get<0>(dpbwtall_ref);
  vector<int> donorid_ref=get<1>(dpbwtall_ref);
  vector<int> startpos_ref=get<2>(dpbwtall_ref);
  vector<int> endpos_ref=get<3>(dpbwtall_ref);
  cout<<"dPBWT works successfully"<<endl;
  for(int ii=0;ii<nref;++ii){
    int count=0;
    vector<vector<int>> matchdata=get_matchdata(queryidall_ref,
                                                donorid_ref,
                                                startpos_ref,
                                                endpos_ref,
                                                ii);
    for(int j=0;j<matchdata.size();++j){
      count=count+matchdata[j][2]-matchdata[j][1];
    }
    sparsity[ii]=static_cast<double>(count)/nref/nsnp;
  }
  return(sparsity);
}





//vector<vector<int>> read_data(string type, int idx) {
//read the data from disc
//  string filename = type + "_match" + to_string(idx) + ".txt";
//  ifstream infile(filename);
//  vector<vector<string>> data;

//  if (infile) {
//    string line;
//    while (getline(infile, line)) {
//      vector<string> row;
//      size_t pos = 0;
//      string token;

//      for (int i = 0; i < 6; i++) {
//        if (i == 4) {
//          token = line.substr(pos + 1, line.find(",", pos) - pos - 1); // Remove most left and most right characters
//          pos += token.length() + 2; // 2 represents the length of the left and right characters to remove
//        }
//        else if (i == 5) {
//          token = line.substr(pos, line.length() - pos - 1); // Remove most right character
//          pos += 4;
//        }
//        else {
//          token = line.substr(pos, line.find(" ", pos) - pos);
//          pos += token.length() + 1;
//        }
//        if (i == 0 || i == 4 || i == 5) {
//          row.push_back(token);
//        }
//      }
//      data.push_back(row);
//    }
//    infile.close();
//  }

//  vector<vector<int>> data_int(data.size(), vector<int>(3));
//  for (int i = 0; i < data.size(); i++) {
//    for (int j = 0; j < 3; j++) {
//      data_int[i][j] = stoi(data[i][j]);
//      if(j==2) data_int[i][j]=data_int[i][j]-1;
//    }
//  }
//  return data_int;
//}



//vector<hMat> paintingall(vector<double>& gd,const vector<int>& refindex, 
//                         const int nind,const string method="Viterbi", 
//                         bool fixrho=true,const int ite_time=10,
//                         const double indfrac=0.1){
//compute painting for all target individuals
//  const int nsnp=gd.size();
//  const int nref=refindex.size();
//  hAnc refidx(refindex);
//  const int npop=refidx.pos.size();
// vector<hMat> painting_all(nind, hMat(npop, nsnp));
//  double rho_use;
//  if(fixrho){
//    cout<<"estimating fixed rho"<<endl;
//    rho_use=est_rho_average(refidx,nref,nsnp,gd,indfrac,ite_time,method);
//  }
//  for(int ii=0;ii<nind;++ii){
//    cout<<"calculating painting for target sample "<<ii<<endl;
//   vector<vector<int>> targetmatchdata=read_data("target",ii);
//   hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
//    if(fixrho){
//      painting_all[ii]=indpainting(mat,gd,rho_use,npop,refindex);
//   }else{
//      vector<int> startpos, endpos;
//      for (const auto& row : targetmatchdata) {
//        startpos.push_back(row[1]);
//        endpos.push_back(row[2]);
//      }
//      rho_use=est_rho_Viterbi(startpos,endpos,nsnp,gd[nsnp-1]-gd[0]);
//      painting_all[ii]=indpainting(mat,gd,rho_use,npop,refindex);
//    }
//  }
//  return(painting_all);
//}



//hMat forwardBackward(vector<vector<bool>> matchmat,
//                                                 vector<double> sameprob,
//                                                 vector<double> otherprob){
//  hMat mat = matchmat2hMatrix(matchmat);
//  hMat forward_prob=forwardProb(mat,sameprob,otherprob);
//  hMat backward_prob=backwardProb(mat,sameprob,otherprob);
//  hMat marginal_prob=marginalProb(forward_prob,backward_prob);
//  return(marginal_prob);
//}


//Rcpp::XPtr<hMat> forwardBackwardR(vector<vector<bool>> matchmat,
//                                  vector<double> sameprob,
//                                  vector<double> otherprob){
//  hMat tmp=forwardBackward(matchmat,sameprob,otherprob);
//  hMat* pd = new hMat()
//}



//}