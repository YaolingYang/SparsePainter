#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include "zlib.h"
#include "pbwt.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
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
                            vector<double>& gd,
                            const int nref){
  //compute the sameprob
  vector<double> sameprob(nsnp);
  for(int j=0;j<nsnp-1;++j){
    sameprob[j]=exp(-rho*(gd[j+1]-gd[j]));
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

hMat indpainting(const hMat& mat,
                 vector<double>& gd, 
                 const double rho,
                 const int npop, 
                 const vector<int>& refindex){
  // return individual painting
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob=cal_sameprob(nsnp,rho,gd,nref);
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
      if(popprob[i]>=0.005){
        marginal_prob_pop.m[j].set(i,popprob[i]);
      }
    }
  }
  return(marginal_prob_pop);
}

void painting(PBWT *p,
              const int L=20,
              bool diff_rho=false,
              const double fixrho=0,
              const double indfrac=0.1,
              bool haploid=false,
              const string mapfile,
              const string popfile,
              const string targetname,
              bool outputpainting=true,
              bool outputaveSNPpainting=true,
              bool outputaveindpainting=true,
              const string paintingfile="painting.txt.gz",
              const string aveSNPpaintingfile="aveSNPpainting.txt.gz",
              const string aveindpaintingfile="painting.txt.gz",
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
  int minmatch=static_cast<int>(ceil(nref*matchfrac));

  // convert pbwt to match_vector
  
  vector<int> queryidall_target;
  vector<int> donorid_target;
  vector<int> startpos_target;
  vector<int> endpos_target;
  
  auto reportMatch = [&](int i, int j, int start, int end) {
    donorid_target.push_back(j);
    startpos_target.push_back(start);
    endpos_target.push_back(end);
    queryidall_target.push_back(i);
  };
  
  matchLongWithin2(p,L, reportMatch); //this stores all matches longer than L

  
  // estimate rho
  
  if(!diff_rho){
    if(fixrho!=0){
      rho_use=fixrho;
      cout << "Using fixed rho "<<rho_use<<endl;
    }else{
      cout<<"Begin estimating fixed rho"<<endl;
      int v_nsamples=static_cast<int>(ceil(queryidx.size()*indfrac));
      vector<int> v_samples=randomsample(queryidx,v_nsamples);
      double rho_sum = 0.0;
      
      // get the matches before the loop
      int v_nhap_left=v_nsamples;
      while(v_nhap_left>0){
        int v_nsamples_use = (ncores < v_nhap_left) ? ncores : v_nhap_left;
        vector<vector<vector<int>>> v_targetmatch_use(v_nsamples_use);
        
        for (int ii = v_nsamples - v_nhap_left; ii < v_nsamples - v_nhap_left + v_nsamples_use; ++ii) {
          // leave one out if the donor file is the same as the target file
          vector<vector<int>> match_data = get_matchdata(queryidall_target,
                                                         donorid_target,
                                                         startpos_target,
                                                         endpos_target,
                                                         v_samples[ii]);
          
          v_targetmatch_use[ii - (v_nsamples - v_nhap_left)] = match_data;
          
#pragma omp parallel for reduction(+:rho_sum)
          for(int ii = v_nsamples - v_nhap_left; ii < v_nsamples - v_nhap_left + v_nsamples_use; ++ii){
            vector<vector<int>> v_targetmatchdata=v_targetmatch_use[ii - (v_nsamples - v_nhap_left)];
            
            vector<int> v_startpos, v_endpos;
            for (const auto& row : v_targetmatchdata) {
              v_startpos.push_back(row[1]);
              v_endpos.push_back(row[2]);
            }
            double rho_estimated=est_rho_Viterbi(v_startpos,v_endpos,nsnp,gdall);
            rho_sum += rho_estimated;
          }
          v_nhap_left=v_nhap_left-v_nsamples_use;
        }
        
        rho_use=rho_sum/v_nsamples;
      }
      cout << "Using fixed rho "<<rho_use<<endl;
    }
  }
  
  // begin painting
  // we only store ncores*2 samples in memory and directly output
  
  omp_set_num_threads(ncores);
  
  int nsamples_use;
  int nhap_left=nhap_use;
  
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
  
  if(outputpainting){
    //output the painting into paintingfile
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
      nsamples_use = (ncores*2 < nhap_left) ? ncores*2 : nhap_left; //ensure both copies are included
      
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
                                                       ii);
        
        targetmatch_use[ii - (nhap_use - nhap_left)] = match_data;
      }
      
      cout<<"Calculating painting for samples "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
#pragma omp parallel for
      for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
        
        
        vector<vector<int>> targetmatchdata=targetmatch_use[ii - (nhap_use - nhap_left)];
        
        hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
        if(!diff_rho){
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
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
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
          
          
          vector<vector<double>> pind_dense=hMatrix2matrix(pind);
          for(int j=0;j<npop;++j){
            for(int k=0;k<nsnp;++k){
              painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
            }
          }
        }
      }
      
      //compute average painting for each SNP
      if(outputaveSNPpainting){
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
      nsamples_use = (ncores*2 < nhap_left) ? ncores*2 : nhap_left; //ensure both copies are included
      
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
                                                       ii);
        
        targetmatch_use[ii - (nhap_use - nhap_left)] = match_data;
      }
      
      cout<<"Calculating painting for samples "<<nhap_use-nhap_left<<"-"<<nhap_use-nhap_left+nsamples_use-1<<endl;
#pragma omp parallel for
      
      for(int ii=nhap_use-nhap_left; ii<nhap_use-nhap_left+nsamples_use; ++ii){
        
        vector<vector<int>> targetmatchdata=targetmatch_use[ii - (nhap_use - nhap_left)];
        
        hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
        
        if(!diff_rho){
          hMat pind=indpainting(mat,gd,rho_use,npop,refindex);
          
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
          
          vector<vector<double>> pind_dense=hMatrix2matrix(pind);
          for(int j=0;j<npop;++j){
            for(int k=0;k<nsnp;++k){
              painting_all[ii-nhap_use+nhap_left][j][k]=pind_dense[j][k];
            }
          }
        }
      }
      
      
      //compute average painting for each SNP
      if(outputaveSNPpainting){
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
      
      vector<vector<vector<double>>>().swap(painting_all);
      nhap_left=nhap_left-nsamples_use;
      looptime++;
    }
    
  }
  
  // get the average painting for each SNP and output
  if(outputaveSNPpainting){
    for(int j=0;j<npop;++j){
#pragma omp parallel for
      for(int k=0;k<nsnp;++k){
        aveSNPpainting[j][k]=aveSNPpainting[j][k]/nhap_use;
      }
    }
    
    if(outputaveSNPpainting){
      //output the average painting for each SNP
      ogzstream outputFile(aveSNPpaintingfile.c_str());
      if (outputFile) {
        outputFile << "population" << " ";
        //the first row is the SNP's physical position
        for (int k = 0; k < nsnp; ++k) {
          outputFile << fixed << setprecision(0) << pd[k];
          if(k != nsnp-1) outputFile << " ";
        }
        outputFile << "\n";
        
        for (int j = 0; j < npop; ++j) {
          outputFile << "pop"<<j << " ";
          for (int k = 0; k < nsnp; ++k) {
            outputFile << fixed << setprecision(4) <<aveSNPpainting[j][k];
            if(k != nsnp-1) outputFile << " ";
          }
          if(j != npop-1) outputFile<< "\n";
        }
        
        outputFile.close();
      } else {
        cerr << "Unable to open file" << aveSNPpaintingfile;
      }
    }
    
  }
  
  
  //output the average painting for each individual
  if(outputaveindpainting){
    //output the average painting for each SNP
    ogzstream outputFile(aveindpaintingfile.c_str());
    if (outputFile) {
      outputFile << "individual" << " ";
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
      cerr << "Unable to open file" << aveindpaintingfile;
    }
  }
  
}


