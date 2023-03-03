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
using namespace Rcpp;
using namespace std;
// When we want to lookup all entries:
// https://stackoverflow.com/questions/28767234/what-container-to-store-unique-values
// How to wrap an object for return to R
// https://stackoverflow.com/questions/12405655/returning-a-custom-object-from-a-wrapped-method-in-rcpp

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
  hVec(int len,double x0){
    // Create a vector of length len filled with x0
    this->len=len;
    this->x0=x0;
  };
  hVec(vector<int> idx,vector<double> val,int len,double x0){
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
  void setnocheck(int p, double val){
    // Set a value, should be known to be in the keys
    v[p]=val;
  };
  void set(int p, double val){
    //Safely set a value
    if(!in(p)) k.push_back(p);
    setnocheck(p,val);
  };
  void setall(vector<int> p, vector<double> val){
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
  hMat(int d1,int d2,double x0=0.0){
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
  void appendColumn(vector<int> idx,vector<double> vals,double x0){
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
    // here the class of it is std::unordered_map<int, vector<int>>::const_iterator
    // we use auto to simplify
    auto it = pos.find(value);
    if (it != pos.end()) {
      return it->second;
    } else { // if the value doesn't exist, return an empty vector
      return vector<int>();
    }
  };
};



double vec_sum(const vector<double>& l){
  // calculate the sum of a vector
  double sum_l=0;
  for(int i=0;i<l.size();++i){
    sum_l=sum_l+l[i];
  }
  return(sum_l);
}

vector<double> vec_multiply(const vector<double>& l1, const vector<double>& l2){
  // calculate the multiply of two vectors
  vector<double> multi_l(l1.size());
  for(int i=0;i<l1.size();++i){
    multi_l[i]=l1[i]*l2[i];
  }
  return(multi_l);
}

vector<int> union_vec(const vector<int>& vec1, const vector<int>& vec2) {
  // get the union of two vectors
  vector<int> result(vec1.size() + vec2.size());
  vector<int>::iterator it;
  it = set_union(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), result.begin());
  result.resize(it - result.begin());
  return(result);
}

vector<int> intesect_vec(const vector<int>& vec1, const vector<int>& vec2) {
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

void update_vec(vector<double>& a, const vector<int>& b, const vector<double>& c) {
  // change the position of b inside a to c
  for (int i = 0; i < b.size(); ++i) {
    a[b[i]] = c[i];
  }
}

double hVecSum(hVec& v,vector<int> idx={}){
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

void hVecScale(hVec& v,const double x){
  //Scale a hash vector by x
  vector<int> idx=v.k;
  vector<double> val=v.getall(idx);
  v.setdefault(x*v.x0);
  vector<double> xvec(idx.size(),x);
  v.setall(idx,vec_multiply(val,xvec));
}


vector<double> hVecdense(hVec& x,const int nrow){
  //Convert a hash vector x into a dense vector with length nrow
  vector<double> vdense(nrow, x.x0);
  vector<int> k=x.k;
  vector<double> vval=x.getall(k);
  update_vec(vdense,k,vval);
  return(vdense);
}


hMat matchmat2hMatrix(const vector<vector<bool>>& matchmat,const double default_val=0.0){
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


std::vector<int> randomsample(const std::vector<int>& popidx, const int number) {
  if(number>popidx.size()) cout<<"number cannot be greater than the size of popidx"<<endl;
  // Initialize the random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  
  // Shuffle the elements of the vector randomly
  std::vector<int> shuffled_popidx = popidx;
  std::shuffle(shuffled_popidx.begin(), shuffled_popidx.end(), gen);
  
  // Take the first number elements of the shuffled vector to create the random sample
  std::vector<int> random_sample(shuffled_popidx.begin(), shuffled_popidx.begin() + number);
  
  return random_sample;
}

vector<double> cal_sameprob(const int nsnp, const double rho, vector<double>& gd){
  //compute the sameprob
  vector<double> sameprob(nsnp);
  for(int j=0;j<nsnp-1;++j){
    sameprob[j]=exp(-rho*(gd[j+1]-gd[j]));
  }
  return(sameprob);
}

vector<double> cal_otherprob(const int nref, const vector<double>& sameprob){
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

hMat forwardBackward(const hMat mat,
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

vector<vector<int>> read_data(string type, int idx) {
  //read the data from disc
  string filename = type + "_match" + to_string(idx) + ".txt";
  ifstream infile(filename);
  vector<vector<string>> data;
  
  if (infile) {
    string line;
    while (getline(infile, line)) {
      vector<string> row;
      size_t pos = 0;
      string token;
      
      for (int i = 0; i < 6; i++) {
        if (i == 4) {
          token = line.substr(pos + 1, line.find(",", pos) - pos - 1); // Remove most left and most right characters
          pos += token.length() + 2; // 2 represents the length of the left and right characters to remove
        }
        else if (i == 5) {
          token = line.substr(pos, line.length() - pos - 1); // Remove most right character
          pos += 4;
        }
        else {
          token = line.substr(pos, line.find(" ", pos) - pos);
          pos += token.length() + 1;
        }
        if (i == 0 || i == 4 || i == 5) {
          row.push_back(token);
        }
      }
      data.push_back(row);
    }
    infile.close();
  }
  
  vector<vector<int>> data_int(data.size(), vector<int>(3));
  for (int i = 0; i < data.size(); i++) {
    for (int j = 0; j < 3; j++) {
      data_int[i][j] = stoi(data[i][j]);
      if(j==2) data_int[i][j]=data_int[i][j]-1;
    }
  }
  return data_int;
}


void removeRowsWithValue(vector<vector<int>>& data, const vector<int>& values) {
  //remove the rows of data which first row is included in values
  auto it = remove_if(data.begin(), data.end(), [&](const vector<int>& row) {
    return find(values.begin(), values.end(), row.front()) != values.end();
  });
  data.erase(it, data.end());
}


hMat matchfiletohMat(const vector<vector<int>>& matchdata, 
                     const int& nref, const int& nsnp){
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

  

double est_rho_EM(hMat& mat, vector<double>& gd,
                  const int ite_time=10){
  //estimate \rho using EM algorithm
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob;
  vector<double> otherprob; 
  double rho_ite=400000/nref;
  vector<double> gl(nsnp-1);
  vector<double> rho_each((nsnp-1));
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
      gl[j]=gd[j+1]-gd[j];
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
    rho_ite=vec_sum(rho_each)/vec_sum(gl);
    cout<<"rho is estimated as "<<rho_ite<<endl;
  }
  //cout<<"rho is estimated as"<<rho_ite<<endl;
  return(rho_ite);
}

// [[Rcpp::export]]
double est_rho_Viterbi(const vector<int>& startpos, const vector<int>& endpos,
                       const int nsnp, const double gdtotal) {
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
      maxend = maxEndMap[j];
      maxkidx = j;
    }
    int maxendnew=maxend;
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
  
  double rho_est = nrec / static_cast<double>(gdtotal);
  cout << "rho is estimated as " << rho_est << endl;
  return rho_est;
}



double est_rho_average(const hAnc& refidx, const int nref, const int nsnp,
                       vector<double>& gd, const double prop=0.1,
                       const int ite_time=20, const string method="Viterbi"){
  // estimate \rho from the reference panel
  int npop=refidx.pos.size();
  vector<double> rho_est;
  
  for(int i=0;i<npop;++i){
    //randomly sample a percentage of prop reference samples
    vector<int> popidx=refidx.findrows(i);
    
    int nref_sample=static_cast<int>(ceil(popidx.size()*prop));
    
    vector<int> samples=randomsample(popidx,nref_sample);
    
    for(int k=0;k<samples.size();k++){
      //leave-one-out
      vector<vector<int>> matchdata=read_data("donor",samples[k]);
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
      cout<<matchdata.size()<<endl;
      
      
      if(method=="Viterbi"){
        vector<int> startpos, endpos;
        for (const auto& row : matchdata) {
          startpos.push_back(row[0]);
          endpos.push_back(row[1]);
        }
        rho_est.push_back(est_rho_Viterbi(startpos,endpos,nsnp,gd[nsnp-1]-gd[0]));
      }else{
        hMat mat=matchfiletohMat(matchdata,nref-npop,nsnp);
        rho_est.push_back(est_rho_EM(mat,gd,ite_time));
      }
    }
  }
  double rho_ave=vec_sum(rho_est)/rho_est.size();
  cout<<"Average estimated rho is "<<rho_ave<<endl;
  return(rho_ave);
}


hMat indpainting(const hMat& mat,vector<double>& gd, const double rho,
                 const hAnc& refidx, const vector<int>& refindex){
  // return individual painting
  int nref=mat.d1;
  int nsnp=mat.d2;
  vector<double> sameprob=cal_sameprob(nsnp,rho,gd);
  vector<double> otherprob=cal_otherprob(nref,sameprob);
  hMat marginal_prob=forwardBackward(mat,sameprob,otherprob);
  int npop=refidx.pos.size();
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


vector<hMat> paintingall(vector<double>& gd,const vector<int>& refindex, 
                         const int nind,const string method="Viterbi", 
                         bool fixrho=true,const int ite_time=10,
                         const double prop=0.1){
  //compute painting for all target individuals
  const int nsnp=gd.size();
  const int nref=refindex.size();
  hAnc refidx(refindex);
  const int npop=refidx.pos.size();
  vector<hMat> painting_all(nind, hMat(npop, nsnp));
  double rho_use;
  if(fixrho){
    cout<<"estimating fixed rho"<<endl;
    rho_use=est_rho_average(refidx,nref,nsnp,gd,prop,ite_time,method);
  }
  for(int ii=0;ii<nind;++ii){
    cout<<"calculating painting for individual "<<ii<<endl;
    vector<vector<int>> targetmatchdata=read_data("target",ii);
    hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
    if(fixrho){
      painting_all[ii]=indpainting(mat,gd,rho_use,refidx,refindex);
    }else{
      vector<int> startpos, endpos;
      for (const auto& row : targetmatchdata) {
        startpos.push_back(row[0]);
        endpos.push_back(row[1]);
      }
      rho_use=est_rho_Viterbi(startpos,endpos,nsnp,gd[nsnp-1]-gd[0]);
      painting_all[ii]=indpainting(mat,gd,rho_use,refidx,refindex);
    }
  }
  return(painting_all);
}


// [[Rcpp::export]]
vector<vector<vector<double>>> paintingalldense(vector<double>& gd,
                                                const vector<int>& refindex,
                                                const int nind,
                                                const string method="Viterbi",
                                                bool fixrho=true,
                                                const int ite_time=10,
                                                const double prop=0.1){
  //compute painting for all target individuals
  const int nsnp=gd.size();
  const int nref=refindex.size();
  hAnc refidx(refindex);
  const int npop=refidx.pos.size();
  cout<<"teststop"<<endl;
  vector<vector<vector<double>>> painting_all(nind, vector<vector<double>>(npop, vector<double>(nsnp)));
  cout<<"teststop2"<<endl;
  double rho_use;
  if(fixrho){
    cout<<"estimating fixed rho"<<endl;
    rho_use=est_rho_average(refidx,nref,nsnp,gd,prop,ite_time,method);
  }
  for(int ii=0;ii<nind;++ii){
    cout<<"calculating painting for individual "<<ii<<endl;
    vector<vector<int>> targetmatchdata=read_data("target",ii);
    hMat mat=matchfiletohMat(targetmatchdata,nref,nsnp);
    if(fixrho){
      hMat pind=indpainting(mat,gd,rho_use,refidx,refindex);
      vector<vector<double>> pind_dense=hMatrix2matrix(pind);
      for(int j=0;j<nref;++j){
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
      rho_use=est_rho_Viterbi(startpos,endpos,nsnp,gd[nsnp-1]-gd[0]);
      hMat pind=indpainting(mat,gd,rho_use,refidx,refindex);
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
double hashmaptest(){
  vector<vector<int>> targetmatchdata=read_data("target",0);
  cout<<targetmatchdata.size()<<endl;
  
  const int length = 1469;
  const double start = 0.0;
  const double end = 0.1468;
  const double increment = (end - start) / (length - 1);
  
  std::vector<double> gd(length);
  
  for (int i = 0; i < length; i++) {
    gd[i] = start + i * increment;
  }
  int nref=20000;
  int nsnp=1469;
  cout<<targetmatchdata.size()<<endl;
  vector<int> removeidx;
  for(int j=0;j<2;++j){
    if(j==0){
      removeidx.push_back(5);
    }else{
      removeidx.push_back(10005);
    }
    //removeidx contains the indices to be removed for leave-one-out
  }

  removeRowsWithValue(targetmatchdata,removeidx);
  cout<<targetmatchdata.size()<<endl;
  //cout<<targetmatchdata.size()<<endl;
  hMat mat=matchfiletohMat(targetmatchdata,nref-2,nsnp);
  double a=est_rho_EM(mat,gd,10);
  return(a);
}

// [[Rcpp::export]]
vector<vector<int>> hashmaptest2(){
  vector<vector<int>> targetmatchdata=read_data("target",0);
  return(targetmatchdata);
}








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