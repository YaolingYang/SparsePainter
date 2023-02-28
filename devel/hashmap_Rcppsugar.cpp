// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
using namespace Rcpp;
// When we want to lookup all entries:
  // https://stackoverflow.com/questions/28767234/what-container-to-store-unique-values
// How to wrap an object for return to R
// https://stackoverflow.com/questions/12405655/returning-a-custom-object-from-a-wrapped-method-in-rcpp

//namespace hMatRcpp {
  
  class hVec { // A sparse vector format
    public:
      IntegerVector k; // the keys that are in the vector
    std::unordered_map<int, double> v; // the values, stored as a map from keys to values
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
    hVec(IntegerVector idx,NumericVector val,int len,double x0){
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
    void setall(IntegerVector p, NumericVector val){
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
    
    NumericVector getall(IntegerVector idx){
      // Get values from the vector: either its set value or the default if not present

      NumericVector values(idx.size());
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
      std::vector<hVec> m; // sparse matrix, i.e. a vector of hVec's
  int d1; // number of rows; currently nominal
  int d2; // Number of columns; should be equal to length(m)
  hMat(int d1){
    // Empty matrix with d1 rows (can append columns)
    this->d1=d1;
    d2=0;
  };
  hMat(int d1,int d2,double x0){
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
  void appendColumn(IntegerVector idx,NumericVector vals,double x0){
    // Append a filled column
    ++d2;
    m.push_back(hVec(idx,vals,d1,x0));
  };

};
  
  //double vec_sum(NumericVector l){
    // calculate the sum of a vector
   // double sum_l=0;
    //for(int i=0;i<l.size();++i){
    //  sum_l=sum_l+l[i];
   // }
   // return(sum_l);
  //}
  
  //NumericVector vec_multiply(NumericVector l1,NumericVector l2){
    // calculate the multiply of two vectors
   // NumericVector multi_l(l1.size());
   // for(int i=0;i<l1.size();++i){
    //  multi_l[i]=l1[i]*l2[i];
    //}
   // return(multi_l);
  //}
  
  //IntegerVector union_vec(const IntegerVector& vec1, const IntegerVector& vec2) {
    // get the union of two vectors
   // IntegerVector result(vec1.size() + vec2.size());
    //IntegerVector::iterator it;
    //it = std::set_union(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), result.begin());
    //result.resize(it - result.begin());
   // return result;
  //}

  //NumericVector update_vec(NumericVector a, IntegerVector b, NumericVector c) {
    // change the position of b inside a to c
    //for (int i = 0; i < b.size(); ++i) {
    //  a[b[i]] = c[i];
   // }
   // return(a);
  //}
  
  double hVecSum(hVec v,IntegerVector idx={}){
    // calculate the sum of the values of a hVec with given indices
    if (idx.size() == 0) {
      idx = v.k;
    }
    return(sum(v.getall(idx)));
  }

  hVec hVecCirc(hVec v1, hVec v2){
    // calculate the circ of two hVec;
    hVec ret=hVec();
    ret.setdefault(v1.x0*v2.x0);
    IntegerVector v1idx=v1.k;
    IntegerVector v2idx=v2.k;
    IntegerVector idx=union_(v1idx,v2idx);
    NumericVector v1val=v1.getall(idx);
    NumericVector v2val=v2.getall(idx);
    ret.setall(idx,v1val*v2val);
    return(ret);
  }

  hVec hVecScale(hVec v,double x){
    //Scale a hash vector by x
    IntegerVector idx=v.k;
    NumericVector val=v.getall(idx);
    v.setdefault(x*v.x0);
    NumericVector xvec(idx.size(),x);
    v.setall(idx,val*xvec);
    return(v);
  }
  
  
  NumericVector hVecdense(hVec x,int nrow){
    //Convert a hash vector x into a dense vector with length nrow
    NumericVector vdense(nrow, x.x0);
    IntegerVector k=x.k;
    NumericVector vval=x.getall(k);
    //vdense=update_vec(vdense,k,vval);
    vdense[k]=vval;
    return(vdense);
  }
  
  
 hMat matchmat2hMatrix(LogicalMatrix matchmat,double default_val=0.0){
    //Convert a logical matrix into a hash matrix
    int nrow=matchmat.nrow();
    int ncol=matchmat.ncol();
    hMat mat(nrow,ncol,default_val);

    for(int j=0;j<ncol;++j){
      for(int i=0;i<nrow;++i){
       if(matchmat(i,j)){
         mat.m[j].set(i,1.0);
       }
      }
    }
    return(mat);
  }
  
  
  NumericMatrix hMatrix2matrix(hMat mat){
    int nrow=mat.d1;
    int ncol=mat.d2;
    NumericMatrix densematrix(nrow,ncol);
    NumericVector vecdence(nrow);
    for(int j=0;j<ncol;++j){
      vecdence=hVecdense(mat.m[j],nrow);
      densematrix(_,j)=vecdence;
    }
    return(densematrix);
  }
  
  
  hMat forwardProb(hMat mat, NumericVector sameprob, NumericVector otherprob){
    //compute normalised forward probability and store in hMat
    int nrow=mat.d1;
    int ncol=mat.d2;
    hMat forward_prob(nrow,ncol,0);
    IntegerVector twj=mat.m[0].k;
    NumericVector l(twj.size(),1.0/twj.size());
    forward_prob.m[0].setall(twj,l);
    for(int j=1;j<ncol;++j){
      twj=mat.m[j].k;
      NumericVector fp=sameprob[j-1]*forward_prob.m[j-1].getall(twj)+otherprob[j-1];
      forward_prob.m[j].setall(twj,fp/sum(fp));
      //NumericVector fprev=forward_prob.m[j-1].getall(twj);
      //for(int i=0;i<twj.size();++i){
        //forward_prob.m[j].set(twj[i],sameprob*fprev[i]+otherprob);
      //}
      //forward_prob.m[j]=hVecScale(forward_prob.m[j],1.0/hVecSum(forward_prob.m[j],twj));
    }
    return(forward_prob);
  }
  
  
  hMat backwardProb(hMat mat, NumericVector sameprob, NumericVector otherprob){
    //compute normalised backward probability and store in hMat
    int nrow=mat.d1;
    int ncol=mat.d2;
    hMat backward_prob(nrow,ncol,1.0/nrow);
    IntegerVector twj;
    NumericVector Bjp1;
    double sumBjp1;
    for(int j=ncol-2;j>=0;--j){
      twj=mat.m[j+1].k;
      Bjp1=backward_prob.m[j+1].getall(twj);
      sumBjp1=sum(Bjp1);
      NumericVector val=otherprob[j]*sumBjp1 + sameprob[j] * Bjp1;
      backward_prob.m[j].setall(twj,val/sumBjp1);
      //for(int i=0;i<twj.size();++i){
      //  backward_prob.m[j].set(twj[i],sameprob*Bjp1[i]+otherprob*sumBjp1);
      //}
      //backward_prob.m[j]=hVecScale(backward_prob.m[j],1.0/sumBjp1);
      backward_prob.m[j].setdefault(otherprob[j]);
    }
    return(backward_prob);
  }

  hMat marginalProb(hMat f, hMat b){
    //compute normalised backward probability and store in hMat
    int nrow=f.d1;
    int ncol=f.d2;
    hMat marginal_prob=hMat(nrow,ncol,1.0/nrow);
    for(int j=0;j<ncol;++j){
      marginal_prob.m[j]=hVecCirc(f.m[j],b.m[j]);
      marginal_prob.m[j]=hVecScale(marginal_prob.m[j],1.0/hVecSum(marginal_prob.m[j]));
    }
    
    return(marginal_prob);
  }
  
  
  
  // [[Rcpp::export]]
  NumericMatrix forwardBackward5(LogicalMatrix matchmat,
                                NumericVector sameprob, NumericVector otherprob){
    hMat mat = matchmat2hMatrix(matchmat);
    hMat forward_prob=forwardProb(mat,sameprob,otherprob);
    hMat backward_prob=backwardProb(mat,sameprob,otherprob);
    hMat marginal_prob=marginalProb(forward_prob,backward_prob);
    return(hMatrix2matrix(marginal_prob));
  }

  
  
  //}