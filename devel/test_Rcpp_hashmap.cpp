// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>

// When we want to lookup all entries:
// https://stackoverflow.com/questions/28767234/what-container-to-store-unique-values
// How to wrap an object for return to R
// https://stackoverflow.com/questions/12405655/returning-a-custom-object-from-a-wrapped-method-in-rcpp

//namespace hMatRcpp {
  
class hVec { // A sparse vector format
public:
  std::vector<int> k; // the keys that are in the vector
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
  hVec(std::vector<int> idx,std::vector<double> val,int len,double x0){
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
    // Safely set a value
    if(!in(p))k.push_back(p);
    setnocheck(p,val);
  };
  bool in(int p){
    // Check if a value has a non-default entry
    if(v.find(p)==v.end()) return(false);
    return(true);
  };
  double get(int p){
    // Get a value from the vector: either its set value or the default if not present
    if(!in(p)) return(x0);
    return(v[p]);
  };
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
      m.push_back(hVec(d1,x0));
    }
  };
  void appendColumn(double x0){
    // Append a default column
    ++d2;
    m.push_back(hVec(d1,x0));
  };
  void appendColumn(std::vector<int> idx,std::vector<double> vals,double x0){
    // Append a filled column
    ++d2;
    m.push_back(hVec(idx,vals,d1,x0));
  };
};


// [[Rcpp::export]]
int hash_test() {

  hMat m(50,10,0.025);

  m.m[4].set(23,25.25);

  Rcpp::Rcout << "Input: " <<4<<","<<23<<" Output: "<< m.m[4].get(23) << std::endl;
  Rcpp::Rcout << "Input: " <<3<<","<<23<<" Output: "<< m.m[3].get(23) << std::endl;
  Rcpp::Rcout << "Input: " <<4<<","<<22<<" Output: "<< m.m[4].get(22) << std::endl;
  return 0;
}

//}

