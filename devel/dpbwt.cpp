#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<random>
#include <algorithm>
#include <utility>

using namespace std;

int L=1000;
int M = 0;
int qM = 0;
int N = 0;
int *dZ;

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