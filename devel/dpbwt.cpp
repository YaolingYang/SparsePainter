#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<random>

using namespace std;

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


tuple<vector<int>,vector<int>,vector<int>,vector<int>> longMatchdpbwt(int L, 
                                                                      int *numMatches,
                                                                      bool **panel, 
                                                                      dpbwt & x){
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
      
      // store the output
      
      dpbwtnode *f, *g, *ftemp, *gtemp;
      int matches = 0;
      f = g = z[0].below;
      for(int k = 0; k<N; ++k){
        ftemp = (panel[Oid][k])? f->u : f->v;
        gtemp = (panel[Oid][k])? g->u : g->v;
        f = (panel[Oid][k])? f->v : f->u;
        g = (panel[Oid][k])? g->v : g->u;
        while (ftemp != gtemp){
          //matchOut << ftemp->originalid << " = q" << i-M+qM << " at [" << dZ[ftemp->id] << ", " << k << ")\n";
          donorid.push_back(ftemp->originalid);
          startpos.push_back(dZ[ftemp->id]);
          endpos.push_back(k-1);
          ++matches;
          ftemp = ftemp->below;
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
        donorid.push_back(f->originalid);
        startpos.push_back(dZ[f->id]);
        endpos.push_back(N-1);
        ++matches;
        f = f->below;
      }
      numMatches[i-M+qM] = matches;
      delete [] z;
      
      //recode the position of the next start position of query haplotype
      queryidall.push_back(startpos.size());
  }
  tuple<vector<int>,vector<int>,vector<int>,vector<int>> results(queryidall,donorid,startpos,endpos);
  return(results);
}

tuple<vector<int>,vector<int>,vector<int>,vector<int>> do_dpbwt(int L=1000, 
                                                                string query="target"){
  string qin = "p_" + query + ".vcf";
  string in = "p_donor.vcf";
  
  bool **panel;
  dpbwt x;
  int *numMatches;
  ReadVCF(in, qin, panel);
  
  numMatches = new int[qM];
  
  return(longMatchdpbwt(L,numMatches, panel, x));
  
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