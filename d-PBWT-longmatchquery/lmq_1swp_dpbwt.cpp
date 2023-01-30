#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<random>

int M = 0;
int qM = 0;
int N = 0;
int L = 1000;
int *dZ;


//-----------------------------------------------------------------------------
//source for time: https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
                    ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif
//-----------------------------------------------------------------------------

struct dpbwtnode{
    dpbwtnode *below, *above, *u, *v;
    int divergence, id, originalid;
};

struct dpbwt{
    dpbwtnode bottomfirstcol;
    dpbwtnode *topfirstcol;
    std::vector<dpbwtnode*> firstcol;
    std::vector<std::vector<bool>> panel;
    int size, count;
};

void ReadVCF(std::string inFile, std::string qinFile, bool ** & panel){
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
            std::cout << "Query file and input file have different numbers of sites. Query has " << qN << ". Panel has " << N << std::endl;
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

void ReadVCF(std::string inFile, bool ** & panel){
    using namespace std;
    {//count M and N
        string line = "1#";
        ifstream in(inFile);
        while (line[1] == '#')
            getline(in, line);
        stringstream linestr;

        linestr.str(line);
        linestr.clear();
        for (int i = 0; i<9; i++)
            linestr >> line;
        N = M = 0;
        while (!linestr.eof()){
            ++M;
            ++M;
            linestr >> line;
        }
        while (getline(in,line)){
            ++N;
        }
        in.close();
    }
    bool *temp = new bool[(long long)M * N];
    panel = new bool*[M];
    for (long long i = 0; i<M; i++){
        panel[i] = &(temp[i*N]);
    }
    string line = "##";
    ifstream in(inFile);
    stringstream linestr;
    int x = 0;
    char y = 0;

    while (line[1] == '#')
        getline(in, line);
    for (int j = 0; j<N; ++j){
        getline(in, line);
        linestr.str(line);
        linestr.clear();
        for (int i = 0; i<9; ++i){
            linestr >> line; 
        }
        for (int i = 0; i<M/2; ++i){
            linestr >> x >> y;
            panel[i*2][j] = (bool)x;
            linestr >> x;
            panel[i*2+1][j] = (bool)x;
        
        }
    }
    in.close();
}

void OutputPanel(const char* outFile, bool** panel, int num){
    std::ofstream out(outFile);
    for (int i = 0; i<num; ++i){
        for (int j = 0; j<N; ++j)
            out << std::setw(3) << (int)panel[i][j];
        out << std::endl;
    }
    out.close();
}

void OutputAll(const char * outFilep, const char * outFiled, const char * outFileu,
        const char * outFilev, dpbwt x){
    std::ofstream outp(outFilep);
    std::ofstream outd(outFiled);
    std::ofstream outu(outFileu);
    std::ofstream outv(outFilev);

    dpbwtnode *top[N+1];
    top[0] = &x.bottomfirstcol;
    while (top[0]->above!=nullptr)
        top[0] = top[0]->above;
    for (int j = 0; j<N; ++j)
        top[j+1] = top[j]->u;
    for (int i=0; i<M-qM; ++i){
        for (int j=1; j<N+1; ++j){
            outp << std::setw(7) << (top[j]->originalid);
            outd << std::setw(7) << top[j]->divergence;
            if (j!=N){
                if (top[j]->u->id != -1)
                    outu << std::setw(7) << top[j]->u->originalid;
                else
                    outu << std::setw(7) << -1;
                if (top[j]->v->id != -1)
                    outv << std::setw(7) << top[j]->v->originalid;
                else
                    outv << std::setw(7) << -1;
            }
            
            top[j] = top[j]->below;
        }
        outp << std::endl;
        outd << std::endl;
        outu << std::endl;
        outv << std::endl;
    }
    outp.close();
    outd.close();
    outu.close();
    outv.close();
}

void singleSweepLongMatch(double *time[], int *numMatches,
        bool **panel, dpbwt & x, int *order, std::ofstream & matchOut, bool debugging){
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
        int Oid = order[i];
        x.panel.push_back(std::vector<bool>());
        x.panel[i].resize(N);
        dpbwtnode *t = &(x.bottomfirstcol),
                  *botk = &(x.bottomfirstcol),
                  *z = new dpbwtnode[N+1];
        x.firstcol.push_back(&z[0]);

        if (i == 0)
            x.topfirstcol = &z[0];

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
            zdtemp = std::min(zdtemp, k);
            bdtemp = std::min(bdtemp, k);
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

    if (debugging)
        OutputAll("1swp_dpbwt_Prefixfull.txt", "1swp_dpbwt_Divergencefull.txt", 
                "1swp_dpbwt_uIDfull.txt", "1swp_dpbwt_vIDfull.txt", x);
    double oldtime = 0, oldcputime = 0;

    for (int i = M-qM; i<M; i++){
        oldtime = get_wall_time();
        oldcputime = get_cpu_time();
        int Oid = order[i];
        dpbwtnode *t = &(x.bottomfirstcol),
                  *botk = &(x.bottomfirstcol),
                  *topk = x.topfirstcol,
                  *z = new dpbwtnode[N+1];

        z[0].id = i;
        z[0].originalid = Oid; //!!!! different from indel benchmark
        z[0].divergence = 0;
        z[0].below = t;
        z[0].above = t->above;

        dpbwtnode *ef, *eg, *f, *g, *ftemp, *gtemp;
        int e = 0;

        ef = eg = f = g = z[0].below;

        int matches = 0;
        for (int k = 0; k<N; ++k){
            t = (panel[Oid][k])? t->v : t->u;
            z[k+1].id = i;
            z[k+1].originalid = Oid;
            z[k+1].below = t;
            z[k+1].above = t->above;

            if (e == k){
                ef = topk;
                eg = botk;
            }
            ef = (panel[Oid][k])? ef->v : ef->u;
            eg = (panel[Oid][k])? eg->v : eg->u;
            if (ef == eg){
                if (ef == topk->u){
                    e = k+1;
                    while (e > 0 && panel[Oid][e-1] == x.panel[eg->id][e-1])
                        e--;
                    if (e != k+1){
                        eg = eg->below;
                        while (eg->divergence <= e)
                            eg = eg->below;
                    }
                }
                else if (ef == botk->v){
                    e = k+1;
                    while (e>0 && panel[Oid][e-1] == x.panel[ef->above->id][e-1])
                        e--;
                    if (e!=k+1){
                        ef = ef->above;
                        while (ef->divergence <= e)
                            ef = ef->above;
                    }
                }
                else{
                    e = ef->divergence - 1;
                    if (e == -1) {
                        e = 0;
                        while (ef != topk->u && ef->divergence <= e)
                            ef = ef->above;
                        while (eg != botk->v && eg->divergence <= e)
                            eg = eg->below;
                    }
                    else {
                        if (!panel[Oid][e]){
                            ef = ef->above;
                            while (e > 0 && panel[Oid][e-1] == x.panel[ef->id][e-1])
                                e--;
                            while (ef != topk->u && ef->divergence <=e)
                                ef = ef->above;
                        }
                        else{
                            while (e > 0 && panel[Oid][e-1] == x.panel[eg->id][e-1])
                                e--;
                            eg = eg->below;
                            while (eg != botk && eg->divergence <= e)
                                eg = eg->below;
                        }
                    }

                }
            }


            ftemp = (panel[Oid][k])? f->u : f->v;
            gtemp = (panel[Oid][k])? g->u : g->v;
            f = (panel[Oid][k])? f->v : f->u;
            g = (panel[Oid][k])? g->v : g->u;
            while (ftemp != gtemp){
                matchOut << ftemp->originalid << " = q" << i-M+qM << " at [" << dZ[ftemp->id] << ", " << k << ")\n";
                ++matches;
                ftemp = ftemp->below;
            }
            if (f==g){
                if (k+1-e == L){
                    if (f == topk->u){
                        dZ[g->id] = k+1 - L;
                        g = g->below;
                    }
                    else if (f == botk->v){
                        f = f->above;
                        dZ[f->id] = k+1 - L;
                    }
                    else{
                        int temp = z[k+1].below->divergence-1;
                        if (temp == -1){
                            temp = 0;
                        }
                        if (!panel[Oid][temp]){
                            f = f->above;
                            dZ[f->id] = k+1 - L;
                        }
                        else {
                            dZ[g->id] = k+1 - L;
                            g = g->below;
                        }
                    }
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
            topk = topk->u;
            botk = botk->v;
        }
        while (f != g){
            //output match
            matchOut << f->originalid << " = q" << i-M+qM << " at [" << dZ[f->id] << ", " << N << ")\n";
            ++matches;
            f = f->below;
        }
        numMatches[i-M+qM] = matches;
        time[1][i-M+qM] = get_cpu_time() - oldcputime;
        time[0][i-M+qM] = get_wall_time() - oldtime;
        delete [] z;
    }
}

void Randomize(int *order, int num){
    std::random_device seed;
    std::default_random_engine generator(seed());
    for (int i = num-1; i>0; --i){
        std::uniform_int_distribution<int> dis(0, i);
        std::swap(order[i], order[dis(generator)]);
    }
}

int main(int argc, char *argv[]){
    std::string in = "panel.vcf", out = "1swp_dpbwt_matches.txt", qin = "query.vcf", tout = "1swp_dpbwt_time.txt", orderFile = "order.txt";
    bool multin = false, readOrder = false, randomWriteOrder = false;
    bool help = false;
    bool debugging = false;
    for (int i = 1; i<argc; ++i){
        switch (argv[i][1]){
            case 'm':
            case 'M':
                multin = true;
                break;
            case 'n':
            case 'N':
                qM = atoi(argv[++i]);
                break;
            case 'i':
            case 'I':
                in = argv[++i];
                break;
            case 'o':
            case 'O':
                out = argv[++i];
                break;
            case 'q':
            case 'Q':
                qin = argv[++i];
                break;
            case 't':
            case 'T':
                tout = argv[++i];
                break;
            case 'l':
            case 'L':
                L = atoi(argv[++i]);
                break;
            case 'd':
            case 'D':
                debugging = true;
                break;
            case 'r':
            case 'R':
                if (randomWriteOrder){
                    std::cout << "Cannot read and generate order at the same time!" << std::endl;
                    throw "Cannot read and generate order at the same time!";
                }
                readOrder = true;
                orderFile = argv[++i];
                break;
            case 'g':
            case 'G':
                if (readOrder){
                    std::cout << "Cannot read and generate order at the same time!" << std::endl;
                    throw "Cannot read and generate order at the same time!";
                }
                randomWriteOrder = true;
                orderFile = argv[++i];
                break;
            case 'h':
            case 'H':
                help = true;
        }
    }
    if (argc < 2 || help){
        std::cout << "This program takes VCF files as input. It outputs matches between query haplotypes and haplotypes in the panel length L or longer. It is recommended for multiallelic sites to be removed from the input VCF. If multiallelic sites are not removed, the code will take 0 as 0 and any other allele as 1. The algorithms used are presented in https://www.biorxiv.org/content/10.1101/2020.01.14.906487v1.\n\noptions:\n-i, input VCF file for the panel, default is panel.vcf\n-m, use separate input files for panel and query haplotypes (default off), if not passed, -n must be specified\n-n, number of query haplotypes (only if query haplotypes are in the input vcf, the last n haplotypes in the vcf will be used), default is 200\n-o, output file for matches, default is matches.txt\n-t, output file for time, default is longMatchTime.txt\n-L, length of match to output. Matches length L or longer will be outputted, default is 1000\n-q, input VCF for query haplotypes to search against the panel, default is query.vcf\n-r, read and input order from file provided.\n-g, generate and use a randomized order, a file name must be provided with this option to write the order used \n-d, use for debugging, this outputs the haplotype, prefix, divergence, u, and v panels. This is meant to only be used with small panels\n-h, help\n-z, run with no parameters\nno parameters, help" << std::endl;
        return 0;
    }
    std::ofstream matchOut(out);
    //even is real, odd is cpu
    double *time[2];
    bool **panel;
    dpbwt x;
    int *numMatches;
    if (multin)
        ReadVCF(in, qin, panel);
    else
        ReadVCF(in, panel);
    std::cout << "M: " << M << "\nN: " << N << "\nqM: " << qM << std::endl;
    if (debugging)
        OutputPanel("1swp_dpbwt_Panel.txt", panel, M);
    time[0] = new double[qM];
    time[1] = new double[qM];
    numMatches = new int[qM];


    int *order = new int[M];
    if (readOrder){
        std::ifstream orderIn(orderFile);
        char c;
        orderIn >> c;
        for (int i = 0; i<M; i++)
            orderIn >> order[i] >> c;
        orderIn.close();
    }
    else {
        for (int i = 0; i<M; ++i){
            order[i] = i;
        }
        if (randomWriteOrder){
            std::ofstream orderOut(orderFile);
            Randomize(order, M);
            orderOut << "{" << order[0];
            for (int i = 1; i<M; ++i)
                orderOut << "," << order[i];
            orderOut << "}";
            orderOut.close();
        }
    }
    
    singleSweepLongMatch(time, numMatches, panel, x, order, matchOut, debugging);

    std::ofstream timeout(tout);
    timeout<< "Real(s)\tCPU(s)\tmatches\n";
    for (int s = 0; s<qM; ++s)
        timeout << time[0][s] << '\t' << time[1][s]
            <<  '\t' << numMatches[s] << '\n';
    timeout.close();
    matchOut.close();
    return 0;
}

