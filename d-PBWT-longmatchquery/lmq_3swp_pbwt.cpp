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

void OutputPanel(const char* outFile, int **panel, int num){
    std::ofstream out(outFile);
    for (int i = 0; i<num; ++i){
        for (int j = 1; j<N+1; ++j)
            out << std::setw(7) << (int)panel[i][j];
        out << std::endl;
    }
    out.close();
}

void OutputPaneluv(const char* outFile, int **panel, int **prefix, int num){
    std::ofstream out(outFile);
    for (int i = 0; i<num; ++i){
        for (int j = 1; j<N; ++j)
            if (panel[i][j]>=num)
                out << std::setw(7) << -1;
            else 
                out << std::setw(7) << (int)prefix[panel[i][j]][j+1];
        out << std::endl;
    }
    out.close();
}

void PBWT(bool **panel, int **prefix, int **divergence,
        int **u, int **v, int order [], int num){
    for (int i = 0; i<num; ++i){
        prefix[i][0] = order[i];
        divergence[i][0] = 0;
    }
    for (int k = 0; k<N; ++k){
        int u2 = 0, v2 = 0, p = k+1, q = k+1;
        std::vector<int> a,b,d,e;
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

void longMatch(double *time[], int *numMatches,
        bool **panel, int **prefix, int **divergence, int **u, 
        int **v, int*order, std::ofstream & matchOut){
    dZ = new int[M];
    for (int i = 0; i<M; ++i){
        dZ[i] = 0;
    }

    int *t = new int[N+1];
    int *zd = new int[N+2], *bd = new int[N+2];
    double oldtime, oldcputime;
    for (int s = 0; s<qM; ++s){
        oldtime = get_wall_time();
        oldcputime = get_cpu_time();
        int Oid = order[M-qM+s];
        t[0] = 0;
        for (int k=0; k<N; ++k)
            if (t[k]!=M-qM)
                if (!panel[Oid][k])
                    t[k+1] = u[t[k]][k];
                else
                    t[k+1] = v[t[k]][k];
            else
                if (!panel[Oid][k])
                    t[k+1] = v[0][k];
                else
                    t[k+1] = M-qM;

        zd[N+1] = bd[N+1] = N;

        for (int k = N; k>=0; --k){
            zd[k] = std::min(zd[k+1],k);
            bd[k] = std::min(bd[k+1],k);
            if (t[k]!=0)
                while(zd[k]>0 && 
                        panel[Oid][zd[k]-1]
                        == panel[prefix[t[k]-1][k]][zd[k]-1])
                    zd[k]--;
            else
                zd[k] = k;
            if (t[k]!=M-qM)
                while(bd[k]>0 && 
                        panel[Oid][bd[k]-1] 
                        == panel[prefix[t[k]][k]][bd[k]-1])
                    bd[k]--;
            else
                bd[k] = k;
        }            
        int f, g, ftemp, gtemp, matches = 0;
        f = g = t[0];
        for (int k = 0; k<N; ++k){
            if (g == M-qM){
                if (f == M-qM){
                    if (!panel[Oid][k]){
                        ftemp = M-qM;
                        f = v[0][k];
                    }
                    else{
                        ftemp = v[0][k];
                        f = M-qM;
                    }
                }
                else{
                    if (!panel[Oid][k]){
                        ftemp = v[f][k];
                        f = u[f][k];
                    }
                    else{
                        ftemp = u[f][k];
                        f = v[f][k];
                    }
                }
                if (!panel[Oid][k]){
                    gtemp = M-qM;
                    g = v[0][k];
                }
                else{
                    gtemp = v[0][k];
                    g = M-qM;
                }
            }
            else
                if (!panel[Oid][k]){
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
                //output match
                matchOut << prefix[ftemp][k+1] << " = q" << s << " at [" << dZ[prefix[ftemp][k+1]] << ", " << k << ")\n";  
                ++matches;
                ++ftemp;
            }
            if (f==g){
                if (k+1-zd[k+1] == L){
                    --f;
                    dZ[prefix[f][k+1]] = k+1-L;
                    //store divergence
                }
                if (k+1-bd[k+1] == L){
                    //store divergence
                    dZ[prefix[g][k+1]] = k+1-L;
                    ++g;
                }
            }
            if (f!=g){
                while (divergence[f][k+1] <= k+1-L){
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
        //TODO output match not even ending at N
        while (f != g){
            //output match
            matchOut << prefix[f][N] << " = q" << s << " at [" << dZ[prefix[f][N]] << ", " << N << ")\n";  
            ++matches;
            ++f;
        }

        numMatches[s] = matches;
        time[1][s] = get_cpu_time() - oldcputime;
        time[0][s] = get_wall_time() - oldtime;
    }
    delete [] t;
    delete [] zd;
    delete [] bd;
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
    std::string in = "panel.vcf", out = "3swp_pbwt_matches.txt", qin = "query.vcf", tout = "3swp_pbwt_time.txt", orderFile = "order.txt";
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
        std::cout << "This program takes VCF files as input. It ignores the multi allelic variations. Then, it takes a query file, this file should have the same number of sites. It outputs matches between query haplotypes and haplotypes in the panel length L or longer. It runs on the PBWT using the triple sweep long match query algorithm presented in https://www.biorxiv.org/content/10.1101/2020.01.14.906487v1.\n\noptions:\n-i, input VCF file for the panel, default is panel.vcf\n-m, use separate input files for panel and query haplotypes (default off), if not passed, -n must be specified\n-n, number of query haplotypes (only if query haplotypes are in the input vcf, the last n haplotypes in the vcf will be used)\n-o, output file for matches, default is matches.txt\n-t, output file for time, default is longMatchTime.txt\n-L, length of match to output. Matches length L or longer will be outputted, default is 1000\n-q, input VCF for query haplotypes to search against the panel, default is query.vcf\n-r, read and input order from file provided.\n-g, generate and use a randomized order, a file name must be provided with this option\n-h, help\n-z, run with no parameters\nno parameters, help" << std::endl;
        return 0;
    }
    std::ofstream matchOut(out);
    //even is real, odd is cpu
    double *time[2];
    bool **panel;
    int **prefix, **divergence, **u, **v;
    int *numMatches;
    if (multin)
        ReadVCF(in, qin, panel);
    else
        ReadVCF(in, panel);
    std::cout << "M: " << M << "\nN: " << N << "\nqM: " << qM << std::endl;
    if (debugging)
        OutputPanel("3swp_pbwt_Panel.txt", panel, M);
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
    
    PBWT(panel, prefix, divergence, u, v, order, M-qM);
    if (debugging) {
        OutputPanel("3swp_pbwt_Prefixfull.txt", prefix, M-qM);
        OutputPanel("3swp_pbwt_Divergencefull.txt", divergence, M-qM);
        OutputPaneluv("3swp_pbwt_uIDfull.txt", u, prefix, M-qM);
        OutputPaneluv("3swp_pbwt_vIDfull.txt", v, prefix, M-qM);
    }
    
    longMatch(time, numMatches, panel, prefix, divergence, u, v, order, matchOut);
    std::ofstream timeout(tout);
    timeout<< "Real(s)\tCPU(s)\tmatches\n";
    for (int s = 0; s<qM; ++s)
        timeout << time[0][s] << '\t' << time[1][s]
            <<  '\t' << numMatches[s] << '\n';
    timeout.close();
    matchOut.close();
    return 0;
}

