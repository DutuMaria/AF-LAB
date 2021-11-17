#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stack>
#include <bitset>
#include <queue>
#include <numeric>
#include <limits>

using namespace std;

const int NMAX = 200001;
const int INF = numeric_limits<int>::max();

//ifstream fin ("apm.in");
//ofstream fout ("apm.out");

ifstream fin("disjoint.in");
ofstream fout("disjoint.out");

auto cmp = [](const pair<int, int>& p1, const pair<int, int>& p2)
{
    return p1.second > p2.second;
};

class DisjointSet{
private:
    int nrMultimi;
    vector<int> tata;
    vector<int> h;
public:
    DisjointSet(int nrMultimi) : nrMultimi(nrMultimi){
        tata.resize(nrMultimi + 1, -1);
        h.resize(nrMultimi + 1, 0);
    }
    void initializare(const int &nr);
    int gasesteReprez(const int &nr);
    void reuneste(const int &x, const int &y);
    bool query(const int &x, const int &y);

};

void DisjointSet::initializare(const int &nr) {
    tata[nr] = 0;
    h[nr] = 0;
}

int DisjointSet::gasesteReprez(const int &nr) {
    if (tata[nr] == 0)
        return nr;
    tata[nr] = gasesteReprez(tata[nr]);
    return tata[nr];
}

void DisjointSet::reuneste(const int &x, const int &y) {
    int repX = gasesteReprez(x), repY = gasesteReprez(y);
    if(h[repX] > h[repY]){
        tata[repY] = repX;
    }
    else{
        tata[repX] = repY;
        if(h[repX] == h[repY])
            h[repY]++;
    }
}

bool DisjointSet::query(const int &x, const int &y) {
    return gasesteReprez(x) == gasesteReprez(y);
}

class Graf{
private:
    int nrNoduri, nrMuchii;
    vector<int> listaAd[NMAX];
    vector<pair<int, int>> listaAdsiCosturi[NMAX];
    vector<tuple<int, int, int>> muchiiCost;
    bitset<NMAX> viz;
    vector<int> distante;
    stack<int> noduriSortTop;
    vector<vector<int>> compBiconexe;
    vector<vector<int>> ctc;
    void dfs(int nod);
    void bfs(int start);
    void dfsSortTop(int nod);
    void dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s);
    void dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, bitset<NMAX> &inStiva);
    void Dijkstra(bitset<NMAX> &viz, vector<int> &distante, priority_queue<pair<int, int>, vector<pair<int,int>>, decltype(cmp)> &minHeap);
    void BellmanFord(queue<int> &coadaNoduri, vector<int> &dist, vector<int> &nrRelaxari, bitset<NMAX> &inCoada, bool &rez);
    int Apm(vector<pair<int, int>> &sol);
public:
    Graf(int nrNoduri, int nrMuchii) : nrNoduri(nrNoduri), nrMuchii(nrMuchii) {};
    void adaugaInListaAdOrientat(int nod1, int nod2);
    void adaugaInListaAdNeorientat(int nod1, int nod2);
    void adaugaMuchie(const int &nod1, const int &nod2, const int &cost);
    int nrCmpConexe();
    void afisDistanteMinBfs(int start);
    void sortareTopologica();
    vector<vector<int>> componenteBiconexe();
    vector<vector<int>> componenteTareConexe();
    void adaugaInListaAdsiCosturiOrientat(const int &nod1, const int &nod2, const int &cost);
    void afisareDijkstra();
    void afisareBellmanford();
    void afisareApm();
};

void Graf::adaugaInListaAdOrientat(int nod1, int nod2) {
    listaAd[nod1].push_back(nod2);
}

void Graf::adaugaInListaAdNeorientat(int nod1, int nod2) {
    listaAd[nod1].push_back(nod2);
    listaAd[nod2].push_back(nod1);
}

void Graf::dfs(int nod) {
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            viz[listaAd[nod][i]] = true;
            dfs(listaAd[nod][i]);
        }
    }
}

int Graf::nrCmpConexe() {
    int nr = 0;
    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i]){
            nr++;
            viz[i] = true;
            dfs(i);
        }
    }
    return nr;
}

void Graf::bfs(int start) {
    distante.resize(nrNoduri + 1, -1);
    queue<int> coada;
    coada.push(start);
    distante[start] = 0;
    while(!coada.empty()){
        int nod = coada.front();
        for (int i = 0;  i < listaAd[nod].size(); i ++){
            if (distante[listaAd[nod][i]] == -1){
                coada.push(listaAd[nod][i]);
                distante[listaAd[nod][i]] = distante[nod] + 1;
            }
        }
        coada.pop();
    }
}

void Graf::afisDistanteMinBfs(int nodStart) {
    bfs(nodStart);
    for (int i = 1; i <= nrNoduri; i ++){
        fout << distante[i] <<" ";
    }
}

void Graf::sortareTopologica() {
    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i]) {
            dfsSortTop(i);
        }
    }
    while (!noduriSortTop.empty()){
        fout << noduriSortTop.top() << " ";
        noduriSortTop.pop();
    }
}

void Graf::dfsSortTop(int nod) {
    viz[nod] = true;
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            dfsSortTop(listaAd[nod][i]);
        }
    }
    noduriSortTop.push(nod);
}

vector<vector<int>> Graf :: componenteBiconexe() {
    vector<int> nivel;
    vector<int> nivelMin;
    nivel.resize(NMAX);
    nivelMin.resize(NMAX);
    stack<int> s;

    dfsCompBiconexe(1,  1, nivel, nivelMin, s);
    return compBiconexe;
}

void Graf::dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s) {
    viz[nod] = true;
    s.push(nod);
    nivel[nod] = nivelCrt;
    nivelMin[nod] = nivelCrt;

    for (int i = 0; i < listaAd[nod].size(); i++){
        int vecin = listaAd[nod][i];
        if (!viz[vecin]){
            dfsCompBiconexe(vecin, nivelCrt + 1, nivel, nivelMin, s);
            nivelMin[nod] = min(nivelMin[nod], nivelMin[vecin]);
            if (nivelMin[vecin] >= nivel[nod]){
                vector<int> compCrt;
                int nodCompCrt;
                do{
                    nodCompCrt = s.top();
                    compCrt.push_back(nodCompCrt);
                    s.pop();
                } while (nodCompCrt != vecin);
                compCrt.push_back(nod);
                compBiconexe.push_back(compCrt);
            }
        }
        else{ // daca e vizitat => muchie de intoarcere => actualizam nivelMin
            nivelMin[nod] = min(nivelMin[nod], nivel[vecin]);
        }
    }
}

vector<vector<int>> Graf :: componenteTareConexe() {
    vector<int> nivel;
    vector<int> nivelMin;
    nivel.resize(NMAX);
    nivelMin.resize(NMAX);
    stack<int> s;
    bitset<NMAX> inStiva;
    int niv = 1;
    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i])
            dfsCompTareConexe(i,  niv, nivel, nivelMin, s, inStiva);
    }
    return ctc;
}

void Graf::dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, bitset<NMAX> &inStiva) {
    viz[nod] = true;
    s.push(nod);
    inStiva[nod] = true;
    nivel[nod] = nivelCrt;
    nivelMin[nod] = nivelCrt;
    nivelCrt++;

    for (int i = 0; i < listaAd[nod].size(); i++) {
        int vecin = listaAd[nod][i];
        if (!viz[vecin]) {
            dfsCompTareConexe(vecin, nivelCrt, nivel, nivelMin, s, inStiva);
            nivelMin[nod] = min(nivelMin[nod], nivelMin[vecin]);
        } else {
            if (inStiva[vecin])
                nivelMin[nod] = min(nivelMin[nod], nivel[vecin]);
        }
    }
    if (nivelMin[nod] == nivel[nod]){
        vector<int> compCrt;
        int nodCompCrt;
        do{
            nodCompCrt = s.top();
            compCrt.push_back(nodCompCrt);
            s.pop();
            inStiva[nodCompCrt] = false;
        } while (nodCompCrt != nod);
        ctc.push_back(compCrt);
    }
}

void HavelHakimi(){
    vector<int> grade;
    int n;
    fin >> n;
    for (int i = 0 ; i < n; i ++){
        int grad;
        fin >> grad;
        grade.push_back(grad);
    }

    if (*max_element(begin(grade), end(grade)) > n - 1 or
        accumulate(grade.begin(), grade.end(), 0) % 2  == 1) {
        // daca maximul din vectorul de grade e mai mare decat n - 1
        // sau suma gradelor este impara => nu se poate reprezenta graful
        fout << "Nu\n";
    }
    else {
        while (true) {
            sort(grade.begin(), grade.end(), greater<>());
            int grad_crt = grade[0];
            if (grad_crt == 0){
                fout << "Da\n";
                return;
            }
            grade.erase(grade.begin());

            for (int i = 0; i < grad_crt; i++){
                grade[i]--;

                if (grade[i] < 0) {
                    fout << "Nu\n";
                    return;
                }
            }
        }
    }
}
/* TEMA 2*/

void Graf::adaugaInListaAdsiCosturiOrientat(const int &nod1, const int &nod2, const int &cost) {
    listaAdsiCosturi[nod1].push_back({nod2, cost});
}

void Graf::afisareDijkstra() {
//    vector<int> tata; //reconstruire drum
//    tata.resize(nrNoduri + 1, 0);
    bitset<NMAX> viz;
    vector<int> dist;
    dist.resize(nrNoduri + 1, INF);
    priority_queue<pair<int, int>, vector<pair<int,int>>, decltype(cmp)> minHeap(cmp);  // primul element e nodul apoi distanta
    dist[1] = 0;
    minHeap.push({1, 0});
    Dijkstra(viz, dist, minHeap);

    for(int i = 2; i <= nrNoduri; i++){
        if(viz[i]){
            fout << dist[i] << " ";
        } else { fout << 0 << " "; }
    }
}

void Graf::Dijkstra(bitset<NMAX> &viz, vector<int> &distante, priority_queue<pair<int, int>, vector<pair<int, int>>, decltype(cmp)> &minHeap) {
    while (!minHeap.empty()){
        auto top = minHeap.top();  // (nod, distanta minima)
        minHeap.pop();
        viz[top.first] = true;
        int nodCrt = top.first;
        int distCrt = top.second;
        int nrVecini = listaAdsiCosturi[nodCrt].size();
        for (int i = 0; i < nrVecini; i++){
            int vecinCrt = listaAdsiCosturi[nodCrt][i].first;
            int costCrt = listaAdsiCosturi[nodCrt][i].second;
            if (!viz[vecinCrt] and (distCrt + costCrt < distante[vecinCrt])){
                distante[vecinCrt] = distCrt + costCrt;
                minHeap.push({vecinCrt, distante[vecinCrt]});
            }
        }
    }
}

void Graf::afisareBellmanford() {
    vector<int> dist;
    queue<int> coadaNoduri;
    coadaNoduri.push(1);
    vector<int> nrRelaxari;
    bitset<NMAX> inCoada;
    inCoada[1] = true;
    nrRelaxari.resize(nrNoduri + 1, 0);
    dist.resize(nrNoduri + 1, INF);
    dist[1] = 0;
    bool rezultat = true;
    BellmanFord(coadaNoduri, dist, nrRelaxari, inCoada, rezultat);

    if (rezultat){
        for (int i = 2; i <= nrNoduri; i++){
            fout << dist[i] << " ";
        }
    } else {fout << "Ciclu negativ!" << "\n";}
}

void Graf::BellmanFord(queue<int> &coadaNoduri, vector<int> &dist, vector<int> &nrRelaxari, bitset<NMAX> &inCoada, bool &rez) {
    while (!coadaNoduri.empty()){
        int front = coadaNoduri.front();
        coadaNoduri.pop();
        inCoada[front] = false;
        int distCrt = dist[front];
        int nrVecini = listaAdsiCosturi[front].size();

        for (int j = 0; j < nrVecini; j++) {
            int vecinCrt = listaAdsiCosturi[front][j].first;
            int costCrt = listaAdsiCosturi[front][j].second;

            if (distCrt + costCrt < dist[vecinCrt]) {
                dist[vecinCrt] = distCrt + costCrt;
                nrRelaxari[vecinCrt] += 1;

                if (nrRelaxari[front] == nrNoduri) {
                    rez = false;
                    return;
                }

                if (!inCoada[vecinCrt]){
                    inCoada[vecinCrt] = true;
                    coadaNoduri.push(vecinCrt);
                }
            }
        }
    }
}

void Graf::adaugaMuchie(const int &nod1, const int &nod2, const int &cost) {
    muchiiCost.push_back(make_tuple(cost, nod1, nod2));
}

void Graf::afisareApm() {
    vector<pair<int, int>> sol;

    fout << Apm(sol) << "\n" << sol.size() << "\n";

    for(int i = 0; i < sol.size(); i ++){
        fout << sol[i].first << " " << sol[i].second << "\n";
    }
}

int Graf::Apm(vector<pair<int, int>> &sol) {
    int costApm = 0;
    DisjointSet disjointSet(nrNoduri);

    for(int i = 1; i <= nrNoduri; i++){
        disjointSet.initializare(i);
    }

    sort(muchiiCost.begin(), muchiiCost.end());

    for(int i = 0; i < nrMuchii; i++){
        int x = get<1>(muchiiCost[i]);
        int y = get<2>(muchiiCost[i]);
        if(disjointSet.gasesteReprez(x) != disjointSet.gasesteReprez(y)){
            disjointSet.reuneste(x, y);
            int costCrt = get<0>(muchiiCost[i]);
            costApm += costCrt;
            sol.push_back({get<1>(muchiiCost[i]), get<2>(muchiiCost[i])});
        }
    }
    return costApm;
}



int main() {
//    int noduri, muchii, s;
//    fin >> noduri >> muchii;
////    //fin >> noduri >> muchii >> s; // pt bfs
//    Graf G(noduri, muchii);
//    for (int i = 0; i < muchii; i++){
//        int n1, n2, cost;
////        fin >> n1 >> n2;
////        G.adaugaInListaAdOrientat(n1, n2);
////        G.adaugaInListaAdNeorientat(n1, n2);
//        fin >> n1 >> n2 >> cost;
//        G.adaugaMuchie(n1, n2, cost);
//    }

//    fout<<G.nrCmpConexe();

//    G.afisDistanteMinBfs(s);

//    G.sortareTopologica();

//    vector<vector<int>> componenteB = G.componenteBiconexe();
//    fout<<componenteB.size()<<"\n";
//    for(int i = 0; i < componenteB.size(); i++) {
//        for(int j = 0; j < componenteB[i].size(); j++) {
//            fout<<componenteB[i][j]<<"  ";
//        }
//        fout<<"\n";
//    }

//    vector<vector<int>> componenteTc = G.componenteTareConexe();
//    fout<<componenteTc.size()<<"\n";
//    for(int i = 0; i < componenteTc.size(); i++) {
//        for(int j = 0; j < componenteTc[i].size(); j++) {
//            fout<<componenteTc[i][j]<<" ";
//        }
//        fout<<"\n";
//    }

//    HavelHakimi();


    int nrM, nrOp;
    fin >> nrM >> nrOp;
    DisjointSet d(nrM);

    for(int i = 1;  i <= nrM; i++){
        d.initializare(i);
    }

    for (int i = 0; i < nrOp;  i++){
        int op, x, y;
        fin >> op >> x >> y;

        if (op == 1){
            d.reuneste(x, y);
        }
        else{
            bool ans = d.query(x, y);
            if (ans)
                fout << "DA\n";
            else
                fout << "NU\n";
        }
    }

//    G.afisareApm();

    return 0;
}

