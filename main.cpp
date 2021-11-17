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

const int NMAX = 100001;
const int INF = numeric_limits<int>::max();

//ifstream fin ("dfs.in");
//ofstream fout ("dfs.out");

//ifstream fin ("bfs.in");
//ofstream fout ("bfs.out");

//ifstream fin ("sortaret.in");
//ofstream fout ("sortaret.out");

//ifstream fin ("hh.in");
//ofstream fout ("hh.out");

//ifstream fin ("biconex.in");
//ofstream fout ("biconex.out");

//ifstream fin ("ctc.in");
//ofstream fout ("ctc.out");

ifstream fin ("apm.in");
ofstream fout ("apm.out");

//ifstream fin ("dijkstra.in");
//ofstream fout ("dijkstra.out");

//ifstream fin ("bellmanford.in");
//ofstream fout ("bellmanford.out");

//ifstream fin("disjoint.in");
//ofstream fout("disjoint.out");


auto cmp = [](const pair<int, int>& p1, const pair<int, int>& p2)
{
    return p1.second > p2.second;
};

class DisjointSet{
private:
    int nrMultimi, nrOperatii;
    vector<int> tata;
    vector<int> h; // vector de inaltimi
public:
    DisjointSet(int nrMultimi) : nrMultimi(nrMultimi){
        tata.resize(nrMultimi + 1, -1);
        h.resize(nrMultimi + 1, 0);
    }
    DisjointSet(int nrMultimi, int nrOperatii) : nrMultimi(nrMultimi), nrOperatii(nrOperatii){
        tata.resize(nrMultimi + 1, -1);
        h.resize(nrMultimi + 1, 0);
    }
    void initializare(const int &nr);
    int gasesteReprez(const int &nr);
    void reuneste(const int &x, const int &y);
    bool query(const int &x, const int &y);
    void afisareDisjointSet();
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

void DisjointSet::afisareDisjointSet() {
    for(int i = 1;  i <= nrMultimi; i++){
        initializare(i);
    }

    for (int i = 0; i < nrOperatii;  i++){
        int op, x, y;
        fin >> op >> x >> y;

        if (op == 1){
            reuneste(x, y);
        }
        else{
            bool ans = query(x, y);
            if (ans)
                fout << "DA\n";
            else
                fout << "NU\n";
        }
    }
}

class Graf{
private:
    int nrNoduri, nrMuchii;
    vector<int> listaAd[NMAX];
    vector<pair<int, int>> listaAdsiCosturi[NMAX];
    vector<tuple<int, int, int>> muchiiCost;
//    bitset<NMAX> viz;

    void dfs(const int &nod, bitset<NMAX> &viz);
    void bfs(const int &start,  vector<int> &distante);
    void dfsSortTop(const int &nod,   stack<int> &noduriSortTop, bitset<NMAX> &viz);
    void dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, bitset<NMAX> &viz, vector<vector<int>> &compBiconexe);
    void dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, bitset<NMAX> &inStiva, bitset<NMAX> &viz, vector<vector<int>> &ctc);
    void Dijkstra(bitset<NMAX> &viz, vector<int> &distante, priority_queue<pair<int, int>, vector<pair<int,int>>, decltype(cmp)> &minHeap);
    void BellmanFord(queue<int> &coadaNoduri, vector<int> &dist, vector<int> &nrRelaxari, bitset<NMAX> &inCoada, bool &rez);
    int Apm(vector<pair<int, int>> &sol);
public:
    Graf(int nrNoduri, int nrMuchii) : nrNoduri(nrNoduri), nrMuchii(nrMuchii) {};
    void adaugaMuchie(const int &nod1, const int &nod2, bool eOrientat);
    void adaugaInListaAdSiCosturi(const int &nod1, const int &nod2, const int &cost);
    void adaugaMuchieCuCost(const int &nod1, const int &nod2, const int &cost);
    int nrCmpConexe();
    void afisareDistanteMinBfs(int start);
    void afisareSortareTopologica();
    void afisareComponenteBiconexe();
    void afisareComponenteTareConexe();
    void afisareDijkstra();
    void afisareBellmanford();
    void afisareApm();
};

void Graf::adaugaMuchie(const int &nod1, const int &nod2, bool eOrientat) {
    listaAd[nod1].push_back(nod2);
    if(!eOrientat){
        listaAd[nod2].push_back(nod1);
    }
}

void Graf::adaugaInListaAdSiCosturi(const int &nod1, const int &nod2, const int &cost) {
    listaAdsiCosturi[nod1].push_back({nod2, cost});
}

void Graf::adaugaMuchieCuCost(const int &nod1, const int &nod2, const int &cost) {
    muchiiCost.push_back(make_tuple(cost, nod1, nod2));
}

void Graf::dfs(const int &nod, bitset<NMAX> &viz) {
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            viz[listaAd[nod][i]] = true;
            dfs(listaAd[nod][i], viz);
        }
    }
}

int Graf::nrCmpConexe() {
    bitset<NMAX> viz;
    int nr = 0;
    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i]){
            nr++;
            viz[i] = true;
            dfs(i, viz);
        }
    }
    return nr;
}

void Graf::afisareDistanteMinBfs(int nodStart) {
    vector<int> distante;
    distante.resize(nrNoduri + 1, -1);
    bfs(nodStart, distante);
    for (int i = 1; i <= nrNoduri; i ++){
        fout << distante[i] <<" ";
    }
}

void Graf::bfs(const int &start, vector<int> &distante) {
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



void Graf::afisareSortareTopologica() {
    bitset<NMAX> viz;
    stack<int> noduriSortTop;
    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i]) {
            dfsSortTop(i, noduriSortTop, viz);
        }
    }
    while (!noduriSortTop.empty()){
        fout << noduriSortTop.top() << " ";
        noduriSortTop.pop();
    }
}

void Graf::dfsSortTop(const int &nod, stack<int> &noduriSortTop, bitset<NMAX> &viz) {
    viz[nod] = true;
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            dfsSortTop(listaAd[nod][i], noduriSortTop, viz);
        }
    }
    noduriSortTop.push(nod);
}

void Graf :: afisareComponenteBiconexe() {
    bitset<NMAX> viz;
    vector<int> nivel;
    vector<int> nivelMin;
    nivel.resize(NMAX);
    nivelMin.resize(NMAX);
    stack<int> s;
    vector<vector<int>> compBiconexe;

    dfsCompBiconexe(1,  1, nivel, nivelMin, s, viz, compBiconexe);

    fout<<compBiconexe.size()<<"\n";
    for(int i = 0; i < compBiconexe.size(); i++) {
        for(int j = 0; j < compBiconexe[i].size(); j++) {
            fout<<compBiconexe[i][j]<<"  ";
        }
        fout<<"\n";
    }
}

void Graf::dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, bitset<NMAX> &viz, vector<vector<int>> &compBiconexe) {
    viz[nod] = true;
    s.push(nod);
    nivel[nod] = nivelCrt;
    nivelMin[nod] = nivelCrt;

    for (int i = 0; i < listaAd[nod].size(); i++){
        int vecin = listaAd[nod][i];
        if (!viz[vecin]){
            dfsCompBiconexe(vecin, nivelCrt + 1, nivel, nivelMin, s, viz, compBiconexe);
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

void Graf :: afisareComponenteTareConexe() {
    bitset<NMAX> viz;
    vector<int> nivel;
    vector<int> nivelMin;
    nivel.resize(NMAX);
    nivelMin.resize(NMAX);
    stack<int> s;
    bitset<NMAX> inStiva;
    vector<vector<int>> ctc;
    int niv = 1;

    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i])
            dfsCompTareConexe(i,  niv, nivel, nivelMin, s, inStiva,viz, ctc);
    }

    fout<<ctc.size()<<"\n";
    for (int i = 0; i < ctc.size(); i++) {
        for(int j = 0; j < ctc[i].size(); j++) {
            fout<<ctc[i][j]<<" ";
        }
        fout<<"\n";
    }
}

void Graf::dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, bitset<NMAX> &inStiva, bitset<NMAX> &viz, vector<vector<int>> &ctc) {
    viz[nod] = true;
    s.push(nod);
    inStiva[nod] = true;
    nivel[nod] = nivelCrt;
    nivelMin[nod] = nivelCrt;
    nivelCrt++;

    for (int i = 0; i < listaAd[nod].size(); i++) {
        int vecin = listaAd[nod][i];
        if (!viz[vecin]) {
            dfsCompTareConexe(vecin, nivelCrt, nivel, nivelMin, s, inStiva, viz, ctc);
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
    int noduri, muchii, s;
    fin >> noduri >> muchii;
//    fin >> noduri >> muchii >> s; // BFS
    Graf G(noduri, muchii);
    for (int i = 0; i < muchii; i++){
        int n1, n2, cost;
//        fin >> n1 >> n2;
//        G.adaugaMuchie(n1, n2, 1);
//        G.adaugaMuchie(n1, n2, 0);
        fin >> n1 >> n2 >> cost;
//        G.adaugaInListaAdSiCosturi(n1, n2, cost); // DIJKSTRA + BELLMAN-FORD
        G.adaugaMuchieCuCost(n1, n2, cost); // APM
    }

//    fout << G.nrCmpConexe(); // DFS, neorientat

//    G.afisareDistanteMinBfs(s); // BFS, orientat

//    G.afisareSortareTopologica(); // Graf orientat

//    G.afisareComponenteBiconexe(); //Graf neorientat

//    G.afisareComponenteTareConexe(); //Graf orientat

//    HavelHakimi(); // Daca se poate construi un graf Neorientat stiind vector de grade

//    G.afisareDijkstra(); // Graf orientat

//    G.afisareBellmanford(); // Graf orientat

//// DisjointSet
//    int nrM, nrOp;
//    fin >> nrM >> nrOp;
//    DisjointSet d(nrM, nrOp);
//    d.afisareDisjointSet();

    G.afisareApm(); // Graf conex neorientat

    return 0;
}

