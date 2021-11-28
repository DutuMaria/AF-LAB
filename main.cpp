#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include <numeric>
#include <limits>
#include <tuple>
#include <utility>

using namespace std;

const int NMAX = 50001;
const int NMAX_SORT_TOP = 50001;
const int NMAX_BFS = 100001;
const int NMAX_DFS = 100001;
const int NMAX_COMP_BICONEXE = 100001;
const int NMAX_COMP_TARE_CONEXE = 100001;
const int NMAX_DIJKSTRA = 50001;
const int NMAX_BELLMAN_FORD = 50001;
const int NMAX_APM = 200001;
const int NMAX_DARB = 100001;

const int INF = numeric_limits<int>::max();

auto cmp = [](const pair<int, int>& p1, const pair<int, int>& p2)
{
    return p1.second > p2.second;
};


class DisjointSets{
private:
    int nrMultimi, nrOperatii;
    vector<int> tata;
    vector<int> h; // vector de inaltimi
public:
    DisjointSets(int nrMultimi) : nrMultimi(nrMultimi){
        tata.resize(nrMultimi + 1, -1);
        h.resize(nrMultimi + 1, 0);
    }
    DisjointSets(int nrMultimi, int nrOperatii) : nrMultimi(nrMultimi), nrOperatii(nrOperatii){
        tata.resize(nrMultimi + 1, -1);
        h.resize(nrMultimi + 1, 0);
    }
    void initializare(const int &nr);
    int gasesteReprez(const int &nr);
    void reuneste(const int &x, const int &y);
    bool query(const int &x, const int &y);
};

void DisjointSets::initializare(const int &nr) {
    tata[nr] = 0;
    h[nr] = 0;
}

int DisjointSets::gasesteReprez(const int &nr) {
    if (tata[nr] == 0)
        return nr;
    tata[nr] = gasesteReprez(tata[nr]);
    return tata[nr];
}

void DisjointSets::reuneste(const int &x, const int &y) {
    int repX = gasesteReprez(x), repY = gasesteReprez(y);
    if (h[repX] > h[repY]){
        tata[repY] = repX;
    } else {
        tata[repX] = repY;
        if (h[repX] == h[repY])
            h[repY]++;
    }
}

bool DisjointSets::query(const int &x, const int &y) {
    return gasesteReprez(x) == gasesteReprez(y);
}


class Graf{
private:
    int nrNoduri, nrMuchii;
    bool eOrientat;
    vector<int> listaAd[NMAX];
    vector<pair<int, int>> listaAdsiCosturi[NMAX]; // pentru Bellman-Ford + Dijkstra
    vector<tuple<int, int, int>> muchiiCost; // pentru APM

    void dfs(const int &nod, vector<bool> &viz);
    void dfsSortTop(const int &nod, stack<int> &noduriSortTop, vector<bool> &viz);
    void dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, vector<bool> &viz, vector<vector<int>> &compBiconexe);
    void dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, vector<bool> &inStiva, vector<bool> &viz, vector<vector<int>> &ctc);
    int nodCuDistantaMaxima(const vector<int> &distante);
public:
    Graf(int nrNoduri, bool eOrientat) : nrNoduri(nrNoduri), eOrientat(eOrientat) {};
    Graf(int nrNoduri, int nrMuchii, bool eOrientat) : nrNoduri(nrNoduri), nrMuchii(nrMuchii), eOrientat(eOrientat) {};
    void adaugaInListaAd(const int &nod1, const int &nod2);
    void adaugaInListaAdSiCosturi(const int &nod1, const int &nod2, const int &cost);
    void adaugaMuchieCuCost(const int &nod1, const int &nod2, const int &cost);
    int nrCmpConexe();
    vector<int> bfs(const int &start = 1);
    stack<int> sortareTopologica();
    vector<vector<int>> componenteBiconexe();
    vector<vector<int>> componenteTareConexe();
    vector<int> Dijkstra(const int &s = 1);
    bool BellmanFord(vector<int> &dist, const int &s = 1);
    int Apm(vector<pair<int, int>> &sol);
    vector<vector<long long>> RoyFloyd(vector<vector<long long>> &matrDrumuriMin);
    int diametruArbore();
};

void Graf::adaugaInListaAd(const int &nod1, const int &nod2) {
    listaAd[nod1].push_back(nod2);
    if (!eOrientat){
        listaAd[nod2].push_back(nod1);
    }
}

void Graf::adaugaInListaAdSiCosturi(const int &nod1, const int &nod2, const int &cost) {
    listaAdsiCosturi[nod1].push_back({nod2, cost});
}

void Graf::adaugaMuchieCuCost(const int &nod1, const int &nod2, const int &cost) {
    muchiiCost.push_back(make_tuple(cost, nod1, nod2));
}

void Graf::dfs(const int &nod, vector<bool> &viz) {
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            viz[listaAd[nod][i]] = true;
            dfs(listaAd[nod][i], viz);
        }
    }
}

int Graf::nrCmpConexe() {
    vector<bool> viz;
    viz.resize(nrNoduri + 1);
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

vector<int> Graf::bfs(const int &start) {
    vector<int> distante;
    distante.resize(nrNoduri + 1, -1);
    queue<int> coada;
    coada.push(start);
    distante[start] = 0;
    while (!coada.empty()){
        int nod = coada.front();
        for (int i = 0;  i < listaAd[nod].size(); i ++){
            if (distante[listaAd[nod][i]] == -1){
                coada.push(listaAd[nod][i]);
                distante[listaAd[nod][i]] = distante[nod] + 1;
            }
        }
        coada.pop();
    }
    return distante;
}

stack<int> Graf::sortareTopologica() {
    vector<bool> viz;
    viz.resize(nrNoduri + 1);
    stack<int> noduriSortTop;

    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i]) {
            dfsSortTop(i, noduriSortTop, viz);
        }
    }
    return noduriSortTop;
}

void Graf::dfsSortTop(const int &nod, stack<int> &noduriSortTop, vector<bool> &viz) {
    viz[nod] = true;
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            dfsSortTop(listaAd[nod][i], noduriSortTop, viz);
        }
    }
    noduriSortTop.push(nod);
}

vector<vector<int>> Graf :: componenteBiconexe() {
    vector<bool> viz;
    viz.resize(nrNoduri + 1);
    vector<int> nivel;
    nivel.resize(nrNoduri + 1);
    vector<int> nivelMin;
    nivelMin.resize(nrNoduri + 1);
    stack<int> s;
    vector<vector<int>> compBiconexe;
    dfsCompBiconexe(1, 1, nivel, nivelMin, s, viz, compBiconexe);

    return compBiconexe;
}

void Graf :: dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, vector<bool> &viz, vector<vector<int>> &compBiconexe) {
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
                do {
                    nodCompCrt = s.top();
                    compCrt.push_back(nodCompCrt);
                    s.pop();
                } while (nodCompCrt != vecin);
                compCrt.push_back(nod);
                compBiconexe.push_back(compCrt);
            }
        } else { // daca e vizitat => muchie de intoarcere => actualizam nivelMin
            nivelMin[nod] = min(nivelMin[nod], nivel[vecin]);
        }
    }
}

vector<vector<int>> Graf :: componenteTareConexe() {
    vector<bool> viz;
    viz.resize(nrNoduri + 1);
    vector<int> nivel;
    nivel.resize(nrNoduri + 1);
    vector<int> nivelMin;
    nivelMin.resize(nrNoduri + 1);
    stack<int> s;
    vector<bool> inStiva;
    inStiva.resize(nrNoduri + 1);
    vector<vector<int>> ctc;
    int niv = 1;

    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i])
            dfsCompTareConexe(i,  niv, nivel, nivelMin, s, inStiva,viz, ctc);
    }
    return ctc;
}

void Graf::dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, vector<bool> &inStiva, vector<bool> &viz, vector<vector<int>> &ctc) {
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
        do {
            nodCompCrt = s.top();
            compCrt.push_back(nodCompCrt);
            s.pop();
            inStiva[nodCompCrt] = false;
        } while (nodCompCrt != nod);
        ctc.push_back(compCrt);
    }
}

vector<int> countingSort(const vector<int> &vect) {
    int maxim = *max_element(vect.begin(), vect.end());
    vector<int> fr1(maxim + 1, 0);

    for (int i : vect) {
        fr1[i]++;
    }

    vector<int> fr2;
    fr2.reserve(vect.size());
    for (int j = 0; j <= maxim; j++) {
        while (fr1[j]) {
            fr2.push_back(j);
            fr1[j]--;
        }
    }
    return fr2;
}

bool HavelHakimi(const int n, const vector<int> &grade){
    vector<int> copieGrade = grade;
    if (*max_element(begin(copieGrade), end(copieGrade)) > n - 1 or
        accumulate(copieGrade.begin(), copieGrade.end(), 0) % 2  == 1) {
        // daca maximul din vectorul de grade e mai mare decat n - 1
        // sau suma gradelor este impara => nu se poate reprezenta graful
       return false;
    } else {
        while (true) {
//            sort(grade.begin(), grade.end(), greater<>());
            copieGrade = countingSort(copieGrade);
            int grad_crt = copieGrade[0];
            if (grad_crt == 0){
                return true;
            }
            copieGrade.erase(copieGrade.begin());

            for (int i = 0; i < grad_crt; i++){
                copieGrade[i]--;

                if (copieGrade[i] < 0) {
                    return false;
                }
            }
        }
    }
}

/* TEMA 2 */

vector<int> Graf::Dijkstra(const int &s) {
    vector<bool> viz;
    viz.resize(nrNoduri + 1);
    vector<int> distante;
    distante.resize(nrNoduri + 1, INF);
    priority_queue<pair<int, int>, vector<pair<int,int>>, decltype(cmp)> minHeap(cmp);  // primul element e nodul apoi distanta
    distante[s] = 0;
    minHeap.push({s, distante[s]});
    while (!minHeap.empty()){
        auto top = minHeap.top();  // (nod, distanta minima)
        minHeap.pop();
        if (viz[top.first]){
            continue;
        }
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
    return distante;
}

bool Graf::BellmanFord(vector<int> &dist, const int &s) {
    queue<int> coadaNoduri;
    vector<int> nrRelaxari;
    nrRelaxari.resize(nrNoduri + 1, 0);
    vector<bool> inCoada;
    inCoada.resize(nrNoduri + 1, false);
    coadaNoduri.push(s);
    inCoada[s] = true;
    dist[s] = 0;
    while (!coadaNoduri.empty()){
        int front = coadaNoduri.front(); // (nod, distanta minima)
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
                    return true;
                }

                if (!inCoada[vecinCrt]){
                    inCoada[vecinCrt] = true;
                    coadaNoduri.push(vecinCrt);
                }
            }
        }
    }
    return false;
}

int Graf::Apm(vector<pair<int, int>> &sol){
    int costApm = 0;
    DisjointSets disjointSet(nrNoduri);

    for (int i = 1; i <= nrNoduri; i++){
        disjointSet.initializare(i);
    }

    sort(muchiiCost.begin(), muchiiCost.end());

    for (int i = 0; i < nrMuchii; i++){
        int x = get<1>(muchiiCost[i]);
        int y = get<2>(muchiiCost[i]);
        if (disjointSet.gasesteReprez(x) != disjointSet.gasesteReprez(y)){
            disjointSet.reuneste(x, y);
            int costCrt = get<0>(muchiiCost[i]);
            costApm += costCrt;
            sol.push_back({get<1>(muchiiCost[i]), get<2>(muchiiCost[i])});
        }
    }
    return costApm;
}

/*
    Floyd-Warshall => Complexitate O(n^3)
=> Graf ORIENTAT cu n noduri, memorat prin matricea ponderilor
 (ponderile pot fi si negative dar NU exista circuite cu cost negativ in G)
=> Cerinta: Pentru oricare doua varfuri x si y ale lui G, sa se determine distanta de la x la y
si un drum minim de la x la y.
=> Idee : Fie d matricea drumurilor (initial egala cu matricea costurilor)
 Parcurgem multimea nodurilor (ca varfuri intermediare) si
 pentru fiecare drum de la un nod i la un nod j, daca varful k este varf intermediar al drumului, atunci
 drumul de la i la j o sa fie minimul dintre drumul gasit anterior si drumul de la i la j prin k.
 d[i][j] = min(d[i][j], d[i][k] + d[k][j]
=> In plus: Afisarea unui drum de la i la j, daca d[i][j] < INF, se face folosind matricea p (matricea predecesorilor)
 (adica cand am reactualizat drumul de la i la j, actualizam si predecesorul, p[i][j] = p[k][j])
=> OBS Pentru un graf neorientat matricea drumurilor minime o sa fie simetrica
*/

vector<vector<long long>> Graf::RoyFloyd(vector<vector<long long>> &matrDrumuriMin) {
    for (int k = 1; k <= nrNoduri; k++){
        for (int i = 1; i <= nrNoduri; i++){
            for (int j = 1; j <= nrNoduri; j++){
                if (matrDrumuriMin[i][j] > matrDrumuriMin[i][k] + matrDrumuriMin[k][j]){
                    matrDrumuriMin[i][j] = matrDrumuriMin[i][k] + matrDrumuriMin[k][j];
                }
            }
        }
    }

    return matrDrumuriMin;
}

int Graf::nodCuDistantaMaxima(const vector<int> &distante) {
    auto  ptNodDistMax= max_element(distante.begin(), distante.end());
    int nod = distance(distante.begin(), ptNodDistMax);
    return nod; //nod cu distanta maxima
}

int Graf::diametruArbore() {
    vector<int> distante = bfs();

    int nod1 = nodCuDistantaMaxima(distante);
    distante = bfs(nod1);
    int nod2 = nodCuDistantaMaxima(distante);
    int diametru = distante[nod2] + 1;

    return diametru;
}

/* -------------------------------------------------------------------------------------------------------------------*/

/*
 DFS, graf Neorientat
*/

void rezultatDfs(){
//    https://www.infoarena.ro/problema/dfs
    ifstream fin ("dfs.in");
    ofstream fout ("dfs.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    Graf G(noduri, muchii, false);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    fout << G.nrCmpConexe() << "\n";
}

/*
 BFS, graf Orientat
*/

void rezultatBfs(){
//    https://www.infoarena.ro/problema/bfs
    ifstream fin ("bfs.in");
    ofstream fout ("bfs.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii >> s;
    Graf G(noduri, muchii, true);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    vector<int> distanteMinBfs = G.bfs(s);
    for (int i = 1; i <= noduri; i++){
        fout << distanteMinBfs[i] <<" ";
    }
}

/*
 Sortare Topologica, graf Orientat
*/

void rezultatSortareTopologica(){
//    https://www.infoarena.ro/problema/sortaret
    ifstream fin ("sortaret.in");
    ofstream fout ("sortaret.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    Graf G(noduri, muchii, true);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    stack<int> noduriSortTop =  G.sortareTopologica();;
    while (!noduriSortTop.empty()){
        fout << noduriSortTop.top() << " ";
        noduriSortTop.pop();
    }
}

/*
 Componente biconexe, graf Neorientat
*/

void rezultatComponenteBiconexe(){
//    https://www.infoarena.ro/problema/biconex
    ifstream fin ("biconex.in");
    ofstream fout ("biconex.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    Graf G(noduri, muchii, false);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    vector<vector<int>> compBiconexe = G.componenteBiconexe();
    fout << compBiconexe.size() << "\n";
    for (int i = 0; i < compBiconexe.size(); i++) {
        for (int j = 0; j < compBiconexe[i].size(); j++) {
            fout << compBiconexe[i][j] << "  ";
        }
        fout << "\n";
    }
}

/*
 Componente Tare Conexe, graf Orientat
*/

void rezultatComponenteTareConexe() {
//    https://www.infoarena.ro/problema/ctc
    ifstream fin("ctc.in");
    ofstream fout("ctc.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    Graf G(noduri, muchii, true);
    for (int i = 0; i < muchii; i++) {
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    vector<vector<int>> compTareConexe = G.componenteTareConexe();
    fout << compTareConexe.size() << "\n";
    for (int i = 0; i < compTareConexe.size(); i++) {
        for (int j = 0; j < compTareConexe[i].size(); j++) {
            fout << compTareConexe[i][j] << " ";
        }
        fout << "\n";
    }
}

/*
 Havel Hakimi => Daca se poate construi un graf Neorientat stiind vectorul de grade
*/

void rezultatHavelHakimi(){
    ifstream fin ("hh.in");
    ofstream fout ("hh.out");
    int n;
    vector<int> grade;
    fin >> n;
    for (int i = 0 ; i < n; i ++){
        int grad;
        fin >> grad;
        grade.push_back(grad);
    }
    bool rezultat = HavelHakimi(n, grade);
    if (rezultat){
        fout << "Da\n";
    } else {
        fout << "Nu\n";
    }
}

/*
 Dijkstra, graf Orientat
*/

void rezultatDijkstra(){
//    https://www.infoarena.ro/problema/dijkstra
    ifstream fin ("dijkstra.in");
    ofstream fout ("dijkstra.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    //    fin >> noduri >> muchii >> s;
    Graf G(noduri, muchii, true);
    for (int i = 0; i < muchii; i++){
        int n1, n2, cost;
        fin >> n1 >> n2 >> cost;
        G.adaugaInListaAdSiCosturi(n1, n2, cost); // pe Orientat
    }
//    G.Dijkstra(s);
    vector<int> dist = G.Dijkstra();
    for (int i = 2; i <= noduri; i++){
        if (dist[i] == INF){
            fout << 0 << " ";
        } else {
            fout << dist[i] << " ";
        }
    }
}

/*
 Bellman-Ford, graf Orientat
*/

void rezultatBellmanFord(){
//    https://www.infoarena.ro/problema/bellmanford
    ifstream fin ("bellmanford.in");
    ofstream fout ("bellmanford.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    //    fin >> noduri >> muchii >> s;
    Graf G(noduri, muchii, true);
    for (int i = 0; i < muchii; i++){
        int n1, n2, cost;
        fin >> n1 >> n2 >> cost;
        G.adaugaInListaAdSiCosturi(n1, n2, cost); //  pe Orientat
    }
    vector<int> distante;
    distante.resize(noduri + 1, INF);
    // bool eCircuit = G.BellmanFord(s);
    bool eCircuit = G.BellmanFord(distante);
    if (!eCircuit){
        for (int i = 2; i <= noduri; i++){
            fout << distante[i] << " ";
        }
    } else {
        fout << "Ciclu negativ!" << "\n";
    }
}

void rezultatDisjoint(){
//    https://www.infoarena.ro/problema/disjoint
    ifstream fin("disjoint.in");
    ofstream fout("disjoint.out");
    int nrMultimi, nrOperatii;
    fin >> nrMultimi >> nrOperatii;
    DisjointSets d(nrMultimi, nrOperatii);
    for (int i = 1;  i <= nrMultimi; i++){
        d.initializare(i);
    }
    for (int i = 0; i < nrOperatii;  i++){
        int op, x, y;
        fin >> op >> x >> y;
        if (op == 1){
            d.reuneste(x, y);
        } else {
            bool ans = d.query(x, y);
            if (ans){
                fout << "DA\n";
            } else {
                fout << "NU\n";
            }
        }
    }
}

/*
 Arbori partiali de cost minim, graf conex Neorientat
*/

void rezultatApm(){
//    https://www.infoarena.ro/problema/apm
    ifstream fin ("apm.in");
    ofstream fout ("apm.out");
    int noduri, muchii;
    fin >> noduri >> muchii;
    Graf G(noduri, muchii, false);
    for (int i = 0; i < muchii; i++){
        int n1, n2, cost;
        fin >> n1 >> n2 >> cost;
        G.adaugaMuchieCuCost(n1, n2, cost);
    }
    vector<pair<int, int>> sol;
    fout << G.Apm(sol) << "\n" << sol.size() << "\n";
    for (int i = 0; i < sol.size(); i ++){
        fout << sol[i].first << " " << sol[i].second << "\n";
    }
}

/*
    Floyd-Warshall, Graf Orientat
*/

void rezultatRoyFloyd(){
//    https://infoarena.ro/problema/royfloyd
    ifstream fin("royfloyd.in");
    ofstream fout("royfloyd.out");
    int nrNoduri;
    fin >> nrNoduri;
    Graf g(nrNoduri, true);
    vector<vector<long long>> matriceDrumuriMin(nrNoduri + 1, vector<long long>(nrNoduri + 1));
    for (int i = 1;  i <= nrNoduri; i++){
        for (int j = 1; j <= nrNoduri; j ++){
            int nr;
            fin >> nr;
            if (nr == 0 and i != j){
                matriceDrumuriMin[i][j] = INF;
            } else {
                matriceDrumuriMin[i][j] = nr;
            }
        }
    }
    matriceDrumuriMin = g.RoyFloyd(matriceDrumuriMin);
    for (int i = 1;  i <= nrNoduri; i++) {
        for (int j = 1; j <= nrNoduri; j++) {
            int nr = matriceDrumuriMin[i][j];
            if (nr == INF){
                fout << 0 << " ";
            } else {
                fout << matriceDrumuriMin[i][j] << " ";
            }
        }
        fout << "\n";
    }
}

/*
    Darb - Diametrul unui arbore
*/

void rezultatDarb(){
//    https://infoarena.ro/problema/darb
    ifstream fin("darb.in");
    ofstream fout("darb.out");
    int nrNoduri;
    fin >> nrNoduri;
    Graf g(nrNoduri, false);
    for(int i = 0; i < nrNoduri; i ++){
        int n1, n2;
        fin >> n1 >> n2;
        g.adaugaInListaAd(n1, n2);
    }
    fout << g.diametruArbore();
}

int main() {
    rezultatDfs();
    rezultatBfs();
    rezultatSortareTopologica();
    rezultatComponenteBiconexe();
    rezultatComponenteTareConexe();
    rezultatHavelHakimi();
    rezultatDijkstra();
    rezultatBellmanFord();
    rezultatDisjoint();
    rezultatApm();
    rezultatRoyFloyd();
    rezultatDarb();

    return 0;
}

