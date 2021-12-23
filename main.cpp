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

const int INF = numeric_limits<int>::max();

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

/*
 Paduri de multimi disjuncte
 COMPLEXITATE => O(M * log N)
    INITIALIZARE: O(1), dar se face de N ori => O(N)
    FIND_REPREZ: O(log N), dar se face de 2M ori =>  O(M * log N)
    UNION: O(log n), dar se face de N ori => O(N * log N)

 Ce am folosit:
    => memoram varfurile fiecarei multimi ca un arbore (vector<int> tata),  avand ca reprezentant radacina
    => pentru inaltimi -> vector<int> h

 Descriere algoritm:
    => INITIALIZARE => fiecare nod reprezinta un arbore (doar el, fiind si radacina)
    => FIND_REPREZ  => pentru un nod, functia returneaza radacina arborelui din care face parte + optimizare (compresie de cale)
        => compresie de cale => tatal varfurilor de pe lantul de la nod (cel pentru care apelam functia) la radacina
                                se va seta ca fiind radacina => pentru ca reprezentantul fiecarui varf de pe acest lant
                                sa fie mai usor de gasit in cautarile ulterioare
                             => inaltimea NU se schimba
    => UNION => calculam reprezentantii pentru nodurile pe care vrem sa le unim
             => daca inaltimea primul nod e mai mare decat inaltimea ceilui de-al doilea atunci
                reprezentatntull delui de-al doilea devine reprezentantul primului, daca nu => invers
             => in cazul in care sunt egale inaltimile alegem oricare din rerezentanti si crestem inaltimea cu 1
*/
void DisjointSets::initializare(const int &nr) {
    tata[nr] = 0;
    h[nr] = 0;
}

int DisjointSets::gasesteReprez(const int &nr) {
    if (tata[nr] == 0)
        return nr;
    tata[nr] = gasesteReprez(tata[nr]); // Optimizare - compresie de cale
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
    bool ePonderat;
    vector<vector<int>> listaAd;
    vector<vector<pair<int, int>>> listaAdsiCosturi; // pentru Bellman-Ford + Dijkstra + MaxFlow
    vector<tuple<int, int, int>> muchiiCost; // pentru APM

    void dfs(const int &nod, vector<bool> &viz);
    void dfsSortTop(const int &nod, stack<int> &noduriSortTop, vector<bool> &viz);
    void dfsCompBiconexe(int nod, int nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, vector<bool> &viz, vector<vector<int>> &compBiconexe);
    void dfsCompTareConexe(int nod, int &nivelCrt, vector<int> &nivel, vector<int> &nivelMin, stack<int> &s, vector<bool> &inStiva, vector<bool> &viz, vector<vector<int>> &ctc);
    int nodCuDistantaMaxima(const vector<int> &distante);
    bool construiesteDrum(const int &start, const int &final, vector<int> &tata, vector<vector<int>> &capacitate, vector<vector<int>> &flux);

public:
    Graf(int nrNoduri, bool eOrientat) : nrNoduri(nrNoduri), eOrientat(eOrientat) {}; // Pentru Floyd-Warshall/Roy-Floyd
    Graf(int nrNoduri, int nrMuchii, bool eOrientat) : nrNoduri(nrNoduri), nrMuchii(nrMuchii), eOrientat(eOrientat){}; // Pentru APM
    Graf(int nrNoduri, int nrMuchii, bool eOrientat, bool ePonderat) : nrNoduri(nrNoduri), nrMuchii(nrMuchii), eOrientat(eOrientat), ePonderat(ePonderat) {
        if(!ePonderat){
            listaAd.resize(nrNoduri + 1, vector<int>());
        } else {
            listaAdsiCosturi.resize(nrNoduri + 1, vector<pair<int, int>>()); // pentru Dijkstra, Bellman-Ford si MaxFlow
        }
    };

    void adaugaInListaAd(const int &nod1, const int &nod2);
    void adaugaInListaAdSiCosturi(const int &nod1, const int &nod2, const int &cost); // pentru Dijkstra, Bellman-Ford si MaxFlow
    void adaugaMuchieCuCost(const int &nod1, const int &nod2, const int &cost); // Pentru APM

    int nrCmpConexe();
    vector<int> bfs(const int &start = 1);
    vector<int> sortareTopologica();
    vector<vector<int>> componenteBiconexe();
    vector<vector<int>> componenteTareConexe();
    vector<int> Dijkstra(const int &s = 1);
    bool BellmanFord(vector<int> &dist, const int &s = 1);
    int Apm(vector<pair<int, int>> &sol);
    vector<vector<long long>> RoyFloyd(vector<vector<long long>> &matrDrumuriMin);
    int diametruArbore();
    int maxFlow(const int &s, const int &t);
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

/*
 DFS, graf Neorientat
 COMPLEXITATE => O(M + N)

 Ce am folosit:
    => vector<bool> viz este vectorul de vizitate
        * viz[i] = false inseamna ca nodul i nu a fost vizitat
        * viz[i] = true inseamna ca nodul i a fost vizitat
    => int nr =  numarul total de componente conexe (rezultatul dorit)

 Descriere algoritm:
   => 2 functii -> dfs => se apeleaza dintr-un nod nevizitat
                       => parcurgem vecinii nodului
                       => daca vecinul curent nu este vizitat atunci il vizitam si apelam dfs pentru vecinul respectiv
                -> nrCmpConexe => parcurgem toate nodurile
                               => cand gasim un nod nevizitat inseamna ca am gasit o noua componenta conexa si apelam dfs pentru nodul respectiv
*/

void Graf::dfs(const int &nod, vector<bool> &viz) {
    for (int i = 0; i < listaAd[nod].size(); i++){
        if (!viz[listaAd[nod][i]]){
            viz[listaAd[nod][i]] = true;
            dfs(listaAd[nod][i], viz);
        }
    }
}

int Graf::nrCmpConexe() {
    vector<bool> viz(nrNoduri + 1);
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

/*
 BFS, graf Orientat
 COMPLEXITATE => O(M + N)

 Ce am folosit:
     => vector<int> distante => distante[i] = x  => inseamna ca distanta de la nodul de start la i este x
                             => nodul de start are distanta 0
     => queue<int> coada => in coada se afla initial nodul de start, ulterior punem pentru fiecare nod vecinii nevizitati in coada

 Descriere algoritm:
    => cat timp coada are elemente
        => luam nodul din front
        => adaugam in coada toti vecinii nevizitati ai nodului
        => crestem distanta
        => dupa ce am procesat toti vecinii nodului il scoatem din coada
*/

vector<int> Graf::bfs(const int &start) {
    vector<int> distante(nrNoduri + 1, -1);
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

/*
 Sortare Topologica, graf Orientat
 COMPLEXITATE => O(M + N)

 Ce am folosit:
    => vector<int> viz - vectorul cu starea nodurilor (vizitat sau nu)
    => stack<int> noduriSortTop - se pun nodurile pe stiva cand cand ne intoarcem din recursivitate
    => vector<int> rezultat - luam de pe stiva nodurile
                            => dau pop in ordinea inversa apelurilor recursive si asa se mentine ordinea nodurilor

 Descriere algoritm:
     => 2 functii -> sortareTopologica => parcurg nodurile si daca nodul curent nu e vizitat apelez pentru nodul respectiv dfsSortTop
                                       => dupa ce le-am parcurs pe toate, pun toate nodurile din stiva in vectorul rezultat
                  -> dfsSortTop => vizitez nodul pentru care fac apelul
                                => parcurg vecinii nodului curent si daca nu sunt vizitati apelam dfsSortTop pentru fiecare vecin
 OBS -> VARIANTA 2 => CU VECTOR DE GRADE INTERIOARE + COADA PENTRU PARCURGERE
*/

vector<int> Graf::sortareTopologica() {
    vector<bool> viz(nrNoduri + 1);
    stack<int> noduriSortTop;
    vector<int> rezultat;

    for (int i = 1; i <= nrNoduri; i++){
        if (!viz[i]) {
            dfsSortTop(i, noduriSortTop, viz);
        }
    }

    while (!noduriSortTop.empty()){
        int nr = noduriSortTop.top();
        rezultat.push_back(nr);
        noduriSortTop.pop();
    }

    return rezultat;
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

/*
 Componente biconexe, graf Neorientat
 COMPLEXITATE => O(M + N)

 Ce am folosit:
    => vector<bool> viz - vector de vizitate
    => vector<int> nivel - nivel[i] = x inseamna ca nodul i are nivelul = x (nivelul pe care se afla nodul)
    => vector<int> nivelMin - nivelMin[i] = x inseamna ca nodul i are nivelul minim = x (nivelul minim la care poate ajunge nodul prin intermediul vecinilor)
    => stack<int> s - pe stiva punem nodurile vizitate, delimitam nodurile din fiecare componenta
    => vector<vector<int>> compBiconexe - vectorul rezultat, tine toate componentele biconexe

 Descriere algoritm:
    => 2 functii -> componenteBiconexe => declararare date
                                       => apel functie dfsCompBiconexe pentru nodul = 1 si nivelCrt = 1
                                       => returnez vectorul compBiconexe
                 -> dfsCompBiconexe => vizitam nodul curent
                                    => il adaugam pe stiva
                                    => nivelul nodului devine nivelul curent
                                    => nivelulMinim devine initial tot nivelul curent
                                    => pentru fiecare vecin al nodului curent
                                        -> daca vecinul nu este vizitat
                                            --> apelam dfsCompBiconexe pentru nodul = vecin si nivelCrt = nivelCrt + 1
                                            --> cand ne intoarcem din recursivitate recalculam nivelul minim al nodului curent
                                            --> daca are un vecin care nu urca mai sus de el in arbore => nod critic => incepe o componenta biconexa
                                            --> luam de pe stiva nodurile pana cand top = vecin + le adaugam in componenta curenta
                                            --> adaugam componenta curenta la vectorul rezultat
                                        -> daca vecinul este vizitat
                                          --> inseamna ca am gasit o muchie de intoarcere => actualizam nivelMin
*/

vector<vector<int>> Graf :: componenteBiconexe() {
    vector<bool> viz(nrNoduri + 1);
    vector<int> nivel(nrNoduri + 1);
    vector<int> nivelMin(nrNoduri + 1);
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

/*
 Componente Tare Conexe, graf Orientat
 COMPLEXITATE => O(M + N)

 Ce am folosit:
    => vector<bool> viz - vector de vizitate
    => vector<int> nivel - nivel[i] = x inseamna ca nodul i are nivelul = x (nivelul pe care se afla nodul)
    => vector<int> nivelMin - nivelMin[i] = x inseamna ca nodul i are nivelul minim = x (nivelul minim la care poate ajunge nodul prin intermediul vecinilor)
    => stack<int> s - pe stiva punem nodurile vizitate, delimitam nodurile din fiecare componenta
    => vector<bool> inStiva -> inStiva[i] = true/false => nodul i este sau nu pe stiva
    => vector<vector<int>> ctc - vectorul rezultat, tine toate componentele tare conexe

 Descriere algoritm:
    => 2 functii -> componenteTareConexe => declararare date
                                         => parcurgem toate nodurile si daca nodul curent nu e vizitat => apelam dfsCompTareConexe pentru nodul respectiv
                                         => returnam vectoul ctc - vector rezultat
                 -> dfsCompTareConexe => parcugem vecinii nodului curent
                                      => daca vecinul nu este vizitat
                                         --> apelam recursiv functia dfsCompTareConexe pentru vecin si cand se intoarce din recursivitate reactualizeaza nivelul minim al nodului
                                      => daca vecinul este vizitat
                                         --> daca nu e pe stiva => apartine altei componenete tare conexe
                                         --> daca e pe stiva => actualizam nivelul minim al nodului curent
                                      => dupa ce parcurgem vecinii, daca nivelul nodului curent este agla cu nivelul minim al sau => radacina componentei tare conexe
                                      => adaugam nodurile de pe stiva in componenta curenta pana cand topul este nodul curent
                                      => adaugam componenta tare conexa la vectorul rezultat
*/

vector<vector<int>> Graf :: componenteTareConexe() {
    vector<bool> viz(nrNoduri + 1);
    vector<int> nivel(nrNoduri + 1);
    vector<int> nivelMin(nrNoduri + 1);
    stack<int> s;
    vector<bool> inStiva(nrNoduri + 1);
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

//Count Sort

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

/*
 Havel Hakimi => Daca se poate construi un graf Neorientat stiind vectorul de grade
 COMPLEXITATE => SORT STL => O(N^2 * log N)
              => CountSort => O(N^2) ( O(N * (N + MAX)), dar MAX -> N pt ca gradul unui nod este < nrNoduri)

 Ce am folosit:
    => vector<int> grade => grade[i] = x (nodul i  are gradul interior = x) (acest vector este constant pentru ca nu vreau sa-l modific)
    => prelucram vectorul copieGrade

 Descriere algoritm:
    => conditii de verificare => daca maximul din vectorul de grade e mai mare decat n - 1 =>  nu se poate reprezenta graful
                              => suma gradelor este impara => nu se poate reprezenta graful
    => intr-un loop infinit la fiecare pas => sortez vector de grade
*/

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

auto cmp = [](const pair<int, int>& p1, const pair<int, int>& p2)
{
    return p1.second > p2.second;
};

/*
 Dijkstra, graf Orientat (MUCHII CU COST POZITIV)
 => lungimea minima a drumului de la nodul 1 la fiecare din nodurile 2, 3, ..., N-1, N

 COMPLEXITATE => CU VECTOR O(N^2)
              => CU HEAP O(M * log N)
              => CU HEAP FIBONACCI O(N * log N + M)

 Ce am folosit:
    =>  vector<bool> viz - vector de vizitate
    =>  vector<int> distante - vectorul rezultat
    =>  priority_queue<pair<int, int>, vector<pair<int,int>>, decltype(cmp)> minHeap(cmp)

 Descriere algoritm:
    => initial -> in minHeap avem perechea (nodul de start = 1, dist = 0)
    => cat timp heap-ul nu este gol, luam front-ul
        -> daca front-ul nu a fost vizitat:
                --> il vizitam
                --> ii parcurgem vecinii
                --> daca gasim un vecin nevizitat caruia ii putem actualiza distanta atunci:
                    => actualizam distanta
                    => dam push in heap
     - daca front-ul a fost vizitat => continue (nu calculcam acelasi lucru de mai multe ori)
    => returnam vectorul de distante
*/

vector<int> Graf::Dijkstra(const int &s) {
    vector<bool> viz(nrNoduri + 1);
    vector<int> distante(nrNoduri + 1, INF);
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

/*
 Bellman-Ford, graf Orientat  => Drumuri minime de la un varf s dat la celelalte varfuri
    => arcele pot sa aiba costuri pozitive sau negative
    => daca exista circuite de cost negativ, algorimul le va detecta => nu are solutie
 COMPLEXITATE => O(M * N)

 Ce am folosit:

 Descriere algoritm:
    => folosim o coada in care punem nodurile actualizate anterior, initial e doar nodul de start in coada
        => cat timp coada nu e goala
            => pentru nodul din front => daca a fost actualizat de N ori => circuit negativ
                                      => daca nu, parcurg vecinii si vedem daca putem actualiza distanta => daca da, actualizam distanta si crestem numarul de relaxari
                                      => daca nodul vecin nu este in coada, il adaugam
    OBS:
    => relaxam muchiile de N - 1 ori
    => daca a N-a oara se tot relaxeaza muchii => circuit de cost negativ
    Optimizare => la un pas este suficient sa relaxam acele arce a caror varfuri au fost actualizate anterior
*/

bool Graf::BellmanFord(vector<int> &dist, const int &s) {
    queue<int> coadaNoduri;
    vector<int> nrRelaxari(nrNoduri + 1, 0);
    vector<bool> inCoada(nrNoduri + 1, false);
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

/*
 Arbori partiali de cost minim, graf conex Neorientat
 COMPLEXITATE:
    PRIM:
               => CU VECTOR DE VIZITATE O(N^2)
               => CU MEMORAREA VARFURILOR INTR-UN MIN HEAP O(N * log N)
               => CU VECTOR DE VIZITATE O(N^2)
    KRUSKAL:
               => SORTARE O(M * log N)
               => CU VECTOR DE REPREZENTANTI O(M * log N + N^2)
               => STRUCTURI PENTRU MULTIMI DISJUNCTE UNION/FIND 0(M * log N)

 Ce am folosit:

 Descriere algoritm:
    => tin muchiile intr-un vector de tuplu  (vector<tuple<int, int, int>> muchiiCost) cost, nod1, nod2
    => sortez crescator vectorul de muchii in functie de cost
    => creez o padure de multimi disjuncte
    => parcurg vectorul de muchii
        => daca nodurile muchiei curente au reprezentanti diferiti => reunim
                                                                   => crestem costul total cu costul muchiei
                                                                   => adaugam muchia la rezultat
    => dupa ce am parcurs toate muchiile => returnez costul total
    (OBS => ne putem opri cand au fost selectate N-1 muchii)

*/

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
 Floyd-Warshall  => Graf ORIENTAT cu n noduri, memorat prin matricea ponderilor
     (ponderile pot fi si negative dar NU exista circuite cu cost negativ in G)
 COMPLEXITATE O(n^3)

 Ce am folosit:

    => Cerinta: Pentru oricare doua varfuri x si y ale lui G, sa se determine distanta minima de la x la y

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

/*
 Darb - Diametrul unui arbore
 COMPLEXITATE => O(N)  OBS => O(N + M), dar fiind arbore M = N - 1

    Diametrul unui arbore reprezintă lungimea drumului (ca numar de noduri) intre cele mai departate două frunze.

 Ce am folosit:

 Descriere algoritm:
    => prima parcurgere BFS ne da un capat al drumului de lungime maxima in arbore (nodul cu distanta maxima)
    => a doua parcurgere BFS (pornind din nodul gasit cu prima parcurgere) ne da al doilea capat al diametrului (nodul cu distanta maxima)
    => diametrul o sa fie egal cu distanta maxima, dupa al doilea BFS, + 1 (+ 1 pentru ca vrem diametrul ca numar de noduri, iar BFS-UL ne da nr de muchii intre 2 noduri)
 */

int Graf::diametruArbore() {
    vector<int> distante = bfs();

    int nod1 = nodCuDistantaMaxima(distante);
    distante = bfs(nod1);
    int nod2 = nodCuDistantaMaxima(distante);
    int diametru = distante[nod2] + 1;

    return diametru;
}

/*
 MaxFlow - Graf Orientat
 COMPLEXITATE => O(N * M^2)

 Ce am folosit:

 Descriere algoritm:
    => Incercăm să trimitem marfă (flux) de la vârful sursă s la destinația t
    => Folosim 2 matrice -> de capacitate (capacitatea maxima care SE POATE TRIMITE pe muchia i j)
                         -> de flux (cantitatea care a fost dusa pana la momentul respectiv)
    => cat timp gasim un drum (cu BFS)
        => adaugam la fluxulMaxim (rezultat) fluxulMinim (gasit pe drumul curent)
        => actualizam fluxul pentru fiecare muchie din drumul curent => daca e muchie directa adaug flux (cat flux am dus)
                                                                     => daca e muchie inversa scad flux (cat flux mai pot duce)
    => daca nu mai exista drumuri => return fluxMaxim
*/

int Graf::maxFlow(const int &s, const int&f) {
    vector<vector<int>> capacitate(nrNoduri + 1, vector<int>(nrNoduri + 1));
    vector<vector<int>> flux(nrNoduri + 1, vector<int>(nrNoduri + 1));
    int fluxMaxim = 0;
    vector<int> tata(nrNoduri + 1, 0);

    for(int i = 1; i <= nrNoduri; i++){
        int nrVecini = listaAdsiCosturi[i].size();
        for(int j = 0; j < nrVecini; j++){
            int nodVecin = listaAdsiCosturi[i][j].first;
            capacitate[i][nodVecin] = listaAdsiCosturi[i][j].second;
        }
    }

    while(construiesteDrum(s, f, tata, capacitate, flux)){
        int fluxDrum = INF;
        int nod = f;
        while (nod != s){
            fluxDrum = min(fluxDrum, capacitate[tata[nod]][nod] - flux[tata[nod]][nod]);
            nod = tata[nod];
        }

        fluxMaxim += fluxDrum;

        nod = f;
        while (nod != s){
            int nodTata = tata[nod];
            flux[nodTata][nod] += fluxDrum;
            flux[nod][nodTata] -= fluxDrum;
            nod = tata[nod];
        }
    }

    return fluxMaxim;

}

bool Graf::construiesteDrum(const int &start, const int &final, vector<int> &tata, vector<vector<int>> &capacitate, vector<vector<int>> &flux) {
   fill(begin(tata), begin(tata) + nrNoduri + 1, 0);
    vector<bool> viz(nrNoduri + 1, false);
    queue<int> coada;

    coada.push(start);
    tata[start] = -1;
    viz[start] = true;

    while (!coada.empty()){
        int nod = coada.front();
        int nrVecini = listaAdsiCosturi[nod].size();
        for (int i = 0; i < nrVecini; i++){
            int vecinCrt = listaAdsiCosturi[nod][i].first;
            if (!viz[vecinCrt] and capacitate[nod][vecinCrt] - flux[nod][vecinCrt] > 0){
                if (vecinCrt == final){
                    tata[vecinCrt] = nod;
                    return true;
                }
                coada.push(vecinCrt);
                viz[vecinCrt] = true;
                tata[vecinCrt] = nod;
            }
        }
        coada.pop();
    }
    return false;
}

/* -------------------------------------------------------------------------------------------------------------------*/

/*
 DFS, graf Neorientat
 Cerinta => Sa se determine numarul componentelor conexe ale grafului.
*/

void rezultatDfs(){
//    https://www.infoarena.ro/problema/dfs
    ifstream fin ("dfs.in");
    ofstream fout ("dfs.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    Graf G(noduri, muchii, false, false);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    fout << G.nrCmpConexe() << "\n";
    fin.close();
    fout.close();
}

/*
 BFS, graf Orientat
 Cerinta => Fiind dat un nod S, sa se determine, pentru fiecare nod X,
            numarul minim de arce ce trebuie parcurse pentru a ajunge din nodul sursa S la nodul X.
*/

void rezultatBfs(){
//    https://www.infoarena.ro/problema/bfs
    ifstream fin ("bfs.in");
    ofstream fout ("bfs.out");
    int noduri, muchii, s;
    fin >> noduri >> muchii >> s;
    Graf G(noduri, muchii, true, false);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    vector<int> distanteMinBfs = G.bfs(s);
    for (int i = 1; i <= noduri; i++){
        fout << distanteMinBfs[i] <<" ";
    }
    fin.close();
    fout.close();
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
    Graf G(noduri, muchii, true, false);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAd(n1, n2);
    }
    vector<int> noduriSortTop =  G.sortareTopologica();
    for(int i = 0; i < noduri; i++){
        fout << noduriSortTop[i] << " ";
    }
    fin.close();
    fout.close();
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
    Graf G(noduri, muchii, false, false);
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
    fin.close();
    fout.close();
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
    Graf G(noduri, muchii, true, false);
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
    fin.close();
    fout.close();
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
    fin.close();
    fout.close();
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
    Graf G(noduri, muchii, true, true);
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
    fin.close();
    fout.close();
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
    Graf G(noduri, muchii, true, true);
    for (int i = 0; i < muchii; i++){
        int n1, n2, cost;
        fin >> n1 >> n2 >> cost;
        G.adaugaInListaAdSiCosturi(n1, n2, cost); //  pe Orientat
    }
    vector<int> distante(noduri + 1, INF);
    // bool eCircuit = G.BellmanFord(s);
    bool eCircuit = G.BellmanFord(distante);
    if (!eCircuit){
        for (int i = 2; i <= noduri; i++){
            fout << distante[i] << " ";
        }
    } else {
        fout << "Ciclu negativ!" << "\n";
    }
    fin.close();
    fout.close();
}

/*
 Paduri de multimi disjuncte
  Cerinta => Se dau N multimi de numere, initial fiecare multime i continand un singur element,
             mai exact elementul i. Asupra acestor multimi se pot face 2 tipuri de operatii, astfel:
    * operatia de tipul 1: se dau doua numere naturale x si y, intre 1 si N.
        Se cere sa reuneasca multimile in care se afla elementul x, respectiv elementul y
        (se garanteaza ca x si y nu se vor afla in aceeasi multime)
    * operatia de tipul 2: se dau doua numere naturale x si y, intre 1 si N.
        Se cere sa afiseze DA daca cele 2 elemente se afla in aceeasi multime, respectiv NU in caz contrar.
*/

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
    fin.close();
    fout.close();
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
    fin.close();
    fout.close();
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
    fin.close();
    fout.close();
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
    Graf g(nrNoduri, nrNoduri - 1, false, false);
    for(int i = 0; i < nrNoduri; i ++){
        int n1, n2;
        fin >> n1 >> n2;
        g.adaugaInListaAd(n1, n2);
    }
    fout << g.diametruArbore();
    fin.close();
    fout.close();
}

/*
    MaxFlow - Graf Orientat
*/

void rezultatMaxFlow(){
//    https://www.infoarena.ro/problema/maxflow
    ifstream fin("maxflow.in");
    ofstream fout("maxflow.out");
    int nrNoduri, nrMuchii, start, final;
    start = 1;
//    fin >>  nrNoduri >> nrMuchii >> start >> final;
    fin >> nrNoduri >> nrMuchii;
    Graf g(nrNoduri, nrMuchii, true, true);
    for(int i = 0; i < nrMuchii; i ++){
        int n1, n2, cost;
        fin >> n1 >> n2 >> cost;
        g.adaugaInListaAdSiCosturi(n1, n2, cost);
    }
    fout << g.maxFlow(start, nrNoduri);
    fin.close();
    fout.close();
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
    rezultatMaxFlow();

    return 0;
}

