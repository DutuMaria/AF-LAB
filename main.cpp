#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stack>
#include <bitset>
#include <queue>
#include <numeric>

using namespace std;

const int NMAX = 100001;
ifstream fin ("f.in");
ofstream fout ("f.out");

class Graf{
private:
    int nrNoduri, nrMuchii;
    vector<int> listaAd[NMAX];
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
public:
    Graf(int nrNoduri, int nrMuchii) : nrNoduri(nrNoduri), nrMuchii(nrMuchii) {};
    void adaugaInListaAdOrientat(int nod1, int nod2);
    void adaugaInListaAdNeorientat(int nod1, int nod2);
    int nrCmpConexe();
    void afisareDistante(int start);
    void sortareTopologica();
    vector<vector<int>> componenteBiconexe();
    vector<vector<int>> componenteTareConexe();
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

void Graf::afisareDistante(int nodStart) {
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
        accumulate(grade.begin(), grade.end(), 0) % 2  == 1)
        fout << "Nu\n";
    else {
        while (true) {
            sort(grade.begin(), grade.end(), greater<>());
            int grad_crt = grade[0];
            if (grad_crt == 0){
                fout << "Da\n";
                return;
            }
            grade.erase(grade.begin());
            if (grad_crt > grade.size()){
                fout << "Nu\n";
                return;
            }
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


int main() {
    int noduri, muchii, s;
    fin >> noduri >> muchii;
    //fin >> noduri >> muchii >> s; // pt bfs
    Graf G(noduri, muchii);
    for (int i = 0; i < muchii; i++){
        int n1, n2;
        fin >> n1 >> n2;
        G.adaugaInListaAdOrientat(n1, n2);
        //G.adaugaInListaAdNeorientat(n1, n2);
    }

    //HavelHakimi();

    return 0;
}

