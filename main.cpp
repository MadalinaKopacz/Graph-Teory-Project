#include <bits/stdc++.h>
#define INF 0x3f3f3f3f
using namespace std;

struct muchie {   //structura pentru a salva toate informatiile legate de muchii
    int x, y; //noduri
    int cost;
    int capacitate;
    int indMuchie;
};

bool operator<(const pair<int, pair<int,int>> &a, const pair<int, pair<int, int>>& b)  //supraincarcare pentru sortare
{
    return a.first > b.first;
}

///Clasa pentru reprezentare de grafuri
class Graf
{
private:
    int nrNoduri;
    int nrMuchii;
    bool orientat; //1 = orientat, 0 = neorientat
    bool cuCost; //1 = cu cost, 0 = fara cost
    bool cuCapacitate; //1 = cu capacitati pt fluxuri, 0 = fara capacitati pt fluxuri
    bool multigraf; //0 = fara muchii care intra si ies in acelasi nod si mai multe muchii intre x si y dist, 1 altfel

    vector<vector<muchie>> lisAdiacenta;

    void DFS(const int nodStart, vector<bool> &vizitat, stack<int> &stiva);    //dfs pentru nr de comp con si pentru sortare topolofica
    void ctcDFS(int x, vector<int> &disc, vector<int> &low, stack<int> &stiva, vector<bool> &gasitStiva, vector<vector<int>> &rez);  //dfs special pt comp tare conexe
    void critConDFS(int x, vector<bool> &viz, vector<int> &disc, vector<int> &low, vector<int> &tata, vector<vector<int>> &rez);     //dfs pentru critical connections
    void biconexeDFS(int x, int tata, vector<int> &disc, vector<int> &low, stack<int> &comp, vector<vector<int>> &rez, int &nrCompBiconexe, vector<bool> &viz);  //dfs pt elemente biconexe

    static bool bfs_HopcroftKarp(vector<vector<int> > &lisAdiacenta, vector<int> &stanga, vector<int> &dreapta, vector<int> &distanta);
    static bool dfs_HopcroftKarp(const int start, vector<vector<int> > &lisAdiacenta, vector<int> &stanga, vector<int> &dreapta, vector<int> &distanta);
public:
    Graf(const int nrNoduriDat, const int nrMuchiiDat, const bool orientat, const bool cuCostDat, const bool cuCapacitate, const bool multigraf);  //constructor
    Graf(const Graf &grafDat);  //constructor de copiere
    ~Graf(); //destructor

    int getNrNoduri(); //getter pt nr de noduri
    int getNrMuchii(); //getter pt nr de muchii

    Graf& operator= (const Graf &grafDat);  //supraincarcare egal
    friend std::ostream& operator<<(std::ostream &out, const Graf &grafDat);  //afisare graf cu toate info despre muchii
    void citireGraf(istream &in);  //functie pentru citire graf
    void citireGrafIndex(istream &in);  //functie pentru citire graf cu nr muchiei ca index


    vector<int> BFS(const int nodStart);   //returneaza dist min de la un nod la celelalte
    int nrCompConexe();                    //returneaza nr de comp con
    vector<vector<int>> Biconex();           //nr si componentele biconexe
    vector<vector<int>> CTC();             //componentel tare conexe
    stack<int> sortareTopologica();        //ret stiva sortarii topologice
    static bool existaGraf(vector<int> &grade);    //0 daca un sir de nr nu poate reprezenta gradele unui graf, 1 altfel
    vector<vector<int>> CriticalConnections();     //muchii critice

    static int findDis(int elem, vector<pair<int,int>> &multimi);    //functie sa vedem daca 2 element sunt in aceeasi multime
    static void unionDis(const int x, const int y, vector<pair<int, int>> &multimi);  //unim doua multimi
    vector<int> Dijkstra(const int start);      //cost minim de la nodul de start pana la celelalte
    vector<int> BellmanFord(const int start);   //cost minim de la nodul de start la celelalte cu much cu costuri de -1
    vector<pair<int, pair<int, int>>> findLisMuchii();  //functie pentru a face lista de muchii din lista de adiacenta
    pair<int, vector<pair<int, int>>> Kruskal(vector<pair<int, pair<int, int>>> muchii);  //cost minim si muchiile din apm

    int maxFlow(int S, int T, const int capacitateMax);  //flux maxim intr un graf
    void BFS_final(const int start, int &ultim, int &distanta);  //diametrul unui arbore
    static vector<vector<int>> royFloyd(const vector<vector<int>>& matrice, const int costMaxim);  //matricea drumurilor minime

    vector<int> cicluEuler(); //ret vector care contine un ciclu eulerian daca exista sau cu -1 in caz contrar
    static pair<int, vector<int>> cuplajMaxim(vector<vector<int>> &lisAdiacenta, const int N, const int M);
    int hamilton();
};

Graf :: Graf(const int nrNoduriDat, const int nrMuchiiDat, const bool orientatDat, const bool cuCostDat, const bool cuCapacitateDat, const bool multigrafDat)  //constructor parametrizat
{
    nrNoduri = nrNoduriDat;
    nrMuchii = nrMuchiiDat;
    orientat = orientatDat;
    cuCost = cuCostDat;
    cuCapacitate = cuCapacitateDat;
    multigraf = multigrafDat;
    vector<muchie> aux(1, {-1, -1, -1, -1, -1});      //initializam prima pozitie cu -1 pentru ca indexarea e de la 1
    for(int i = 0; i <= nrNoduri+1; ++i) {
        lisAdiacenta.push_back(aux);
    }
}

Graf :: Graf(const Graf & grafDat)   //constructor copiere
{
    nrNoduri = grafDat.nrNoduri;
    nrMuchii = grafDat.nrMuchii;
    orientat = grafDat.orientat;
    cuCost = grafDat.cuCost;
    cuCapacitate = grafDat.cuCapacitate;
    multigraf = grafDat.multigraf;
    lisAdiacenta = grafDat.lisAdiacenta;
}

Graf :: ~Graf()
{
    lisAdiacenta.clear();
}

int Graf :: getNrNoduri()    //getter numar noduri
{
    return nrNoduri;
}

int Graf :: getNrMuchii()    //getter numar muchii
{
    return nrMuchii;
}

int Graf :: findDis(int elem, vector<pair<int,int>> &multimi)   //gasim reprezentantul multimii respectice
                                                              //care e de fapt radacina arborelui
{
    int radacina = elem;
    while(multimi[radacina].first != radacina) {
        radacina= multimi[radacina].first;
    }
    return radacina;
}

void Graf :: unionDis(const int x, const int y, vector<pair<int, int>> &multimi)
{
    int radX = findDis(x, multimi), radY = findDis(y, multimi);     //vedem care sunt reprezentantii fiecarei multimi

    if(multimi[radX].second > multimi[radY].second) {     //legam arborede mic de cel mai mare
        multimi[radX].second = multimi[radX].second + multimi[radY].second;
        multimi[radY].first = radX;
        multimi[radY].second = multimi[radX].second;
    } else {
        multimi[radY].second = multimi[radY].second + multimi[radX].second;
        multimi[radX].first = radY;
        multimi[radX].second = multimi[radY].second;
    }
}

///SUPRAINCARCARI OPERATORI
Graf& Graf :: operator= (const Graf &grafDat)//supraincarcare operator egal
{
    if(this != &grafDat) {
        this->nrNoduri = grafDat.nrNoduri;
        this->nrMuchii = grafDat.nrMuchii;

        this->orientat = grafDat.orientat;
        this->cuCost = grafDat.cuCost;
        this->cuCapacitate = grafDat.cuCapacitate;
        if(!this->lisAdiacenta.empty())
            lisAdiacenta.clear();

        this->lisAdiacenta = grafDat.lisAdiacenta;
    }
    return *this;
}

ostream& operator<< (std::ostream& out, const Graf& grafDat)  //supraincarcare operator afisare
{
    out << "Numarul de noduri ale grafului este: " << grafDat.nrNoduri << "\n";
    out << "Numarul de muchii ale grafului este: " << grafDat.nrMuchii << "\n";

    for (int i = 1; i <= grafDat.nrNoduri; ++i) {
        out << i << ": ";
        for (int j = 1; j < grafDat.lisAdiacenta[i].size(); ++j) {
            out <<"{" << grafDat.lisAdiacenta[i][j].x << " " << grafDat.lisAdiacenta[i][j].y<<" "<< grafDat.lisAdiacenta[i][j].cost<<" " << grafDat.lisAdiacenta[i][j].capacitate;
            out << " " << grafDat.lisAdiacenta[i][j].indMuchie << "}";
        }
        out<<"\n";
    }

    return out;
}

///METODE CITIRI

void Graf :: citireGraf(istream &in)  //citim de la tastatura sau din fisier un graf neorientat
{
    int x, y, cost, capacitate, multigraf;

    for(int i = 1; i <= nrMuchii; ++i) {
        in >> x >> y;
        muchie m;
        m.x = x;
        m.y = y;
        if (cuCost) {
            in >> cost;
            m.cost = cost;
        } else {
            m.cost = -1;
        }

        if (cuCapacitate) {
            in >> capacitate;
            m.capacitate = capacitate;
        } else {
            m.capacitate = -1;
        }

        if (multigraf) {
            m.indMuchie = i;
        }
        else {
            m.indMuchie = -1;
        }
        lisAdiacenta[m.x].push_back(m);
        if (!orientat) {
            swap(m.x, m.y);
            lisAdiacenta[m.x].push_back(m);
        }

    }
}

///Afisare distante BFS
vector<int> Graf :: BFS(const int nodStart)     //parcurgere in latime
{
    queue<int> Q;
    int x, y;
    vector<int> distanta(nrNoduri+1, -1);

    distanta[nodStart] = 0;
    Q.push(nodStart);

    while(!Q.empty()) {
        x = Q.front();   //eliminam nodul curent din coada
        Q.pop();

        for (int i = 1; i < lisAdiacenta[x].size(); ++i)
            if (distanta[lisAdiacenta[x][i].y] == -1) {              //daca nu a fost vizitat vecinul
                distanta[lisAdiacenta[x][i].y] = distanta[x] + 1;
                Q.push(lisAdiacenta[x][i].y);
            }
    }

    return distanta;
}

///Aflare nr componente conexe
void Graf:: DFS(const int x, vector<bool> &vizitat, stack<int> &stiva)   //parcurgere in adancime
{
    vizitat[x] = 1;

    for(int i = 1; i < lisAdiacenta[x].size(); i++)
        if(vizitat[lisAdiacenta[x][i].y] == 0)    //daca vecinul sau nu a fost vizitat, mergem la el
            DFS(lisAdiacenta[x][i].y, vizitat, stiva);
    stiva.push(x);
}

int Graf :: nrCompConexe()   //numar de componente conexe facut cu DFS
{
    stack<int> stiva;       //implementata formal, pentru a putea utiliza DFS ul si la sortareTopologica()
    int nrCompCon = 0;
    vector<bool> vizitat(nrNoduri+1, 0);

    for(int i = 1; i <= nrNoduri; ++i )
        if(!vizitat[i]) {
            nrCompCon++;
            DFS(i, vizitat, stiva);
        }

    return nrCompCon;
}

///Biconex
//https://www.geeksforgeeks.org/biconnected-components/
void Graf :: biconexeDFS(int x, int tata, vector<int> &disc, vector<int> &low, stack<int> &comp, vector<vector<int>> &rez, int &nrCompBiconexe, vector<bool> &viz)
{
    viz[x] = 1;
    disc[x] = disc[tata] + 1;
    low[x] = disc[x];

    for (int i = 1; i < lisAdiacenta[x].size(); ++i)
        if (lisAdiacenta[x][i].y != tata) {
            if (!viz[lisAdiacenta[x][i].y]) {
                comp.push(lisAdiacenta[x][i].y);
                biconexeDFS(lisAdiacenta[x][i].y, x, disc, low, comp, rez, nrCompBiconexe, viz);

                low[x] = min(low[x], low[lisAdiacenta[x][i].y]);   //conexiune a fiului cu un stramos al lui x
                if (disc[x] <= low[lisAdiacenta[x][i].y]) {   //am gasit o muchie critica
                    nrCompBiconexe++;
                    rez.push_back(vector<int>(1));
                    comp.push(x);
                    while (!comp.empty() && comp.top() != lisAdiacenta[x][i].y) {
                        rez[nrCompBiconexe-1].push_back(comp.top());
                        comp.pop();
                    }
                    if (!comp.empty()) {   //adaugam si radacina componentei biconexe
                        rez[nrCompBiconexe-1].push_back(comp.top());
                        comp.pop();
                    }
                }
            } else if (disc[lisAdiacenta[x][i].y] < low[x]) //daca e vizitat si are timpul de disc mai mic, atunci modificam timp min x
                low[x] = disc[lisAdiacenta[x][i].y];
        }
}

vector<vector<int>> Graf :: Biconex()
{
    vector<int> disc(nrNoduri+1, -1);
    vector<int> low(nrNoduri+1, -1);
    vector<bool> viz(nrNoduri + 1, 0);
    stack<int> comp;
    vector<vector<int>> rez;

    int nrCompBiconexe = 0;

    disc[1] = 1;
    biconexeDFS(2, 1, disc, low, comp, rez, nrCompBiconexe, viz);

    return rez;
}

///Componente Tare Conexe - Algoritmul lui Tarjan
//https://www.geeksforgeeks.org/tarjan-algorithm-find-strongly-connected-components/
void Graf::ctcDFS(int x, vector<int> &disc, vector<int> &low, stack<int> &stiva, vector<bool> &gasitStiva, vector<vector<int>> &rez)
{
    static int timp = 0;


    //low-> cel mai mic timp de descoperire al unui nod dintr o comp con
    //disc-> timpul de descoperire al nodului
    //disc  = low pt radacina unui arb de comp con
    disc[x] = low[x] = ++timp;
    stiva.push(x);
    gasitStiva[x] = 1;


    for (int i = 1; i < lisAdiacenta[x].size(); i++) {
        if (disc[lisAdiacenta[x][i].y] == -1)    {
            ctcDFS(lisAdiacenta[x][i].y, disc, low, stiva, gasitStiva, rez);
            low[x] = min(low[x], low[lisAdiacenta[x][i].y]);   //conexiune a fiului cu un stramos al lui x
        } else if (gasitStiva[lisAdiacenta[x][i].y] == 1)
            low[x] = min(low[x], disc[lisAdiacenta[x][i].y]);  //muchie de intoarcere
    }

    vector<int> compNoua;
    if (low[x] == disc[x]) {       //daca am gasit nod de start pt componenta tare conexa scoatem nodurile din stiva
        while (stiva.top() != x) {
            compNoua.push_back(stiva.top());
            gasitStiva[stiva.top()] = 0;
            stiva.pop();
        }

        compNoua.push_back(stiva.top());
        gasitStiva[stiva.top()] = 0;
        rez.push_back(compNoua);
        stiva.pop();
    }
}

vector<vector<int>> Graf :: CTC()
{
    vector<int> disc(nrNoduri+1, -1);
    vector<int> low(nrNoduri+1, -1);
    stack<int> stiva;
    vector<bool> gasitStiva(nrNoduri+1, 0);
    vector<vector<int>> rez;

    for (int i = 1; i <= nrNoduri; ++i) {   //aflam daca sunt mai multe
        if (disc[i] == -1)
            ctcDFS(i, disc, low, stiva, gasitStiva, rez);
    }

    return rez;
}

stack<int> Graf :: sortareTopologica()  //afisam o sortare topologica posibila
{
    stack<int> stiva;
    vector<bool> vizitat(nrNoduri+1, 0);

    for(int i = 1; i <= nrNoduri; ++i )
        if(!vizitat[i])
            DFS(i, vizitat, stiva);

    return stiva;

}

///algoritmul Havel-Hakimi pt a afla daca o secventa de numere poate forma graf simplu
//https://www.geeksforgeeks.org/find-if-a-degree-sequence-can-form-a-simple-graph-havel-hakimi-algorithm/
//O(n^2)
bool Graf :: existaGraf(vector<int> &grade)
{
    vector<int> copieGrade;
    copieGrade = grade;

    int S = 0;      //suma grade
    for (int i = 0; i < copieGrade.size(); ++i)
        S += copieGrade[i];
    if (S % 2 == 1) {    //suma gradelor intr-un graf neorientat e mereu para
        return 0;
    }

    sort(grade.begin(), grade.end(), greater<int>());  //luam gradele in ordine descrescatoare
    while (true) {

        if (grade[0] == 0) {
            return 1;
        }

        int fst = grade[0];
        grade.erase(grade.begin() + 0);
        if (fst > grade.size()) {  //daca avem un grad mai mare decat numarul de noduri ramase
            return 0;
        }

        for (int i = 0; i < fst; ++i) {
            grade[i]--;
            if (grade[i] < 0) {
                return 0;
            }
        }

        int j = fst;                            //reordonam descrescator
        for (int i = 0; i < fst; ++i)   {       //fiind deja ordonate descrescator
            if(j >= grade.size())               //atunci cand scadem 1 valorile nu vor mai fi ordonate decat daca sunt egale initial
                break;                           //si fst e mai mic decat numarul de valori egale intial
            if (grade[i] < grade[j]) {
                swap(grade[i], grade[j]);
                j++;
            }
        }
    }

}
///Muchii Critice
//https://www.geeksforgeeks.org/bridge-in-a-graph/
void Graf::critConDFS(int x, vector<bool> &viz, vector<int> &disc, vector<int> &low, vector<int> &tata, vector<vector<int>> &rez)
{
    static int timp = 0;

    viz[x] = 1;
    disc[x] = low[x] = ++timp;

    for (int i = 1; i < lisAdiacenta[x].size(); ++i) {
        if (!viz[lisAdiacenta[x][i].y])    {
            tata[lisAdiacenta[x][i].y] = x;
            critConDFS(lisAdiacenta[x][i].y, viz, disc, low, tata, rez);
            low[x] = min(low[x], low[lisAdiacenta[x][i].y]);    //conexiune a fiului cu un stramos al lui x

            if (low[lisAdiacenta[x][i].y] > disc[x]) {           //daca nu putem ajunge in fiu altfel (are timpul minim mai mare decat timpul de desc)
                vector<int> aux;
                aux.push_back(x);
                aux.push_back(lisAdiacenta[x][i].y);
                rez.push_back(aux);
            }
        } else if (lisAdiacenta[x][i].y != tata[x])  //daca e vizitat deja si nu e parinte lui x, modificam valoarea lui low[x] pt urmatoarele functii
            low[x] = min(low[x], disc[lisAdiacenta[x][i].y]);
    }
}

vector<vector<int>> Graf :: CriticalConnections()
{
    vector<bool> viz(nrNoduri+1, 0);
    vector<int> disc(nrNoduri+1);
    vector<int> low(nrNoduri+1);
    vector<int> tata(nrNoduri+1, -1);
    vector<vector<int>> rez;

    for (int i = 1; i <= nrNoduri; ++i) {   //pt fiecare comp
        if (!viz[i])
            critConDFS(i, viz, disc, low, tata, rez);
    }
    return rez;
}

vector<int> Graf :: Dijkstra(const int start)
{
    //minHeap ordonat dupa cost modelat cu priority queue
    priority_queue< pair<int, int>, vector <pair<int, int> >, greater<pair<int, int>>> pQueue;
    vector<int> distDij(nrNoduri + 1, INF);
    vector<bool> apartPQ(nrNoduri+ 1, 0);
    pQueue.push(make_pair(0, start)), distDij[start] = 0;

    while (!pQueue.empty()) {
        int x = pQueue.top().second;
        pQueue.pop();
        if (!apartPQ[x]) {     //sa nu vizitam acelasi nod de mai multe ori
            apartPQ[x] = 1;
            for (int i = 1; i < lisAdiacenta[x].size(); i++)
                if (distDij[lisAdiacenta[x][i].y] > distDij[x] + lisAdiacenta[x][i].cost) {     //verificam daca am gasit cost mai mic
                    distDij[lisAdiacenta[x][i].y] = distDij[x] + lisAdiacenta[x][i].cost;
                    pQueue.push(make_pair(distDij[lisAdiacenta[x][i].y], lisAdiacenta[x][i].y));
                }
        }
    }
    return distDij;
}

vector<int> Graf :: BellmanFord(const int start)
{
    vector<int> distBMF(nrNoduri + 1, INF);
    vector<int> viz(nrNoduri + 1, 0);
    vector<bool> apartCoada(nrNoduri + 1, false);
    queue<int> coada;
    int faraCiclNeg = 1;

    coada.push(start), apartCoada[start] = 1;
    distBMF[start] = 0;

    while (!coada.empty() && faraCiclNeg) {
        int x = coada.front();
        coada.pop();
        apartCoada[x] = 0;

        for (int i = 1; i < lisAdiacenta[x].size(); i++)
            if (distBMF[x] + lisAdiacenta[x][i].cost < distBMF[lisAdiacenta[x][i].y]) {
                distBMF[lisAdiacenta[x][i].y] = distBMF[x] + lisAdiacenta[x][i].cost;      //relaxam
                viz[lisAdiacenta[x][i].y]++;

                if(!apartCoada[lisAdiacenta[x][i].y]) {   //daca nu am mai parcurs nodul resp
                    coada.push(lisAdiacenta[x][i].y);
                    apartCoada[lisAdiacenta[x][i].y] = 1;
                }
                if(viz[lisAdiacenta[x][i].y] >= nrNoduri) {   //verificam daca avem ciclu
                    faraCiclNeg = 0;
                }
            }
    }
    if(!faraCiclNeg)
        distBMF.clear();
    return distBMF;
}

vector<pair<int, pair<int, int>>> Graf :: findLisMuchii()
{
    vector<pair<int, pair<int, int>>> lisMuchii(1);
    for(int i = 1; i < lisAdiacenta.size(); ++i)
        for(auto nod : lisAdiacenta[i])
            if(orientat || (!orientat && i < nod.y)) {
                pair<int, pair<int, int>> much;
                much.first = nod.cost;     /// put the cost first so we can sort by cost
                much.second.first = nod.x;
                much.second.second = nod.y;
                lisMuchii.push_back(much);
            }
    return lisMuchii;
}
pair<int, vector<pair<int, int>>> Graf :: Kruskal(vector<pair<int, pair<int, int>>> muchii)
{
    int N = nrNoduri, M = nrMuchii;
    vector<pair<int, int>> multimi(N + 1, {0, 0});
    vector<pair<int,int>> rezMuchii(1, {-1, -1});
    int costTotal = 0;

    for(int i = 1; i <= N; ++i) {
        multimi[i] = {i, 1};
    }

    sort(muchii.begin() + 1, muchii.end());    //sortam muchiile

    for(int i = 1; i <= M && rezMuchii.size() < N; ++i) {
        if( findDis(muchii[i].second.first, multimi) != findDis(muchii[i].second.second, multimi) ) {
            unionDis(muchii[i].second.first, muchii[i].second.second, multimi);
            costTotal = costTotal + muchii[i].first;
            rezMuchii.push_back({muchii[i].second.first, muchii[i].second.second});
        }
    }

    return {costTotal, rezMuchii};
}
int BFS_EK(vector<vector<int>> &capacitati, int S,int T, vector<int> &tati, vector<vector<int>> &flux,vector<bool> &viz, vector<vector<int>> &rezidual)
{
    tati.assign(capacitati.size(), 0);
    queue<int> coada;
    coada.push(S);

    tati[S] = -1;

    viz.clear();
    viz.resize(capacitati.size(), 0);
    viz[S] = 1;

    while(!coada.empty() && tati[T] == 0) {
        int x = coada.front();
        coada.pop();
        for(int v : rezidual[x]) {
            if(!viz[v] && capacitati[x][v] > flux[x][v]) {
                coada.push(v);
                tati[v] = x;
                viz[v] = 1;
            }
        }
    }
    return tati[T];
}

int Graf :: maxFlow( int S, int T, const int capacitateMax)
{

    vector<vector<int>> capacitati(lisAdiacenta.size(), vector<int>(lisAdiacenta.size(), 0));
    vector<vector<int>> flux(lisAdiacenta.size(), vector<int>(lisAdiacenta.size(), 0));
    vector<int> tati(lisAdiacenta.size(), 0);
    vector<bool> viz(lisAdiacenta.size(), 0);
    vector<vector<int>> rezidual(lisAdiacenta.size());

    for(int i = 1; i < lisAdiacenta.size(); ++i) {
        for(int j = 1; j < lisAdiacenta[i].size(); ++j) {
            capacitati[i][lisAdiacenta[i][j].y] = lisAdiacenta[i][j].capacitate;
            rezidual[i].push_back(lisAdiacenta[i][j].y);
            rezidual[lisAdiacenta[i][j].y].push_back(i);
        }
    }

    int fluxMaxim = 0;
    while(BFS_EK(capacitati, S, T, tati, flux, viz, rezidual) != 0) {
        for(int i : rezidual[T])
            if(viz[i] && flux[i][T] < capacitati[i][T]) {
                tati[T] = i;
                int fluxDrum  = capacitateMax;
                for(int j = T; j != S; j = tati[j]) {
                    int x = tati[j];
                    fluxDrum = min(fluxDrum, capacitati[x][j]-flux[x][j]);
                }
                if(fluxDrum > 0) {
                    for(int j = T; j != S; j = tati[j]) {
                        int x = tati[j];
                        flux[x][j] += fluxDrum;
                        flux[j][x] -= fluxDrum;

                    }
                    fluxMaxim += fluxDrum;
                }
            }
    }
    return fluxMaxim;
}

void Graf :: BFS_final(const int start, int &ultim, int &distanta)
{
    vector<bool> viz(nrNoduri+1, 0);
    queue<int> coada;

    vector<int> dist(nrNoduri+1, 0);
    coada.push(start);
    viz[start] = 1;
    dist[start] = 1;
    distanta = 0;
    while(!coada.empty()) {
        int x = coada.front();
        coada.pop();
        for(int i = 1; i < lisAdiacenta[x].size(); ++i) {
            if(!viz[lisAdiacenta[x][i].y]) {
                dist[lisAdiacenta[x][i].y] = dist[x] + 1;
                viz[lisAdiacenta[x][i].y] = 1;
                coada.push(lisAdiacenta[x][i].y);
            }
        }
    }
    for (int i = 1; i < lisAdiacenta.size(); ++i) {
        if(dist[i] > distanta) {
            distanta = dist[i];
            ultim = i;
        }
    }
}

vector<vector<int>> Graf :: royFloyd(const vector<vector<int>>& matrice, const int costMaxim)
{
    vector<vector<int>> distante = matrice;
    //initializare matrice distante
    for (int i = 1; i <  distante.size(); ++i) {
        for (int j = 1; j < distante.size(); j++) {
            if(!matrice[i][j] && i != j) {
                distante[i][j] = costMaxim;
            }
        }
    }

    //calculare distante
    for (int k = 1; k < distante.size(); k++) {
        for (int i = 1; i < distante.size(); i++) {
            for (int j = 1; j < distante.size(); j++) {
                if (distante[i][j] > distante[i][k] + distante[k][j] && i!=j) {
                    distante[i][j] =   distante[i][k] + distante[k][j];
                }
            }
        }
    }

    return distante;
}


vector<int> Graf :: cicluEuler()
{

    vector<int> ciclu;
    vector<bool> eliminat(nrMuchii+1, 0);

    for(int i = 1; i < lisAdiacenta.size(); ++i) {
        if((lisAdiacenta[i].size() - 1) %2 == 1) {
            ciclu.push_back(-1);
            ciclu.push_back(-1);
            return ciclu;
        }
    }

    stack<int> drumCurent;
    //incepem din primul nod
    drumCurent.push(1);


    while (!drumCurent.empty()) {
        int nodCurent = drumCurent.top();
        if (lisAdiacenta[nodCurent].size() > 1) {
            muchie nodUrm = lisAdiacenta[nodCurent].back();
            lisAdiacenta[nodCurent].pop_back();

            if(!eliminat[nodUrm.indMuchie])    {
                drumCurent.push(nodUrm.y);
                eliminat[nodUrm.indMuchie]=1;
            }
        } else {
            ciclu.push_back(nodCurent);
            drumCurent.pop();
        }
    }

    return ciclu;
}


bool Graf :: bfs_HopcroftKarp(vector<vector<int> > &lisAdiacenta, vector<int> &stanga, vector<int> &dreapta, vector<int> &distanta)
{
    queue<int> coada;

    for(int i = 1; i < stanga.size(); ++i) {
        if(!stanga[i]) {
            distanta[i] = 0;
            coada.push(i);
        }
        else {
             distanta[i] = INF;
        }
    }
    distanta[0] = INF;

    while(!coada.empty()) {
        int x = coada.front();
        coada.pop();

        if(distanta[x] < distanta[0]) {
            for(int i = 0; i < lisAdiacenta[x].size(); ++i) {
                if(distanta[dreapta[lisAdiacenta[x][i]]] == INF) {
                    distanta[dreapta[lisAdiacenta[x][i]]] = distanta[x] + 1;
                    coada.push(dreapta[lisAdiacenta[x][i]]);
                }
            }
        }
    }

    return (distanta[0] != INF);
}

bool Graf :: dfs_HopcroftKarp(const int start, vector<vector<int> > &lisAdiacenta, vector<int> &stanga, vector<int> &dreapta, vector<int> &distanta)
{
    if(start) {
        for(int i = 0; i < lisAdiacenta[start].size(); ++i){
            if(distanta[dreapta[lisAdiacenta[start][i]]] == distanta[start]+1) {
                if(dfs_HopcroftKarp(dreapta[lisAdiacenta[start][i]], lisAdiacenta, stanga, dreapta, distanta)){
                    dreapta[lisAdiacenta[start][i]] = start;
                    stanga[start] = lisAdiacenta[start][i];
                    return 1;
                }
            }
        }
        distanta[start] = INF;
        return 0;
    }
    return 1;
}

pair<int, vector<int>> Graf :: cuplajMaxim(vector<vector<int>> &lisAdiacenta, const int N, const int M)
{
    vector<int> stanga(N + 1, 0), dreapta(M + 1, 0);
    vector<int> distanta(N + 1);
    int rezultat = 0;

    while(bfs_HopcroftKarp(lisAdiacenta, stanga, dreapta, distanta)) {
        for(int i = 1; i <= N; ++i) {
            if(!stanga[i] && dfs_HopcroftKarp(i, lisAdiacenta, stanga, dreapta, distanta)) {
                rezultat++;
            }
        }
    }
    return {rezultat, stanga};
}

int Graf :: hamilton()
{
    vector<vector<int> > dpCost(1 << nrNoduri, vector<int>(nrNoduri, INF));
    int costTotal;

    dpCost[1][0] = 0;
    for (int i = 0; i < 1 << nrNoduri; ++i) { //toate ciclurile
		for (int j = 0; j < nrNoduri; ++j) {
			if (i & (1<<j)) {
				for (int k = 1; k < lisAdiacenta[j].size(); ++k){
					if (i & (1<<lisAdiacenta[j][k].y)) {
                        dpCost[i][j] = min(dpCost[i][j], dpCost[i ^ (1<<j)][lisAdiacenta[j][k].y] + lisAdiacenta[j][k].cost);
					}
				}
			}
        }
    }
    costTotal = INF;
	for (int i = 1; i < lisAdiacenta[0].size(); ++i) {
		costTotal = min(costTotal, dpCost[(1<<nrNoduri) - 1][lisAdiacenta[0][i].y] + lisAdiacenta[0][i].cost);
    }
    if(costTotal != INF)
        return costTotal;
    else return -1;
}

//corpuri functii pt rezolvare probleme ajutatoare

Graf CitireGraf(string numeFisier, bool orientat, bool cuCost, bool cuCapacitate, bool multigraf)
{
    ifstream fin(numeFisier);

    int N, M;
    fin >> N >> M;

    Graf graf(N, M, orientat, cuCost, cuCapacitate, multigraf);
    graf.citireGraf(fin);
    fin.close();
    return graf;
}

void Rezolva_BFS(string fisierIntrare, string fisierIesire)  //https://www.infoarena.ro/problema/bfs - 100pc
{
    ifstream fin(fisierIntrare);
    ofstream fout(fisierIesire);

    int N, M, S;
    fin >> N >> M >> S;

    Graf graf(N, M, 1, 0, 0, 0);
    graf.citireGraf(fin);

    vector<int> rez = graf.BFS(S);
    for(int i = 1; i < rez.size() ; ++i) {
        fout << rez[i] << " ";
    }
    fout.close();
}

void Rezolva_DFSComponenteConexe(string fisierIntrare, string fisierIesire) //https://infoarena.ro/problema/dfs - 100pc
{
    ofstream fout(fisierIesire);
    fout << CitireGraf(fisierIntrare, 0, 0, 0, 0).nrCompConexe();
    fout.close();
}

void Rezolva_Biconex(string fisierIntrare, string fisierIesire)    // https://infoarena.ro/problema/biconex  - 100 pc
{
   vector<vector<int>> rez = CitireGraf(fisierIntrare, 0, 0, 0, 0).Biconex();
    ofstream fout(fisierIesire);
    fout << rez.size() << "\n";
    for(int i = 0; i < rez.size(); ++i) {
        for(int j = 1; j < rez[i].size(); ++j)
            fout << rez[i][j] << " ";
        fout << "\n";
    }
    fout.close();
}

void Rezolva_CTC(string fisierIntrare, string fisierIesire)      //https://infoarena.ro/problema/ctc - 100 pc
{
    vector<vector<int>> rez = CitireGraf(fisierIntrare, 1, 0, 0, 0).CTC();
    ofstream fout(fisierIesire);
    fout << rez.size() << "\n";
    for (int i = 0; i < rez.size(); ++i) {
        for (int j = 0; j < rez[i].size(); ++j) {
            fout << rez[i][j] << " ";
        }
        fout << "\n";
    }
    fout.close();
}

void Rezolva_SortareTopologica(string fisierIntrare, string fisierIesire)   //https://infoarena.ro/problema/sortaret - 100pc
{
    ofstream fout(fisierIesire);
    stack<int> rez = CitireGraf(fisierIntrare, 1, 0, 0, 0).sortareTopologica();
    while (!rez.empty()) {
        fout << rez.top() << " ";
        rez.pop();
    }
    fout.close();
}

//exemple: (DA:(5 5 5 5 5 5),(5 5 5 5 4 4),(3 3 3 3),(4 4 4 4 4),(2 2 1 1),(2 2 1 1 0),(2 0 1 2 1),(2 1 2 1 2 1 2 3)
//(NU:(5 5 5 5 5 4),(3 2 1 0), (3 3 3 3 3)
void Rezolva_Havel_Hakimi(string fisierIntrare, string fisierIesire)
{
    ifstream fin(fisierIntrare);
    ofstream fout(fisierIesire);
    vector<int> grade;
    int N, x;
    fin >> N;

    for (int i = 0; i < N; i++) {
        fin >> x;
        grade.push_back(x);
    }
    fin.close();
    if(Graf :: existaGraf(grade)) fout << "Da";
    else fout << "Nu";
    fout.close();
}

void Rezolva_Critical_Connection()    //https://leetcode.com/problems/critical-connections-in-a-network/  - Accepted
{
    int N, M;
    cin >> N >> M;

    Graf graf(N, M, 0, 0, 0, 0);
    graf.citireGraf(cin);
    vector<vector<int>> rez =  graf.CriticalConnections();

    cout << "[";
    for(int i = 0; i < rez.size() ; ++i) {
        cout << "[" << rez[i][0] << "," << rez[i][1] << "]";
    }
    cout << "]";
}

void Rezolva_Disjoint(string fisierIntrare, string fisierIesire)  //https://www.infoarena.ro/problema/disjoint - 100 pc
{
    ifstream fin(fisierIntrare);
    ofstream fout(fisierIesire);

    int N, M;
    fin >> N >> M;
    vector<bool> rez;
    vector<pair<int, int>> multimi(N + 1, {0,0}); //indexare de la 1

    for(int i = 1; i <= N; ++i)
        multimi[i].first = i, multimi[i].second = 1;  //fiecare numar e singur intr o multime

    for(int i = 1; i <= M; ++i) {
        int x, y, task;
        fin >> task >> x >> y;

        if (task == 1) {
            if(Graf::findDis(x, multimi) != Graf::findDis(y, multimi)) {
                Graf::unionDis(x, y, multimi);
            }
        } else {
            if(Graf::findDis(x, multimi) == Graf::findDis(y, multimi))
                fout << "DA\n";
            else
                fout << "NU\n";
        }
    }
    fin.close();
    fout.close();
}

void Rezolva_Dijkstra(string fisierIntrare, string fisierIesire)    // https://www.infoarena.ro/problema/dijkstra  - 100 pc
{
    ofstream fout(fisierIesire);
    vector<int> distDij = CitireGraf(fisierIntrare, 1, 1, 0, 0).Dijkstra(1);

    for(int i = 2; i < distDij.size(); i++)
        if(distDij[i] == INF) fout << 0 << " ";
        else fout << distDij[i] << " ";
}

void Rezolva_BMF(string fisierIntrare, string fisierIesire)    // https://www.infoarena.ro/problema/bellmanford  - 100 pc
{
    ofstream fout(fisierIesire);
    vector<int> distBMF = CitireGraf(fisierIntrare, 1, 1, 0, 0).BellmanFord(1);

    if(distBMF.size()) {
        for(int i = 2; i < distBMF.size(); i++) {
            fout << distBMF[i] << " ";
        }
    } else
        fout << "Ciclu negativ!";
}

void Rezolva_APM(string fisierIntrare, string fisierIesire)    //https://www.infoarena.ro/problema/apm  - 100 pc
{
    ofstream fout(fisierIesire);
    Graf graf = CitireGraf(fisierIntrare, 0, 1, 0, 0);
    vector<pair<int, pair<int, int>>> muchii =  graf.findLisMuchii();
    pair<int, vector<pair<int, int>>> rez = graf.Kruskal(muchii);

    fout << rez.first << "\n";
    fout << rez.second.size() - 1 << "\n";
    for(int i = 1; i < rez.second.size(); ++i)
        fout << rez.second[i].first << " " << rez.second[i].second << "\n";
}

void Rezolva_Max_Flow(string fisierIntrare, string fisierIesire)    //https://www.infoarena.ro/problema/maxflow  - 100 pc
{
    ofstream fout(fisierIesire);
    Graf graf = CitireGraf(fisierIntrare, 1, 0, 1, 0);

    fout << graf.maxFlow(1, graf.getNrNoduri(), 110005);
}


void Rezolva_Darb(string fisierIntrare, string fisierIesire)    //https://www.infoarena.ro/problema/darb  - 100 pc
{
    ifstream fin(fisierIntrare);
    int N;
    fin >> N;
    Graf graf(N, N-1, 0, 0, 0, 0);
    graf.citireGraf(fin);

    int u1, u2, distanta;
    graf.BFS_final(1, u1, distanta);
    graf.BFS_final(u1, u2, distanta);

    ofstream fout(fisierIesire);
    fout << distanta;
}

void Rezolva_RoyFloyd(string fisierIntrare, string fisierIesire)    //https://www.infoarena.ro/problema/royfloyd  - 100 pc
{
    ifstream fin(fisierIntrare);
    int N;
    fin >> N;

    vector<vector<int>> matrice;
    matrice.resize(N+1);
    for (int i = 1; i <= N; ++i) {
        matrice[i].resize(N+1);
        for (int j = 1; j <= N; j++) {
            int c;
            fin >> c;
            matrice[i][j] = c;
        }
    }

    ofstream fout(fisierIesire);
    vector<vector<int>> distante = Graf :: royFloyd(matrice, 1005);
    for(int i = 1; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            fout << distante[i][j] << " ";
        }
        fout << "\n";
    }
}

void Rezolva_Eulerian(string fisierIntrare, string fisierIesire)
{
    ifstream fin(fisierIntrare);
    int N;
    Graf graf = CitireGraf(fisierIntrare, 0, 0, 0, 1);

//    cout << graf;
    ofstream fout(fisierIesire);
    vector<int> ciclu = graf.cicluEuler();
    for(int i = 0; i < ciclu.size()-1; ++i) {
        fout << ciclu[i] << " ";
    }
    fout.close();
}

void Rezolva_Cuplaj(string fisierIntrare, string fisierIesire)
{
    ifstream fin(fisierIntrare);
    int N, M, E;
    fin >> N >> M >> E;

    vector<vector<int>> lisAdiacenta(N+1);
    for(int i = 1; i <= E; i ++){
        int x, y;
        fin >> x >> y;
        lisAdiacenta[x].push_back(y);
    }

    fin.close();

    pair<int, vector<int>> cuplaj = Graf :: cuplajMaxim(lisAdiacenta, N, M);

    ofstream fout(fisierIesire);
    fout<< cuplaj.first <<"\n";
    for(int i = 1; i <= N; ++i)
        if(cuplaj.second[i])
            fout << i<< " " << cuplaj.second[i] << "\n";
    fout.close();
}

void Rezolva_Hamilton(string fisierIntrare, string fisierIesire)
{
    ifstream fin(fisierIntrare);
    Graf graf = CitireGraf(fisierIntrare,1,1,0,0);
    int rezultat = graf.hamilton();

    ofstream fout(fisierIesire);
    if(rezultat != -1) {
        fout << rezultat;
    }
    else {
        fout << "Nu exista solutie";
    }
}

///Driver code
int main()
{
    Rezolva_Hamilton("hamilton.in", "hamilton.out");
    return 0;
}
