#include <cstdio>
#include <iostream>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include "cpxmacro.h"


using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

// error status and message buffer
int status;
char errmsg[BUF_SIZE];

//per registrare il tempo
timeval start, stop;
double elapsedTime;


const int NAME_SIZE = 512;
char name[NAME_SIZE];

int N = 6; // |N| numero di nodi

//price matrix
double** C;

void inizialize_matrix(ifstream &myfile) {

    if (myfile.is_open())
        myfile >> N;
        else cerr << "Unable to open file!!!!!";

    //cout << "numero di nodi letto" << N << endl;

    C= new double* [N];
        for(int i = 0; i < N; ++i)
            C[i] = new double[N];


    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
           myfile >> C[i][j];
            //cout << "numero letto" << C[i][j] << endl;
        }
    }

}



void setupLP(CEnv env, Prob lp) {

    int xCounter = 0;

    // hmap per le variabili x, trova l'indice di x_i_j
    std::vector<std::vector<int> > xMap;
    xMap.resize(N);
    for (int i = 0; i < N; ++i) {
        xMap[i].resize(N);
        for (int j = 0; j < N; ++j) {
            xMap[i][j] = -1;
        }
    }

    // x_i_j = numero di unità di flusso trasportate dal nodo i al nodo j, ∀ (i, j) ∈ A
    {

    char varType = 'I';
    double obj = 0.0;
    double lb = 0.0;
    double ub = CPX_INFBOUND;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            snprintf(name, NAME_SIZE, "x_%i_%i", i, j);
            char* varName = (char*) (&name[0]);
            CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, &varType, &varName);

            xMap[i][j] = xCounter;
            ++xCounter;
        }
    }

    }

    // hmap per le variabili y, trova l'indice di y_i_j
    std::vector<std::vector<int> > yMap;
    yMap.resize(N);
    for (int i = 0; i < N; ++i) {
        yMap[i].resize(N);
        for (int j = 0; j < N; ++j) {
            yMap[i][j] = -1;
        }
    }

    //inserimentoi y in hmap
    // y_i_j = 1 se l'arco (i, j) viene utilizzato, 0 altrimenti, ∀ (i, j) ∈ A
    {
    char varType = 'B';
    double lb = 0.0;
    double ub = 1.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {

            double obj = C[i][j];

            snprintf(name, NAME_SIZE, "y_%i_%i", i, j);
            char* varName = (char*) (&name[0]);
            CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, &varType, &varName);

            yMap[i][j] = xCounter;
            ++xCounter;
        }
    }
    }

    // primo vincolo

    {
        char sense = 'E';
        int matbeg = 0; // 1 solo vincolo, al primo posto
        double rhs = N;
        std::vector<int> idx(N); // N variabili
        std::vector<double> coef(N, 1.0);

        for (int j = 0; j < N; j++) {
            idx[j] = xMap[0][j];
        }

        CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], 0, 0);
    }

    // secondo vincolo
    for (int k = 1; k < N; ++k) {
        std::vector<int> idx(N * 2);
        std::vector<double> coef(N * 2);

        for (int i = 0; i < N; ++i) {
            idx[i] = xMap[i][k];
            coef[i] = 1.0;
        }

        for (int j = 0, ind = N; j < N; ++j, ++ind) {
            idx[ind] = xMap[k][j];
            coef[ind] = -1.0;
        }

        //std::cout << std::endl;
        char sense = 'E';
        int matbeg = 0;
        double rhs = 1.0;

        CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], 0, 0);
    }

    // terzo vincolo

    {
    char sense = 'E';
    int matbeg = 0;
    double rhs = 1.0;
    std::vector<int> idx(N);
    std::vector<double> coef(N, 1.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            idx[j] = yMap[i][j];
        }

        CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], 0, 0);
    }

    }

   {
    // quarto vincolo
    char sense = 'E';
    int matbeg = 0;
    double rhs = 1.0;
    std::vector<int> idx(N);
    std::vector<double> coef(N, 1.0);

    for (int j = 0; j < N; ++j) {

        for (int i = 0; i < N; ++i) {
            idx[i] = yMap[i][j];
        }

        CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], 0, 0);
    }

    }

    // vincolo quinto
    char sense = 'L';
    int matbeg = 0;
    double rhs = 0.0;
    std::vector<int> idx(2);
    std::vector<double> coef(2);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {

            idx[0] = xMap[i][j];
            idx[1] = yMap[i][j];

            coef[0] = 1.0;
            coef[1] = -N;

            CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], 0, 0);
        }
    }

    // vieto l'utilizzo dei cappi
   {
    double coef = 1.0;
    char sense = 'E';
    int matbeg = 0;
    double rhs = 0.0;

    for (int i = 0, j = 0; i < N; ++i, ++j) {
        int idx = yMap[i][j];
        CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx, &coef, 0, 0);
    }
    }
}

int main(int argc, char const *argv[]) {

    if (argc < 2) { // per la ricerca locale servono 2 parametri
            throw std::runtime_error("usage: ./main filename.dat");}

	ifstream myfile (argv[1], std::ios::in);

	inizialize_matrix(myfile);

	myfile.close();

  try {
        // init


        DECL_ENV(env);
        DECL_PROB(env, lp, "risolutore quadro elettrico"); // macro per dare il nome

        setupLP(env, lp); // setup LP

        CPXsetdblparam(env, CPX_PARAM_TILIM, 1800.0); //upper bound of time for Cplex resolution in seconds, return the best solution it has so far

        // print problem description to a file
        CHECKED_CPX_CALL(CPXwriteprob, env, lp, "quadro_elettrico.lp", 0);


        struct timeval  tv1, tv2;
		gettimeofday(&tv1, NULL);
        gettimeofday(&tv1, NULL);
        clock_t start = clock();
        CHECKED_CPX_CALL(CPXmipopt, env, lp); // optimize
        clock_t end = clock();
        gettimeofday(&tv2, NULL);



        // stampa il valore della funzione obiettivo
        double objval;
        CHECKED_CPX_CALL(CPXgetobjval, env, lp, &objval);
        std::cout << "Obj val: " << objval << std::endl;

        // get the number and value of all variables
        int cur_numcols = CPXgetnumcols(env, lp);
        std::vector<double> varVals;
        varVals.resize(cur_numcols);
        CHECKED_CPX_CALL(CPXgetx, env, lp, &varVals[0], 0, cur_numcols - 1);

        int surplus; // will contain the space missing to save the column names

        // se a CPXgetcolname passiamo una array troppo piccolo per salvarci
        // tutti i nomi delle variabili, lui ritorna un codice d'errore e scrive
        // in surplus quanto spazio manca (un valore negativo)
        // NULL arguments to obtain the total memory required
        status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0, cur_numcols - 1);
        int cur_colnamespace = -surplus; // the space needed to save the names

        // allocate memory
        char** cur_colname = (char **) malloc(sizeof (char *)*cur_numcols);
        char* cur_colnamestore = (char *) malloc(cur_colnamespace);

        // get the names
        CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace, &surplus, 0, cur_numcols - 1);

        // print index, name and value of each column
        //  for (int i = 0; i < cur_numcols; i++) {
        //    std::cout << "Column " << i << ", " << cur_colname[i] << " = " << varVals[i] << std::endl;
        //}

        // stampa la soluzione in formato standard
        CHECKED_CPX_CALL(CPXsolwrite, env, lp, "quadro_elettrico.sol");

        std::cout << "TEMPO: " << (double)(tv2.tv_sec+tv2.tv_usec*1e-6 - (tv1.tv_sec+tv1.tv_usec*1e-6)) << " seconds (user time)\n";
        cout << " TEMPO: " << ((float)(end - start))/ CLOCKS_PER_SEC << " seconds (CPU time)\n";

        // free
        free(cur_colname);
        free(cur_colnamestore);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);

        gettimeofday(&stop, NULL);

        for(int i = 0; i < N; ++i) {
                delete [] C[i];
            }
            delete [] C;


    } catch (std::exception& e) {
        std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
    }

    return 0;
}
