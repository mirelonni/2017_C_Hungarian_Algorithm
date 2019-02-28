#ifndef __SOLVE_H__
#define __SOLVE_H__

#include <stdlib.h>

typedef struct {
    int nn; // numar noduri
    int **mc; // matrice de costuri
    int **ml; // matrice de lucru
    int *lin; // vector pentru identificarea zerourilor pe linii
    int *col; // vector pentru identificarea zerourilor pe coloane
} TGraphM;

void alloc_matrix(TGraphM * g, int n)
{
    g->nn = n;
    // ma alocat spatiu corespunzator pentru cate elemente folosesc
    // o matrice de costuri
    // o matrice pe care lucrez (aplic algoritmul unguresc)
    // 2 vectori pentru alegerea elementelor folosite la suma de final
    g->mc = (int **)malloc(n * sizeof(int*));
    g->ml = (int **)malloc(n * sizeof(int*));
    g->lin = (int*)malloc(n * sizeof(int));
    g->col = (int*)malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        g->mc[i] = (int *)calloc(n, sizeof(int));
        g->ml[i] = (int *)calloc(n, sizeof(int));
    }
}

void insert_edge_matrix(TGraphM *g, int v1, int v2, int cost)
{
    // am inserat si in ambele matrici aceleasi valori
    g->mc[v1][v2] = cost;
    g->ml[v1][v2] = cost;
}

void make_zeros_row(TGraphM *g)
{
    // functia care scade minimul de pe linii
    int i, j, min_r;

    for (i = 0; i < g->nn; i++) {
        min_r = g->ml[i][0];

        for (j = 0; j < g->nn; j++) {
            if (g->ml[i][j] < min_r) {
                min_r = g->ml[i][j];
            }
        }

        for (j = 0; j < g->nn; j++) {
            g->ml[i][j] -= min_r;
        }
    }
}

void make_zeros_col(TGraphM *g)
{
    // functia care scade minimul de pe coloane
    int i, j, min_c;

    for (i = 0; i < g->nn; i++) {
        min_c = g->ml[0][i];

        for (j = 0; j < g->nn; j++) {
            if (g->ml[j][i] < min_c) {
                min_c = g->ml[j][i];
            }
        }

        for (j = 0; j < g->nn; j++) {
            g->ml[j][i] -= min_c;
        }
    }
}

int fix(TGraphM *g)
{
    // for care gaseste minimul de pe elementele nemarcate
    int min = 1000, i, j;

    for (i = 0; i < g->nn; i++) {
        for (j = 0; j < g->nn; j++) {
            if (g->lin[i] == -1 && g->col[j] == -1) {
                if (g->ml[i][j] < min) {
                    min = g->ml[i][j];
                }
            }
        }
    }

    // for care scade minimul de pe elementele nemarcate
    // si il aduna la cele marcate de 2 ori
    for (i = 0; i < g->nn; i++) {
        for (j = 0; j < g->nn; j++) {
            if (g->lin[i] == -1 && g->col[j] == -1) {
                g->ml[i][j] -= min;
            }

            if (g->lin[i] != -1 && g->col[j] != -1) {
                g->ml[i][j] += min;
            }
        }
    }

    // functia returneaza minimul deoarece avem nevoie de el in implementarea
    // algoritmului unguresc
    return min;
}

void chose(TGraphM *g)
{
    int cond = 0, min = 1000, k, nr_0_col, nr_0_lin, i, j;

    while (cond != g->nn) {
        // while pana cand marcam un nr necesar de elemente
        // (egale cu dimensiunea matricii)
        if (min != 0) {
            // daca minimul intors de fix nu este 0 adica daca facem modificari
            // in matrice
            cond = 0;      // conditia de while este readusa la 0

            for (i = 0; i < g->nn; i++) {
                g->lin[i] = -1; // vectorii de elemente marcate sunt adusi
                g->col[i] = -1; // la conditia de fals
            }
        }

        for (i = 0; i < g->nn; i++) {
            k = 0;
            nr_0_col = 0; // contorul de zerouri este initializat la 0

            for (j = 0; j < g->nn; j++) {
                if (g->ml[i][j] == 0) {
                    if (g->col[j] == -1 && g->lin[i] == -1) {
                        // daca elementul este 0 si nu au mai fost alte elemente
                        // marcate pe linia sau coloana respectiva
                        nr_0_col++;
                        k = j;
                    }
                }

            }

            if (nr_0_col == 1) {
                // daca este un singur 0 se atribuie indicele coloanei
                // in pozitia corespunzatoare a vectorului col (k)
                g->col[k] = i;
                cond++;
            }
        }

        for (i = 0; i < g->nn; i++) {
            k = 0;
            nr_0_lin = 0;

            for (j = 0; j < g->nn; j++) {
                if (g->ml[j][i] == 0) {
                    if (g->col[i] == -1 && g->lin[j] == -1) {
                        nr_0_lin++;
                        k = j;
                    }
                }
            }

            if (nr_0_lin == 1) {
                g->lin[k] = i;
                cond++;
            }
        }

        if (cond != g->nn) {
            // daca dupa parcurgerea pe linii si pe coloane sunt mai putine
            // taieturi decat dimensiunea matricii se apeleaza fix pentru a
            // se reusi pasul de marcare a zerourilor din matrice cu un lumar
            // egal de taieturi pentru finalizarea algoritmului unguresc
            min = fix(g);
        }
    }
}

int get_sum(TGraphM *g)
{
    // se calculeaza suma de final ce corespunde indicilor din lin si col
    int sum = 0;

    for (int i = 0; i < g->nn; i++) {
        if (g->lin[i] != -1) {
            sum += g->mc[i][g->lin[i]];
        }

        if (g->col[i] != -1) {
            sum += g->mc[g->col[i]][i];
        }
    }

    return sum;
}

void free_the_people(TGraphM *g)
{
    // functie de eliberare a memoriei ocupate de structura
    for (int i = 0; i < g->nn; i++) {
        free(g->ml[i]);
        free(g->mc[i]);
    }

    free(g->ml);
    free(g->mc);
    free(g->lin);
    free(g->col);
}

int solve(char *testInputFileName) {
    FILE * file = fopen(testInputFileName,"r");

    if (file == NULL) return -1;

    int n;
    fscanf(file, "%d", &n);
    TGraphM g;
    alloc_matrix(&g, n);
    int i, j, cost;

    // citirea matriciei din fisier si introducerea ei in structura
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fscanf(file, "%d", &cost);
            insert_edge_matrix(&g, i, j, cost);
        }
    }

    fclose(file);
    // primul pas al algoritmului unguresc
    make_zeros_row(&g);
    make_zeros_col(&g);
    // ceilalti pasi sunt inclusi in functia chose
    chose(&g);
    // ultimul pas din algoritmul unguresc
    int sum = get_sum(&g);
    // eliberarea memoriei alocate
    free_the_people(&g);
    return sum;
}

#endif
