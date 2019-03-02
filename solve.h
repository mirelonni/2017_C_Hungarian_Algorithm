#ifndef __SOLVE_H__
#define __SOLVE_H__

#include <stdlib.h>

/*
  struct for saving the bipartite graph
*/
typedef struct {
    int nn; // number of nodes
    int **mc; // cost matrix
    int **ml; // work matrix
    int *lin; // vector for identifiing zeros on rows
    int *col; // vector for identifiing zeros on cols
} TGraphM;


/*
  graph alloc function
*/
void alloc_matrix(TGraphM * g, int n)
{
    g->nn = n;

    g->mc = (int **)malloc(n * sizeof(int*));
    g->ml = (int **)malloc(n * sizeof(int*));
    g->lin = (int*)malloc(n * sizeof(int));
    g->col = (int*)malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        g->mc[i] = (int *)calloc(n, sizeof(int));
        g->ml[i] = (int *)calloc(n, sizeof(int));
    }
}

/*
  graph insert function
*/
void insert_edge_matrix(TGraphM *g, int v1, int v2, int cost)
{

    g->mc[v1][v2] = cost;
    g->ml[v1][v2] = cost;
}

/*
  function that decreases each row with the lowest element
*/
void make_zeros_row(TGraphM *g)
{
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

/*
  function that decreases each col with the lowest element
*/
void make_zeros_col(TGraphM *g)
{

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
    // for loop that finds the min of unmarked elements
    int min = 99999999, i, j;

    for (i = 0; i < g->nn; i++) {
        for (j = 0; j < g->nn; j++) {
            if (g->lin[i] == -1 && g->col[j] == -1) {
                if (g->ml[i][j] < min) {
                    min = g->ml[i][j];
                }
            }
        }
    }

    // for loop that decreases unmarked ele with min
    // and adds min to the marked ones
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

    return min;
}

void chose(TGraphM *g)
{
    int cond = 0, min = 1000, k, nr_0_col, nr_0_lin, i, j;

    while (cond != g->nn) {
        // while loop until there are marked enough elements
        // (equal to the size of the matrix)
        if (min != 0) {
            // if fix doesn't return 0 it means that we have modified the matrix

            cond = 0;      // the while condition is brougt back to 0

            for (i = 0; i < g->nn; i++) {
                g->lin[i] = -1; // the elements vectors are brought back
                g->col[i] = -1; // to the false condition
            }
        }

        for (i = 0; i < g->nn; i++) {
            k = 0;
            nr_0_col = 0; // the zeros counter is initialized to 0

            for (j = 0; j < g->nn; j++) {
                if (g->ml[i][j] == 0) {
                    if (g->col[j] == -1 && g->lin[i] == -1) {
                        // if the element is 0 ant there are none other
                        // marked elements on that row or column
                        nr_0_col++;
                        k = j;
                    }
                }

            }

            if (nr_0_col == 1) {
                // if there is a single 0 the index of the column is saved in
                // the corespondent columns vector position(k)
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
            // if there are less cuts than the matrix size fix is called to
            // succed in the cutting step of the hungarian algorithm,
            // so the cuts will eventualy be equal to the matrix size
            min = fix(g);
        }
    }
}

/*
  returns the sum of the lowest costs
*/
int get_sum(TGraphM *g)
{

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

/*
  Function to free the allocated memory
*/
void free_the_people(TGraphM *g)
{

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

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fscanf(file, "%d", &cost);
            insert_edge_matrix(&g, i, j, cost);
        }
    }

    fclose(file);

    make_zeros_row(&g);
    make_zeros_col(&g);
    chose(&g);
    int sum = get_sum(&g);

    free_the_people(&g);
    return sum;
}

#endif
