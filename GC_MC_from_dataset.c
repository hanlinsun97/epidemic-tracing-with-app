/*
//-------------------- #GC_from_dataset.c ------------------------------- #

This code can be redistributed and / or modified
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed ny the authors in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


If you use this code please cite

Ginestra Bianconi, Hanlin Sun, Giacomo Rapisardi, and Alex Arenas.
"Message-passing approach to epidemic tracing and mitigation with apps."
Physical Review Research 3, no.1(2021) : L012014.

(c)
Hanlin Sun(hanlin.sun @qmul.ac.uk)
Ginestra Bianconi(g.bianconi @qmul.ac.uk)

*/

//  This program implements the Monte Carlo simulation on real dataset
//  (Friendship networks from the music streaming site Deezer in the countries of Romania,
//  Hungary and Croatia)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define Nrunmax 5      // Max iteration of Monte Carlo
#define kc 10            // degree threshold of adopting an app. See detailed definition in the paper.

int *vis1, *size_cluster1, **knn1, *k1, c1, c2, c3;
int **y, z;

/**************************************************************************
 Recurrence is a subrutine for calulating the giant component in link percolation
 **************************************************************************/

int Recurrence(int i, int cluster_size, int ncluster)
{
    int j, n3, aus_cluster_size;
    cluster_size++;
    vis1[i] = ncluster;
    for (n3 = 0; n3 < k1[i]; n3++)
    {
        j = knn1[i][n3];
        if ((y[i][j] == 1) && (vis1[j] == 0))
        {
            aus_cluster_size = Recurrence(j, cluster_size, ncluster);
            cluster_size = aus_cluster_size;
        }
    }
    return cluster_size;
}

int main(int argc, char **argv)
{
    int i, j, n, it, ncluster1, ncluster2, ncluster3, GMCC, cluster_size, m1, m2, m3, m1_aus, m2_aus, m3_aus, c1_aus, c2_aus, c3_aus, *sigma, nrun, nc, np, nc2;
    int ii, jj, N;
    int s1, s2, Nc1, Nc2, Nc3, aus, aus3, **adj1, **adj2, **adj3, N0, **s, nrun2, **denR, *T, in;
    double p, f, **x, *Sm, **n1, **n2, MGC1, MGC2, MGC3, MGC4, nsum1, nsumold1, nsum2, nsumold2, aus1, aus2, sigma11, sigma10, *sigma11m, *sigma10m, **GCA, rho, r, *Tx;
    int ausi, GC;
    double pc, kav, kappa, kappaT, *kx, Norm1, pc_sum, *pc_list, sum_degree;
    char filename[60], string1[50], string2[50];

    FILE *gp2, *gp, *ffp, *gp3;
    FILE *fp = 0;



    //&&&&&&&&&CHANGE INPUT FILENAME%%%%%%%%%%%%%%%%%%%%%
    sprintf(filename, "REF9_RO_edges.txt");

    //%%%%%%%%%%%CHANGE OUTPUT FILE NAME%%%%%%%%%%%%%%%%
    sprintf(string2, "MonteCarlo_RO.txt");

    // The txt file contains 3 column:
    //          p: probability of retaining a link,
    //          rho: probability of adopting the app. See detailed definition in the paper.
    //          R: fraction of nodes in the giant component.





    N = 0;

    fp = fopen(filename, "r");

    // Get the number of nodes
    while (!feof(fp))
    {
        fscanf(fp, "%d %d ", &i, &j);
        if (i > N)
            N = i;
        if (j > N)
            N = j;
    }
    fclose(fp);
    N = N + 1; // The node ID starts from zero.

    gp = fopen(string2, "w");

    Tx = (double *)calloc(N, sizeof(double));
    vis1 = (int *)calloc(N, sizeof(int));
    kx = (double *)calloc(N, sizeof(double));
    k1 = (int *)calloc(N, sizeof(int));
    x = (double **)calloc(N, sizeof(double *));
    y = (int **)calloc(N, sizeof(int *));
    T = (int *)calloc(N, sizeof(int));

    pc_list = (double *)calloc(N, sizeof(double));
    knn1 = (int **)calloc(N, sizeof(int *));
    adj1 = (int **)calloc(N, sizeof(int *));
    GCA = (double **)calloc(101, sizeof(double *));

    for (i = 0; i < N; i++)
    {
        knn1[i] = (int *)calloc(N, sizeof(int));
        adj1[i] = (int *)calloc(N, sizeof(int));
        x[i] = (double *)calloc(N, sizeof(double));
        y[i] = (int *)calloc(N, sizeof(int));
        if (i < 101)
        {
            GCA[i] = (double *)calloc(101, sizeof(double));
        }
    }

    size_cluster1 = (int *)calloc(N, sizeof(int));

    for (i = 0; i < N; i++)
    {
        k1[i] = 0;
    }

    /* Obtain the degree sequence k1 and network structure knn1 from dataset or random generated network */
    /* knn1 is a matrix. The (i,j)-th entry is the j-th neighbor of node i. */

    ffp = fopen(filename, "r");

    while (!feof(ffp))
    {
        fscanf(ffp, "%d %d ", &i, &j);
        if (i != j)
        {
            k1[i]++;
            k1[j]++;
            knn1[i][k1[i] - 1] = j;
            knn1[j][k1[j] - 1] = i;
        }
    }
    fclose(ffp);


 
    for (nrun = 0; nrun < Nrunmax; nrun++)
    {
        printf("nrun: %d\n", nrun);

        // Average of uncorrelated network ensemble starts. (Theoretical analysis)

        for (nc = 0; nc < 1; nc++) // Only care about rho=0 case. Calculate for all rho from 0 to 1 by setting nc = 100
        {
            rho = nc * 0.01;

            kappa = 0.;
            kappaT = 0.;
            kav = 0.;
            for (i = 0; i < N; i++)
            {
                T[i] = 0;
                if ((Tx[i] < rho) || ((k1[i] > kc)))
                {
                    T[i] = 1;
                }
                kappa += k1[i] * (k1[i] - 1.);
                kappaT += k1[i] * (k1[i] - 1.) * T[i];
                kav += k1[i];
            }

            kappa = kappa / kav;
            kappaT = kappaT / kav;
            pc = 0.5 * (-1. + sqrt(1. + 4. * kappaT / (kappa - kappaT))) / kappaT;
            if (pc > 1.)
                pc = 1;

            pc_list[nc] += pc / Nrunmax;

            // Average of uncorrelated network ensemble ends. (Theoretical analysis)

            // Monte Carlo simulation starts.

            for (i = 0; i < N; i++)
            {
                Tx[i] = drand48();
                for (n = 0; n < k1[i]; n++)
                {
                    j = knn1[i][n];
                    x[i][j] = drand48();
                    x[j][i] = x[i][j];
                }
            }

            for (nc2 = 0; nc2 < 101; nc2++)
            {
                p = nc2 * 0.01;
                GC = 0.;
                for (i = 0; i < N; i++)
                {
                    for (n = 0; n < k1[i]; n++)
                    {
                        j = knn1[i][n];
                        y[i][j] = 0;
                        y[j][i] = 0;
                        if ((x[i][j] < p) && ((T[i] * T[j]) != 1))
                        {
                            y[i][j] = 1;
                            y[j][i] = 1;
                        }
                    }
                }

                ncluster1 = 0;
                for (i = 0; i < N; i++)
                {
                    vis1[i] = 0;
                }
                m1 = 0;
                ncluster1 = 0;
                c1 = 0;
                for (n = 0; n < N; n++)
                {
                    if (vis1[n] == 0)
                    {
                        cluster_size = 0;
                        ncluster1++;
                        cluster_size = Recurrence(n, cluster_size, ncluster1);
                        size_cluster1[ncluster1] = cluster_size;
                        if (cluster_size > m1)
                        {
                            m1 = cluster_size;
                            c1 = ncluster1;
                        }
                    }
                }

                Nc1 = c1;
                GC = 0.;
                for (i = 0; i < N; i++)
                {

                    if (vis1[i] == Nc1)
                    {
                        GC++;
                    }

                    if ((T[i] == 1) && (vis1[i] != Nc1))
                    {
                        ausi = 0;
                        for (n = 0; n < k1[i]; n++)
                        {
                            j = knn1[i][n];
                            if ((vis1[j] == Nc1) && (T[j] == 1) && (x[i][j] < p))
                            {
                                ausi = 1;
                            }
                        }
                        GC += ausi;
                    }
                }

                GCA[nc][nc2] += (double)GC;
            }
        }
    }

    for (nc = 0; nc < 1; nc++) // Only care about rho=0 case. Calculate for all rho from 0 to 1 by setting nc = 100
    {
        rho = (double)nc * 0.01;
        for (nc2 = 0; nc2 < 101; nc2++)
        {
            p = (double)nc2 * 0.01;
            fprintf(gp, "%f %lf  %f\n", (double)p, (double)rho, ((double)GCA[nc][nc2]) / ((double)(N * Nrunmax)));
        }
    }

    fclose(gp);
    return 0;
}