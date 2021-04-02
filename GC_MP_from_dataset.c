/*
//-------------------- #GC_from_generator.c ------------------------------- #

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
// This program implements the(averaged) message passing algorithm on real datasets (Friendship networks from the 
// music streaming site Deezer in the countries of Romania, Hungary and Croatia)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#define kc 10
#define err_tol 0.01

int main(int argc, char **argv)
{
    int i, j, n, k, iter, np, N, nc;
    int **knn1, *k1;
    double **sigma_N, **sigma_T, *sigma, *T;
    double rho, kappa, kappaT, kav, pcp, f, p, nsum, nsum_new, prod_N, prod_TN, prod, MGC, pc, nsumold;
    char filename[60], string1[50], string2[50];
    FILE *gp2, *gp, *ffp, *fp;

    srand48(time(NULL));

    //&&&&&&&&&CHANGE INPUT FILENAME%%%%%%%%%%%%%%%%%%%%%
    sprintf(filename, "REF9_RO_edges.txt");

    //%%%%%%%%%%%CHANGE OUTPUT FILE NAME%%%%%%%%%%%%%%%%
    sprintf(string2, "MessagePassing_RO.txt");

    /* Count the number of nodes of the dataset */
    N = 0;
    fp = fopen(filename, "r");
    while (!feof(fp))
    {
        fscanf(fp, "%d %d ", &i, &j);
        if (i > N)
            N = i;
        if (j > N)
            N = j;
    }
    fclose(fp);
    N = N + 1;      // Node ID starts from zero.
    
    /* Initialize all the variable needed. Allocate the memery. */

    gp = fopen(string2, "w");

    T = (double *)calloc(N, sizeof(double));
    k1 = (int *)calloc(N, sizeof(int));
    sigma = (double *)calloc(N, sizeof(double));
    knn1 = (int **)calloc(N, sizeof(int *));
    sigma_N = (double **)calloc(N, sizeof(double *));   // Messages
    sigma_T = (double **)calloc(N, sizeof(double *));

    for (i = 0; i < N; i++)
    {
        knn1[i] = (int *)calloc(N, sizeof(int));
        sigma_N[i] = (double *)calloc(N, sizeof(double));
        sigma_T[i] = (double *)calloc(N, sizeof(double));
    }

    for (i = 0; i < N; i++)
    {
        k1[i] = 0;
    }

    nsum = 1;

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

            sigma_N[i][j] = (drand48()); //    Initialize the messages needed for message passing.
            sigma_N[j][i] = (drand48());
            sigma_T[i][j] = (drand48());
            sigma_T[j][i] = (drand48());
        }
    }
    fclose(ffp);
    

    for(nc = 0; nc < 1; nc++){
        rho = nc * 0.01;
        printf("rho: %lf\n", rho);
        rho = 0.;
        kappa = 0.;
        kappaT = 0.;
        kav = 0.;
        for (i = 0; i < N; i++)
        {
            sigma[i] = 0.;

            /* Analytical results */

            kav += k1[i];
            T[i] = rho;
            if (k1[i] > kc)
            {
                T[i] = 1;
            }
            kappa += k1[i] * (k1[i] - 1.);
            kappaT += k1[i] * (k1[i] - 1.) * T[i];
        }
        kappa = kappa / kav;
        kappaT = kappaT / kav;
        pc = 0.5 * (-1. + sqrt(1. + 4. * kappaT / (kappa - kappaT))) / kappaT;
        if (pc > 1.)
        {
            pc = 1.;
        }
        // fprintf(gp2, "%lf %lf\n", pc, rho);


        /* Message passing averaged over the ensemble */
        for (f = 0.; f < 1.01; f += 0.01)
        {
            p = 1 - f;
            printf("p = %lf\n", p);
            nsumold = 10000;

            /* Initializaing messages */

            while (fabs(nsum - nsumold) > err_tol)
            {

                // printf("Err:%lf\n", fabs(nsum - nsumold));
                nsumold = nsum;
                nsum = 1;

                for (i = 0; i < N; i++)
                {
                    for (n = 0; n < k1[i]; n++)
                    {
                        j = knn1[i][n];
                        prod_N = 1.;
                        prod_TN = 1.;
                        for (np = 0; np < k1[i]; np++)
                        {
                            if (np != n)
                            {
                                prod_TN = prod_TN * (1. - sigma_N[knn1[i][np]][i] - sigma_T[knn1[i][np]][i]);
                                prod_N = prod_N * (1. - sigma_N[knn1[i][np]][i]);
                            }
                        }
                        nsum -= (sigma_N[i][j] + sigma_T[i][j]);
                        sigma_T[i][j] = p * T[i] * (1. - prod_N);
                        sigma_N[i][j] = p * (1 - T[i]) * (1. - prod_TN);
                        nsum += (sigma_N[i][j] + sigma_T[i][j]);
                    }
                }
            }
            
        
            /* After converge */
            MGC = 0.;
            for (i = 0; i < N; i++)
            {
                if (k1[i] > 0)
                {
                    prod = 1.;
                    for (n = 0; n < k1[i]; n++)
                    {
                        prod = prod * (1 - sigma_N[knn1[i][n]][i] - sigma_T[knn1[i][n]][i]);
                    }
                    MGC += (1 - prod);
                }
            }
            fprintf(gp, "%lf %lf %lf\n", p, rho, (double)MGC / (double)N);
            
        }
    }
    fclose(gp);
    // fclose(gp2);

    return 0;
}
