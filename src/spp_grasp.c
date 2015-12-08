/* spp_grasp.c
 * Authors : DELAVERNHE Florian, LEGRU Guillaume
 * Date : 17th of Nov. 2015
 *
 * Implementation of GRASP metaheuristic applied to SPP.
 *
 * use : argv[0] datafile alpha
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>

struct timeval start_utime, stop_utime;

void crono_start()
{
	struct rusage rusage;
	
	getrusage(RUSAGE_SELF, &rusage);
	start_utime = rusage.ru_utime;
}

void crono_stop()
{
	struct rusage rusage;
	
	getrusage(RUSAGE_SELF, &rusage);
	stop_utime = rusage.ru_utime;
}

double crono_ms()
{
	return (stop_utime.tv_sec - start_utime.tv_sec) * 1000 +
    (stop_utime.tv_usec - start_utime.tv_usec) / 1000 ;
}

typedef struct
{
	int nbctr;
	int nbvar;
	int** matrix;
	int* coef;
} data;

void display (const int* sol, const int size)
{
	int i;
	for (i = 0; i < size; i++)
	{
		printf("%d ", sol[i]);
	}
	printf("\n");
}

double cout(int* sol, data* dat)
{
	int i;
	double val = 0;
	for (i = 0; i < dat->nbvar; ++i)
	{
		val += sol[i] *(double) dat->coef[i];
		
	}
	return val;
}

double localSearch(data* dat, int * sol,int* actCtr)
{
	int i,j,k;
	int present[dat->nbvar];
	double currentCost = cout(sol, dat);
	int fin = 0;
	int currentSol[dat->nbvar];
	int solProp[dat->nbvar];
	for (i = 0; i < dat->nbvar; i++) currentSol[i] = sol[i];
		
	while(!fin)
	{
		fin = 1;
		for(i = 0; i < dat->nbvar; ++i)
		{
			if(sol[i])
			{
				for(j = 0; j < dat->nbvar;++j)
				{
					present[j] = 1;
				}

				for(j=0;j< dat->nbctr;++j)
				{
					for(k=0;k<dat->nbvar;++k)
					{
						if(actCtr[j] == 0)
						{
							if(dat->matrix[j][k] != dat->matrix[j][i])
							{
								present[k] = 0;
							}
						}	
					}
				}

				int currentOne = i;
				

				for(k=0;k<dat->nbvar;++k)  
				{
					if(present[k] && k!=currentOne)
					{
						for (i = 0; i < dat->nbvar; i++) solProp[i] = currentSol[i];
						solProp[k] = 1;
						solProp[currentOne] = 0;

						double c = cout(solProp,dat);
						if (c > currentCost)
						{
							currentCost = c;
							for (i = 0; i < dat->nbvar; i++) currentSol[i] = solProp[i];
							fin=0;
						}
					}
				}
			}
		}
	}
	for (i = 0; i < dat->nbvar; i++) sol[i] = currentSol[i];
}

int readfile(data* dat, char* datafile)
{
	FILE *fin;
	int val; 
	int i,j,k;
	int ctrSize;
	
	fin = fopen(datafile, "r");

	//nbctr nbvar
	fscanf(fin, "%d", &val);
	dat->nbctr = val;
	fscanf(fin, "%d", &val);
	dat->nbvar = val;
	
	dat->coef = (int *) malloc (dat->nbvar * sizeof(int));
	for (i = 0; i < dat->nbvar; ++i)
	{
		fscanf(fin, "%d", &val);
		dat->coef[i] = val;
	}
	
	dat->matrix = (int **) malloc (dat->nbctr * sizeof(int*));
	for (i = 0; i <= dat->nbctr; ++i)
	{
		dat->matrix[i] = (int *) malloc (dat->nbvar * sizeof(int));
		for (k = 0; k < dat->nbvar; ++k)
		{
			dat->matrix[i][k] = 0;
		}
		
		fscanf(fin, "%d", &ctrSize);
		for (j=0; j < ctrSize; ++j)
		{
			fscanf(fin, "%d", &val);
			dat->matrix[i][val] = 1;

						
		}
	}
	return 1; 
}

void setUtility(double* utility, data* dat, int* fixVar, int*actCtr)
{
	int k;
	int i,j;	
	for (j = 0; j < dat->nbvar; ++j)
	{
		k = 0;

		if(!fixVar[j])
		{
			for (i = 0; i < dat->nbctr; ++i)
			{
				if (actCtr[i])
				{
					k += dat->matrix[i][j];
				}
			}
			utility[j] = (k == 0) ? 0 : (double) dat->coef[j] / (double) k;
		}
		else
		{
			utility[j] = 0;
		}
	}
}

const int allDisabled(const int* actCtr, const data* dat)
{
	int res = 1;
	int i;
	for (i = 0; i < dat->nbctr; ++i)
	{
		if(actCtr[i])
		{
			res = 0;
		}
	}
	return res;
}

int main (int argc, char** argv)
{
	data dat;
	double* utility;
	int* fixVar;
	int* actCtr;
	int* sol;
	double Umax, bound;
	double alpha[] = {0.0, 0.25, 0.5, 0.7, 0.8, 0.85, 0.90, 0.95, 1.0};
	int* aux;
	int choice;
	int nbIter;

	double temps;
	
	int i,j,k,l;
	int cpt;
	
	printf("Vérification du nombre d'argument ... \n");
	
	if (argc != 3)
	{
		printf("error : usage : %s <datafile> <nbIterations>\n", argv[0]);
		exit(1);
	}
	printf("argv[0] : %s\n", argv[0]);
	printf("argv[1] : %s\n", argv[1]);
	printf("argv[2] : %s\n", argv[2]);
	
	printf("OK!\n");
	
	nbIter = atoi(argv[2]);

	printf("Lecture du fichier de données... \n"); 
	if (readfile(&dat, argv[1]))
	{
		crono_start();
		
		printf("OK!\n");
		printf("Allocation mémoire...\n");
		
		sol = (int*) malloc (dat.nbvar* sizeof(int));
		fixVar = (int*) malloc (dat.nbvar* sizeof(int));
		actCtr = (int*) malloc (dat.nbctr* sizeof(int));
		utility = (double *) malloc (dat.nbvar * sizeof(double));

		for (i = 0; i < dat.nbctr; ++i)
		{
			actCtr[i] = 1;
		}

		for (j = 0; j < dat.nbvar; ++j)
		{
			fixVar[j] = 0;
		}
		
		for (j = 0; j < dat.nbvar; ++j)
		{
			sol[j] = 0;
		}

		printf("OK!\n");

		for (l = 1; l < 9; ++l) {
			double valmin = 9999999;
			double valmax =0;
			double valmoy = 0;

			for (k = 0; k < nbIter; ++k)
			{
				//printf("Construction d'une solution initiale ...\n");

				for (i = 0; i < dat.nbctr; ++i)
				{
					actCtr[i] = 1;
				}

				for (j = 0; j < dat.nbvar; ++j)
				{
					fixVar[j] = 0;
				}
				for (j = 0; j < dat.nbvar; ++j)
				{
					sol[j] = 0;
				}
				
				
				double Umax=-1;
				while (!allDisabled(actCtr, &dat) && Umax!=0.00)
				{
					// printf("Calcul des utilitées ...\n");
					setUtility(utility, &dat, fixVar, actCtr);
					// printf(" ... OK!\n");

			
					// printf("Calcul de Umax ...\n");
					Umax = -1;
					for (j = 0; j < dat.nbvar; ++j)
					{
						if ( utility[j] > Umax)
						{
							Umax  = utility[j];
						}
					}
					// printf("OK : Umax = %f!\n", Umax);
					assert (Umax != -1); 
					if(Umax != 0.0)
					{
						// printf("Determination des candidats en fonction de Umax et de alpha ...\n");
						bound = Umax * alpha[l];
						aux = (int *) malloc (dat.nbvar*sizeof(int)); 
						cpt = 0;
						for (j = 0; j < dat.nbvar; ++j)
						{
							if(utility[j] >= bound)
							{
								aux[cpt] = j;
								++cpt;
							}

						}
						// printf("OK : Seuil  = %f pour un Umax de %f\n", bound, Umax);

						srand((unsigned int) time(0));
						choice  = (rand() % cpt) ;
						// printf("Choix aléatoire de la valeur fixé : %d d'utilité %f et de valeur %d \n", aux[choice], utility[aux[choice]],dat.coef[aux[choice]]);
						
						if(fixVar[aux[choice]] == 0)
						{
							sol[aux[choice]] = 1;
						
							fixVar[aux[choice]] = 1;
	
							for (i = 0; i < dat.nbctr; ++i)
							{
								if (dat.matrix[i][aux[choice]] == 1)
								{
									actCtr[i] = 0;
									for (j = 0; j < dat.nbvar; ++j)
									{
										if(dat.matrix[i][j] == 1)
										{
											fixVar[j] = 1;
										}
									}
								}
							}
						}
					}
				}
				// printf("Solution initiale : OK!\n");
				
				for (i=0;i<dat.nbctr;i++)
				{
					cpt = 0;
					for (j=0;j<dat.nbvar;j++)
					{
						if(sol[j] == 1 && dat.matrix[i][j] == 1)
						{
							cpt++;	
						
						}
					
					}
					if(cpt>=2)
					{
						printf("----------------------ERROR---------------------");
					}
					
				}
				
				localSearch(&dat, sol, actCtr);

				for (i=0;i<dat.nbctr;i++)
				{
					cpt = 0;
					for (j=0;j<dat.nbvar;j++)
					{
						if(sol[j] == 1 && dat.matrix[i][j] == 1)
						{
							cpt++;	
						
						}
					
					}
					if(cpt>=2)
					{
						printf("----------------------ERROR---------------------");
					}
					
				}
								
				double res  =cout(sol, &dat);
				if(res >= valmax)
				{
					valmax = res;
				}
				if(res <= valmin)
				{
					valmin = res;
				}
				valmoy += res;
			}

			
			valmoy /= (double) nbIter;
			printf("Alpha = %f\n", alpha[l]);
			printf("val min  : %f!\n",valmin);
			printf("val max  : %f!\n",valmax);
			printf("val moy  : %f!\n",valmoy);
			printf("-------------------------\n");
		}
		crono_stop();
		temps = crono_ms()/1000,0;

		printf("Temps de résolution : %f\n", temps);
		
		for (i = 0; i < dat.nbctr; ++i)
		{
			free(dat.matrix[i]);
		}
		free(dat.matrix);
		free(dat.coef);
		free(utility);
		free(actCtr);
		free(fixVar);
		free(sol);
		free(aux);
	}
	else
	{
		printf("impossible to read : wrong data file !");
		exit(1);
	}

	

	return 0;
}
