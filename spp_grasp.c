/* spp_grasp.c
 * Authors : DELAVERNHE Florian, LEGRU Guillaume
 * Date : 17th of Nov. 2015
 *
 * Implementation of GRASP metaheuristic applied to SPP.
 *
 * use : argv[0] datafile alpha
 */


/* TODO relire bien la matrix*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

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
		val += (sol[i]) ? dat->coef[i] : 0;  
	}
	return val;
}

double localSearch(data* dat, int * sol)
{
	int i,j,k;
	int present[dat->nbvar];
	double currentCost = cout(sol, dat);
	int fin = 0;
	int currentSol[dat->nbvar];
	for (i = 0; i < dat->nbvar; i++) currentSol[i] = sol[i];
		
	while(!fin)
	{
		for(i = 0; i < dat->nbvar; ++i)
		{

			if(sol[i])
			{
				for(j = 0; j <dat->nbvar;++j)
				{
					present[j] = 1;
				}

				for(j=0;j< dat->nbctr;++j)
				{
					for(k=0;k<dat->nbvar;++k)
					{
						if(dat->matrix[j][k] != dat->matrix[j][i])
						{
							present[k] = 0;
						}
							
					}
				}

				int currentOne = i;
				
				fin = 1;
				for(k=0;k<dat->nbvar;++k)  
				{
					
					if(present[k] && k!=currentOne)
					{
						currentSol[k] = 1;
						currentSol[currentOne] = 0;
						

						double c = cout(currentSol,dat);
						if (c > currentCost)
						{
							currentCost = c;
							for (i = 0; i < dat->nbvar; i++) sol[i] = currentSol[i];
							currentOne  = k;
							fin=0;
						}
						else
						{
							currentSol[k] = 0;
							currentSol[currentOne] = 1;
						}
						// display(currentSol, dat->nbvar);
					}
				}
			}
		}
	}
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
	
	int i,j,k,l;
	int cpt;
	
	printf("Vérification du nombre d'argument ... \n");
	
	if (argc != 3)
	{
		printf("error : usage : %s <datafile> <alpha>\n", argv[0]);
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
		/* for (i = 0; i < dat.nbctr; i++) */
		/* { */
		/* 	for (j = 0; j < dat.nbvar; j++) */
		/* 	{ */
		/* 		printf(" %d", dat.matrix[i][j]); */
		/* 	} */
		/* 	printf("\n"); */
		/* } */

		//exit(0);
		
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
		printf("Construction d'une solution initiale ...\n");

		for (l = 1; l < 9; ++l) {
			for (k = 1; k <= nbIter; ++k)
			{
				while (!allDisabled(actCtr, &dat))
				{
					printf("Calcul des utilitées ...\n");
					setUtility(utility, &dat, fixVar, actCtr);
					printf("Utilitées : ");
					/* for (j = 0; j < dat.nbvar; j++) */
					/* { */
					/* 	printf("%f ", utility[j]); */
					/* } */
					printf(" ... OK!\n");

			
					printf("Calcul de Umax ...\n");
					Umax = -1;
					for (j = 0; j < dat.nbvar; ++j)
					{
						if ( utility[j] > Umax)
						{
							Umax  = utility[j];
						}
					}
					printf("OK : Umax = %f!\n", Umax);
					assert (Umax != -1); 

					printf("Determination des candidats en fonction de Umax et de alpha ...\n");
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
					printf("OK : Seuil  = %f pour un Umax de %f\n", bound, Umax);

					srand((unsigned int) time(0));
					choice  = rand() % cpt;
					printf("Choix aléatoire de la valeur fixé : %d d'utilité %f\n", aux[choice], utility[aux[choice]]);

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
					//exit(0);
				}
				printf("Solution initiale : OK!\n");
				display(sol, dat.nbvar);
			   
				localSearch(&dat, sol);
				display(sol, dat.nbvar);
				printf("Cout : %f\n", cout(sol, &dat));
			}
		}
		
	}
	else
	{
		printf("impossible to read : wrong data file !");
		exit(1);
	}

	

	return 0;
}
