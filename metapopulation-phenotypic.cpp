#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>
#include<tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <unistd.h>

using namespace std;

#define PI 3.1415927

#define NIL (0)    
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long idum;
gsl_rng *gerador;


struct populacao{
  int *N;
  int N0;
  double mut;
  int **step;
  int *nstep;
  int *contind;
  int ntraits;
  double sigma_r;
  double sigma_trait;
  double sigma_trait_mut;
  int n_subpopul;
  double alpha_r;
  int K;
  double migration;
  double lambda;
};

struct caminhada{

  int nstep;
  double fitness_final;
  int sequence_final;
};

struct otimo{
  
  double **traits;
  int tau;
  double *v;
  double **trait_shift;
  
};



struct sequencia{
  int *cont_mut;
  int **ind;
  int L;
  int ind_min;
  int *global_phase;
  int ind_max;
  int ind_max_ant;
  int ind_max_landscape;
  int dham;
  int Nmax;
  int *max;
  double roughness;
  double rough_local;
  //  int *path;
  int path_size;
  double **fitness;
  double ***traits;
  double **traits_new;
  double W_max;
  double *trait_init;
};

void fitnessinit(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
void dynamics(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
void migration(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
void fitness_evaluation(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);
void density(struct sequencia *sequence, struct populacao *popul, struct otimo *optimum);

int main(int ac, char **av)
{
  FILE *ptt1, *ptt2, *ptt3;

  otimo optimum;

  sequencia sequence;

  populacao popul;

  caminhada walk;

  int nwalks, max1, cont_path_total, cont_path, cont_newpath, i, j, k, cont_walk, verif, max, conf, config, t, tmax, model_dominant[2], *cont_seq, int_conf, min_step, cont_max, n_landscapes, popul_max, indice, ind_sub, ind_max, cont_extinction, cont_nonextinction, deme, sum_cont_extinction, cont_ind[2], cont_rescue_total, cont_rescue_partial;

  double fit_max, fit_med, soma1, cont_conf, sum_trait_2, *drop_fitness, *trait, sum_t_2, distance_phen, divergence, f_max;

  long seed;  //define a seed

  char arq1[1000], arq2[1000], arq3[1000];

  if (ac < 12 + atoi(av[1]) - 1) {  // Update the condition to reflect the new count
        std::cout << "Start the program like this:\n" << av[0]
                  << " <N> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> <K1> <K2> ... <Kn>\n"
                  << std::endl;
        return -1;
    }

  /*  if (ac!=12)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << " <N> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> \n"
            << endl;
      exit (-1);
      }*/

  j=0;
  
  popul.n_subpopul = atoi (av[++j]);
  popul.K = atoi (av[++j]);
  popul.ntraits = atoi (av[++j]);
  popul.mut = atof (av[++j]);
  popul.migration = atof (av[++j]);
  n_landscapes = atoi (av[++j]);
  popul.alpha_r = atof (av[++j]);
  popul.lambda = atof (av[++j]);
  drop_fitness = new double[popul.n_subpopul];
  for (int i = 0; i < popul.n_subpopul; ++i) {
        drop_fitness[i] = atof(av[++j]);  // Assign each K value from the arguments
  }
  sequence.W_max = atof (av[++j]);
  tmax = atoi (av[++j]);
    
  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set(gerador, seed); //give the seed to random generator
  
  srand(time(NULL));

  max1 = (int) pow(2.,(double)(sequence.L));
  if( sequence.L==32 )
    max1 = RAND_MAX;

  popul_max = 100000;

  popul.N = new int[popul.n_subpopul];

  sequence.traits = new double**[popul.n_subpopul];
  for( i=0; i<popul.n_subpopul; i++ )
    sequence.traits[i] = new double*[2*popul.K];
  for( i=0; i<popul.n_subpopul; i++ )
    for( j=0; j<(2*popul.K); j++ )
      sequence.traits[i][j] = new double[popul.ntraits];

  sequence.traits_new = new double*[2*popul.K];
  for( i=0; i<(2*popul.K); i++ )
    sequence.traits_new[i] = new double[popul.ntraits];

  sequence.trait_init = new double[popul.ntraits];

  optimum.traits = new double*[popul.n_subpopul];
  for( i=0; i<popul.n_subpopul; i++ )
    optimum.traits[i] = new double[popul.ntraits];

  optimum.trait_shift = new double*[popul.n_subpopul];
  for( i=0; i<popul.n_subpopul; i++ )
    optimum.trait_shift[i] = new double[popul.ntraits];

  sequence.fitness = new double*[popul.n_subpopul];
  for( i=0; i<popul.n_subpopul; i++ )
    sequence.fitness[i] = new double[2*popul.K];

  trait = new double[popul.ntraits];

  config = n_landscapes;
  conf = 0;
  cont_conf = 0;

  sum_cont_extinction = 0;

  cont_extinction = 0;
  cont_nonextinction = 0;
  cont_rescue_total = cont_rescue_partial = 0;
  while( (++conf)<=config )
    {

      for( k=0; k<popul.n_subpopul; k++ )
	popul.N[k] = popul.K;    
 
       for( deme=0; deme<popul.n_subpopul; deme++ )
	{
	  sum_trait_2 = -2*log((1-drop_fitness[deme])/sequence.W_max);
	  sum_t_2 = 0;
	  for( k=0; k<popul.ntraits; k++ )
	    {
	      trait[k] = gsl_ran_gaussian(gerador, 1.);
	      sum_t_2 += (trait[k]*trait[k]);
	    }
	  
	  for( k=0; k<popul.ntraits; k++ )
	    trait[k] = trait[k]*sqrt(sum_trait_2)/sqrt(sum_t_2);
      
	  for( k=0; k<popul.ntraits; k++ )
	    optimum.trait_shift[deme][k] = trait[k];
      
      
	  for( k=0; k<popul.ntraits; k++ )
	    optimum.traits[deme][k] = optimum.trait_shift[deme][k];
	}
     

      fitnessinit(&sequence,&popul,&optimum);

      t = 0;
      verif = 0;
      while( (t<=tmax) && ((popul.N[0]>0)||(popul.N[1]>0)) )
	{
	  t++;

	  fitness_evaluation(&sequence,&popul,&optimum);

	  dynamics(&sequence,&popul,&optimum);

 	  migration(&sequence,&popul,&optimum);

          density(&sequence,&popul,&optimum);

	  cont_ind[0] = cont_ind[1] = 0;

	  for( deme=0; deme<popul.n_subpopul; deme++ )
	    for( k=0; k<popul.N[deme]; k++ )
	      if( sequence.fitness[deme][k]>1 )
		cont_ind[deme]++;
	  
	  if( (cont_ind[0]>100) && (cont_ind[1]>100) )
	    {
	      cont_rescue_total++;
	      break;
	    }

	  if( (popul.N[0]==0) && (popul.N[1]==0) )
	    {
	      cont_extinction++;
	      break;
	    }	  
	  
	}

      if( ((cont_ind[0]>100) && (cont_ind[1]<100)) || ((cont_ind[0]<100) && (cont_ind[1]>100)) )
	  cont_rescue_partial++;

    }

 
 
  sprintf(arq2,"Metapopulation-K%d-Lambda%g-ntraits%d-U%g-mig%g-drop%g-Wmax%g.dat",popul.K,popul.lambda,popul.ntraits,popul.mut,popul.migration,drop_fitness[1],sequence.W_max);
  ptt2 = fopen(arq2,"a");
  fprintf(ptt2,"%g \t %g  \t %g  \n",drop_fitness[0],((double)cont_extinction/n_landscapes),((double)cont_rescue_total/n_landscapes));
  fclose(ptt2);

}



void fitnessinit(sequencia *sequence, populacao *popul, otimo *optimum)
{
  int i, max, max1, ind, j, k;

  max1 = (int) pow(2.,(double)(sequence->L));
  if( sequence->L==32 )
    max1 = RAND_MAX;

  for( i=0; i<popul->ntraits; i++ )
    sequence->trait_init[i] = 0.;

  for( k=0; k<popul->n_subpopul; k++ )
    for( i=0; i<popul->N[k]; i++ )
      {
	for( j=0; j<popul->ntraits; j++ )
	  sequence->traits[k][i][j] = sequence->trait_init[j];
	//	sequence->nmut[k][i] = 0;
      }  
  
}


void dynamics(sequencia *sequence, populacao *popul, otimo *optimum)
{
  int i, k, m, contbin, kinf, aleat, dig, oper, j, bit, verif, k1, cont, nlabelaux, size, indtau, *ndel, auxind, *ind1, pos, oper1, ind_pilha[2], max1, nmut, n_offspring, **pilha, seq_ind, deme;

  unsigned int *N_newborn;
     
  double x, soma, aux, soma1, aux1, mut, Fit, *auxF, malp, mutb, mb, r, Fmax, fitness, saux, F_tot, distance_2, med, a, factor, sum_2, y[100], norm, r_mut;

  double *prob;

  max1 = (int) pow(2.,(double)(sequence->L));
  if( sequence->L==32 )
    max1 = RAND_MAX;

  for( deme=0; deme<popul->n_subpopul; deme++ )
    {
      cont = 0;
      for( i=0; i<popul->N[deme]; i++ )
	{
	  fitness = sequence->fitness[deme][i];
	  n_offspring = gsl_ran_poisson(gerador,fitness);

	  for( j=0; j<n_offspring; j++ )
	    {
	      for( k=0; k<popul->ntraits; k++ )
		sequence->traits_new[cont][k] = sequence->traits[deme][i][k];
	      // sequence->nmut1[cont] = sequence->nmut[deme][i];
	      r = gsl_ran_flat(gerador, 0., 1.);
	      if( (r < popul->mut) )
		{
		  r_mut = gsl_ran_exponential(gerador, popul->lambda);
		  sum_2 = 0;
		  for( k=0; k<popul->ntraits; k++ )
		    {
		      y[k] = gsl_ran_gaussian(gerador, 1.);
		      sum_2 +=  y[k]*y[k];
		    }
		  norm = pow(sum_2,0.5);
		  for( k=0; k<popul->ntraits; k++ )
		    sequence->traits_new[cont][k] += r_mut*y[k]/norm;
		  //  sequence->nmut1[cont]++;
		}
	      cont++;
	    }

	}
 
      popul->N[deme] = cont;
      for( i=0; i<popul->N[deme]; i++ )
	{
	  // sequence->nmut[deme][i] = sequence->nmut1[i];
	  for( j=0; j<popul->ntraits; j++ )
	    sequence->traits[deme][i][j] = sequence->traits_new[i][j];
	}

    }
  

}




void migration(sequencia *sequence, populacao *popul, otimo *optimum)
{
  int NDinc[2], nmig[2], k, i, j, ind_group, k_aux, C[20], indice, n_ind, n_migration;

  double med, prob_survival, r, ***trait_mig;

  trait_mig = new double**[2];
  for( i=0; i<2; i++ )
    trait_mig[i] = new double*[1000];
  for( i=0; i<2; i++ )
    for( j=0; j<(1000); j++ )
      trait_mig[i][j] = new double[popul->ntraits];

  for( i=0; i<2; i++ )
    NDinc[i] = nmig[i] = 0;

  C[0] = 1;
  C[1] = 0;

  for( i=0; i<2; i++ )
    {
      n_migration = gsl_ran_poisson(gerador,(popul->migration*popul->N[i]));

      if( n_migration>popul->N[i] )
	n_migration = popul->N[i];
      
      n_ind = 0;

      while( n_ind<n_migration )
	{
	  indice = gsl_ran_flat(gerador, 0., 1.)*popul->N[i];
	  for( j=0; j<popul->ntraits; j++ )
	    trait_mig[C[i]][NDinc[C[i]]][j] = sequence->traits[i][indice][j];
	  NDinc[C[i]]++;
	  n_ind++;
	  
	  for( k=indice; k<(popul->N[i]-1); k++ )
	    for( j=0; j<popul->ntraits; j++ )
	      sequence->traits[i][k][j] = sequence->traits[i][k+1][j];
	  popul->N[i]--;
	}
      
    }
      
 
  for( i=0; i<2; i++ )
    {
      for( j=0; j<NDinc[i]; j++ )
	for( k=0; k<popul->ntraits; k++ )
	  sequence->traits[i][popul->N[i]+j][k] = trait_mig[i][j][k];

      popul->N[i] += NDinc[i];
    }

 
  for( i=0; i<2; i++ )
    for( j=0; j<1000; j++ )
      delete[] trait_mig[i][j];
  for( i=0; i<2; i++ )
    delete[] trait_mig[i];
  delete[] trait_mig;
  
  
}


void fitness_evaluation(sequencia *sequence, populacao *popul, otimo *optimum)
{
  int k, indice, max1, deme;
  double delta_r[100], distance, factor;

  factor = sequence->W_max;
  for( deme=0; deme<popul->n_subpopul; deme++ )
    for( indice=0; indice<popul->N[deme]; indice++ )
      {
	distance = 0;
	for( k=0; k<popul->ntraits; k++ )
	  distance += (optimum->traits[deme][k]-sequence->traits[deme][indice][k])*(optimum->traits[deme][k]-sequence->traits[deme][indice][k]);

	sequence->fitness[deme][indice] = factor*exp(-distance/(2*popul->alpha_r*popul->alpha_r));
      }  
  
}



#include <vector>
#include <algorithm>
#include <ctime>

// Assuming appropriate definitions for sequencia, populacao, and otimo

void density(sequencia *sequence, populacao *popul, otimo *optimum)
{
    int j, deme;

    for (deme = 0; deme < popul->n_subpopul; deme++)
      {
        if (popul->N[deme] > popul->K)
	  {
            // Create an index vector for shuffling
            std::vector<int> indices(popul->N[deme]);
            for (int i = 0; i < popul->N[deme]; ++i) indices[i] = i;

            // Initialize random seed
            std::srand(unsigned(std::time(0)));

            // Random shuffle the indices
            std::random_shuffle(indices.begin(), indices.end());

            // Create new fitness and traits arrays with resized population
            std::vector<double> new_fitness(popul->K);
            std::vector<std::vector<double>> new_traits(popul->K, std::vector<double>(popul->ntraits));

	    //  for (int i = 0; i < popul->N[deme]; ++i)
	    //  cout << i << "\t" << sequence->fitness[deme][i] << endl;

            for (int i = 0; i < popul->K; ++i)
	      {
                new_fitness[i] = sequence->fitness[deme][indices[i]];
                for (j = 0; j < popul->ntraits; j++)
		  {
		    new_traits[i][j] = sequence->traits[deme][indices[i]][j];
		  }
	      }
	    
	    
            // Copy back the trimmed population
            for (int i = 0; i < popul->K; ++i)
	      {
                sequence->fitness[deme][i] = new_fitness[i];
                for (j = 0; j < popul->ntraits; j++)
		  {
		    sequence->traits[deme][i][j] = new_traits[i][j];
		  }
	      }
	    	    
            popul->N[deme] = popul->K;
	    //	    for (int i = 0; i < popul->N[deme]; ++i)
	    //   cout << i << "\t" << sequence->fitness[deme][i] << endl;
	    
	    
	  }
      }
}
