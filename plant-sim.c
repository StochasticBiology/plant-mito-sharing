#include <stdio.h>
#include <stdlib.h>

#define RND drand48()

#define OUTPUTSTATE 0

// max number of mitos
#define _NM 200
// max number of mt-encoded genes
#define _NG 100
// timesteps for simulation (currently minutes)
#define MAXT 300    
// number of samples for each parameter set
#define NEXPT 4

int NM, NG;

// structure describing a mitochondrion's DNA/RNA/protein complement
typedef struct {
  int DNA[_NG];
  int mRNA[_NG];
  int protein[_NG];
} Mito;

// structure containing parameters of model
// these are "probabilities of this event per timestep"
typedef struct {
  double transcription;
  double translation;
  double mrnadegradation;
  double burstsize;
  double loneproteindegradation;
  double jointproteindegradation;
  double socialencounter;
  double dnaexchange;
  double mrnaexchange;
  double proteinexchange;
} Params;

// set parameter values
// current estimate with minute timesteps
void SetParams(Params *P)
{
  // average copy number around 2 (given decay rate)
  P->transcription = 2./60;
  // not sure  
  P->translation = 10./60;
  P->burstsize = 1;
  // 4hr half life
  P->mrnadegradation = 0.25/60;
  // 1hr half life
  P->loneproteindegradation = 1./60;
  // >>1hr half life
  P->jointproteindegradation = 0.01/60;
  // mixing 1-2hr
  P->socialencounter = 30./60;
  
  // probabilities of exchange given that an encounter has occurred
  P->dnaexchange = 0.1;
  P->mrnaexchange = 1.;
  P->proteinexchange = 1.;
}

// initialise mitochondrial population
// everyone gets one gene's worth of DNA, and no mRNA or protein
void Initialise(Mito *M)
{
  int i, j;
  for(i = 0; i < NM; i++)
    {
      for(j = 0; j < NG; j++)
	{
	  if(NG > NM) 
	    M[i].DNA[j] = (i == j%NM);
	  else
	    M[i].DNA[j] = (i%NG == j);
	  M[i].mRNA[j] = 0;
	  M[i].protein[j] = 0;
	}
    }
}

// convenience function for printing aligned integers of different sizes
void myprint(int i)
{
  if(i < 10) printf("%i   ", i);
  else if(i < 100) printf("%i  ", i);
  else if(i < 1000) printf("%i ", i);
  else printf("%i", i);
}

int main(void)
{
  Mito M[_NM], newM[_NM];
  Params P;
  int t;
  int i, j, k;
  int min, ref, partner;
  int socials[_NM];
  FILE *fp;
  int expt;
  
  double propzero, meancomplex, meanprotein, meanmrna;
  
  // set parameters and initialise population
  SetParams(&P);
  
  fp = fopen("k2-out.csv", "w");
  printf("social,n.mito,n.gene,expt,mean.complex,mean.protein,mean.rna,prop.zero\n");
  fprintf(fp, "social,n.mito,n.gene,expt,mean.complex,mean.protein,mean.rna,prop.zero\n");

  for(P.socialencounter = 0; P.socialencounter < 1.1/60; P.socialencounter += 0.25/60)
    {
      for(NM = 10; NM < 200; NM += 10)
	{
	  for(NG = 1; NG < 20; NG ++)
	    {
	      for(expt = 0; expt < NEXPT; expt++)
		{
		  Initialise(M);

		  // loop through timescale of simulation
		  for(t = 0; t < MAXT; t++)
		    {
		      // output current state of simulation
		      // printf("==== %i, %i, Time %i --\n", NM, NG, t);
		      if(OUTPUTSTATE)
			{
			  for(i = 0; i < NM; i++)
			    {
			      printf("Mito %i --\n", i);
			      for(j = 0; j < NG; j++)
				printf("%i   ", M[i].DNA[j]);
			      printf("\n");
			      for(j = 0; j < NG; j++)
				myprint(M[i].mRNA[j]);
			      printf("\n");
			      for(j = 0; j < NG; j++)
				myprint(M[i].protein[j]);
			      printf("\n\n");
	    
			    }
			}
	      
		      // first we'll do gene expression and degradation dynamics internal to each mitochondrion
		      // loop over mitochondria
		      for(i = 0; i < NM; i++)
			{
			  // set the mito in the next timestep to initially be equal to this one now
			  newM[i] = M[i];

			  // we want to find how many unbound (out-of-complex) protein subunits we have
			  // so we identify the gene that has the minimal protein count, and say there are this many intact complexes
			  // e.g. for protein counts 2, 5, 3, 10 we say there are 2 intact complexes and the remaining 0, 3, 1, 8 are unbound
			  min = -1;
			  for(j = 0; j < NG; j++)
			    {
			      if(min == -1 || M[i].protein[j] < min)
				min = M[i].protein[j];
			    }
			  // now do expression and degradation gene-by-gene
			  for(j = 0; j < NG; j++)
			    {
			      // transcribe if we have DNA for this gene
			      if(M[i].DNA[j])
				{
				  if(RND < P.transcription)
				    newM[i].mRNA[j]++;
				}
			      // translate and/or degrade this mRNA
			      for(k = 0; k < M[i].mRNA[j]; k++)
				{
				  if(RND < P.translation)
				    newM[i].protein[j] += P.burstsize;
				  if(RND < P.mrnadegradation)
				    newM[i].mRNA[j]--;
				}
			      // degrade this protein, with rate depending on whether it's unbound or not
			      for(k = 0; k < M[i].protein[j]; k++)
				{
				  if(M[i].protein[j] > min)
				    {
				      // this protein has subunits floating around that aren't in one of the "min" complexes -- degrade these preferentially
				      if(RND < P.loneproteindegradation)
					newM[i].protein[j]--;
				    }
				  else
				    {
				      // this protein only has subunits embedded in complexes, so apply this rate
				      if(RND < P.jointproteindegradation)
					newM[i].protein[j]--;
				    }
				}
			    }
			}
		      // update system state and initialise social interactions across mitos
		      for(i = 0; i < NM; i++)
			{
			  M[i] = newM[i];
			  socials[i] = 0;
			}

		      // now we consider exchange between mitos
		      // loop through each mito
		      for(i = 0; i < NM; i++)
			{
			  // if this one hasn't had a social exchange in this timestep, it can do so
			  if(socials[i] == 0 && RND < P.socialencounter)
			    {
			      // find a partner that also hasn't exchanged yet (and isn't this mito)
			      do{
				partner = RND*NM;
			      }while(i == partner || socials[partner] == 1);
			      //printf("%i meets %i\n", i, partner);
			      // record this interaction
			      socials[i] = socials[partner] = 1;
			      // loop over genes
			      for(j = 0; j < NG; j++)
				{
				  // exchange DNA 
				    if(M[i].DNA[j] && RND < P.dnaexchange)
				      {
					ref = 0;
					for(k = 0; k < NG; k++)
					  ref = (M[partner].DNA[k] == 1 ? k : ref);
					M[i].DNA[ref] = 1; M[i].DNA[j] = 0;
					M[partner].DNA[j] = 1; M[partner].DNA[ref] = 0;
				      }
				  // aggregate and partition mRNA
				  if(RND < P.mrnaexchange)
				    {
				      ref = (M[i].mRNA[j]+M[partner].mRNA[j]);
				      M[i].mRNA[j] = RND*ref;
				      M[partner].mRNA[j] = ref-M[i].mRNA[j];
				    }
				  // aggregate and partition proteins
				  if(RND < P.proteinexchange)
				    {
				      ref = (M[i].protein[j]+M[partner].protein[j]);
				      M[i].protein[j] = RND*ref;
				      M[partner].protein[j] = ref-M[i].protein[j];
				    }
				}
			    }
			}
		    }
		  propzero = meancomplex = meanprotein = meanmrna = 0;
		  for(i = 0; i < NM; i++)
		    {
		      min = -1;
		      for(j = 0; j < NG; j++)
			{
			  if(min == -1 || M[i].protein[j] < min)
			    min = M[i].protein[j];
			  meanprotein += M[i].protein[j];
			  meanmrna += M[i].mRNA[j];
			}
		      if(min == 0) propzero++;
		      meancomplex += min;
		    }
		  meanprotein /= NM*NG;
		  meancomplex /= NM*NG;
		  meanmrna /= NM*NG;
		  propzero /= NM;
		  printf("%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f\n", P.socialencounter, NM, NG, expt, meancomplex, meanprotein, meanmrna, propzero);
		  fprintf(fp, "%.3f,%i,%i,%i,%.3f,%.3f,%.3f,%.3f\n", P.socialencounter, NM, NG, expt, meancomplex, meanprotein, meanmrna, propzero);
    
		}
	    }
	}
    }
  //  fclose(fp);
  
  return 0;
}
			
