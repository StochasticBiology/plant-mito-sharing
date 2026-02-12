#include <stdio.h>
#include <stdlib.h>

#define RND drand48()

#define NMITO 200    // number of mitochondria
#define MAXT 1000   // timescale
#define NSAMP 100   // number of samples per expt
#define NSUBS 5
#define DNA NSUBS
#define NCHEMS (NSUBS+1)
#define MAXCHEMS 100
#define NREP 10

int NFUSE;       // number of fusion events per timestep
int LIFE;        // subunit lifespan
int EXPRESSION;  // expression rate
int IMPORT;      // nucleus-mito import rate

// different output statistics for a mitochondrion
#define QUERY_COMPLEX 0
#define QUERY_NUCLEOPROTEIN 1
#define QUERY_MULTIPLE_DNA_COMPLEX 2
#define QUERY_DNA 3

// different rules for targetted protein import
#define TARGET_RANDOM 0
#define TARGET_DNA_BEARING 1
#define TARGET_NO_NUCLEOPROTEIN 2
#define TARGET_DNA_NO_COMPLEX 3
#define TARGET_RANDOM_BATCH 4

// different rules for sharing content upon fusion
#define SWAP_NONE 0
#define SWAP_SUBUNITS 1
#define SWAP_SUBUNITS_DNA 2
#define SWAP_NUCLEOPROTEIN 3
#define SWAP_COMPLEXES 4

// social rules
#define SOCIAL_MTDNA 0
#define SOCIAL_ALL_MITOS 1

// structure for enzyme and metabolite content of a compartment (cytosol or mito)
typedef struct tagCompartment
{
  int copies[NCHEMS];
  int birthdates[NCHEMS][MAXCHEMS];
} Compartment;

// dynamic parameters
typedef struct tagParams
{
  double cytomito;        // cytosol-mito exchange
  double bind, unbind;    // PQR complex binding/unbinding rates
  double rate;            // metabolic reaction rate
  int complexneeded;      // is PQR complex (as opposed to colocalised constituents) needed for reaction?
  int enzymes;            // count of each enzyme across cell
  double nfuse;              // number of fusion events per timestep
} Params;

void Empty(Compartment *C)
{
  int i;
  for(i = 0; i < NCHEMS; i++)
    C->copies[i] = 0;
}

// qtype == 0: number of complexes, == 1: number of DNA-complexes, == 2: number of DNAs if at least one complex
int Query(Compartment C, int qtype)
{
  int i;
  int min = MAXCHEMS;

  if(qtype == QUERY_DNA) return C.copies[DNA];
  for(i = 0; i < (qtype == QUERY_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
    {
      if(C.copies[i] < min) min = C.copies[i];
    }
  if(qtype == QUERY_MULTIPLE_DNA_COMPLEX) return (min > 0 && C.copies[DNA] > 1 ? C.copies[DNA] : 0);
  return min;
}

void Output(Compartment C)
{
  int i, j;

  for(i = 0; i < NCHEMS; i++)
    {
      printf("%i (", C.copies[i]);
      for(j = 0; j < C.copies[i]; j++) printf("%i ", C.birthdates[i][j]);
      printf("), ");
    }
  printf("--> %i %i %i\n", Query(C, QUERY_COMPLEX), Query(C, QUERY_NUCLEOPROTEIN), Query(C, QUERY_MULTIPLE_DNA_COMPLEX));
}

void Pop(Compartment *C, int chem, int ref)
{
  if(ref > C->copies[chem]) return;
  int i;
  for(i = ref; i < C->copies[chem]-1; i++)
    {
      C->birthdates[chem][i] = C->birthdates[chem][i+1];
    }
  C->copies[chem] -= 1;
}

void Push(Compartment *C, int chem, int birthdate)
{
  if(C->copies[chem] > MAXCHEMS-1) return;
  int i = C->copies[chem];
  C->birthdates[chem][i] = birthdate;
  C->copies[chem] += 1;
}

void Transfer(Compartment *C1, Compartment *C2, int chem, int ref)
{
  if(ref >= C1->copies[chem]) return;
  int birthdate = C1->birthdates[chem][ref];
  Pop(C1, chem, ref);
  Push(C2, chem, birthdate);
}

void Mix(Compartment *C1, Compartment *C2, int mtype)
{
  // case 0: mix all subunits and DNA randomly
  // case 1: mix subunits randomly, freeze DNA
  // case 2: mix DNA-complexes randomly, leftover subunits and DNA too
  // case 3: mix complexes and leftovers randomly, freeze DNA
  
  Compartment t1, t2;
  int i, j;
  int q1, q2;
  int dnain, dnaout;

  if(mtype == SWAP_NONE) return;
  
  Empty(&t1);
  Empty(&t2);

  dnain = C1->copies[DNA]+C2->copies[DNA];
  //Output(*C1);
  //Output(*C2);
  // if we're mixing complexes or DNA complexes, first compute how many we have and partition them intact
  if(mtype == SWAP_COMPLEXES || mtype == SWAP_NUCLEOPROTEIN)
    {
      q1 = Query(*C1, (mtype == SWAP_NUCLEOPROTEIN ? QUERY_NUCLEOPROTEIN : QUERY_COMPLEX));
      q2 = Query(*C2, (mtype == SWAP_NUCLEOPROTEIN ? QUERY_NUCLEOPROTEIN : QUERY_COMPLEX));

      //   printf("%i %i\n", q1, q2);
      
      // first partition complexes from C1
      for(j = 0; j < q1; j++)
	{
	  if(RND < 0.5)
	    {
	      for(i = 0; i < (mtype == SWAP_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
		Push(&t1, i, C1->birthdates[i][j]);
	    }
	  else
	    {
	      for(i = 0; i < (mtype == SWAP_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
		Push(&t2, i, C1->birthdates[i][j]);
	    }
	}
      // then partition complexes from C2
      for(j = 0; j < q2; j++)
	{
	  if(RND < 0.5)
	    {
	      for(i = 0; i < (mtype == SWAP_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
		Push(&t1, i, C2->birthdates[i][j]);
	    }
	  else
	    {
	      for(i = 0; i < (mtype == SWAP_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
		Push(&t2, i, C2->birthdates[i][j]);
	    }
	}
      // remove the elements correspond to those complexes we moved, so we can use the followup code to partition the leftovers
      for(i = 0; i < (mtype == SWAP_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
	{
	  for(j = 0; j < q1; j++)
	    Pop(C1, i, 0);
	  for(j = 0; j < q2; j++)
	    Pop(C2, i, 0);
	}

      /*      printf("Here q1 was %i q2 was %i\n", q1, q2);
	      Output(*C1);
	      Output(*C2);
	      Output(t1);
	      Output(t2);
	      printf("end\n");*/
    }

  // go through biomolecules independently
  for(i = 0; i < (mtype == SWAP_SUBUNITS_DNA || mtype == SWAP_NUCLEOPROTEIN ? NCHEMS : NSUBS); i++)
    {
      // for each one, choose which daughter mito to put it in
      for(j = 0; j < C1->copies[i]; j++)
	{
	  if(RND < 0.5) Push(&t1, i, C1->birthdates[i][j]);
	  else Push(&t2, i, C1->birthdates[i][j]);

	  // if(RND < 0.5) { printf("1 to 1\n"); Push(&t1, i, C1->birthdates[i][j]);}
	  //else {printf("1 to 2\n"); Push(&t2, i, C1->birthdates[i][j]); }
	}
      // printf("- %i, %i in 1, %i in 2\n", j, t1.copies[NSUBS], t2.copies[NSUBS]);
      for(j = 0; j < C2->copies[i]; j++)
	{
	  if(RND < 0.5) Push(&t1, i, C2->birthdates[i][j]);
	  else Push(&t2, i, C2->birthdates[i][j]);
	  //	        	  if(RND < 0.5) { printf("2 to 1\n"); Push(&t1, i, C2->birthdates[i][j]);}
	  // else {printf("2 to 2\n"); Push(&t2, i, C2->birthdates[i][j]); }

	}
      //      printf("- %i, %i in 1, %i in 2\n", j, t1.copies[NSUBS], t2.copies[NSUBS]);
    }
  // if we shared DNA, we've already done so. otherwise reconstruct the original DNA profiles
  if(mtype == SWAP_SUBUNITS || mtype == SWAP_COMPLEXES)
    {
      i = DNA;
      for(j = 0; j < C1->copies[i]; j++)
	Push(&t1, i, C1->birthdates[i][j]);
      for(j = 0; j < C2->copies[i]; j++)
	Push(&t2, i, C2->birthdates[i][j]);
    }

  *C1 = t1; *C2 = t2;

  // printf("---, %i in 1, %i in 2\n", C1->copies[NSUBS], C1->copies[NSUBS]);
	       	       
  //  Output(*C1);
  //Output(*C2);
  dnaout = C1->copies[DNA]+C2->copies[DNA];
  if(dnain != dnaout)
    {
      printf("Warning: changed DNA!\n");
      exit(0);
    }
}

void Decay(Compartment *C, int t)
{
  int i, j;
  for(i = 0; i < NSUBS; i++)
    {
      for(j = 0; j < C->copies[i]; j++)
	{
	  if(C->birthdates[i][j] < t - LIFE)
	    {
	      Pop(C, i, j);
	      j--;
	    }
	}
    }
}

// binomial sampler with probability r
int Count(int n, double r)
{
  int sign, count, i;
  
  sign = 1; count = 0;
  if(n < 0) { sign = -1; n = -n; }
  for(i = 0; i < n; i++)
    {
      if(RND < r) count++;
    }
  return count*sign;
}

// create a compartment from a definition -- just for testing
void Create(Compartment *C, int c1, int c2, int c3, int c4, int c5, int c6, int birthdate)
{
  int i;
  int j;

  C->copies[0] = c1; for(i = 0; i < c1; i++) C->birthdates[0][i] = birthdate;
  C->copies[1] = c2; for(i = 0; i < c2; i++) C->birthdates[1][i] = birthdate;
  if(NSUBS == 5) {
  C->copies[2] = c3; for(i = 0; i < c3; i++) C->birthdates[2][i] = birthdate;
  C->copies[3] = c4; for(i = 0; i < c4; i++) C->birthdates[3][i] = birthdate;
  C->copies[4] = c5; for(i = 0; i < c5; i++) C->birthdates[4][i] = birthdate;
  C->copies[5] = c6; for(i = 0; i < c6; i++) C->birthdates[5][i] = birthdate;
  }
}

void RunTest(void)
{
  Compartment C1, C2;
  int i;

  printf("Testing transfer...\n");
  Create(&C1, 0,0,1,0,0,0, 5);
  Create(&C2, 3,3,3,3,3,3, 7);
  Transfer(&C1, &C2, 2, 0);
  Output(C1); Output(C2);
  printf("\n");
 
  
  printf("Testing mixing 1...\n");
  for(i = 0; i <= 10; i++)
    {
      Create(&C1, 0,0,0,0,0,1, 5);
      Create(&C2, 3,3,3,3,3,3, 7);

      Mix(&C2, &C1, SWAP_NUCLEOPROTEIN);
      printf("-- %i\n", i);
      Output(C1); Output(C2);
      printf("\n");
    }

  printf("Testing mixing 2...\n");
  for(i = 0; i <= 3; i++)
    {
      Create(&C1, 0,2,0,3,0,1, 5);
      Create(&C2, 3,3,4,3,3,4, 7);

      Mix(&C1, &C2, i);
      printf("-- %i\n", i);
      Output(C1); Output(C2);
      printf("\n");
    }

}

int main(void)
{
  Compartment N, *M;
  int i;
  int t;
  int r, k, m1, m2;
  int completes, completes2;
  FILE *fp;
  int TARGET, SWAP, SOCIAL, POISSON;
  int rep;
  int avprot, avdna, avprotempty, avprotfull, maxdna, onedna;
  char fstr[100];
  
  // target: 0 random, 1 only DNA, 2 only without DNA-complex, 3 only with DNA without complex, 4 all in one
  // swap: 0 random subunits, 1 random subunits and DNA, 2 subunit sets, 3 random subunits and DNA only for mitos with DNA; 4 no action
  // questions:
  // does targetting help?
  // does targetting mitochondria with DNA help? (target == 1)
  // does targetting mitochondria without a functioning DNA-complex help? (target == 2)
  // does targetting mitochondria with DNA but without a functioning complex help? (target == 3)
  // does fission-fusion help?
  // does exchanging DNA help? (swap == 1)
  // does exchanging subunit clusters help? (swap == 2)

  // control params: expression rate, import rate, protein lifetime, fusion rate, number of mitos, number of DNAs
  // compare to: average proteins per mito (2-20), average DNAs per mito (0.5), timescale of spread through chondriome (hours), protein lifespan (days)
  // and distributions thereof!

  srand48(112);
  RunTest();
  //return 0;

  // say we have 10min as a time unit
  // NFUSE fusions -> each mito undergoes 2*NFUSE/NMITO fusions per unit time, so connected in NMITO/(2*NFUSE) timesteps
  // -> require NMITO/(2*NFUSE) = 1 hour -> NFUSE = NMITO/2 / 1hr -> NFUSE = NMITO/2 / 6units = NMITO / 12
  // i.e. each mito has 1/6 fusion events per 10min unit
  // LIFE = 30h = 180 units


  // say we have 1hr as a time unit
  // NFUSE fusions -> each mito undergoes 2*NFUSE/NMITO fusions per unit time, so connected in NMITO/(2*NFUSE) timesteps
  // -> require NMITO/(2*NFUSE) = 1 hour -> NFUSE = NMITO/2 / 1hr -> NFUSE = NMITO/2 / 1units = NMITO / 2
  // i.e. each mito has 1 fusion events per 1hr unit
  // LIFE = 30h = 30 units

  M = (Compartment*)malloc(sizeof(Compartment)*NMITO);

  sprintf(fstr, "sim-out-new-%i.csv", NSUBS);
  fp = fopen(fstr, "w");
  fprintf(fp, "poisson,target,swap,social,nfuse,life,expression,import,rep,t,completes,completes2,avprot,avdna,avprotempty,avprotfull,maxdna,onedna\n");

  for(POISSON = 0; POISSON <= 1; POISSON++)
    {
      for(TARGET = 0; TARGET <= 4; TARGET++)
	{
	  for(SOCIAL = 0; SOCIAL <= 1; SOCIAL++)
	    {
	      for(SWAP = 0; SWAP <= 4; SWAP++)
		{
		  NFUSE = NMITO/2;
		  // for(NFUSE = 0; NFUSE <= NMITO/2; NFUSE += NMITO/2)
		  {
		    LIFE = 30;
		    //	      for(LIFE = 1; LIFE < 100; LIFE *= 2)
		    {
		      for(EXPRESSION = 100; EXPRESSION <= 1000; EXPRESSION *= 10)
			{
			  IMPORT = EXPRESSION;
			  //		      for(IMPORT = 100; IMPORT <= 1000; IMPORT *= 2)
			  {
			    printf("%i,%i,%i,%i,%i,%i,%i\n", TARGET, SWAP, SOCIAL, NFUSE, LIFE, EXPRESSION, IMPORT);
			    for(rep = 0; rep < NREP; rep++)
			      {
				// empty nucleus
				for(i = 0; i < NMITO; i++)
				  Empty(&N);
				// empty all mitos
				for(i = 0; i < NMITO; i++)
				  Empty(&(M[i]));
				// push mtDNA into some mitos
				for(i = 0; i < NMITO/2; i++)
				  {
				    r = RND*NMITO;
				    if(POISSON == 0)
				      Push(&(M[i]), DNA, 0);
				    else
				      Push(&(M[r]), DNA, 0);
				  }
				for(t = 0; t <= MAXT; t++)
				  {
				    // produce new subunits
				    for(i = 0; i < EXPRESSION; i++)
				      {
					r = RND*NSUBS;
					Push(&N, r, t);
				      }
				    //printf("  import\n");
				    // transfer random nuclear content to random mito
				    m1 = RND*NMITO;
				    for(i = 0; i < IMPORT; i++)
				      {
					r = RND*NSUBS;
					k = RND*N.copies[r];
					completes = 0;
					switch(TARGET)
					  {
					    // just pick a random mito
					  case 0: m1 = RND*NMITO; break;
					    // pick a random mito with DNA
					  case 1: 
					    do{
					      m1 = RND*NMITO;
					    }while(!(M[m1].copies[DNA] != 0)); break;
					    // pick a random mito without a DNA-complex
					  case 2:
					    do{
					      m1 = RND*NMITO;
					      completes++;
					    }while(!(Query(M[m1], QUERY_NUCLEOPROTEIN) == 0) && completes < 10); break;
					    // pick a random mito with DNA and without a DNA-complex
					  case 3:
					    do{
					      m1 = RND*NMITO;
					      completes++;
					    }while(!(Query(M[m1], QUERY_NUCLEOPROTEIN) == 0 && Query(M[m1], QUERY_DNA) != 0) && completes < 10); break;
					    // same mito gets all content
					  case 4: break;
					  }
					if(completes < 10)
					  Transfer(&N, &(M[m1]), r, k);
				      }
				    // fuse and exchange
				    //printf("  fusion\n");

				    for(i = 0; i < NFUSE; i++)
				      {
					/// first choose the mitos
					if(SOCIAL == SOCIAL_ALL_MITOS)
					  {
					    // just choose random mitos
					    do{ 
					      m1 = RND*NMITO;
					      m2 = RND*NMITO;
					    }while(m1 == m2);
					  }
					else
					  {
					    // choose random mitos bearing mtDNA
					    do{
					      m1 = RND*NMITO;
					      m2 = RND*NMITO;
					    }while(m1 == m2 || M[m1].copies[DNA] == 0 || M[m2].copies[DNA] == 0);
					  }
					//// then choose what to exchange
					//    Output(M[m1]);
					//Output(M[m2]);
					//printf("%i %i: %i %i\n", i, SWAP, m1, m2);
					// running out of mitos with DNA?
					Mix(&(M[m1]), &(M[m2]), SWAP);
					//		      Output(M[m1]);
					// Output(M[m2]);
		
				      }
				    //printf("  decay\n");
				    // decay old subunits
				    Decay(&N, t);
				    for(i = 0; i < NMITO; i++)
				      Decay(&(M[i]), t);
      
				    // output state
				    //Query(N);
				    if(t == 100 || t == 900 || t == 1000)
				      {
					completes = completes2 = avprot = avdna = avprotempty = avprotfull = maxdna = onedna = 0;
					for(i = 0; i < NMITO; i++)
					  {
					    //	  printf("  ");
					    completes += Query(M[i], QUERY_NUCLEOPROTEIN);
					    completes2 += Query(M[i], QUERY_MULTIPLE_DNA_COMPLEX);
					    avprot += M[i].copies[0];
					    avdna += M[i].copies[DNA];
					    avprotempty += (M[i].copies[DNA] == 0 ? M[i].copies[0] : 0);
					    avprotfull  += (M[i].copies[DNA] != 0 ? M[i].copies[0] : 0);
					    onedna += (M[i].copies[DNA] == 1);
					    if(M[i].copies[DNA] > maxdna) maxdna = M[i].copies[DNA];
					  }
				    
					fprintf(fp, "%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%i,%.3f\n", POISSON, TARGET, SWAP, SOCIAL, NFUSE, LIFE, EXPRESSION, IMPORT, rep, t, completes, completes2, (float)avprot/NMITO, (float)avdna/NMITO, (float)avprotempty/NMITO, (float)avprotfull/NMITO, maxdna, (float)onedna/NMITO);
				      }
				  }
			      }
			   }
			}
		     } 
		  }
	       }
	    }
	}
    }
  
  fclose(fp);
  
  return 0;
}
