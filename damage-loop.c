#include <stdio.h>
#include <stdlib.h>

#define RND drand48()

// experiment types
// 0 -- BER pathway where we don't require all proteins together simultaneously
// 1 -- TR pathway
// 2 -- BER pathway where we do require all proteins together simultaneously
#define EXPT_NORMAL 0
#define EXPT_TEMPLATE 1
#define EXPT_COMPLEX 2

// ----- change this to run different experiments
// -- for EXPT_NORMAL and EXPT_COMPLEX we have NPROTS = 5, requiring monomers
// -- for EXPT_TEMPLATE we have NPROTS = 1, requiring dimers
int EXPT, NPROTS;
int CHOOSE_RANDOM;

//#define EXPT EXPT_COMPLEX
//#define NPROTS (EXPT == EXPT_TEMPLATE ? 1 : 5)

#define NMITO 300    // number of mitochondria
#define MAXT 1000   // timescale
#define NSAMP 100   // number of samples per expt
#define MAXTYPES 5
#define MAXPROTS 100
#define MAXDNA 200
#define NREP 30

// different output statistics for a mitochondrion
#define QUERY_COMPLEX 0
#define QUERY_NUCLEOPROTEIN 1
#define QUERY_TEMPLATE_CAPACITY 2
#define QUERY_DNA 3
#define QUERY_FREE_SUBUNITS 4
#define QUERY_DAMAGE 5

// different rules for targetted protein import
#define TARGET_RANDOM 0
#define TARGET_DNA_BEARING 1
#define TARGET_NO_NUCLEOPROTEIN 2
#define TARGET_DNA_NO_COMPLEX 3
#define TARGET_RANDOM_BATCH 4
#define TARGET_DAMAGE 5

// different rules for sharing content upon fusion
#define SWAP_NONE 0
#define SWAP_SUBUNITS 1
#define SWAP_SUBUNITS_DNA 2
#define SWAP_NUCLEOPROTEIN 3
#define SWAP_COMPLEXES 4
#define SWAP_NUCLEOPROTEINS_AND_COMPLEXES 5

// requirement types for molecules
#define REQ_PROTEIN 0
#define REQ_DNA 1

// social rules
#define SOCIAL_MTDNA 0
#define SOCIAL_ALL_MITOS 1

// structure containing mitochondrial (or nuclear) protein and DNA content
// with birthdates (for lifespan) and damage states respectively
// 0 = no damage, nonzero = some step on the repair pathway
typedef struct tagCompartment {
  int proteins[MAXTYPES];
  int birthdates[MAXTYPES][MAXPROTS];
  int DNA;
  int damage[MAXDNA];
} Compartment;

// parameters for experiments
typedef struct tagParams {
  int POISSON;
  float MUT;
  int TARGET;
  int SOCIAL;
  int SWAP;
  int NFUSE;
  float LIFE;
  float EXPRESSION;
  float IMPORT;
  float SCANRATE;
} Params;

// useful statistics of a compartment
typedef struct tagOutput {
  int completes, completes2;
  float avprot, avdna, avprotempty, avprotfull;
  int maxdna;
  float onedna, freesubs;
  int damage;
} Outputs;

// pop a protein from a compartment's set
void Pop(Compartment *C, int chem, int ref)
{
  if(ref >= C->proteins[chem]) {
    printf("Oops, attempt to pop absent protein?\n");
    return;
  }
  int i;
  for(i = ref; i < C->proteins[chem]-1; i++)
    {
      C->birthdates[chem][i] = C->birthdates[chem][i+1];
    }
  C->proteins[chem] -= 1;
}

// push a protein, with birthdate, into a compartment's set
void Push(Compartment *C, int chem, int birthdate)
{
  if(C->proteins[chem] > MAXPROTS-1) {
    // allow quiet fail here
    //printf("Oops, too much protein?\n");
    return;
  }
  int i = C->proteins[chem];
  C->birthdates[chem][i] = birthdate;
  C->proteins[chem] += 1;
}

// push DNA, with damage ref, into a compartment
void PushDNA(Compartment *C, int damage)
{
  if(C->DNA >= MAXDNA)
    {
      printf("Somehow too much DNA!\n");
      exit(0);
    }
  C->damage[C->DNA] = damage;
  (C->DNA)++;
}

// pop DNA from a compartment
void PopDNA(Compartment *C, int ref)
{
  if(ref >= C->DNA) {
    printf("Oops, attempt to pop absent DNA?\n");
    return;
  }
  int i;
  for(i = ref; i < C->DNA-1; i++)
    {
      C->damage[i] = C->damage[i+1];
    }
  C->DNA -= 1;
}

// empty a compartment
void Empty(Compartment *C)
{
  int i;
  for(i = 0; i < NPROTS; i++)
    C->proteins[i] = 0;
  C->DNA = 0;
}

// how many, of which molecule type, are required for a repair step in this experiment
int required(int moltype)
{
  if(EXPT == EXPT_TEMPLATE) {
    if(moltype == REQ_PROTEIN) return 2;
    if(moltype == REQ_DNA) return 2;
  }
  if(moltype == REQ_PROTEIN) return 1;
  if(moltype == REQ_DNA) return 0;
  return 0;
}

// see if this compartment has a protein complement (and DNA?) that helps to fix DNA damage
void Process(Compartment *C, Params P)
{
  int i, j, k;
  
  // loop through DNA in compartment
  for(i = 0; i < C->DNA; i++)
    {
      // loop through steps in damage repair pathway
      for(j = 0; j < NPROTS; j++)
	{
	  // if we are at this step and have sufficient proteins, repair
	  if(C->damage[i] == j+1 && C->proteins[j] >= required(REQ_PROTEIN) && C->DNA >= required(REQ_DNA) )
	    {
	      // if this is the first damage step, give every protein a chance at finding the damage and making the first step
	      if(j+1 == 1) {
		for(k = 0; k < C->proteins[j]; k++)
		  {
		    if(RND < P.SCANRATE)
		      C->damage[i] = j+1 + 1;
		  }
	      } else {
		// otherwise just keep fixing if the protein is there
	        C->damage[i]++;
	      }
	    }
	}
      // if we've reached the end of the pathway, we are done
      if(C->damage[i] == NPROTS+1) C->damage[i] = 0;
      // if we need all steps to occur simultaneously, and this hasn't happened, reset damage state
      // (effectively modelling no activity for this mtDNA)
      if(EXPT == EXPT_COMPLEX && C->damage[i] > 0)
	C->damage[i] = 1;
    }
}

// cause damage to DNA with characteristic rate
void DNADamage(Compartment *C, Params P)
{
  int i;
  
  // randomly mutate DNA
  for(i = 0; i < C->DNA; i++)
    {
      if(C->damage[i] == 0 && RND < P.MUT)
	C->damage[i] = 1;
    }
}

// move a protein from one compartment to another 
void Transfer(Compartment *C1, Compartment *C2, int chem, int ref)
{
  if(ref >= C1->proteins[chem]) return;
  int birthdate = C1->birthdates[chem][ref];
  Pop(C1, chem, ref);
  Push(C2, chem, birthdate);
}

// remove proteins over a given age
void ProteinDecay(Compartment *C, int t, Params P)
{
  int i, j;
  for(i = 0; i < NPROTS; i++)
    {
      for(j = 0; j < C->proteins[i]; j++)
	{
	  if(C->birthdates[i][j] < t - P.LIFE)
	    {
	      Pop(C, i, j);
	      j--;
	    }
	}
    }
}

// ask questions about a compartment's content
/*#define QUERY_COMPLEX 0 - number of complexes
  #define QUERY_NUCLEOPROTEIN 1 - number of nucleoprotein complexes
  #define QUERY_TEMPLATE_CAPACITY 2 - number of DNAs with at least one other DNA and a complex
  #define QUERY_DNA 3 - number of DNAs
  #define QUERY_FREE_SUBUNITS 4 - number of proteins outside complexes
  #define QUERY_DAMAGE 5 - amount of damage */
int Query(Compartment C, int qtype)
{
  int count = 0;
  int i;
  int complexes;
  
  if(qtype == QUERY_DAMAGE)
    {
      for(i = 0; i < C.DNA; i++)
	count += (C.damage[i] > 0);
    }
  if(qtype == QUERY_DNA)
    {
      count = C.DNA;
    }
  if(qtype == QUERY_COMPLEX)
    {
      count = MAXPROTS;
      for(i = 0; i < NPROTS; i++)
	{
	  if(C.proteins[i] < count)
	    count = C.proteins[i];
	}
    }
  if(qtype == QUERY_FREE_SUBUNITS)
    {
      int complexes;
      complexes = MAXPROTS;
      for(i = 0; i < NPROTS; i++)
	{
	  if(C.proteins[i] < complexes)
	    complexes = C.proteins[i];
	}
      for(i = 0; i < NPROTS; i++)
	{
	  count += C.proteins[i] - complexes;
	}
    }

  if(qtype == QUERY_NUCLEOPROTEIN)
    {
      count = MAXPROTS;
      for(i = 0; i < NPROTS; i++)
	{
	  if(C.proteins[i] < count)
	    count = C.proteins[i];
	}
      if(C.DNA < count) count = C.DNA;
    }
  if(qtype == QUERY_TEMPLATE_CAPACITY)
    {
      count = MAXPROTS;
      for(i = 0; i < NPROTS; i++)
	{
	  if(C.proteins[i] < count)
	    count = C.proteins[i];
	}
      if(count >= 2 && C.DNA >= 2) count = C.DNA;
      else count = 0;
    }

  return count;
}

// output some properties of a compartment
void Output(Compartment C, int showtimes)
{
  int i, j;

  for(i = 0; i < NPROTS; i++)
    {
      printf("%i ", C.proteins[i]);
      if(showtimes == 1)
	{
	  printf("(");
	  for(j = 0; j < C.proteins[i]; j++) printf("%i ", C.birthdates[i][j]);
	  printf("), ");
	}
    }
  printf("+ %i: ", C.DNA);
  for(j = 0; j < C.DNA; j++)
    printf("%i, ", C.damage[j]);
  printf("--> %i comp %i DNAcomp %i temps\n", Query(C, QUERY_COMPLEX), Query(C, QUERY_NUCLEOPROTEIN), Query(C, QUERY_TEMPLATE_CAPACITY));
}

// mix contents of two compartments according to a mix rule
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

  dnain = C1->DNA+C2->DNA;
  
  // if we're  DNA complexes, first compute how many we have and partition them intact
  if(mtype == SWAP_NUCLEOPROTEIN || mtype == SWAP_NUCLEOPROTEINS_AND_COMPLEXES)
    {
      q1 = Query(*C1, QUERY_NUCLEOPROTEIN);
      q2 = Query(*C2, QUERY_NUCLEOPROTEIN);

      // first partition complexes from C1
      for(j = 0; j < q1; j++)
	{
	  if(RND < 0.5)
	    {
	      PushDNA(&t1, C1->damage[j]);
	      for(i = 0; i < NPROTS; i++)
		Push(&t1, i, C1->birthdates[i][j]);
	    }
	  else
	    {
	      PushDNA(&t2, C1->damage[j]);
	      for(i = 0; i < NPROTS; i++)
		Push(&t2, i, C1->birthdates[i][j]);
	    }
	}
      // then partition complexes from C2
      for(j = 0; j < q2; j++)
	{
	  if(RND < 0.5)
	    {
	      PushDNA(&t1, C2->damage[j]);
	      for(i = 0; i < NPROTS; i++)
		Push(&t1, i, C2->birthdates[i][j]);
	    }
	  else
	    {
	      PushDNA(&t2, C2->damage[j]);
	      for(i = 0; i < NPROTS; i++)
		Push(&t2, i, C2->birthdates[i][j]);
	    }
	}
      // remove the elements correspond to those complexes we moved, so we can use the followup code to partition the leftovers
      for(j = 0; j < q1; j++)
	{

	  for(i = 0; i < NPROTS; i++)
	    {
	      Pop(C1, i, 0);
	    }
	  PopDNA(C1, 0);
	}
      for(j = 0; j < q2; j++)
	{
	  for(i = 0; i < NPROTS; i++)
	    {
	      Pop(C2, i, 0);
	    }
	  PopDNA(C2, 0);
	}
    }

  // if we're also partitioning complexes, do the same for those that remain
  if(mtype == SWAP_COMPLEXES || mtype == SWAP_NUCLEOPROTEINS_AND_COMPLEXES)
    {
      q1 = Query(*C1, QUERY_COMPLEX);
      q2 = Query(*C2, QUERY_COMPLEX);

      //   printf("%i %i\n", q1, q2);
      
      // first partition complexes from C1
      for(j = 0; j < q1; j++)
	{
	  if(RND < 0.5)
	    {
	      for(i = 0; i < NPROTS; i++)
		Push(&t1, i, C1->birthdates[i][j]);
	    }
	  else
	    {
	      for(i = 0; i < NPROTS; i++)
		Push(&t2, i, C1->birthdates[i][j]);
	    }
	}
      // then partition complexes from C2
      for(j = 0; j < q2; j++)
	{
	  if(RND < 0.5)
	    {
	      for(i = 0; i < NPROTS; i++)
		Push(&t1, i, C2->birthdates[i][j]);
	    }
	  else
	    {
	      for(i = 0; i < NPROTS; i++)
		Push(&t2, i, C2->birthdates[i][j]);
	    }
	}
      // remove the elements correspond to those complexes we moved, so we can use the followup code to partition the leftovers
      for(i = 0; i < NPROTS; i++)
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

  // go through biomolecules independently -- first proteins
  for(i = 0; i < NPROTS; i++) //(mtype == SWAP_SUBUNITS_DNA || mtype == SWAP_NUCLEOPROTEIN || mtype == SWAP_NUCLEOPROTEINS_AND_COMPLEXES ? NPROTS : NPROTS); i++)
    {
      // for each one, choose which daughter mito to put it in
      for(j = 0; j < C1->proteins[i]; j++)
	{
	  if(RND < 0.5) Push(&t1, i, C1->birthdates[i][j]);
	  else Push(&t2, i, C1->birthdates[i][j]);

	  // if(RND < 0.5) { printf("1 to 1\n"); Push(&t1, i, C1->birthdates[i][j]);}
	  //else {printf("1 to 2\n"); Push(&t2, i, C1->birthdates[i][j]); }
	}
      // printf("- %i, %i in 1, %i in 2\n", j, t1.proteins[NPROTS], t2.proteins[NPROTS]);
      for(j = 0; j < C2->proteins[i]; j++)
	{
	  if(RND < 0.5) Push(&t1, i, C2->birthdates[i][j]);
	  else Push(&t2, i, C2->birthdates[i][j]);
	  //	        	  if(RND < 0.5) { printf("2 to 1\n"); Push(&t1, i, C2->birthdates[i][j]);}
	  // else {printf("2 to 2\n"); Push(&t2, i, C2->birthdates[i][j]); }

	}
      //      printf("- %i, %i in 1, %i in 2\n", j, t1.proteins[NPROTS], t2.proteins[NPROTS]);
    }
  // if we're sharing DNA, parititon what's left
  if(!(mtype == SWAP_SUBUNITS || mtype == SWAP_COMPLEXES))
    {
      for(j = 0; j < C1->DNA; j++)
	{
	  if(RND < 0.5) PushDNA(&t1, C1->damage[j]);
	  else PushDNA(&t2, C1->damage[j]);
	}
      for(j = 0; j < C2->DNA; j++)
	{
	  if(RND < 0.5) PushDNA(&t1, C2->damage[j]);
	  else PushDNA(&t2, C2->damage[j]);
	}
    } else  // reconstruct original DNA profiles
    {
      for(j = 0; j < C1->DNA; j++)
	PushDNA(&t1, C1->damage[j]);
      for(j = 0; j < C2->DNA; j++)
	PushDNA(&t2, C2->damage[j]);
    }

  *C1 = t1; *C2 = t2;

  //printf("---, %i in 1, %i in 2\n", C1->DNA, C2->DNA);
	       	       
  //Output(*C1);
  //Output(*C2);
  dnaout = C1->DNA+C2->DNA;
  if(dnain != dnaout)
    {
      printf("Warning: changed DNA!\n");
      exit(0);
    }
}


// create a compartment from a definition -- just for testing
void Create(Compartment *C, int c1, int c2, int c3, int c4, int c5, int c6, int birthdate)
{
  int i;
  int j;

  if(NPROTS == 1)
    {
      C->proteins[0] = c1; for(i = 0; i < c1; i++) C->birthdates[0][i] = birthdate;
      C->DNA = c2; for(i = 0; i < c2; i++) C->damage[i] = birthdate;
    }
  if(NPROTS == 5) {
    C->proteins[0] = c1; for(i = 0; i < c1; i++) C->birthdates[0][i] = birthdate;
    C->proteins[1] = c2; for(i = 0; i < c2; i++) C->birthdates[1][i] = birthdate;
    C->proteins[2] = c3; for(i = 0; i < c3; i++) C->birthdates[2][i] = birthdate;
    C->proteins[3] = c4; for(i = 0; i < c4; i++) C->birthdates[3][i] = birthdate;
    C->proteins[4] = c5; for(i = 0; i < c5; i++) C->birthdates[4][i] = birthdate;
    C->DNA = c6; for(i = 0; i < c6; i++) C->damage[i] = birthdate;

  }
}

// test mixing process
void RunTest(void)
{
  Compartment C1, C2;
  int i, j;
  Params P;

  P.SCANRATE = 1;
  P.MUT = 1;
  
  printf("Testing transfer...\n");
  Create(&C1, 0,0,1,0,0,0, 5);
  Create(&C2, 3,3,3,3,3,3, 7);
  Transfer(&C2, &C1, 0, 0);
  Output(C1, 0); Output(C2, 0);
  printf("\n");
 
  
  /* #define SWAP_NONE 0
     #define SWAP_SUBUNITS 1
     #define SWAP_SUBUNITS_DNA 2
     #define SWAP_NUCLEOPROTEIN 3
     #define SWAP_COMPLEXES 4
     #define SWAP_NUCLEOPROTEINS_AND_COMPLEXES 5
  */

  for(j = 0; j <= 5; j++)
    {

      printf("Testing mixing ");
      switch(j) {
      case 0: printf("none"); break;
      case 1: printf("subs"); break;
      case 2: printf("subs + DNA"); break;
      case 3: printf("nucleoprot"); break;
      case 4: printf("complexes"); break;
      case 5: printf("nucleoprot + complexes"); break;
      }
      printf("...\n");
      for(i = 0; i <= 3; i++)
	{
	  Create(&C1, 0,2,0,3,0,1, 5);
	  Create(&C2, 3,3,4,3,3,4, 7);

	  Mix(&C1, &C2, j);
	  printf("-- %i\n", i);
	  Output(C1, 0); Output(C2, 0);
	  printf("\n");
	}
    }
  printf("Testing repair...\n");
  Create(&C1, 2,2,0,3,0,1, 0);
  Create(&C2, 3,3,4,3,3,1, 0);
  DNADamage(&C1, P);
  DNADamage(&C2, P);
  Output(C1, 0); Output(C2, 0);
  Process(&C1, P); Process(&C2, P);
  Output(C1, 0); Output(C2, 0);
}

// print out some useful statistics
void OutputStats(Outputs O)
{
  printf("  %.3f average protein\n", O.avprot);
  printf("  %i nucleoprot %i templaters\n", O.completes, O.completes2);
  if(O.avdna > 0)
    printf("  %i damage %.3f av DNA %.3f damage percent\n", O.damage, O.avdna, (float)O.damage/(O.avdna*NMITO));
  else
    printf("Somehow no average DNA!\n");
}

// get some useful statistics from a compartment
void GetStats(Compartment *M, Outputs *O)
{
  int i;
  O->completes = O->completes2 = O->avprot = O->avdna = O->avprotempty = O->avprotfull = O->maxdna = O->onedna = O->freesubs = O->damage = 0;
  for(i = 0; i < NMITO; i++)
    {
      O->completes += Query(M[i], QUERY_NUCLEOPROTEIN);
      O->completes2 += Query(M[i], QUERY_TEMPLATE_CAPACITY);
      O->avprot += M[i].proteins[0];
      O->avdna += Query(M[i], QUERY_DNA);
      O->onedna += (Query(M[i], QUERY_DNA) == 1);
      O->damage += Query(M[i], QUERY_DAMAGE);
      if(Query(M[i], QUERY_DNA) > O->maxdna) O->maxdna = Query(M[i], QUERY_DNA);
      O->freesubs += Query(M[i], QUERY_FREE_SUBUNITS);
    }

  O->avprot /= NMITO; O->avdna /= NMITO; O->avprotempty /= NMITO; O->avprotfull /= NMITO; O->onedna /= NMITO; O->freesubs /= NMITO;
}

// run a simulation of mitochondrial sharing
void Simulate(Params P, FILE *fp, int rep, Compartment *Mret, Outputs *Oret)
{
  Compartment N, *M;
  int i;
  int t;
  int avprot, avdna, avprotempty, avprotfull, maxdna, onedna, freesubs, avdamage;
  int r;
  int m1, m2;
  int k;
  int completes, completes2;
  Outputs O;
  int check;
  int posslist[NMITO];
  int nposs;
  
  M = (Compartment*)malloc(sizeof(Compartment)*NMITO);
  
  // empty nucleus
  Empty(&N);
  // empty all mitos
  for(i = 0; i < NMITO; i++)
    Empty(&(M[i]));
  // push mtDNA into some mitos
  for(i = 0; i < NMITO/3; i++)
    {
      r = RND*NMITO;
      if(P.POISSON == 0)
	PushDNA(&(M[i]), 0);
      else
	PushDNA(&(M[r]), 0);
    }
  for(t = 0; t <= MAXT; t++)
    {
      // produce new subunits
      for(i = 0; i < P.EXPRESSION; i++)
	{
	  r = RND*NPROTS;
	  Push(&N, r, t);
	}
      
      // construct a list of mitos for possible import according to import criterion
      nposs = 0;
      switch(P.TARGET)
	{
	  // just pick a random mito
	case TARGET_RANDOM:
	  for(check = 0; check < NMITO; check++)
	    posslist[nposs++] = check;
	  break;
	      
	  // pick a random mito with DNA
	case TARGET_DNA_BEARING:
	  nposs = 0;
	  for(check = 0; check < NMITO; check++)
	    {
	      if(Query(M[check], QUERY_DNA) != 0)
		{
		  posslist[nposs++] = check;
		}
	    }
	  break;
	      
	  // pick a random mito without a DNA-complex
	case TARGET_NO_NUCLEOPROTEIN:
	  nposs = 0;
	  for(check = 0; check < NMITO; check++)
	    {
	      if(Query(M[check], QUERY_NUCLEOPROTEIN) == 0)
		{
		  posslist[nposs++] = check;
		}
	    }
	  break;
	      
	  // pick a random mito with DNA and without a DNA-complex
	case TARGET_DNA_NO_COMPLEX:
	  nposs = 0;
	  for(check = 0; check < NMITO; check++)
	    {
	      if(Query(M[check], QUERY_NUCLEOPROTEIN) == 0 && Query(M[check], QUERY_DNA) != 0)
		{
		  posslist[nposs++] = check;
		}
	    }
	  break;
	      
	  // same mito gets all content
	case TARGET_RANDOM_BATCH:
	  posslist[0] = RND*NMITO;
	  nposs = 1;
	  break;

	  // pick a random mito with damage
	case TARGET_DAMAGE:
	  nposs = 0;
	  for(check = 0; check < NMITO; check++)
	    {
	      if(Query(M[check], QUERY_DAMAGE) != 0)
		{
		  posslist[nposs++] = check;
		}
	    }
	  break;
	}

      //  transfer random nuclear content to mitos randomly chosen from this list
      for(i = 0; i < P.IMPORT; i++)
	{
	  r = RND*NPROTS;
	  if(N.proteins[r] == 0) continue;
	  k = RND*N.proteins[r];

	  if(nposs == 0)
	    m1 = -1;
	  else
	    m1 = posslist[(int)(RND*nposs)];
	  
	  if(m1 == -1 && CHOOSE_RANDOM) 
	    m1 = RND*NMITO;
	  
	  if(m1 != -1)
  	    Transfer(&N, &(M[m1]), r, k);
	}
      // fuse and exchange
      //printf("  fusion\n");

      // make list of fuseable mitos
      nposs = 0;
      if(P.SOCIAL == SOCIAL_ALL_MITOS)
	{
	  for(check = 0; check < NMITO; check++)
	    posslist[nposs++] = check;
	}
      else {
	for(check = 0; check < NMITO; check++)
	  {
	    if(Query(M[check], QUERY_DNA) != 0)
	      {
		posslist[nposs++] = check;
	      }
	  }
      }
			      
      // fuse mitos randomly chosen from this list
      for(i = 0; i < P.NFUSE; i++)
	{
	  if(nposs == 0) {
	    printf("Somehow no mitos with DNA!\n");
	    exit(0);
	  }
	  if(nposs > 1) {
	    do{
	      m1 = posslist[(int)(RND*nposs)];
	      m2 = posslist[(int)(RND*nposs)];
	    }while(m1 == m2);
	  } else m1 = m2 = posslist[0];
	    
	  Mix(&(M[m1]), &(M[m2]), P.SWAP);
    		
	}
      //printf("  decay\n");
      // decay old subunits
      ProteinDecay(&N, t, P);
      for(i = 0; i < NMITO; i++) {
	ProteinDecay(&(M[i]), t, P);
	DNADamage(&(M[i]), P);
	Process(&(M[i]), P);
      }
				    
      // output state
      //Query(N);
      if(t == 100 || t == 900 || t == 1000)
	{ 
	  GetStats(M, &O);
	  if(fp != NULL) {
	    fprintf(fp, "%.3e,%.3e,%i,%i,%i,%i,%.3e,", P.SCANRATE, P.MUT, P.TARGET, P.SWAP, P.SOCIAL, P.NFUSE, P.LIFE);
	    fprintf(fp, "%.3e,%.3e,%i,%i,", P.EXPRESSION, P.IMPORT, rep, t);
	    fprintf(fp, "%i,%i,%.3e,%.3e,%.3e,%.3e,",O.completes, O.completes2, O.avprot, O.avdna, O.avprotempty, O.avprotfull);
	    fprintf(fp, "%i,%.3e,%.3e,%i\n", O.maxdna, O.onedna, O.freesubs, O.damage);
	  }
	}
    }
  for(i = 0; i < NMITO; i++)
    Mret[i] = M[i];

  *Oret = O;
  
  free(M);
}

// routine for testing behaviour under dynamic simulation
void RunSimTest(void)
{
  Compartment *M;
  FILE *fp;
  int i;
  Params P;
  Outputs O;
  
  printf("Testing simulation...\n");

  M = (Compartment*)malloc(sizeof(Compartment)*NMITO);
  fp = NULL;

  P.POISSON = 1; P.MUT = 0.1; P.TARGET = 0; P.SOCIAL = 0; P.SWAP = 0; P.LIFE = 24*7; P.EXPRESSION = 15; P.IMPORT = 15;
  P.NFUSE = NMITO/2;
  Simulate(P, fp, 0, M, &O);
  for(i = 0; i < NMITO; i++)
    Output(M[i], 0);
  OutputStats(O);
  free(M);
  
}

int main(int argc, char *argv[])
{
  Compartment N, *M;
  int i;
  int t;
  FILE *fp;
  int rep;
  char fstr[100];
  int minEXPRESSION, maxEXPRESSION;
  double minMUT, maxMUT;
  Params P;
  Outputs O;
  int EXPTlabel;
  
  if(argc != 3) {
    printf("Which experiment should I run? 0-5 and should I choose random mitos when I can't choose principled? 0-1\n");
    return 0;
  }
  EXPTlabel = atoi(argv[1]);
  if(EXPTlabel < 0 || EXPTlabel > 5) {
    printf("Experiment not recognised. 0-5\n");
    return 0;
  }
  CHOOSE_RANDOM = atoi(argv[2]);
  if(CHOOSE_RANDOM !=  0 && CHOOSE_RANDOM != 1) {
    printf("Choose-random not recognised. 0-1\n");
    return 0;
  }

  // choose expression levels to obtain correct scale of per-mito copy number
  // (4 for 5-protein pathway; 2 for MSH1-like templater)

  minMUT = 0.02; maxMUT = 0.081;
  switch(EXPTlabel) {
  case 0: EXPT = EXPT_NORMAL; NPROTS = 5; minEXPRESSION = 40; maxEXPRESSION = 90; break;
  case 1: EXPT = EXPT_TEMPLATE; NPROTS = 1; minEXPRESSION = 4; maxEXPRESSION = 20; break;
  case 2: EXPT = EXPT_COMPLEX; NPROTS = 5; minEXPRESSION = 40; maxEXPRESSION = 90; break;
  case 3: EXPT = EXPT_NORMAL; NPROTS = 5; minEXPRESSION = 4; maxEXPRESSION = 20; break;
  case 4: EXPT = EXPT_COMPLEX; NPROTS = 5; minEXPRESSION = 4; maxEXPRESSION = 20; break;
  case 5: EXPT = EXPT_TEMPLATE; NPROTS = 1; minEXPRESSION = 4; maxEXPRESSION = 20; minMUT = 0.002; maxMUT = 0.0081; break;
  }
  
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

  M = (Compartment*)malloc(sizeof(Compartment)*NMITO);

  RunSimTest();
    
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

  sprintf(fstr, "sim-damage-update-%i-%i.csv", EXPTlabel, CHOOSE_RANDOM);
  fp = fopen(fstr, "w");
  fprintf(fp, "scanrate,mut,target,swap,social,nfuse,life,expression,import,rep,t,completes,completes2,avprot,avdna,avprotempty,avprotfull,maxdna,onedna,freesubs,avdamage\n");

  P.POISSON = 1;
  for(P.SCANRATE = 0.1; P.SCANRATE <= 1.1; P.SCANRATE += 0.45)
    {
      for(P.MUT = minMUT; P.MUT <= maxMUT; P.MUT *= 2)
	{
	  for(P.TARGET = 0; P.TARGET <= 5; P.TARGET++)
	    {
	      for(P.SOCIAL = 0; P.SOCIAL <= 1; P.SOCIAL++)
		{
		  for(P.SWAP = 0; P.SWAP <= 5; P.SWAP++)
		    {
		      P.NFUSE = NMITO/2;
		      // for(NFUSE = 0; NFUSE <= NMITO/2; NFUSE += NMITO/2)
		      {
			P.LIFE = 24*7;
			//	      for(LIFE = 1; LIFE < 100; LIFE *= 2)
			{
			  for(P.EXPRESSION = minEXPRESSION; P.EXPRESSION <= maxEXPRESSION; P.EXPRESSION += (maxEXPRESSION-minEXPRESSION) )
			    {
			      P.IMPORT = P.EXPRESSION;
			      //		      for(IMPORT = 100; IMPORT <= 1000; IMPORT *= 2)
			      {
				printf("%i,%i,%i,%i,%.3e,%.3e,%.3e\n", P.TARGET, P.SWAP, P.SOCIAL, P.NFUSE, P.LIFE, P.EXPRESSION, P.IMPORT);
				for(rep = 0; rep < NREP; rep++)
				  {
				    Simulate(P, fp, rep, M, &O); 
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

