#define THRESHOLD 5
#define sizeSEED 10
#define sizeREAD 96
#define sizeIDX2 (1<<20)
#define sizeIDX1 (1<<20)

#include "nmmintrin.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iobank.h"
#include <time.h>
int key (char *s, int j)
{
  int k=0;
  int i, c;
  for (i=0; i<sizeSEED; i++)
    {
      c = s[i+j];
      k = k<<2;
      k = k+((c>>1)&3);
    }
  return k;
}

int __attribute__ ((noinline)) distance (char *reads, char *genome, int k)
{
  int d=0;
  unsigned long int r0=0, r1=0, r2=0;
  unsigned long int g0=0, g1=0, g2=0;
  int i;
  for (i=0; i<32; i++)
    {
      r0 = r0<<2;
      r0 = r0|((reads[i]>>1)&3);
      r1 = r1<<2;
      r1 = r1|((reads[i+32]>>1)&3);
      r2 = r2<<2;
      r2 = r2|((reads[i+64]>>1)&3);

      g0 = g0<<2;
      g0 = g0|((genome[k+i]>>1)&3);
      g1 = g1<<2;
      g1 = g1|((genome[k+i+32]>>1)&3);
      g2 = g2<<2;
      g2 = g2|((genome[k+i+64]>>1)&3);
    }
/*	printf("\nreads =");
        for(i=0; i<96; i++){
                printf("%c", reads[i]);
        }

	printf("\ngenome=");
	for(i=0; i<96; i++){
		printf("%c", genome[i+k]);
	}
        printf("\nr0 =%lx\n",r0);
        printf("g0 =%lx\n",g0);
        printf("r1 =%lx\n",r1);
        printf("g1 =%lx\n",g1);
        printf("r2 =%lx\n",r2);
        printf("g2 =%lx\n",g2);
*/
        unsigned long int s00 =r0^g0;
  //      printf("s00=%lx\n", s00);
        unsigned long int s10 = 0xAAAAAAAAAAAAAAAA & s00;               // prise en compte des parties gauches
    //    printf("s10=%lx\n", s10);
        unsigned long int s20 = 0x5555555555555555 & s00;               // prise en compte des parties droites
      //  printf("s20=%lx\n", s20);
        s20 = s20<<1;
       // printf("s20 decal=%lx\n", s20);
        s10 = s10|s20;
       // printf("s10 ou =%lx\n", s10);

        unsigned long int s01 =r1^g1;
       // printf("s01=%lx\n", s01);
        unsigned long int s11 = 0xAAAAAAAAAAAAAAAA & s01;
       // printf("s11=%lx\n", s11);
        unsigned long int s21 = 0x5555555555555555 & s01;
       // printf("s21=%lx\n", s21);
        s21 = s21<<1;
       // printf("s21 decal=%lx\n", s21);
        s11 = s11|s21;
       // printf("s11 ou =%lx\n", s11);

        unsigned long int s02 = r2^g2;
       // printf("s02=%lx\n", s02);
        unsigned long int s12 = 0xAAAAAAAAAAAAAAAA & s02;
       // printf("s12=%lx\n", s12);
        unsigned long int s22 = 0x5555555555555555 & s02;
       // printf("s22=%lx\n", s22);
        s22 = s22<<1;
       // printf("s22 decal=%lx\n", s22);
        s12 = s12|s22;
       // printf("s12 ou =%lx\n", s12);

        d = _mm_popcnt_u64(s10)+_mm_popcnt_u64(s11)+_mm_popcnt_u64(s12);
        //printf("score0=%d\n",_mm_popcnt_u32(s10));
       // printf("score1=%d\n",_mm_popcnt_u32(s11));
       // printf("score2=%d\n",_mm_popcnt_u32(s12));
       // printf("scoret=%d\n",d);                  
        return d;

}

void __attribute__ ((noinline)) process_read(struct s_bank *query, int nq, char *genome, int *IDX1, int*IDX2, int *IDX)
{
  int r=0;
  int iq, x, d, j, k;
  int RES[10];

  for (iq = 0; iq < query->len[nq]-sizeSEED+1; iq++) // for each words of sizeSEED characters
    {
      k = key(query->seq[nq],iq); // get the key
      for (j = IDX2[k]; j<IDX2[k]+IDX1[k]; j++) // get all position in the genome where such key is
	{
	  d = distance(query->seq[nq],genome,IDX[j]-iq); // compute a distance
 if (d <= THRESHOLD)
	    {
	      for (x=0; x<r; x++) if (RES[x] == IDX[j]-iq) break; // test if the position has already been marked as a good one
	      if (x == r) // no
		{
		  RES[r] = IDX[j]-iq; // mark this position as a visited position
		  r++;
		 printf ("%d %d\n",IDX[j]-iq,d); // print results if distance <= threshold
		}
	    }
	}
    }
}

int main (int argc, char *argv[])
{
  struct s_bank *ref;    // reference genome
  struct s_bank *query;  // reads
  int *IDX, *IDX1, *IDX2;
  int nq, ir, i, k;
  if (argc <= 2) { fprintf (stderr,"Usage : %s query bank\n",argv[0]); exit (0); }

  // read the query bank and the ref bank
  query   = readBank(argv[1]);
  fprintf (stderr,"query:   %d seq %d nt  (max=%d)\n",query->nb_seq,query->nb_res,query->max_len);
  ref = readBank(argv[2]);
  fprintf (stderr,"bank:    %d seq %d nt  (max=%d)\n",ref->nb_seq,ref->nb_res,ref->max_len);

  // genome indexing
  IDX  = (int *)malloc(ref->nb_res*sizeof(int));
  IDX1 = (int *)malloc(sizeIDX1*sizeof(int));
  IDX2 = (int *)malloc(sizeIDX2*sizeof(int));
  // initialization of IDX1
  for (ir=0; ir<ref->len[0]-sizeSEED+1; ir++)
    {
      k = key(ref->seq[0],ir);
      IDX1[k]++;
    }
  // initialization of IDX2
  IDX2[0] = 0;
  for (i=1; i<sizeIDX2; i++) IDX2[i] = IDX2[i-1]+IDX1[i-1];
  // initialization of IDX;
  for (ir=0; ir<ref->len[0]-sizeSEED+1; ir++)
    {
      k = key(ref->seq[0],ir);
      IDX[IDX2[k]] = ir;
      IDX2[k]++;
    }
  // initialization of IDX2
  IDX2[0] = 0;
  for (i=1; i<sizeIDX2; i++) IDX2[i] = IDX2[i-1]+IDX1[i-1];

  // process all the reads sequentially
  for (nq=0; nq<query->nb_seq; nq++) 
    {
      printf ("%s\n",query->com[nq]);
      process_read(query,nq,ref->seq[0],IDX1,IDX2,IDX);
    }
  free(IDX);
  free(IDX1);
  free(IDX2);
  return 1;
}
  
