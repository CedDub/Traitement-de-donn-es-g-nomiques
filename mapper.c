#define THRESHOLD 5
#define sizeSEED 10
#define sizeREAD 96
#define sizeIDX2 (1<<20)
#define sizeIDX1 (1<<20)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iobank.h"
#include "immintrin.h"
#include "nmmintrin.h"

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

/*int __attribute__ ((noinline)) distance (char *reads, char *genome, int k)
{
  int d=0;
  int i;
  for (i=0; i<sizeREAD; i++)
    {
      if (reads[i]!=genome[i+k]) d++;
    }
  return d;
}*/

int __attribute__ ((noinline)) distance (char *reads, char *genome, int k)
{
  int d=0;
  long int r0, r1, r2;
  long int g0, g1, g2;
  int i;
  for (i=0; i<64; i++)
    {
      r0 = r0<<2;
      r0 = r0+((reads[i]>>1)&3);
      r1 = r1<<2;
      r1 = r1+((reads[i+64]>>1)&3);
      r2 = r2<<2;
      r2 = r2+((reads[i+128]>>1)&3);

      g0 = g0<<2;
      g0 = g0+((genome[i]>>1)&3);
      g1 = g1<<2;
      g1 = g1+((genome[i+64]>>1)&3);
      g2 = g2<<2;
      g2 = g2+((genome[i+128]>>1)&3);
    }
  
  // XOR
/*
  long int s00 = _mm_xor_si64(r0, g0);
  long int s01 = _mm_xor_si64(r1, g1);
  long int s02 = _mm_xor_si64(r2, g2);

  // Séléction des chaines de gauche et de droite par des masques
  long int s10 = _mm_and_si64(2863311530, s00);
  long int s20 = _mm_and_si64(1431655765, s00);
  long int s11 = _mm_and_si64(2863311530, s01);
  long int s21 = _mm_and_si64(1431655765, s01);
  long int s12 = _mm_and_si64(2863311530, s02);
  long int s22 = _mm_and_si64(1431655765, s02);

  // Décalage à gauche
  s10 = _m_pslldi(s10, 1);
  s11 = _m_pslldi(s10, 1);
  s12 = _m_pslldi(s10, 1);

  // OR
  s10 = _mm_or_si64(s10, s20);
  s11 = _mm_or_si64(s11, s21);
  s12 = _mm_or_si64(s12, s22);
*/

	long int s00 =r0^g0;
	long int s10 = 2863311530 & s00;
	long int s20 = 1431655765 & s00;
	s10 = s10>>1;
	s10 = s10|s20;

	long int s01 =r1^g1;
	long int s11 = 2863311530 & s01;
	long int s21 = 1431655765 & s01;
	s11 = s11>>1;
	s11 = s11|s21;

	long int s02 =r2^g2;
	long int s12 = 2863311530 & s02;
	long int s22 = 1431655765 & s02;
	s12 = s12>>1;
	s12 = s12|s22;

	d = _popcnt64(s10) + _popcnt64(s11) + _popcnt64(s12);

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
  
