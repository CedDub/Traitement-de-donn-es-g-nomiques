#define THRESHOLD 5
#define sizeSEED 10
#define sizeREAD 96
#define sizeIDX2 (1<<20)
#define sizeIDX1 (1<<20)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iobank.h"

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
  int i;
  for (i=0; i<sizeREAD; i++)
    {
      if (reads[i]!=genome[i+k]) d++;
    }
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
  
