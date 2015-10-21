#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ExitError(char *msg)
{
  fprintf (stderr,"%s\nexit...\n",msg);
  exit (0);
}

struct s_bank {
  int nb_seq;               // nombre de sequences dans la banque
  int nb_res;               // nombre de residus
  int max_len;              // taille de la plus grande sequence
  char *data;               // image memoire de la banque
  char **seq;               // tableau de sequences
  char **com;               // tableau de commentaires
  int  *len;                // tableau de taille de sequence
};


struct s_bank * readBank(char *fname)
{
  FILE *fseq;                              // descripteur du fichier de la banque
  int i,j,k,l;
  struct s_bank *B = (struct s_bank *) malloc (sizeof(struct s_bank));

  if ((fseq=fopen(fname,"r"))==NULL) 
    { fprintf (stderr,"cannot open %s\n",fname); exit (0); }

  fseek (fseq,0,SEEK_END);
  k = ftell(fseq);
  if ((B->data = (char *) malloc(k*sizeof(char)))==NULL) 
    { fprintf (stderr,"Memory Overflow (readBank)\n"); exit (0); }
  rewind(fseq);
  fread (B->data,1,k,fseq);
  fclose(fseq);


  // 1ere passe de la banque
  // on compte le nombre de séquences

  B->nb_seq  = 0;
  i = 0;
  while (i<k) 
    {
      if (B->data[i]=='>')
	{
	  while (B->data[i]!='\n') i++;
	  B->nb_seq++;
	}
      else
	{
	  while ((i<k)&&(B->data[i]!='>')) i++;
	}
    }

  if ((B->seq = (char **)  malloc ((B->nb_seq)*sizeof(char *)))==NULL) 
    { fprintf (stderr,"Memory Overflow (readBank)\n"); exit (0); }
  if ((B->com = (char **)  malloc ((B->nb_seq)*sizeof(char *)))==NULL) 
    { fprintf (stderr,"Memory Overflow (readBank)\n"); exit (0); }
  if ((B->len = (int *)  malloc ((B->nb_seq)*sizeof(int)))==NULL) 
    { fprintf (stderr,"Memory Overflow (readBank)\n"); exit (0); }


  // 2eme passe sur la banque
  // on mémorise les séquences et les commentaires

  B->max_len = 0;
  B->nb_seq=0;
  B->nb_res=0;
  i=0;
  while (i<k) 
    {
      if (B->data[i]=='>')
	{
	  B->com[B->nb_seq] = &B->data[i];
	  while (B->data[i]!='\n') i++;
	  B->data[i]='\0';
	  i++;
	  B->nb_seq++;
	}
      else
	{
	  B->seq[B->nb_seq-1] = &B->data[i];
	  j=i; l=0;
	  while ((i<k)&&(B->data[i]!='>'))
	    {
	      if (B->data[i]>='A') { B->data[j]=B->data[i]; j++; l++; B->nb_res++; }
	      i++;
	    }
	  B->len[B->nb_seq-1] = l;
	  if (l > B->max_len) B->max_len = l;
	}
    }
  return B;
}


void writeBank(struct s_bank *B, char *fname)
{
 FILE *fseq;                              // descripteur du fichier de la banque

  int i,k;

  if ((fseq=fopen(fname,"w"))==NULL) 
    { fprintf (stderr,"cannot open %s\n",fname); exit (0); }


  for (i=0; i<B->nb_seq; i++)
    {
      fprintf (fseq,"%s",B->com[i]);
      for (k=0; k<B->len[i]; k++)
	{
	  fprintf (fseq,"%c",B->seq[i][k]);
	  if (k%80 == 0) fprintf (fseq,"\n");
	}
      if ((k-1)%80 != 0) fprintf (fseq,"\n");
    }

  fclose (fseq);
}
int CODE_AA [256] =
  { 
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,           /*   0 -  15 */
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,           /*  16 -  31 */
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,           /*  32 -  47 */
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,           /*  48 -  63 */
    23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,           /*  64 -  79 */
    14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,           /*  80 -  95 */
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,           /*  96 - 111 */
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,           /* 112 - 127 */
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23
  };


