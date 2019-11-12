/* This program calculates residence time autocorrelation functions to evaluate exchange kinetics
   in and out of a defined region. */

/* Inspired from the methodology described in:
S. A. Kislenko - R. H. Amirov - I. S. Samoylov 
Phys. Chem. Chem. Phys., 12, 11245 (2010) */

/* 
Intermittent function = <hI(t).hI(0)>/<hI(0)^2>
Continuous function = <hC(t).hC(0)>/<hC(0)^2>
*/

/* To execute the program, you just have to launch it with an input file giving
the following information (program < inputfile)
first line: name of the xyz file with the positions (> name)
second line: number of configurations in this file (> nconfigs)
third line: time between two configurations (> dtime)
fourth line: number of species for which the calculation will be done (> nspecies) !MAX 1000!
next lines: atom type for each species (> atomtype)
next line: limits in the z direction (> zmin,zmax)
next line: number of correlations steps (> ncorr)
*/

/* !ONLY WORKS FOR A CONSTANT NUMBER OF ATOMS WRITTEN IN THE SAME ORDER FOR ALL CONFIGURATIONS! */
/* Note that no conversion is done so the time unit in the output is the same as the input. */ 


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define L 1000


char name[50];
/* Possibility to have up to 1000 atom types */
char atomtype[1001][5],type[5];
int **normI, **normC, nconfigs, nspecies, nionstot, ncorr;
double zmin, zmax;
double **storeC, **cfC, dstoreC, rvacfC, dtime, prodstoreC;
double **storeI, **cfI, dstoreI, rvacfI;


void read();
void calculate();
void write();

/* For dynamic allocation */
int **imatrix(int nl, int nc);
double **dmatrix(int nl, int nc);


int main ( void )
{
 printf("\nHello!\n\n");
 read();
 calculate();
 write();
 printf("\nDone!\n\n");

 return 0;
}


void read()
{
 int i,j,count;
 FILE *in;
 char ligne[L];

 /* Reads the information in the input file */
 scanf("%s",name);		/* name of the file with the positions */
 scanf("%d",&nconfigs);		/* number of configurations in the file */
 scanf("%lf",&dtime);		/* time between two configurations in ps */
 scanf("%d",&nspecies);		/* number of species for which the calculation will be done */
 for(i=1;i<=nspecies;i++)
    {
     scanf("%s",atomtype[i]);   /* atom type for each species */
    }
 scanf("%lf %lf",&zmin,&zmax);	/* limits in the z direction */
 scanf("%d",&ncorr);		/* number of correlations steps */

 /* Read the total number of atoms in the xyz file */
 in=fopen(name,"r");
 fgets(ligne,L,in);
 sscanf(ligne,"%d",&nionstot);
 fclose(in);
  
 printf("Total number of atoms: %d.\n",nionstot);

 /* Allocates the matrices */
 storeI=dmatrix(nionstot+1,ncorr+1); 
 storeC=dmatrix(nionstot+1,ncorr+1); 
 cfI=dmatrix(ncorr+1,nspecies+1); normI=imatrix(ncorr+1,nspecies+1);
 cfC=dmatrix(ncorr+1,nspecies+1); normC=imatrix(ncorr+1,nspecies+1);
}


void calculate()
{
 int i,ii,j,k,l,n,itype;
 double z;
 FILE *in;
 char ligne[L];

 in=fopen(name,"r"); 

 /* Calculation of the autocorrelation function */
 for(i=1;i<=nconfigs;i++)
    {
     /* Skip the line with the number of atoms */
     fgets(ligne,L,in);
     /* Skip the comment line */
     fgets(ligne,L,in);
     if((i%100)==0) printf("%d\n",i); 
     if(i<=ncorr) 
       {
        for(j=1;j<=nionstot;j++) 
	   {	
	    /* Get the atom type and z position for each atom */
	    fgets(ligne,L,in);
	    sscanf(ligne,"%s %*lf %*lf %lf",type,&z);
	    /* Check if the atom is considered in this calculation */
	    itype=0;
	    for(l=1;l<=nspecies;l++)
	       {
	        if(strcmp(atomtype[l],type)==0) itype=l;
	       }
	    if(itype!=0)
	      {
	       /* Main calculation: heaviside and autocorrelation */
	       if((z<=zmax)&&(z>=zmin)) 
	         {
		  storeI[j][i]=1; storeC[j][i]=1;
		 }
	       if((z>zmax)||(z<zmin)) 
	         {
	          storeI[j][i]=0; storeC[j][i]=0;
	         }
	       dstoreI=storeI[j][i]; 
	       dstoreC=storeC[j][i];
	       for(k=i;k>=1;k--) 
	          {
	           rvacfI=dstoreI*storeI[j][k]; cfI[i-k][itype]+=rvacfI; 
                   normI[i-k][itype]++;
		   prodstoreC=storeC[j][k];
                   n=k+1;
		   while((n<i)&&(prodstoreC==1))
		        {
		         prodstoreC*=storeC[j][n];
		         n++;
		        }
		   rvacfC=dstoreC*prodstoreC; cfC[i-k][itype]+=rvacfC; 
                   normC[i-k][itype]++;
		  }
	      }
      	   }	
	 }
    if(i>ncorr) 
      {
       ii=i%ncorr;
       for(j=1;j<=nionstot;j++) 
	  {
	   /* Get the atom type and z position for each atom */
	   fgets(ligne,L,in);
	   sscanf(ligne,"%s %*lf %*lf %lf",type,&z);
	   /* Check if the atom is considered in this calculation */
	   itype=0;
	   for(l=1;l<=nspecies;l++)
	      {
	       if(strcmp(atomtype[l],type)==0) itype=l;
	      }
	   if(itype!=0)
	     {
	      /* Main calculation: heaviside and autocorrelation */
	      if((z<=zmax)&&(z>=zmin)) 
	        {
	         storeI[j][ii]=1; storeC[j][ii]=1;
	        }
	      if((z>zmax)||(z<zmin)) 
	        {
	         storeI[j][ii]=0; storeC[j][ii]=0;
	        }
	      dstoreI=storeI[j][ii];
              dstoreC=storeC[j][ii];
	      for(k=ncorr;k>=ii+1;k--) 
	         {
	          rvacfI=dstoreI*storeI[j][k]; cfI[ii-k+ncorr][itype]+=rvacfI; 
                  normI[ii-k+ncorr][itype]++;
	          prodstoreC=storeC[j][k];
	          n=k+1;
	          while((n<=ncorr)&&(prodstoreC==1))
	    	       {
			prodstoreC*=storeC[j][n];
			n++;
		       }
		  n=1;
		  while((n<ii)&&(prodstoreC==1))
		       {
		 	prodstoreC*=storeC[j][n];
			n++;
		       }
		  rvacfC=dstoreC*prodstoreC; cfC[ii-k+ncorr][itype]+=rvacfC; 
                  normC[ii-k+ncorr][itype]++;
		 }
	      for(k=ii;k>=1;k--) 
	         {
	          rvacfI=dstoreI*storeI[j][k]; cfI[ii-k][itype]+=rvacfI; 
		  normI[ii-k][itype]++;
		  prodstoreC=storeC[j][k];
                  for(n=k+1;n<ii;n++) 
		     {
		      if(n>ii) printf("Problem\n"); 
		      prodstoreC*=storeC[j][n];
		     }
		  rvacfC=dstoreC*prodstoreC; cfC[ii-k][itype]+=rvacfC; 
 		  normC[ii-k][itype]++;
		 }		    
	     }
	   } /* end loop nionstot */
      } /* end if i/ncorr */
    } /* end loop on configs */

 /* Normalisation of the autocorrelation function */
 for(i=0;i<ncorr;i++)
  {
    for(j=1;j<=nspecies;j++)
       {
        cfI[i][j]/=(normI[i][j]*1.0);
        cfC[i][j]/=(normC[i][j]*1.0);
       }
  }

 fclose(in);


}

void write()
{
 int i,j;
 FILE *out;

 for(j=1;j<=nspecies;j++)
    {
     sprintf(name,"%s%d.dat","resid-time-acf",j);
     out=fopen(name,"w");
     fprintf(out,"# Time - Intermittent function - Continuous function\n");
     for(i=0;i<ncorr;i++) 
	{
	 fprintf(out,"%lf	%lf	%lf\n",i*dtime,cfI[i][j]/cfI[0][j],cfC[i][j]/cfC[0][j]);
	}
     fclose(out); 
    }

}


/********************************************************************************************/
/*Allocation dynamique*/

int **imatrix(int nl, int nc)
{
  int i;
  int **m;
  m=(int **) malloc(nl*sizeof(int*));
  if (m) { m[0]=(int *) malloc(nl*nc*sizeof(int));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

double **dmatrix(int nl, int nc)
{
  int i;
  double **m;
  m=(double **) malloc(nl*sizeof(double*));
  if (m) { m[0]=(double *) malloc(nl*nc*sizeof(double));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}




