/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.   
        usage: ms nsam howmany  -t  theta  [options]
		or
	       ms nsam howmany -s segsites  [options] 

	   nsam is the number of gametes per sample.
	   howmany is the number of samples to produce.
	   With -t the numbers of segregating sites will randomly vary 
		from one sample to the next.
	   with -s segsites,  the number of segregating sites will be
		segsites in each sample.

           Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be 
		specified.) The output will have the gametes in order such that
		the first n1 gametes are from the first island, the next n2 are
		from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate 
           rho: recombination rate between ends of segment times 4No
	   f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
	   track_len:  mean length of conversion track in units of sites.  The 
		       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No. 
           howmany: howmany samples to generate.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation. 
	A seed file called "seedms" will be created  if it doesn't exist. The
		seed(s) in this file will be modified by the program. 
		So subsequent runs
		will produce new output.  The initial contents of seedms will be
		printed on the second line of the output.
        Output consists of one line with the command line arguments and one
	 	line with the seed(s).
		The samples appear sequentially following that line.
		Each sample begins with "//", then the number of segregating sites, the positions
		of the segregating sites (on a scale of 0.0 - 1.0). On the following
		lines are the sampled gametes, with mutants alleles represented as
		ones and ancestral alleles as zeros.
	To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
		or:  cc -o ms ms.c streec.c rand2.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
		some other optimization switches might be usefully employed with some 
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
	Modifications to add gene conversion 5 Mar 2001
	Added demographic options -d  13 Mar 2001
	Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
	Changed cleftr() to check for zero rand() .13 June 2001
	Move seed stuff to subroutine seedit()  11 July 2001
	Modified streec.c to handle zero length demographic intervals 9 Aug 2001
	Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
	Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
	Changed many command line options, allowing arbitrary migration matrix, and subpopulation
	   sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
	   split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
	!!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
	Fixed bug which resulted in incorrect results for the case where
             rho = 0.0 and gene conversion rate > 0.0. This case was not handled
	    correctly in early versions of the program. 5 Apr 2004.  (Thanks to
	    Vincent Plagnol for pointing out this problem.) 
	Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
		 1 Jul 2004.
	Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
	Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
	Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
	Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
	    to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
	Added -seeds option 4 Nov 2006
	Added "tbs" arguments feature 4 Nov 2006
	Added -L option.  10 May 2007
	Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  May 19 2007.
	fixed bug with memory allocation when using -I option. This caused problems expecially on Windows
          machines.  Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one.
          Oct. 17, 2007.
     Modified pickb() and pickbmf() to eliminate rare occurrence of fixed mutations Thanks to J. Liechty and K. Thornton. 4 Feb 2010.
	Added -p n  switch to allow position output to be higher precision.  10 Nov 2012.
     Changed spot from int to long in the function re().  29 July 2013.  ( Thanks to Yuri D'Elia for this suggestion.)
      Changed function definitions, and misc other things to comply with c99
     requirements.  4 Mar 2014.
***************************************************************************/

#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "getnums.h"

extern "C" { 
    #include "ms.h"
}

using namespace Rcpp;

unsigned maxsites = SITESINC ;

struct segl {
	int beg;
	struct node *ptree;
	int next;
	};

double *posit ;
double *agevec ;
double alleleage ;
double segfac ;
int count, ntbs, nseeds ;


// [[Rcpp::export]]
Rcpp::StringVector ms_main(NumericVector nsam, NumericVector nreps, NumericVector t, NumericVector variable_list_rcpp, IntegerVector I_rcpp, NumericVector migration, NumericMatrix en, NumericMatrix ej)
{

    int i, k, howmany, segsites , afreq;
    char **list, **cmatrix(int gsam, int len), **tbsparamstrs ;
    FILE *pf; //*fopen() ;
    double probss, tmrca, ttot ;
    struct params pars;
    //deleted maxsites code here

    ntbs = 0 ;

    count=0;

    
    pars = getnums(&howmany, nsam, nreps, t, variable_list_rcpp, I_rcpp, migration, en, ej, pars) ;   /* results are stored in global variable, pars */

    if( !pars.commandlineseedflag ) seedit( "s");
    pf = stdout ;

    if( pars.mp.segsitesin ==  0 ) {
        list = cmatrix(pars.cp.nsam,maxsites+1);
        posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
        agevec = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
    }
    else {
        list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
        posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
        agevec = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
        if( pars.mp.theta > 0.0 ){
            segfac = 1.0 ;
            for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
        }
    }

    while( howmany-count++ ) {
        
        fprintf(pf,"\n//");
        segsites = gensam( list, &probss, &tmrca, &ttot, maxsites, pars) ;
        if( pars.mp.timeflag ) fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
            if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
                if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))
                    fprintf(pf,"prob: %g\n", probss ) ;
                fprintf(pf,"segsites: %d\n",segsites);
                if( segsites > 0 )    fprintf(pf,"positions: ");
                for( i=0; i<segsites; i++)
                    fprintf(pf,"%6.*lf ", pars.output_precision,posit[i] );
                fprintf(pf,"\n");
                if( (segsites > 0) && pars.mp.ageflag ){
                    fprintf(pf,"allele ages: ");
                    for( i=0; i<segsites; i++)
                        fprintf(pf,"%6.*lf ", pars.output_precision,agevec[i] );
                    fprintf(pf,"\n");
                    fprintf(pf,"allele freqs: ");
                    for( i=0; i<segsites; i++){
                        for(k=afreq=0; k<pars.cp.nsam; k++) afreq += ( (list[k][i] == '1') ? 1: 0 ) ;
                        fprintf(pf,"%d ", afreq );
                    }
                    fprintf(pf,"\n");
                }
                if( segsites > 0 )
                    for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }
            }
    }
    if( !pars.commandlineseedflag ) seedit( "end" );

}
