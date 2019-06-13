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

struct params getnums(int *phowmany, NumericVector nsam, NumericVector nreps, NumericVector t,
        NumericVector variable_list_rcpp, IntegerVector I_rcpp, NumericVector migration,
        NumericMatrix en, NumericMatrix ej, struct params pars){
    int i, j, m, n, r, sum , npop , npop2;
    double migr, mij, psize, palpha;
    FILE *pf;
    
    //void caseI(), caseen(), caseej(), casema();
    int commandlineseed( char ** ) ;
    void free_eventlist( struct devent *pt, int npop );
    struct devent *ptemp , *pt ;
    char ch3;
    NumericVector rcpp_row_en, rcpp_row_ej;
    
    int variable_list[variable_list_rcpp.length()];
    for(i = 0; i < variable_list_rcpp.length(); i++){
        variable_list[i] = variable_list_rcpp[i];
    }
    
    pars.cp.nsamin = nsam[0]; //*
    pars.cp.nsam = pars.cp.nsamin;
    *phowmany = nreps[0]; //*
    
    pars.commandlineseedflag =  variable_list[0]; //*
    pars.output_precision =  variable_list[1];
    pars.cp.r = pars.mp.theta =  pars.cp.f = variable_list[2];
    pars.cp.track_len = variable_list[3] ;
    pars.cp.npop = npop = variable_list[4] ;
    pars.mp.segsitesin = variable_list[5] ;
    pars.mp.treeflag = variable_list[6] ;
    pars.mp.timeflag =  variable_list[7];
    pars.mp.ageflag =  variable_list[8];
    pars.mp.mfreq = variable_list[9];
    pars.cp.nsites = variable_list[10] ; //*
    
    pars.cp.deventlist = NULL;
    
    pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
    pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
    
    pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
    (pars.cp.config)[0] = pars.cp.nsamin ;
    pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
    (pars.cp.size)[0] = 1.0  ;
    pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
    (pars.cp.alphag)[0] = 0.0  ;
    /*End of intitialization*/
    
    /* Case t */
    pars.mp.theta = t[0];
    
    /* Case T - Default = 0 */
    pars.mp.treeflag = 1 ;
    
    /* Case I */
    migr = 0.0;
    int I[I_rcpp.length()];
    for(m = 0; m < I_rcpp.length(); m++){
        I[i] = I_rcpp[i];
    }
    pars = caseI(I, migr, pars);
    
    /* Case ma - did not code other m cases without the a */

    
    double migmat_array[migration.length()];
    for(m = 0; m < migration.length(); m++){
        migmat_array[m] = migration[m];
    }
    
    pars = casema(migmat_array, pars);
    
    double row[3];
    for(i = 0; i <= en.nrow(); i++){
        Rcpp::NumericVector rcpp_row_en = en(i,_);
        for(r = 0; r < 3; r++){ //hard coded as three because the en and ej should only have 3 numbers
            row[r] = rcpp_row_en[r];
        }
        pars = caseen(row, pars);
    }
    
    for(i = 0; i <= ej.nrow(); i++){
        Rcpp::NumericVector rcpp_row_ej = ej(i,_);
        for(r = 0; r < 3; r++){
            row[r] = rcpp_row_ej[r];
        }
        pars = caseej(row, pars);
    }
    return pars;
}

struct params caseen(double *en, struct params pars){
    void addtoelist( struct devent *pt, struct devent *elist );
    struct devent *ptemp , *pt ;
    
    pt = (struct devent *)malloc( sizeof( struct devent) ) ;
    pt->detype = 'n';
    
    pt->time = en[0];
    pt->nextde = NULL ;
    if( pars.cp.deventlist == NULL )
        pars.cp.deventlist = pt ;
    else if ( pt->time < pars.cp.deventlist->time ) {
        ptemp = pars.cp.deventlist ;
        pars.cp.deventlist = pt ;
        pt->nextde = ptemp ;
    }
    else {
        addtoelist( pt, pars.cp.deventlist ) ;
        pt->popi =  en[1] -1 ; /* keep this as - 1? */
        pt->paramv = en[2] ; }
    
    return pars;
}

struct params caseej(double *ej, struct params pars){
    void addtoelist( struct devent *pt, struct devent *elist );
    struct devent *ptemp , *pt ;
    
    pt = (struct devent *)malloc( sizeof( struct devent) ) ;
    pt->detype = 'j';
    
    pt->time = ej[0];
    pt->nextde = NULL ;
    if( pars.cp.deventlist == NULL )
        pars.cp.deventlist = pt ;
    else if ( pt->time < pars.cp.deventlist->time ) {
        ptemp = pars.cp.deventlist ;
        pars.cp.deventlist = pt ;
        pt->nextde = ptemp ;
    }
    else{
        addtoelist( pt, pars.cp.deventlist ) ;
        pt ->popi = ej[1] - 1;
        pt ->popj = ej[2] - 1; }
    
    return pars;
}

struct params caseI(int *I, double migr, struct params pars){
    int i, j, npop;
    pars.cp.npop = I[0];
    pars.cp.config = (int *) realloc( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
    npop = pars.cp.npop ;
    
    for( i=1; i<= pars.cp.npop; i++) {
        pars.cp.config[i-1] = I[i];
    }
    pars.cp.mig_mat = (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
    pars.cp.mig_mat[0] = (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
    for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] = (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
    pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
    pars.cp.alphag = (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
    for( i=1; i< pars.cp.npop ; i++) {
        (pars.cp.size)[i] = (pars.cp.size)[0]  ;
        (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
    }
    for( i=0; i<pars.cp.npop; i++)
        for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
    for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
    
    return pars;
}

struct params casema(double *migmat_array, struct params pars){
    int pop, pop2, tempArg, npop, npop2;
    
    tempArg = -1;
    for( pop = 0; pop <npop; pop++)
        for( pop2 = 0; pop2 <npop; pop2++){
            pars.cp.mig_mat[pop][pop2]= migmat_array[tempArg++];}
    for( pop = 0; pop < npop; pop++) {
        pars.cp.mig_mat[pop][pop] = 0.0 ;
        for( pop2 = 0; pop2 < npop; pop2++){
            if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
        }
    }
    return pars;
}

