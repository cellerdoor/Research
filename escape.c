//
//  ratio.c
//  
//
//  Created by Lavisha Jindal on 2016-01-06.
//
//

#include "PrA_plsmd_mc.h"
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)



extern double ran2(long);
int main(int argc, char *argv[])
{
    FILE *output1;
    FILE *output2;
    FILE *output3;

    //INTRODUCTION OF VARIABLES
    int grid_max; //number of ParA molecules
    int bndry;
    int pl_max; //number of molecules in the plasmid chain
    int num; //that pl_mol that is chosen to shift
    int tmax; //maximum steps for the MC simulation
    int i,j,k,l,z; //counters
    int the;
    int ent_pro; //calculate the number of ParA molecules 'burned'
    double xdist,ydist;
    double randno;
    double wall_dist;
    double dist; // stores distance between parA and plasmid
    double pot_h;
    int *ratio; //this is the amount of ParB on plasmid mols
    int *removed; //stores information whether the parA has been removed
    double *Apos;  // stores position of the one-d parA molecules
    double *xpos;  //stores the position of the plasmid chain molecules
    double *ypos;  //stores the position of the plasmid chain molecules
    double rand_no; //a random number for probability calculation
    double rad_pl; //radius of the plasmid molecules
    double dx,dy,dA; //changes in the plasmid distance , ParA separation
    double cod_pl; //stores the center of mass of the plasmid
    double E_change; //stores the total change in energy due to the move
    double E_ch_parA; // change in energy due to the parA
    double r_0;
    double v_teether;
    double d_c;
    double frac; //ratio of beads which are ParB in the plasmid
    double d_old, d_new;
    double cod_x,cod_y; //the x and y of the centre of mass for the plasmid
    int grid_start,grid_end; // the numbers close to the plasmid that have an effect to begin with
    char filename1[256];	/*used for filename*/
    char filename2[256];
    double low, high;
    double p_r;
    
    
    //rough variables:
    double rads;
    int plasmid_nos;
    double f_good;
    int n_good,ntot;
    
    //INITIALIZATIONS
    grid_max = 150;
    bndry = 15;
    pl_max = atol(argv[1]);   //no of plasmid beads
    plasmid_nos = pl_max + 2;
    frac = 100*1.0/100.;          //fraction of parB beads in plasmid
    tmax = 3000000;   //number of MC runs, reduce to 0 if testing potentials
    dA = 1.0;        //diameter of the parA beads
    dx = .1;       //distance moved in one step
    r_0 = 1.5;  // maximum extention for the plasmid beads
    rad_pl = 1.0;  //diameter of the plasmid beads
    pot_h = -1.25;   //height of the potential felt bet parA and plasmid
    d_c = 1.5;    // that distance at which parA and plasmid interact has to be more than sum of radius
    wall_dist = 10.0;  //the constraining mag force distance, cut
    v_teether = 10.0;  //the teethering potential bw two plasmid neighbours
    srand(time(NULL));
    p_r = atol(argv[2])*1.0/1000000;
    
    //Filenames decided here:
    sprintf(filename1, "cod_h%2f_fb%2ld_pr_%0.3f.dat",wall_dist,atol(argv[1]),p_r);
    
    long seed;
    seed = (long)time(NULL);
    ran2(-seed);
    
    Apos = (double *) malloc(grid_max * sizeof(double));
    xpos = (double *) malloc(plasmid_nos * sizeof(double));
    ypos = (double *) malloc(plasmid_nos * sizeof(double));
    ratio = (int *) malloc(pl_max * sizeof(int));
    removed = (int *) malloc(grid_max * sizeof(int));
    
    for (i=1; i<=grid_max; i++) {
        Apos[i] = i*dA ;
        removed[i] = 1;
    }
    rads = pl_max*rad_pl/2./M_PI;
    for (i=1; i<=pl_max/2; i++) {
        xpos[i] = grid_max/2.-pl_max*1.0/4.*rad_pl+rad_pl*i;
        xpos[pl_max+1-i] = grid_max/2.-pl_max*1.0/4.*rad_pl+rad_pl*i;
        ypos[i] = wall_dist - 1.0 ;
        ypos[pl_max+1-i] = wall_dist - 2.0;
    }
    for (i=1; i<=pl_max; i++) {
    	  ratio[i] = i  <= frac*pl_max ? 1 : 0;
    }
    xpos[pl_max+1] = xpos[1];
    xpos[0] = xpos[pl_max];
    ypos[pl_max+1] = ypos[1];
    ypos[0] = ypos[pl_max];


    //Equilibration period*/
    for (k=1;k<=50000;k++){
        for (z=1;z<=1;z++){
            for (j=1; j <=pl_max; j++){
                E_change = 0.0;
                ntot++;
                
                xpos[pl_max+1] = xpos[1];
                xpos[0] = xpos[pl_max];
    			ypos[pl_max+1] = ypos[1];
    			ypos[0] = ypos[pl_max];
    				 
                //MONTE CARLO INTERACTIONS
                //pick a random particle and move it
                num = (int)(ran2(seed)/(1./pl_max)) + 1;
                xdist = xpos[num] + dx*2.*(.5-ran2(seed));
                ydist = ypos[num] + dx*2.*(.5-ran2(seed));
                if (ydist >= wall_dist || ydist <= 0.0){
                    continue;
                }
                
                //start calculating energy changes- neighbour
                d_old = pow(xpos[num]-xpos[num+1],2.) + pow(ypos[num]-ypos[num+1],2.0);
                d_new = pow(xdist-xpos[num+1],2.0) + pow(ydist-ypos[num+1],2.);
                E_change = E_p(d_new,r_0,v_teether) - E_p(d_old,r_0,v_teether);
                d_old = pow(xpos[num]-xpos[num-1],2.) + pow(ypos[num]-ypos[num-1],2.0);
                d_new = pow(xdist-xpos[num-1],2.0) + pow(ydist-ypos[num-1],2.);
                E_change = E_change + E_p(d_new,r_0,v_teether) - E_p(d_old,r_0,v_teether);
                
                //energy changes- other plasmid beads
                for (i=num+1;i<=pl_max;i++){
                    dist = pow(xdist-xpos[i],2.)+pow(ydist-ypos[i],2.);
                    if (dist >= 1.5*rad_pl){
                        continue;
                    }
                    E_change = E_change + E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25);
                    dist = pow(xpos[num]-xpos[i],2.)+pow(ypos[num]- ypos[i],2.);
                    E_change = E_change - E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25);
                }
                for (i=1;i<num;i++){
                    dist = pow(xdist-xpos[i],2.)+pow(ydist- ypos[i],2.);
                    if (dist >= 1.5*rad_pl){
                        continue;
                    }
                    E_change = E_change + E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25) ;
                    dist = pow(xpos[num]-xpos[i],2.)+pow(ypos[num]-ypos[i],2.);
                    E_change = E_change - E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25) ;
                }
                
                //energy calculation - parA layer - if farther than 2.0*rad_pl then no calculation of anything
                if (ydist >= 2.5*rad_pl){
                    randno = ran2(seed);
                    if (E_change <= 0.0 || randno <= exp(-E_change)) {
                        xpos[num] = xdist;
                        ypos[num] = ydist;
                        n_good++;
                    }
                    continue;
                }
                
                //plasmid-ParA hardwall calculation
                // caluculate the lowest grid number it needs to evaluate within
                low = xpos[1];
                high = xpos[1];
                for (i=2 ;i <= pl_max;i++){
                	  if (xpos[i] < low){
                	  		low = xpos[i];	
                	  }
                	  if (xpos[i] > high){
                	  		high = xpos[i];
                	  }
                }
                grid_start = (int)low - 5;
                grid_end = (int)high + 5; 
                
                for (i=grid_start;i<=grid_end;i++){
                    dist = pow(xdist-Apos[i],2.)+pow(ydist,2.);
                    E_change = E_change + E_lj(dist,dA,1.122*dA,1.0,-0.25);
                    dist = pow(xpos[num]-Apos[i],2.)+ pow(ypos[num],2.);
                    E_change = E_change - E_lj(dist,dA,1.122*dA,1.0,-0.25);
                }
                
                //plasmid-ParA: if bead is not ParB then no attraction
                if (ratio[num] == 0){
                    randno = ran2(seed);
                    if (E_change <= 0.0 || randno <= exp(-E_change)) {
                        xpos[num] = xdist;
                        ypos[num] = ydist;
                        n_good++;
                    }
                    continue;
                }
                
                //energy calculation - parA layer- attractive- only if that particular grid number is still a parA
                for (i=grid_start;i<=grid_end;i++){
                    dist = pow(xdist-Apos[i],2.)+pow(ydist,2.);
                    E_change = E_change + removed[i]*E_sq(dist,d_c,pot_h);
                    dist = pow(xpos[num]-Apos[i],2.)+pow(ypos[num],2.);
                    E_change = E_change - removed[i]*E_sq(dist,d_c,pot_h);
                }
                
                randno = ran2(seed);
                if (E_change <= 0.0 || randno <= exp(-E_change)) {
                    xpos[num] = xdist;
                    ypos[num] = ydist;
                    n_good++;
                }
            }
        }
        
    }
    ent_pro = 0;
    output1 = fopen(filename1,"w");
    for (k=1;k<=tmax;k++){
        for (z=1;z<=1000;z++){
            for (j=1; j <=pl_max; j++){
            
                E_change = 0.0;
                ntot++;
                
                xpos[pl_max+1] = xpos[1];
                xpos[0] = xpos[pl_max];
                ypos[pl_max+1] = ypos[1];
                ypos[0] = ypos[pl_max];
                    
                //MONTE CARLO INTERACTIONS
                //pick a random particle and move it
                num = (int)(ran2(seed)/(1./pl_max)) + 1;
                xdist = xpos[num] + dx*2.*(.5-ran2(seed));
                ydist = ypos[num] + dx*2.*(.5-ran2(seed));
                if (ydist >= wall_dist || ydist <= 0.0){
                    continue;
                }

                //start calculating energy changes- neighbour
                d_old = pow(xpos[num]-xpos[num+1],2.) + pow(ypos[num]-ypos[num+1],2.0);
                d_new = pow(xdist-xpos[num+1],2.0) + pow(ydist-ypos[num+1],2.);
                E_change = E_p(d_new,r_0,v_teether) - E_p(d_old,r_0,v_teether);
                d_old = pow(xpos[num]-xpos[num-1],2.) + pow(ypos[num]-ypos[num-1],2.0);
                d_new = pow(xdist-xpos[num-1],2.0) + pow(ydist-ypos[num-1],2.);
                E_change = E_change + E_p(d_new,r_0,v_teether) - E_p(d_old,r_0,v_teether);
                
                //energy changes- other plasmid beads
                for (i=num+1;i<=pl_max;i++){
                    dist = pow(xdist-xpos[i],2.)+pow(ydist-ypos[i],2.);
                    if (dist >= 1.5*rad_pl){
                        continue;
                    }
                    E_change = E_change + E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25);
                    dist = pow(xpos[num]-xpos[i],2.)+pow(ypos[num]- ypos[i],2.);
                    E_change = E_change - E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25);
                }
                for (i=1;i<num;i++){
                    dist = pow(xdist-xpos[i],2.)+pow(ydist- ypos[i],2.);
                    if (dist >= 1.5*rad_pl){
                        continue;
                    }
                    E_change = E_change + E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25) ;
                    dist = pow(xpos[num]-xpos[i],2.)+pow(ypos[num]-ypos[i],2.);
                    E_change = E_change - E_lj(dist,rad_pl,1.122*rad_pl,1.0,-0.25) ;
                }
                    
                //energy calculation - parA layer - if farther than 2.0*rad_pl then no calculation of anything
                if (ydist >= 2.5*rad_pl){
                    randno = ran2(seed);
                    if (E_change <= 0.0 || randno <= exp(-E_change)) {
                        xpos[num] = xdist;
                        ypos[num] = ydist;
                        n_good++;
                    } 
                    continue;
                }
                
                //plasmid-ParA hardwall calculation
                low = xpos[1];
                high = xpos[1];
                for (i=2 ;i <= pl_max;i++){
                	  if (xpos[i] < low){
                	  		low = xpos[i];	
                	  }
                	  if (xpos[i] > high){
                	  		high = xpos[i];
                	  }
                }

                grid_start = (int)low - 5;
                grid_end = (int)high + 5; 

                
                
                for (i=grid_start;i<=grid_end;i++){
                    dist = pow(xdist-Apos[i],2.)+pow(ydist,2.);
                    E_change = E_change + E_lj(dist,dA,1.122*dA,1.0,-0.25);
                    dist = pow(xpos[num]-Apos[i],2.)+ pow(ypos[num],2.);
                    E_change = E_change - E_lj(dist,dA,1.122*dA,1.0,-0.25);
                }
                
                //plasmid-ParA: if bead is not ParB then no attraction
                if (ratio[num] == 0){
                    randno = ran2(seed);
                    if (E_change <= 0.0 || randno <= exp(-E_change)) {
                        xpos[num] = xdist;
                        ypos[num] = ydist;
                        n_good++;
                    }
                    continue;
                }
                
                //energy calculation - parA layer- attractive
                for (i=grid_start;i<=grid_end;i++){
                    dist = pow(xdist-Apos[i],2.)+pow(ydist,2.);
                    E_change = E_change + removed[i]*E_sq(dist,d_c,pot_h);
                    dist = pow(xpos[num]-Apos[i],2.)+pow(ypos[num],2.);
                    E_change = E_change - removed[i]*E_sq(dist,d_c,pot_h);
                }
                
                randno = ran2(seed);
                if (E_change <= 0.0 || randno <= exp(-E_change)) {
                    xpos[num] = xdist;
                    ypos[num] = ydist;
                    n_good++;
                }
            }
            low = xpos[1];
            high = xpos[1];
            for (i=2 ;i <= pl_max;i++){
                if (xpos[i] < low){
                    low = xpos[i];
                }
                if (xpos[i] > high){
                    high = xpos[i];
                }
            }
            
            grid_start = (int)low - 5;
            grid_end = (int)high + 5;
            for (i=grid_start; i<grid_end; i++) {   //counter for parA
                E_change = 0.0;   //calculating total energy due to parA bead
                for (l=1; l<=pl_max; l++){   //counter for plasmid ring
                    dist = sqrt(pow(xpos[l]-Apos[i],2.)+pow(ypos[l],2.));
                    if (dist <= d_c){
                        E_change = E_change + pot_h*ratio[l];
                    }
                }
                //printf("%d\t%f\t%f\t%f\n",i,cod_x - Apos[i],cod_y,E_change);
                if (E_change < 0){
                    if (ran2(seed) <= p_r){
                        removed[i] = 0;
                        ent_pro = ent_pro +1;
                    }
                }
            }
        }
        cod_y = cod(ypos,pl_max);
        cod_x = cod(xpos,pl_max);
        //printf("%f\t%f\n",cod_x,cod_y);
        if (fabs(cod_x-75.0) >= bndry){
	    printf("%f\t%d\t %d\n",cod_x,k,ent_pro);
            break;
        }
        
    }
    fclose(output1);
    
    
    /* free matrices and vectors */
    free(xpos);
    free(ypos);
    free(ratio);
    free(Apos);

    return(0);
}


double ran2(long idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (idum <= 0) {
        if (-(idum) < 1) idum=1;
        else idum = -(idum);
        idum2=(idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(idum)/IQ1;
            idum=IA1*(idum-k*IQ1)-k*IR1;
            if (idum < 0) idum += IM1;
            if (j < NTAB) iv[j] = idum;
        }
        iy=iv[0];
    }
    k=(idum)/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

    
