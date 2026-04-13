#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "all_barrier.h"
#include <vector>
#include <boost/dynamic_bitset.hpp>	/* for dynamic_bitset objects (www.boost.org)*/
#include <fstream>
#include <string>
#include <array>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>

using namespace std;

//UNIQUE
bool lireFichier(const char *param_file, unsigned int &seed, int &n_ind, int &nb_loci, int &dim, int &gen_Max, double &mu, double &mu_mod, double &mig, double &rho, double &rho_LB, double &rho_CL, double &sig, double &amp, double &epistasis, vector<double> &opt1, vector<double> &opt2, double &var_pref, double &var_cost, int &nb_dmi, double &eps, int &dmi_neutral, double &ini_pref, double & ini_cost, int & skip_mate_choice, int & common_origin, int & export_metrics, int &cycle_mig_sel, int &load_mutation_file)
{
    FILE * fichierE;
    fichierE = fopen(param_file,"r");
    if (!fichierE)
    {
        cout << "Parameter file does not exist"<< endl;
        return 1;
    }
    int x;
    double temp;
    bool term;
    do
    {
        x = fgetc(fichierE);
    }
    while (!((x == '*') || (x == EOF)));
    if (x == EOF)
        term = true;
    else
    {
        term = false;
        if(fscanf(fichierE,"%*s %d ",&seed)!=1)
            term=true;
        //cout <<term << " " <<seed<< " "<<fgetc(fichierE) <<endl;
        if(fscanf(fichierE,"%*s %d ",&n_ind)!=1)
            term=true;
        //cout <<term << " "<<n_ind<<endl;
        if(fscanf(fichierE,"%*s %d ",&nb_loci)!=1)
            term=true;
        //cout <<term << " "<<nb_loci<<endl;
        if(fscanf(fichierE,"%*s %d ",&dim)!=1)
            term=true;
        //cout <<term << " "<<dim<<endl;
        if(fscanf(fichierE,"%*s %d ",&gen_Max)!=1)
            term=true;
        //cout <<term << " "<<gen_Max<<endl;
        if(fscanf(fichierE,"%*s %lf ",&mu)!=1)
            term=true;
        //cout <<term << " "<<mu<<endl;
        if(fscanf(fichierE,"%*s %lf ",&mu_mod)!=1)
            term=true;
        //cout <<term << " "<<mu<<endl;
        if(fscanf(fichierE,"%*s %lf ",&mig)!=1)
            term=true;
        //cout <<term << " "<<mig<<endl;
        if(fscanf(fichierE,"%*s %lf ",&rho)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %lf ",&rho_LB)!=1)
            term=true;
        //cout <<term << " "<<rho<<endl;
        if(fscanf(fichierE,"%*s %lf ",&rho_CL)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %lf ",&sig)!=1)
            term=true;
        //cout <<term << " "<<sig<<endl;
        if(fscanf(fichierE,"%*s %lf ",&amp)!=1)
            term=true;
        //cout <<term << " "<<amp<<endl;
        if(fscanf(fichierE,"%*s %lf ",&epistasis)!=1)
            term=true;
        //cout <<term << " "<<epistasis<<endl;
        for (int i=0; i<dim; i++)
        {
            if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
                term=true;
            //cout <<temp <<endl;
            opt1.push_back(temp);
        }
        for (int i=0; i<dim; i++)
        {
            if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
                term=true;
            opt2.push_back(temp);
        }
        if(fscanf(fichierE,"%*s %lf ",&var_pref)!=1)
            term=true;
        //cout <<term << " "<<var_pref<<endl;
        if(fscanf(fichierE,"%*s %lf ",&var_cost)!=1)
            term=true;
        //cout <<term << " "<<var_cost<<endl;
        if(fscanf(fichierE,"%*s %d ",&nb_dmi)!=1)
            term=true;
        //cout <<term << " "<<nb_dmi<<endl;
        if(fscanf(fichierE,"%*s %lf ",&eps)!=1)
            term=true;
        //cout <<term << " "<<eps<<endl;
        if(fscanf(fichierE,"%*s %d ",&dmi_neutral)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %lf ",&ini_pref)!=1)
            term=true;
        //cout <<term << " "<<ini_pref<<endl;
        if(fscanf(fichierE,"%*s %lf ",&ini_cost)!=1)
            term=true;
        //cout <<term << " "<<ini_cost<<endl;
        if(fscanf(fichierE,"%*s %d ",&skip_mate_choice)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %d ",&common_origin)!=1)
            term=true;
        //cout <<term << " "<<common_origin<<endl;
        if(fscanf(fichierE,"%*s %d ",&export_metrics)!=1)
            term=true;
        if (fscanf(fichierE,"%*s %d ",&cycle_mig_sel)!=1)
            term=true;
        if (fscanf(fichierE,"%*s %d ",&load_mutation_file)!=1)
            term=true;
        x = fgetc(fichierE);
        if (x != '*')
            term=true;
    }
    fclose(fichierE);
    return term;
}

//UNIQUE
void migration(gsl_rng* r,Ind *pop1, Ind *pop2, Ind *temp1, const double &mig, int &nb_ind_1, int &nb_ind_2)
{
    int i, nb_mig_1to2, nb_mig_2to1;
    if (mig==0) // skip if migration rate is 0
        return;

    nb_mig_1to2=gsl_ran_poisson(r,mig*nb_ind_1); //number of migrants from pop1 to 2
    nb_mig_2to1=gsl_ran_poisson(r,mig*nb_ind_2); //number of migrants from pop2 to 1
    if (nb_mig_1to2+nb_mig_2to1==0)
        return;

    if (nb_mig_1to2>nb_ind_1)
        nb_mig_1to2=nb_ind_1;//nb migrants capped by nb of individuals

    if (nb_mig_2to1>nb_ind_2)
        nb_mig_2to1=nb_ind_2;//nb migrants capped by nb of individuals

    nb_ind_1=nb_ind_1-nb_mig_1to2;//number of individual staying in the pop1
    for (i=0; i<nb_mig_1to2; i++)// move the individual from pop 1 in a temporary vector; starting from individual in position nb_ind_1
        temp1[i]=pop1[i+nb_ind_1];

    nb_ind_2=nb_ind_2-nb_mig_2to1;//number of individual staying in the pop2
    for (i=0; i<nb_mig_2to1; i++)// transfer directly individuals from pop2 to pop1
        pop1[i+nb_ind_1]=pop2[i+nb_ind_2];// move the migrant at the end of pop1, starting in position nb_ind_1, directly from pop2

    nb_ind_1+=nb_mig_2to1;// add migrants to counter
    for (i=0; i<nb_mig_1to2; i++)// transfer from the temporary vector them to the pop2
        pop2[i+nb_ind_2]=temp1[i];// move the migrant at the end of pop2, starting in position nb_ind_2,  from the temp placeholder

    nb_ind_2+=nb_mig_1to2;// add migrants to counter

}




void initialize_pop(Ind * pop,const int &n_ind, const int &nb_Lblocks, const int & nb_loci_per_LG, const int &dim, const double &pref, const double &cost, const int &common_origin, const vector<double> &opt, const int &locus, const int &nb_dmi, const int *dmi_pos, const int * target_chrom, const double * mutations)
{
    int index_ind, index_chrom, index_dim;
    for (index_ind=0; index_ind<n_ind; index_ind++)
    {
        for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)
        {
            pop[index_ind].chr_a[index_chrom].sel.resize(nb_loci_per_LG);
            pop[index_ind].chr_b[index_chrom].sel.resize(nb_loci_per_LG);
            pop[index_ind].chr_a[index_chrom].pheno.resize(dim);
            pop[index_ind].chr_b[index_chrom].pheno.resize(dim);

        }
        pop[index_ind].pheno.resize(dim);
        pop[index_ind].pref_a=pref;
        pop[index_ind].pref_b=pref;
        pop[index_ind].cost_a=cost;
        pop[index_ind].cost_b=cost;
    }
    //ancestral state, no LA, no DMI
    // we determine it only for individual 0, at it may be modified downstream

    for (index_dim=0; index_dim<dim; index_dim++)
    {
        for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)
        {
            pop[0].chr_a[index_chrom].pheno[index_dim]=0;
            pop[0].chr_b[index_chrom].pheno[index_dim]=0;
        }
        pop[0].pheno[index_dim]=0;
    }

    switch (common_origin)
    {
    case 0:
        break;
    case 1:// with LA, pop close to optimum
    {
        move_to_optimum(pop[0],dim, nb_Lblocks, nb_loci_per_LG, opt, target_chrom, mutations);
        break;
    }
    case 2://with DMI, no LA
    {
        add_dmi(pop[0],dim,nb_dmi,nb_loci_per_LG, locus, dmi_pos, target_chrom, mutations); // only chr_a has been edited to accurately reflect what we want; chr_b has been edited to keep track of the loci to not edit for changing the phenotype

        vector <double> opt_ancestral(dim,0); // of an haplotype

        move_to_optimum(pop[0],dim, nb_Lblocks, nb_loci_per_LG, opt_ancestral, target_chrom, mutations);

        break;
    }
    case 3:  //with DMI, and LA
    {
        add_dmi(pop[0],dim,nb_dmi, nb_loci_per_LG, locus, dmi_pos, target_chrom, mutations);// only chr_a has been edited to accurately reflect what we want; chr_b has been edited to keep track of the loci to not edit for changing the phenotype

        move_to_optimum(pop[0],dim, nb_Lblocks, nb_loci_per_LG, opt, target_chrom, mutations);

        break;
    }
    case 4:
    {
        break;
    }
    default:
    {
        cout << "wrong values for parameter common origin. it has to be 0 (no LA and no DMI), 1 (LA), 2 (DMI), 3 (LA and DMI) or 4 (loading from a file)"<<endl;
        exit(0);
    }
    }
    for (index_ind=1; index_ind<n_ind; index_ind++)
    {
        // apply to all other individuals
        for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)
        {
            pop[index_ind].chr_a[index_chrom].sel=pop[0].chr_a[index_chrom].sel;
            pop[index_ind].chr_b[index_chrom].sel=pop[0].chr_b[index_chrom].sel;
            pop[index_ind].chr_a[index_chrom].pheno=pop[0].chr_a[index_chrom].pheno;
            pop[index_ind].chr_b[index_chrom].pheno=pop[0].chr_b[index_chrom].pheno;
        }
        pop[index_ind].pheno=pop[0].pheno;
    }
    for (index_dim=0; index_dim<dim; index_dim++)
        cout << pop[0].pheno[index_dim]<< " ";
    cout <<endl;
}


// function to edit a single individual to be relatively close to a fix point. If DMI are initially fixed, there are sites to not edit, respectively tracked by chr_a, with the derived sites we want to keep as such, and by chr_b with site we want to keep ancestral

void move_to_optimum(Ind &ind, const int &dim, const int &nb_Lblocks, const int & nb_loci_per_LG, const vector<double> &opt, const int * target_chrom, const double * mutations)
{
    int index_dim, index_chrom, index_locus;
    double d;
    vector <double> starting_pheno(dim,0); // of an haplotype
    for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)// calculate the current phenotype due to possible preexisting mutations at the haplotype level
    {
        for (index_dim=0; index_dim<dim; index_dim++)
            starting_pheno[index_dim]+=ind.chr_a[index_chrom].pheno[index_dim];
    }
    vector <double> temp_pheno(dim,0);
    double new_pheno_distance=0;
    double current_best_pheno_distance=0;// calculate current distance
    for (index_dim=0; index_dim<dim; index_dim++)
    {
        current_best_pheno_distance+=(2*starting_pheno[index_dim]-opt[index_dim])*(2*starting_pheno[index_dim]- opt[index_dim]); // this is the distance between the fully incompatible diploid individual and the local optimum
        //cout <<"dim "<< index_dim << " optimum " << opt[index_dim] << " starting pheno " << starting_pheno[index_dim] << endl;
    }
    //cout << "starting distance " << current_best_pheno_distance <<endl;
    for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)
    {
        for (index_locus=0; index_locus<nb_loci_per_LG; index_locus++) // together with previous line, this loops over all loci
        {
            if((!ind.chr_a[index_chrom].sel.test(index_locus))&(!ind.chr_b[index_chrom].sel.test(index_locus))) // we only want cases where the locus does not harbor a DMI, be it A or B
            {
                new_pheno_distance=0;
                for (index_dim=0; index_dim<dim; index_dim++)
                {
                    //recalculate square phenotypic distance for diploid individual and with one extra mutations as homozygote
                    d=(2*starting_pheno[index_dim]+ 2*mutations[(index_chrom*nb_loci_per_LG+index_locus)*dim+index_dim] -opt[index_dim]);
                    new_pheno_distance+= d*d;
                    //    cout << mutations[(index_chrom*nb_loci_per_LG+index_locus)*dim+index_dim] << " ";
                }
                //  cout << "new_pheno_distance "<<new_pheno_distance<< endl;
                //cout << new_pheno_distance <<"--"<<current_best_pheno_distance<<" ";
                if (new_pheno_distance<current_best_pheno_distance)// we only test smaller beacsue we are using the square distance
                {
                    // given that this is initialization,we only modify the genome "a" of the first individual and will then paste it on all at the end.
                    ind.chr_a[index_chrom].sel[index_locus].flip();
                    current_best_pheno_distance=new_pheno_distance;// update best pheno square distance
                    //  cout << "starting pheno " ;
                    for (index_dim=0; index_dim<dim; index_dim++)// update new starting_pheno
                    {
                        starting_pheno[index_dim]+=mutations[(index_chrom*nb_loci_per_LG+index_locus)*dim+index_dim];
                        ind.chr_a[index_chrom].pheno[index_dim]+=mutations[(index_chrom*nb_loci_per_LG+index_locus)*dim+index_dim];
                        //  cout << starting_pheno[index_dim] << " ";
                    }
                    // cout << endl;
                }
            }
            //      else
            //   {
            //       cout << "mutated site"<<endl;
            //   }
        }
        // finally we are updating chr_b,
        ind.chr_b[index_chrom].sel=ind.chr_a[index_chrom].sel;
        ind.chr_b[index_chrom].pheno=ind.chr_a[index_chrom].pheno;
        //compute the individual fitness
        for (index_dim=0; index_dim<dim; index_dim++)
        {
            ind.pheno[index_dim]+=ind.chr_a[index_chrom].pheno[index_dim]+ind.chr_b[index_chrom].pheno[index_dim];
        }
    }
    // cout << endl;
    //cout << "ending distance " << current_best_pheno_distance <<  " ind pheno "<< ind.pheno[0]<<endl;

}

//function add DMI to the genome - only target chr_a of the target - chr_b will be later updated by the moved optimum function - here chr_b is used to track the DMI sites - the locus variable is used to indicate whether we want the A locus to be derived (0) or the B locus (1). This is done for a single individual.
void add_dmi(Ind &ind, const int &dim, const int &nb_dmi, const int &nb_loci_per_LG, const int &locus, const int *dmi_pos, const int *target_chrom, const double *mutations)
{
    int start=0;
    int alt=0;
    int index_chrom, index_locus, index_dim;
    if (locus==0) //decide whether we want the first locus of the DMI or the second
    {
        start=0;
        alt=1;
    }
    else
    {
        start=1;
        alt=0; //we need to do this because, we need to know which positions are involved in the DMIs
    }
    for (int i=0; i<nb_dmi; i++)
    {
        index_chrom=target_chrom[2*i+start];
        index_locus=dmi_pos[2*i+start];
        //only chr_A is edited to track the DMI we want; chr_B has a different role; see below
        ind.chr_a[index_chrom].sel.flip(index_locus);

        for (index_dim=0; index_dim<dim; index_dim++)
            ind.chr_a[index_chrom].pheno[index_dim]+=mutations[(index_chrom*nb_loci_per_LG+index_locus)*dim+index_dim];// the phenotype is updated accordingly

        //this is a temporary editing of chr_b to track where the alternative locus are tracked for the DMI. Chr_b will be overwritten in the move_optimum function.
        index_chrom=target_chrom[2*i+alt];
        index_locus=dmi_pos[2*i+alt];
        ind.chr_b[index_chrom].sel.flip(index_locus);
    }
}



void new_gamete(gsl_rng* r, double * mutations, const int &n_ind,const double &mu_GW, const double &mu_mod, const int *parent1_ID, const int *parent2_ID, Ind *parent1, Ind *parent2, Ind *offspring, const int &nb_loci, const int &dim, const int &nb_LD_blocks,const double &rho, const double &rho_LG, const double &rho_LB, const double &rho_CL, const double & var_pref, const double & var_cost, int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
{
    // int mut, ns, j, site, co_nb,index_dim,inherit_mkL;
    int  index_offsprg, index_chrom;
    double mut_pref, mut_cost,has_rec;
    bool strand_a_p1=false, strand_a_p2=false;

    for (index_offsprg = 0; index_offsprg < n_ind; index_offsprg++) // loop over all offspring chromosomes inherited from a given parent; origin decide if we are looping through odd (not choosy parent) or even (choosy parent) indices; move by steps of 2 to get from one individual to the next.
        // this means all element of offspring are being overwritten at least from position 0 to n_ind-1
    {
        //cout <<"ind_nb "<<index_offsprg<<endl;
        if (gsl_rng_uniform(r) > 0.5)
            offspring[index_offsprg].sex=1;
        else
            offspring[index_offsprg].sex=0;

        // inherit the freely recombining loci
        //inherit_FR_loci(r,offspring[index_offsprg].pref_a,parent1[parent1_ID[index_offsprg]].pref_a,parent1[parent1_ID[index_offsprg]].pref_b);
        //inherit_FR_loci(r,offspring[index_offsprg].pref_b,parent2[parent2_ID[index_offsprg]].pref_a,parent2[parent2_ID[index_offsprg]].pref_b);
        inherit_FR_loci(r,offspring[index_offsprg].cost_a,parent1[parent1_ID[index_offsprg]].cost_a,parent1[parent1_ID[index_offsprg]].cost_b);
        inherit_FR_loci(r,offspring[index_offsprg].cost_b,parent2[parent2_ID[index_offsprg]].cost_a,parent2[parent2_ID[index_offsprg]].cost_b);
        inherit_FR_loci(r,offspring[index_offsprg].neutral_a,parent1[parent1_ID[index_offsprg]].neutral_a,parent1[parent1_ID[index_offsprg]].neutral_b);
        inherit_FR_loci(r,offspring[index_offsprg].neutral_b,parent2[parent2_ID[index_offsprg]].neutral_a,parent2[parent2_ID[index_offsprg]].neutral_b);

        // randomly choosing which strand to start from
        if (gsl_rng_uniform(r) > 0.5)
        {
            strand_a_p1=true;
        }
        else
        {
            strand_a_p1=false;
        }
        if (gsl_rng_uniform(r) > 0.5)
        {
            strand_a_p2=true;
        }
        else
        {
            strand_a_p2=false;
        }

        // inheritance of the prefe locus if it is freely recombining
        if (rho_CL==0.5)
        {
            inherit_FR_loci(r,offspring[index_offsprg].pref_a,parent1[parent1_ID[index_offsprg]].pref_a,parent1[parent1_ID[index_offsprg]].pref_b);
            inherit_FR_loci(r,offspring[index_offsprg].pref_b,parent2[parent2_ID[index_offsprg]].pref_a,parent2[parent2_ID[index_offsprg]].pref_b);
        }
        else
        {
            // inheritance of the pref locus in case of non FR case - strand A for the offspring
            // determine whether the mate choice locus recombines away from the first linkage block
            has_rec=gsl_rng_uniform(r);
            if (strand_a_p1 & (has_rec>=rho_CL))
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_a; // inherit the strand from the A chromosome from parent 1; no rec
            }
            else if (strand_a_p1 & (has_rec<rho_CL))
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_b; // inherit the strand from the B chromosome from parent 1; rec
            }
            else if ((!strand_a_p1) & (has_rec>=rho_CL))
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_b; // inherit the strand from the B chromosome from parent 1; no rec
            }
            else if ((!strand_a_p1) & (has_rec<rho_CL))
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_a; // inherit the strand from the A chromosome from parent 1; rec
            }
            // inheritance of the pref locus in case of non FR case - strand B for the offspring
            // determine whether the mate choice locus recombines away from the first linkage block
            has_rec=gsl_rng_uniform(r);
            if (strand_a_p2 & (has_rec>=rho_CL))
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_a; // inherit the strand from the A chromosome from parent 2; no rec
            }
            else if (strand_a_p2 & (has_rec<rho_CL))
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_b; // inherit the strand from the B chromosome from parent 2; rec
            }
            else if ((!strand_a_p2) & (has_rec>=rho_CL))
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_b; // inherit the strand from the B chromosome from parent 2; no rec
            }
            else if ((!strand_a_p2) & (has_rec<rho_CL))
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_a; // inherit the strand from the A chromosome from parent 2; rec
            }
        }

        //cout << "nb CO " << co_nb<<endl;
        //cout << (*(genotypes[chr_ID])).sel<< endl;
        //cout << (*(genotypes[chr_ID_alt])).sel<< endl;

        if(gsl_rng_uniform(r)<mu_mod) // mutation at the preference locus, strand A
        {
            mut_pref=gsl_ran_gaussian(r,var_pref);
            offspring[index_offsprg].pref_a += mut_pref;
            if (offspring[index_offsprg].pref_a<0)
                offspring[index_offsprg].pref_a=0;
            if (offspring[index_offsprg].pref_a>0.5)
                offspring[index_offsprg].pref_a=0.5;
        }
        if(gsl_rng_uniform(r)<mu_mod) // mutation at the preference locus, strand B
        {
            mut_pref=gsl_ran_gaussian(r,var_pref);
            offspring[index_offsprg].pref_b += mut_pref;
            if (offspring[index_offsprg].pref_b<0)
                offspring[index_offsprg].pref_b=0;
            if (offspring[index_offsprg].pref_b>0.5)
                offspring[index_offsprg].pref_b=0.5;
        }


        if(gsl_rng_uniform(r)<mu_mod)
        {
            mut_cost=gsl_ran_gaussian(r,var_cost);
            offspring[index_offsprg].cost_a += mut_cost;
            if (offspring[index_offsprg].cost_a<0.5)
                offspring[index_offsprg].cost_a=0.5;
        }
        if(gsl_rng_uniform(r)<mu_mod)
        {
            mut_cost=gsl_ran_gaussian(r,var_cost);
            offspring[index_offsprg].cost_b += mut_cost;
            if (offspring[index_offsprg].cost_b<0.5)
                offspring[index_offsprg].cost_b=0.5;
        }


        for (index_chrom=0; index_chrom<nb_LD_blocks; index_chrom++)
        {
            gamete_linkage_block(r, rho_LG, rho_LB, dim, nb_loci, mu_GW, parent1[parent1_ID[index_offsprg]], index_chrom, offspring[index_offsprg].chr_a[index_chrom], mutations, strand_a_p1, Co_pos, off1, off2, rec, pos);
            gamete_linkage_block(r, rho_LG, rho_LB, dim, nb_loci, mu_GW, parent2[parent2_ID[index_offsprg]], index_chrom, offspring[index_offsprg].chr_b[index_chrom], mutations, strand_a_p2, Co_pos, off1, off2, rec, pos);

        }

    }
}




//
void update_mean_metrics(const int & nb_ind, const int &dim, double &mean_fit, double &mean_fit_abs, double &mean_fit_dmi, vector <double> &mean_pheno, double &mean_pref, double &mean_cost, const Ind *pop)
{
    mean_fit_abs=0;
    mean_fit_dmi=0;
    fill(mean_pheno.begin(), mean_pheno.end(), 0);
    mean_pref=0;
    mean_cost=0;
    mean_fit=0;
    int i, j;
    for (i=0; i <nb_ind; i++)
    {
        mean_fit+=(pop[i].abs_fit*pop[i].dmi_fit);
        mean_fit_abs+=pop[i].abs_fit;
        mean_fit_dmi+=pop[i].dmi_fit;
        mean_pref+=pop[i].pref;
        mean_cost+=pop[i].cost;
        for (j=0; j<dim; j++)
        {
            mean_pheno[j]+=pop[i].pheno[j];
        }
        //cout << pop[i].pheno[0] <<" ";
    }
    mean_fit_abs/=nb_ind;
    mean_fit_dmi/=nb_ind;
    mean_pref/=nb_ind;
    mean_cost/=nb_ind;
    mean_fit/=nb_ind;
    for (j=0; j<dim; j++)
    {
        mean_pheno[j]/=nb_ind;
    }
}


void export_pop(ofstream &fout, const int &nb_ind, const int &dim, const int &nb_Lblocks, Ind *pop)
{
    for (int i=0; i< nb_ind; i++)
    {
        fout << pop[i].neutral_a << " " << pop[i].pref_a << " " << pop[i].cost_a << " " ;
        for (int j=0; j<nb_Lblocks; j++)
        {
            for (int k =0; k<dim; k++)
                fout << pop[i].chr_a[j].pheno[k]<< " ";
            fout << pop[i].chr_a[j].sel <<" ";
        }
        fout <<endl;
        fout << pop[i].neutral_b << " " << pop[i].pref_b << " " << pop[i].cost_b << " " ;
        for (int j=0; j<nb_Lblocks; j++)
        {
            for (int k =0; k<dim; k++)
                fout << pop[i].chr_b[j].pheno[k]<< " ";
            fout << pop[i].chr_b[j].sel <<" ";
        }
        fout <<endl;
    }
}


//
//void inherit_FR_loci(gsl_rng* r, bool &new_value, const bool &chr_1, const bool &chr_2)
//{
//    if (chr_1==chr_2)
//    {
//        new_value=chr_1;
//    }
//    else
//    {
//        if (gsl_ran_bernoulli(r,0.5))
//            new_value = chr_1;
//        else
//            new_value = chr_2;
//    }
//}


//void more_statistics(const int &nb_ind, Chr **hap, const vector<double> &opt, const double &amp, const double &half_epistasis, double &mean_hom_fit, double &h_0, const double &nb_locus, const int &dim, double &mean_origin)
//{
//    int i,j;
//    double d1, d2;
//    h_0=0;
//    mean_hom_fit=0;
//    mean_origin=0;
//    for (i=0; i <nb_ind; i++)
//    {
//        d1=0;
//        d2=0;
//        for (j=0; j<dim; j++)
//        {
//            d1+=((*(hap[2*i])).pheno[j]*2-opt[j]) * ((*(hap[2*i])).pheno[j]*2-opt[j]);
//            d2+=((*(hap[2*i+1])).pheno[j]*2-opt[j]) * ((*(hap[2*i+1])).pheno[j]*2-opt[j]);
//        }
//        mean_origin+=(*(hap[2*i])).neutral+(*(hap[2*i+1])).neutral;
//        mean_hom_fit+=exp(-amp * pow(d1, half_epistasis));
//        mean_hom_fit+=exp(-amp * pow(d2, half_epistasis));
//        for (j=0; j<nb_locus; j++)
//        {
//            h_0+=double((*(hap[2*i])).sel[j] ^ (*(hap[2*i+1])).sel[j]);
//        }
//    }
//    mean_hom_fit/=(2*nb_ind);
//    mean_origin/=(2*nb_ind);
//    h_0/=(nb_ind*nb_locus);
//}

//void parent_mean_metrics(const int & nb_ind, double &mean_fit_abs, double &mean_fit_dmi, double &mean_pheno, double &mean_pref, double &mean_cost, const Ind *pop, const int* parent)
//{
//    mean_fit_abs=0;
//    mean_fit_dmi=0;
//    mean_pheno=0;
//    mean_pref=0;
//    mean_cost=0;
//    int i;
//    for (i=0; i <nb_ind; i++)
//    {
//        mean_fit_abs+=pop[parent[i]].abs_fit;
//        mean_fit_dmi+=pop[parent[i]].dmi_fit;
//        mean_pheno+=pop[parent[i]].pheno[0];
//        mean_pref+=pop[parent[i]].pref;
//        mean_cost+=pop[parent[i]].cost;
//        // if(pop[parent[i]].dmi_fit<1){
//        //     cout <<"met "<<pop[i].dmi_fit<<endl;
//        // }
//    }
//    mean_fit_abs/=nb_ind;
//    mean_fit_dmi/=nb_ind;
//    mean_pheno/=nb_ind;
//    mean_pref/=nb_ind;
//    mean_cost/=nb_ind;
//}
//
//void parent_stat(const int & nb_ind, const int* parent, const int & nb_parent, double &var)
//{
//    double mean=nb_ind/nb_parent;
//    var=0;
//    int i;
//    int * parent_count = new int [nb_parent];
//
//    for (i=0; i<nb_parent; i++)
//        parent_count[i]=0;
//
//    for (i=0; i<nb_ind; i++)
//        parent_count[parent[i]]++;
//
//    for (i=0; i<nb_parent; i++)
//        var+= (parent_count[i]-mean)*(parent_count[i]-mean);
//    var/=nb_parent;
//    delete [] parent_count;
//
//}


//void choose_parents(gsl_rng* r, int *parent, const int &n_off, const int &nb_choosy, const double *fitness_choosy, const int &nb_other, const double *fitness_other)
//{
//	int candidate;
//	for (int j = 0; j < n_off; j++){
//	//cout << j<< endl;
//	//cout <<nb_choosy << " " << nb_other<<endl;
//		do
//		{
//			candidate = int(gsl_rng_uniform(r) * nb_choosy);// choose random individual
//		} while (gsl_rng_uniform(r) > fitness_choosy[candidate]); // if it passes test
//
//		parent[2*j] = candidate; //assign this individual as parent 1; the entry here corresponds to the individual j, funf in position 2*j and 2j+1 in the pop vactor.
//		do
//		{
//			candidate = int(gsl_rng_uniform(r) * nb_other);
//			//cout << candidate<<endl;
//
//		} while (gsl_rng_uniform(r) > fitness_other[candidate]);
//
//		parent[2*j+1] = candidate;
//		//cout<< parent[2*j]<< " " <<parent[2*j+1] <<"-";
//	}
//} // modify the parent vector; a 2N vector, with parents given in 2j and 2j+1 positions
//
//void choose_parents(gsl_rng* r, int *parent, const int &n_off, const int &nb_choosy, const double *fitness_choosy, const int &nb_other, const double *fitness_other, const double *pref, const double *cost)
//{
//	int candidate, cumul_cost;
//	for (int j = 0; j < n_off; j++){
//	//cout << j<< endl;
//	//cout <<nb_choosy << " " << nb_other<<endl;
//        cumul_cost=0;
//		do
//		{
//			candidate = int(gsl_rng_uniform(r) * nb_choosy);// choose random individual
//		} while (gsl_rng_uniform(r) > fitness_choosy[candidate]); // if it passes test
//
//		parent[2*j] = candidate; //assign this individual as parent 1; the entry here corresponds to the individual j, funf in position 2*j and 2j+1 in the pop vactor.
//		do
//		{
//			candidate = int(gsl_rng_uniform(r) * nb_other); //NEED SMARTER WAY HERE. MULTINOM?
//			//cout << candidate<<endl;
//
//		} while ((gsl_rng_uniform(r) > fitness_other[candidate]) & (gsl_rng_uniform(r)>cumul_cost*cost[j]));
//
//		parent[2*j+1] = candidate;
//		//cout<< parent[2*j]<< " " <<parent[2*j+1] <<"-";
//	}
//} // modify the parent vector; a 2N vector, with parents given in 2j and 2j+1 positions


//void compress_chr()
//
//void update_fitness(const int &nb_ind, const int &dim, Chr **pop, const double &amp, const double &half_epistasis, const vector<double> &opt, double *fitness)
//{
//	double w_max = 0;
//	double d, x;
//	int i;
//	for (i = 0; i < nb_ind; i++){
//		d = 0;
//		for (int j = 0; j < dim; j++){
//			x = (*(pop[2*i])).pheno[j]+ (*(pop[2*i+1])).pheno[j];
//			d += (x - opt[j]) * (x - opt[j]); // "d" is the square of the distance to the optimum
//		}
//		// fitness of individual i:
//		fitness[i] = exp(-amp * pow(d, half_epistasis));
//		if (w_max < fitness[i])
//			w_max = fitness[i];
//	}
//	for (i =0; i < nb_ind; i++){
//		fitness[i] /= w_max;
////		cout <<fitness[i]<<" ";
//		}
//}
//
//void update_fitness(const int &nb_ind, const int &dim, Chr **pop, const double &amp, const double &half_epistasis, const vector<double> &opt, double *fitness, double *pref, double *cost)
//{
//	double w_max = 0;
//	double d, x;
//	int i;
//	for (i = 0; i < nb_ind; i++){
//		d = 0;
//		pref[i]=(*(pop[2*i])).pref+ (*(pop[2*i+1])).pref;
//		cost[i]=(*(pop[2*i])).cost+ (*(pop[2*i+1])).cost;
//		for (int j = 0; j < dim; j++){
//			x = (*(pop[2*i])).pheno[j]+ (*(pop[2*i+1])).pheno[j];
//			d += (x - opt[j]) * (x - opt[j]); // "d" is the square of the distance to the optimum
//		}
//		// fitness of individual i:
//		fitness[i] = exp(-amp * pow(d, half_epistasis));
//		if (w_max < fitness[i])
//			w_max = fitness[i];
//	}
//	for (i =0; i < nb_ind; i++){
//		fitness[i] /= w_max;
////		cout <<fitness[i]<<" ";
//		}
//}

//
//void new_gamete(gsl_rng* r, double * mutations, const int &n_ind,const double &mu_GW, const double &mu_mod, const int *parent1_ID, const int *parent2_ID, Ind *parent1, Ind *parent2, Ind *offspring, const int &nb_loci, const int &dim, const int &nb_LD_blocks,const double &rho, const double &rho_LG, const double &rho_LB, const double & var_pref, const double & var_cost, int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
//{
//
//    int mut, ns, j, site, co_nb, index_offsprg, index_chrom, index_dim;
//    double mut_pref, mut_cost, pheno;
//    bool strand_a_p1=false, strand_a_p2=false;
//
//    for (index_offsprg = 0; index_offsprg < n_ind; index_offsprg++) // loop over all offsprings chromosome inherited from a given parent; origin decide if we are looping through odd (not choosy parent) or even (choosy parent) indices; move by steps of 2 to get from ind to the next.
//    {
//        //cout <<"ind_nb "<<index_offsprg<<endl;
//        if (gsl_rng_uniform(r)>0.5)
//            offspring[index_offsprg].sex=1;
//        else
//            offspring[index_offsprg].sex=0;
//
//        // inherit the freely recombining loci
//        inherit_FR_loci(r,offspring[index_offsprg].pref_a,parent1[parent1_ID[index_offsprg]].pref_a,parent1[parent1_ID[index_offsprg]].pref_b);
//        inherit_FR_loci(r,offspring[index_offsprg].pref_b,parent2[parent2_ID[index_offsprg]].pref_b,parent2[parent2_ID[index_offsprg]].pref_b);
//        inherit_FR_loci(r,offspring[index_offsprg].cost_a,parent1[parent1_ID[index_offsprg]].cost_a,parent1[parent1_ID[index_offsprg]].cost_b);
//        inherit_FR_loci(r,offspring[index_offsprg].cost_b,parent2[parent2_ID[index_offsprg]].cost_a,parent2[parent2_ID[index_offsprg]].cost_b);
//        inherit_FR_loci(r,offspring[index_offsprg].neutral_a,parent1[parent1_ID[index_offsprg]].neutral_a,parent1[parent1_ID[index_offsprg]].neutral_b);
//        inherit_FR_loci(r,offspring[index_offsprg].neutral_b,parent2[parent2_ID[index_offsprg]].neutral_a,parent2[parent2_ID[index_offsprg]].neutral_b);
//
//
//        if (gsl_rng_uniform(r)>0.5)
//            strand_a_p1=true;
//        if (gsl_rng_uniform(r)>0.5)
//            strand_a_p2=true;
//
//
//        //cout << "nb CO " << co_nb<<endl;
//        //cout << (*(genotypes[chr_ID])).sel<< endl;
//        //cout << (*(genotypes[chr_ID_alt])).sel<< endl;
//        if(gsl_rng_uniform(r)<mu_mod)
//        {
//            mut_pref=gsl_ran_gaussian(r,var_pref);
//            offspring[index_offsprg].pref_a += mut_pref;
//            if (offspring[index_offsprg].pref_a<0)
//                offspring[index_offsprg].pref_a=0;
//            if (offspring[index_offsprg].pref_a>0.5)
//                offspring[index_offsprg].pref_a=0.5;
//        }
//        if(gsl_rng_uniform(r)<mu_mod)
//        {
//            mut_pref=gsl_ran_gaussian(r,var_pref);
//            offspring[index_offsprg].pref_b += mut_pref;
//            if (offspring[index_offsprg].pref_b<0)
//                offspring[index_offsprg].pref_b=0;
//            if (offspring[index_offsprg].pref_b>0.5)
//                offspring[index_offsprg].pref_b=0.5;
//        }
//
//
//        if(gsl_rng_uniform(r)<mu_mod)
//        {
//            mut_cost=gsl_ran_gaussian(r,var_cost);
//            offspring[index_offsprg].cost_a += mut_cost;
//            if (offspring[index_offsprg].cost_a<0.5)
//                offspring[index_offsprg].cost_a=0.5;
//        }
//
//
//        if(gsl_rng_uniform(r)<mu_mod)
//        {
//            mut_cost=gsl_ran_gaussian(r,var_cost);
//            offspring[index_offsprg].cost_b += mut_cost;
//            if (offspring[index_offsprg].cost_b<0.5)
//                offspring[index_offsprg].cost_b=0.5;
//        }
//
//
//        for (index_chrom=0; index_chrom<nb_LD_blocks; index_chrom++)
//        {
//            co_nb = gsl_ran_poisson(r,rho_LG);
//            if (co_nb>(nb_loci-1))
//                co_nb=(nb_loci-1);
//            mut = gsl_ran_poisson(r,mu_GW); // number of mutations at selected loci on genome segment
//            if (((co_nb==0)|(parent1[parent1_ID[index_offsprg]].chr_a[index_chrom].sel==parent1[parent1_ID[index_offsprg]].chr_b[index_chrom].sel)))
//            {
//                //cout << "no recomb"<<endl;
//                offspring[index_offsprg].chr_a[index_chrom].sel=parent1[parent1_ID[index_offsprg]].chr_a[index_chrom].sel;
//
//                for (index_dim = 0; index_dim < dim; index_dim++)
//                    offspring[index_offsprg].chr_a[index_chrom].pheno[index_dim] = parent1[parent1_ID[index_offsprg]].chr_a[index_chrom].pheno[index_dim];
//                if (mut>0)
//                {
//                    for (j = 0; j < mut; j++)
//                    {
//                        site = int(gsl_rng_uniform(r)*nb_loci);// note that this allows for double mutation at the same site (with probability 1/nb_loci²)
//                        //cout << "site " << site<< endl;
//
//                        ns = (nb_loci*index_chrom+site)*dim;
//                        // cout <<"site "<<site <<endl;
//                        // cout << chrom[free].sel[site]<< endl;
//                        if (offspring[index_offsprg].chr_a[index_chrom].sel.test(site))
//                        {
//                            for (index_dim = 0; index_dim < dim; index_dim++)
//                                offspring[index_offsprg].chr_a[index_chrom].pheno[index_dim] -= mutations[ns + index_dim];
//                        }
//                        else
//                        {
//                            for (index_dim = 0; index_dim < dim; index_dim++)
//                            {
//                                //    cout << k<< " "<<chrom[free].pheno[k]<< endl;
//                                offspring[index_offsprg].chr_a[index_chrom].pheno[index_dim] += mutations[ns + index_dim];
//                                //      cout <<"done"<<endl;
//                            }
//                            //    cout <<"loop done ";
//                        }
//                        //cout <<"flip "<< site<< endl;
//                        // cout << chrom[free].sel<<endl;
//                        offspring[index_offsprg].chr_a[index_chrom].sel.flip(site);
//                        //   cout <<"flip done"<<endl;
//                    }
//                    //cout << "new chr " << chrom[free].sel<<endl;
//                }
//                //    cout <<nb << " "<<free<<"\n";
//                //  cout <<nb << " "<<free<<"\n";
//            }
//            else
//            {
//                //cout << "recomb"<<endl;
//
//                //cout << (*(pop[chr_ID])).nbchr<<endl;
//                //cout << "parent"<< chr_ID << "-"<< chr_ID_alt << " all_chr_ent"<< (*(genotypes[chr_ID])).ID<< " "<< (*(genotypes[(chr_ID_alt)])).ID<<endl;
//                //rec(r, chrom[free], (*(genotypes[chr_ID])), (*(genotypes[(chr_ID_alt)])), co_nb, nb_loci);
//                //cout << parent1[parent1_ID[nb]].chr_a[i].sel<<" "<< parent1[parent1_ID[nb]].chr_b[i].sel<<endl;
//                //cout <<index_offsprg <<" "<<index_chrom<< " "<<offspring[index_offsprg].chr_a[index_chrom].sel<<endl;
//                if(strand_a_p1)
//                    rec_v2(r, offspring[index_offsprg].chr_a[index_chrom], parent1[parent1_ID[index_offsprg]].chr_a[index_chrom], parent1[parent1_ID[index_offsprg]].chr_b[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec, pos);
//                else
//                    rec_v2(r, offspring[index_offsprg].chr_a[index_chrom], parent1[parent1_ID[index_offsprg]].chr_b[index_chrom], parent1[parent1_ID[index_offsprg]].chr_a[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec,pos);
////cout << "new gamete done"<< endl;
//                if (mut>0)
//                {
//                    for (j = 0; j < mut; j++)
//                    {
//                        site = int(gsl_rng_uniform(r)*nb_loci);
//                        offspring[index_offsprg].chr_a[index_chrom].sel.flip(site);
//                        //cout << "site " << site<< endl;
//                    }
//                }
//                for (index_dim=0; index_dim<dim; index_dim++)
//                {
//                    //chrom[free].pheno[i]=0;
//                    pheno=0;
//                    for (j = 0; j < nb_loci; j++)
//                    {
//                        if (offspring[index_offsprg].chr_a[index_chrom].sel.test(j))
//                            pheno+= mutations[(index_chrom*nb_loci+j)*dim+index_dim];
//                        //chrom[free].pheno[i] += mutations[j*dim+i];
//                    }
//                    offspring[index_offsprg].chr_a[index_chrom].pheno[index_dim]=pheno;
//                }
//            }
//
//
//
//            if ((co_nb%2!=0)^(gsl_rng_uniform(r)>rho_LB))
//                strand_a_p1=!strand_a_p1;
//
//            co_nb = gsl_ran_poisson(r,rho_LG);
//            mut = gsl_ran_poisson(r,mu_GW); // number of mutations at selected loci on genome segment
//            if (((co_nb==0)|(parent2[parent2_ID[index_offsprg]].chr_a[index_chrom].sel==parent2[parent2_ID[index_offsprg]].chr_b[index_chrom].sel)))
//            {
//                //cout << "no recomb"<<endl;
//                offspring[index_offsprg].chr_b[index_chrom].sel=parent2[parent2_ID[index_offsprg]].chr_a[index_chrom].sel;
//
//                for (index_dim = 0; index_dim < dim; index_dim++)
//                    offspring[index_offsprg].chr_b[index_chrom].pheno[index_dim] = parent2[parent2_ID[index_offsprg]].chr_a[index_chrom].pheno[index_dim];
//                if (mut>0)
//                {
//                    for (j = 0; j < mut; j++)
//                    {
//                        site = int(gsl_rng_uniform(r)*nb_loci);// note that this allows for double mutation at the same site (with probability 1/nb_loci²)
//                        //cout << "site " << site<< endl;
//                        ns = (nb_loci*index_chrom+site)*dim;
//                        // cout <<"site "<<site <<endl;
//                        // cout << chrom[free].sel[site]<< endl;
//                        if (offspring[index_offsprg].chr_b[index_chrom].sel.test(site))
//                        {
//                            for (index_dim = 0; index_dim < dim; index_dim++)
//                                offspring[index_offsprg].chr_b[index_chrom].pheno[index_dim] -= mutations[ns + index_dim];
//                        }
//                        else
//                        {
//                            for (index_dim = 0; index_dim < dim; index_dim++)
//                            {
//                                //    cout << k<< " "<<chrom[free].pheno[k]<< endl;
//                                offspring[index_offsprg].chr_b[index_chrom].pheno[index_dim] += mutations[ns + index_dim];
//                                //      cout <<"done"<<endl;
//                            }
//                            //    cout <<"loop done ";
//                        }
//                        //cout <<"flip "<< site<< endl;
//                        // cout << chrom[free].sel<<endl;
//
//                        offspring[index_offsprg].chr_b[index_chrom].sel.flip(site);
//                        //   cout <<"flip done"<<endl;
//                    }
//                    //cout << "new chr " << chrom[free].sel<<endl;
//                }
//                //    cout <<nb << " "<<free<<"\n";
//                //  cout <<nb << " "<<free<<"\n";
//            }
//            else
//            {
//                //cout << "recomb"<<endl;
//
//                //cout << (*(pop[chr_ID])).nbchr<<endl;
//                //cout << "parent"<< chr_ID << "-"<< chr_ID_alt << " all_chr_ent"<< (*(genotypes[chr_ID])).ID<< " "<< (*(genotypes[(chr_ID_alt)])).ID<<endl;
//                //rec(r, chrom[free], (*(genotypes[chr_ID])), (*(genotypes[(chr_ID_alt)])), co_nb, nb_loci);
//                //cout << parent2[parent2_ID[index_offsprg]].chr_a[i].sel<<" "<< parent2[parent2_ID[index_offsprg]].chr_b[i].sel<<endl;
//                if(strand_a_p2)
//                    rec_v2(r, offspring[index_offsprg].chr_b[index_chrom], parent2[parent2_ID[index_offsprg]].chr_a[index_chrom], parent2[parent2_ID[index_offsprg]].chr_b[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec,pos);
//                else
//                    rec_v2(r, offspring[index_offsprg].chr_b[index_chrom], parent2[parent2_ID[index_offsprg]].chr_b[index_chrom], parent2[parent2_ID[index_offsprg]].chr_a[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec,pos);
////cout << "new gamete done"<< endl;
//
//                if (mut>0)
//                {
//                    for (j = 0; j < mut; j++)
//                    {
//                        site = int(gsl_rng_uniform(r)*nb_loci);
//                        offspring[index_offsprg].chr_b[index_chrom].sel.flip(site);
//                        //cout << "site " << site<< endl;
//                    }
//                }
//                for (index_dim=0; index_dim<dim; index_dim++)
//                {
//                    //chrom[free].pheno[i]=0;
//                    pheno=0;
//                    for (j = 0; j < nb_loci; j++)
//                    {
//                        if (offspring[index_offsprg].chr_b[index_chrom].sel.test(j))
//                            pheno+= mutations[(index_chrom*nb_loci+j)*dim+index_dim];
//                        //chrom[free].pheno[i] += mutations[j*dim+i];
//                    }
//                    offspring[index_offsprg].chr_b[index_chrom].pheno[index_dim]=pheno;
//                }
//            }
//
//
//
//            if ((co_nb%2!=0)^(gsl_rng_uniform(r)>rho_LB))
//                strand_a_p2=!strand_a_p2;
//        }
//
//
//    }
//}
