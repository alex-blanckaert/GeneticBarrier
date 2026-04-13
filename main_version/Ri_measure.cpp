#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "Ri_measure.h"
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>	/* for dynamic_bitset objects (www.boost.org)*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <array>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>




using namespace std;

// array are fixed size (compared to vector)


int main(int argc, char *argv[])
{
    //***** Definition of parameters *****
    // check that the program is given 4 or 5 arguments
    if (!(argc==5))
    {
        cout <<"only "<< argc << " arguments given to the program instead of 5 "<<endl;
        return 1;
    }

    //Define parameter variables.
    unsigned int seed; // seed
    const int nb_Lblocks=5; // number of likage blocks
    int n_ind; // population size (diploid)
    int gen_Max;//Max number of generation
    int nb_loci;//Number of locus
    double mu;//Mutation rate
    double mu_mod;//Mutation rate modifier locus
    double mig;//Migration rate
    int dim;//Dimension of Fisher landscape
    double sig; //Variance in mutation effect
    double fit_alt; // fitness of the optimal phenotype for the alternative environment
    double epistasis; //Q parameter in Gros 2009
    double rho; // recombination rate between adjacent loci
    double rho_LB; //recombination rate between linkage blocks
    double rho_CL; // recombnation between choice locus and 1 locus of the "main" chromosome
    double var_pref; // variance for the deviation of mutations of the preference trait
    double var_cost; // variance for the deviation of mutations of  the cost trait
    int nb_dmi; // nb of DMIs
    double eps; // epistasis for a DMI
    int nb_trials; // nb of invasion trial
    int n_boot; // nb of individuals used to measure F1 mrtrics
    int skip_mate_choice; // variable to control whether mate choice is skipped (1) or not (0)
    int cycle_mig_sel; // variable to control whether selection happens before (1) or after selection (0).
    int load_mutation_file;// variable to indicate whether a mutational file should be created(0) or loaded (1)
    //const vector<double> opt_1= {1,1,1,1,1}; // Position optimum 1
    //const vector<double> opt_2= {-1,-1,-1,-1,-1}; // Position optimum 2

    vector<double> opt_1; // Position optimum 1
    vector<double> opt_2; // Position optimum 2
    vector<int> pos; // placeover to store crossover position
    vector<double> mean_pheno; // variable to store the mean phenotype

    // read parameter file
    if (lireFichier_RI(argv[1], seed, n_ind, nb_loci, dim, gen_Max, mu, mu_mod, mig, rho, rho_LB, rho_CL, sig, fit_alt, epistasis, opt_1, opt_2, var_pref, var_cost, nb_dmi, eps, skip_mate_choice, nb_trials, n_boot, cycle_mig_sel, load_mutation_file))
    {
        cout <<"wrong parameter file"<<endl;
        return 1;
    }
    // check that initial preference and initial cost are not out of bounds
    if (nb_loci % nb_Lblocks !=0)
    {
        cout << "number of loci must be a multiple of "<<nb_Lblocks <<endl;
        return 1;
    }

    // initalize RNG
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    //printf ("generator type: %s\n", gsl_rng_name (r));

    // define variable that are not parameters
    int i;
    const int two_n=2*n_ind;
    //const int four_n=4*n_ind;
    //const int eight_n=8*n_ind;
    const double mu_GW=mu*nb_loci;
    //const double mu_mod=mu*10;
    const double half_epistasis=epistasis/2;
    const int nS=nb_loci*dim;
    const int two_nb_dmi=2*nb_dmi;
    double dist_opt_sq=0.0;
    int nb_loci_per_LG=nb_loci/nb_Lblocks;
    double rho_LG=rho*(nb_loci_per_LG-1);

    // distance between the two phenotypes that are optimal for each environment
    for (i =0; i<dim; i++)
        dist_opt_sq+=(opt_1[i]-opt_2[i])*(opt_1[i]-opt_2[i]);
    double amp; // average effect of mutations (see Gros 2009)

    Ind * pop1_males= new Ind [n_ind];// define table  of 2*N individuals - males of population 1
    Ind * pop1_females= new Ind [n_ind];
    Ind * pop2_males= new Ind [n_ind];
    Ind * pop2_females= new Ind [n_ind];

    Ind * temp_1= new Ind [n_ind]; //define placeholder table for offspring
    Ind * temp_2= new Ind [n_ind];

    Ind * pop1_males_ini= new Ind [n_ind];// define table  of 2*N individuals - conserve the initial state of the population to reset it for each trial
    Ind * pop1_females_ini= new Ind [n_ind];
    Ind * pop2_males_ini= new Ind [n_ind];
    Ind * pop2_females_ini= new Ind [n_ind];

    int * parent_pop1_choosy =new int [n_ind]; // define table to hold ID of the parents (N) of the next generation
    int * parent_pop1_other =new int [n_ind];
    int * parent_pop2_choosy =new int [n_ind];
    int * parent_pop2_other =new int [n_ind];
    // counter of number of males and females in each population and for the number of migrating individuals
    int nb_male_1, nb_female_1, nb_male_2, nb_female_2, nb_male_1_ini, nb_female_1_ini, nb_male_2_ini, nb_female_2_ini, mig1to2_fem, mig2to1_fem, mig1to2_mal, mig2to1_mal, gen, pop1_mk, pop2_mk, pop1_tf, pop2_tf, pop1_mkL, pop2_mkL, pop1_tfL, pop2_tfL, nb_F1, nb_mating;
    double mean_fitness, mean_abs_fit, mean_dmi, mean_pref, mean_cost;

    // "mutations" will hold the effect of the 1 allele at each locus
    // on each phenotypic axis:

    double * mutations = new double [nS];
    int * pos_dmi = new int [two_nb_dmi]; // vector holding position of the DMI pair with the even number corresponding to A_i and the odd one to B_i
    int * target_chrom= new int [two_nb_dmi];
    int * Co_pos= new int [(nb_loci_per_LG-1)];//vector to hold futur position for crossover
    for (i=0; i<(nb_loci_per_LG-1); i++)
        Co_pos[i]=i;
    // if the file given in argument exist, read the file otherwise generate new mutations
    if (lireMutationFile(argv[2],r, nb_Lblocks, nb_loci_per_LG,mutations, pos_dmi, target_chrom, dim, nb_loci, nS, sig, nb_dmi,0,load_mutation_file))
    {
        cout <<"could not read or write mutation file "<<endl;
        return 1;
    }
    //placehodlers for recombination function - will be used when creating new gametes through recombination
    boost::dynamic_bitset<> recomb;
    boost::dynamic_bitset<> off1;
    boost::dynamic_bitset<> off2;
    // give the rght size to the placeholder
    off1.resize(nb_loci_per_LG);
    off2.resize(nb_loci_per_LG);
    recomb.resize(nb_loci_per_LG);

    // define output stremas
    string nomFichier(argv[4]);
    ofstream fout;
    fout.open(nomFichier);



    //***** Initialization of the simulation *****




    // initialize all individuals
    initialize_pop_anc(temp_1, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(temp_2, n_ind, nb_Lblocks, nb_loci_per_LG, dim);

    initialize_pop_anc(pop1_males, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(pop1_females, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(pop2_males, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(pop2_females, n_ind, nb_Lblocks, nb_loci_per_LG, dim);

    initialize_pop_anc(pop1_males_ini, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(pop1_females_ini, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(pop2_males_ini, n_ind, nb_Lblocks, nb_loci_per_LG, dim);
    initialize_pop_anc(pop2_females_ini, n_ind, nb_Lblocks, nb_loci_per_LG, dim);

    // load the population file into the vector_ini
    if(reload_population(argv[3],n_ind,nb_loci, dim, mu, mu_mod, mig, rho, rho_LB, rho_CL, sig, fit_alt, epistasis, opt_1, opt_2, var_pref, var_cost, nb_dmi, eps, nb_female_1_ini, nb_male_1_ini, nb_female_2_ini, nb_male_2_ini, pop1_females_ini, pop1_males_ini, pop2_females_ini, pop2_males_ini, nb_Lblocks, nb_loci_per_LG, cycle_mig_sel))
    {
        cout <<"could not load population file "<<endl;
        return 1;
    }
    amp=-log(fit_alt)/pow(dist_opt_sq,half_epistasis);

    //update fitness in case the fitness landscape has change
    if (epistasis==2)
    {
        update_fitness_v2(nb_female_1_ini, nb_Lblocks, dim, amp, opt_1, pop1_females_ini,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness_v2(nb_female_2_ini, nb_Lblocks, dim, amp, opt_2, pop2_females_ini,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness_v2(nb_male_1_ini,   nb_Lblocks, dim, amp, opt_1, pop1_males_ini,  eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness_v2(nb_male_2_ini,   nb_Lblocks, dim, amp, opt_2, pop2_males_ini,  eps,nb_dmi,pos_dmi,target_chrom);
    }
    else
    {
        update_fitness(nb_female_1_ini, nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females_ini,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness(nb_female_2_ini, nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_females_ini,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness(nb_male_1_ini,   nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_males_ini,  eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness(nb_male_2_ini,   nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males_ini,  eps,nb_dmi,pos_dmi,target_chrom);
    }



    // export parameters to output file and metrics of the starting population


    ifstream fichierE(argv[3]);
    string line;
    getline(fichierE, line);
    fout << line<<endl;
    fichierE.close();

    if (n_boot>n_ind) // we used the temp variable to hold the fake F1 gnerated - this prevent overflow
        n_boot = n_ind;

    fout  << "seed " << seed << " pop_size " <<n_ind <<" nb_loci "<< nb_loci << " nb_dim "<< dim << " mut_rate "<< mu << " mut_rate_mod " << mu_mod <<" mig_rate "<<mig<< " rec_rate "<<rho<< " rec_rate_between_LB "<< rho_LB << " rec_rate_LC "<< rho_CL << " var_mut "<<sig<< " fit_alt "<< fit_alt<< " epistasis "<<epistasis<< " var_pref "<<var_pref<< " var_cost "<<var_cost << " nb_dmi " << nb_dmi<< " dmi_e "<<eps<< " skip_mate_choice " << skip_mate_choice << " life_cycle_is_mig_sel " << cycle_mig_sel << " opt_one ";
    for (i=0; i<dim; i++)
    {
        fout << opt_1[i]<<" ";
    }
    fout <<"opt_two ";
    for (i=0; i<dim; i++)
    {
        fout << opt_2[i]<<" ";
    }
    fout<<endl;
    // recompute mean metrics of the population
    fout << "nb_trials "<< nb_trials<< " nb_paring_for_F1_metrics "<< n_boot;
    update_mean_metrics(nb_female_1_ini, dim, mean_abs_fit, mean_dmi, mean_fitness, mean_pref, mean_cost, mean_pheno, pop1_females_ini);
    fout << " mean_fit_eco_f1 " << mean_abs_fit << " mean_dmi_f1 " << mean_dmi << " mean_fit_f1 " << mean_fitness<< " mean_pref_f1 " << mean_pref << "  mean_cost_f1 "<< mean_cost;
    for (i =0; i<dim; i++)
    {
        fout << " mean_pheno_f1_t" << i << " " << mean_pheno[i];
    }
    update_mean_metrics(nb_female_2_ini, dim, mean_abs_fit, mean_dmi, mean_fitness, mean_pref, mean_cost, mean_pheno, pop2_females_ini);
    fout << " mean_fit_eco_f2 " << mean_abs_fit << " mean_dmi_f2 " << mean_dmi << " mean_fit_f2 " << mean_fitness<< " mean_pref_f2 " << mean_pref << "  mean_cost_f2 "<< mean_cost;
    for (i =0; i<dim; i++)
    {
        fout << " mean_pheno_f2_t" << i << " " << mean_pheno[i];
    }
    update_mean_metrics(nb_male_1_ini, dim, mean_abs_fit, mean_dmi, mean_fitness, mean_pref, mean_cost, mean_pheno, pop1_males_ini);
    fout << " mean_fit_eco_m1 " << mean_abs_fit << " mean_dmi_m1 " << mean_dmi << " mean_fit_m1 " << mean_fitness<< " mean_pref_m1 " << mean_pref << "  mean_cost_m1 "<< mean_cost;
    for (i =0; i<dim; i++)
    {
        fout << " mean_pheno_m1_t" << i << " " << mean_pheno[i];
    }
    update_mean_metrics(nb_male_2_ini, dim, mean_abs_fit, mean_dmi, mean_fitness, mean_pref, mean_cost, mean_pheno, pop2_males_ini);
    fout << " mean_fit_eco_m2 " << mean_abs_fit << " mean_dmi_m2 " << mean_dmi << " mean_fit_m2 " << mean_fitness<< " mean_pref_m2 " << mean_pref << "  mean_cost_m2 "<< mean_cost;
    for (i =0; i<dim; i++)
    {
        fout << " mean_pheno_m2_t" << i << " " << mean_pheno[i];
    }
    fout << endl;

//***** Main loop for female*****
    for (int rep=0; rep<nb_trials; rep++)
    {
        //cout <<rep<<endl;
        //initialize counters with 2 copy of the marker allele and titme of loss to 0. 0 is a default (and impossible to obtain value). Observing it in the exported file means the marker was never loss not fixed.
        pop1_mk=2;
        pop2_mk=2;
        pop1_tf=0;
        pop2_tf=0;
        pop1_mkL=2;
        pop2_mkL=2;
        pop1_tfL=0;
        pop2_tfL=0;

        //copy the populations from the initial templates using https://cplusplus.com/reference/algorithm/copy/
        copy(pop1_females_ini,pop1_females_ini+n_ind,pop1_females);
        copy(pop2_females_ini,pop2_females_ini+n_ind,pop2_females);
        copy(pop1_males_ini,pop1_males_ini+n_ind,pop1_males);
        copy(pop2_males_ini,pop2_males_ini+n_ind,pop2_males);

        nb_female_1=nb_female_1_ini;
        nb_female_2=nb_female_2_ini;
        nb_male_1=nb_male_1_ini;
        nb_male_2=nb_male_2_ini;

        // switch one individual from pop 1 to pop 2 and vice versa, deciding whether it is male or female and computing F1 statistics

        //choose female migrants
        mig1to2_fem=int(gsl_rng_uniform(r) * nb_female_1_ini);
        mig2to1_fem=int(gsl_rng_uniform(r) * nb_female_2_ini);

        //cout << mig1to2_fem << "\t" << mig2to1_fem << "\n";

        // exchange two focal individual here
        swap_one_individual(pop1_females, pop2_females, temp_1, mig1to2_fem, mig2to1_fem);

        //measure F1 fitness in pop 1 by making n_boot attemped crosses
        if (skip_mate_choice==0)
            choose_parents_test_other(r, dim, parent_pop1_choosy, parent_pop1_other, n_boot, nb_male_1, nb_female_1, pop1_males, pop1_females, dist_opt_sq, mig1to2_fem, nb_F1, nb_mating); //choosy is male - used therefore a special function for the female; one male is picked randomly ans tres to mate mate n time, n being determined by the cost
        else
            choose_parents_test_migrant(r, dim, parent_pop1_other, parent_pop1_choosy, n_boot, nb_female_1, nb_male_1, pop1_females, pop1_males, mig1to2_fem, nb_F1, nb_mating); // mate choice is skipped here

        new_gamete(r, mutations, nb_F1, mu_GW, 0, parent_pop1_choosy, parent_pop1_other, pop1_males, pop1_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, 0, 0, Co_pos,off1, off2, recomb, pos);// generates the gamtes; choosy has to be male
        //upudate fitness of F1 and export their fitness
        if (epistasis==2)
            update_fitness_v2(nb_F1, nb_Lblocks, dim, amp, opt_1, temp_1, eps, nb_dmi, pos_dmi,target_chrom);
        else
            update_fitness(nb_F1, nb_Lblocks, dim, amp, half_epistasis, opt_1, temp_1, eps, nb_dmi, pos_dmi, target_chrom);

        update_mean_metrics(nb_F1, dim, mean_abs_fit, mean_dmi, mean_fitness, temp_1); // this is only fitness of the F1 cross
        fout << "f1m\t";
            for (i =0; i<dim; i++)
        {
            fout << pop1_females[mig1to2_fem].pheno[i] <<"\t";
        }
        fout << pop1_females[mig1to2_fem].abs_fit << "\t" << pop1_females[mig1to2_fem].dmi_fit << "\t" << pop1_females[mig1to2_fem].pref <<"\t" << nb_F1 << "\t" << nb_mating << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t"<<mean_fitness<<"\t";

        //measure F1 fitness in pop 2 by making n_boot attemped crosses - same as above except all 1 indices are replaced by 2
        if (skip_mate_choice==0)
            choose_parents_test_other(r, dim, parent_pop2_choosy, parent_pop2_other, n_boot, nb_male_2, nb_female_2, pop2_males, pop2_females, dist_opt_sq, mig2to1_fem, nb_F1, nb_mating);
        else
            choose_parents_test_migrant(r, dim, parent_pop2_other, parent_pop2_choosy, n_boot, nb_female_2, nb_male_2, pop2_females, pop2_males, mig2to1_fem, nb_F1, nb_mating);

        new_gamete(r, mutations, nb_F1, mu_GW, 0, parent_pop2_choosy, parent_pop2_other, pop2_males, pop2_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, 0, 0, Co_pos, off1, off2, recomb, pos);
        if (epistasis==2)
            update_fitness_v2(nb_F1, nb_Lblocks, dim, amp, opt_2, temp_1, eps,nb_dmi,pos_dmi,target_chrom);
        else
            update_fitness(nb_F1, nb_Lblocks, dim, amp, half_epistasis, opt_2, temp_1, eps, nb_dmi, pos_dmi, target_chrom);

        update_mean_metrics(nb_F1, dim, mean_abs_fit, mean_dmi, mean_fitness, temp_1); // this is only fitness of the F1 cross
        fout << "f2m\t";
        for (i =0; i<dim; i++)
        {
            fout << pop2_females[mig2to1_fem].pheno[i] <<"\t";
        }
        fout << pop2_females[mig2to1_fem].abs_fit << "\t" << pop2_females[mig2to1_fem].dmi_fit << "\t" << pop2_females[mig2to1_fem].pref <<"\t" << nb_F1 << "\t" << nb_mating << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t"<<mean_fitness<<"\t";



        // cout <<"begin loop\n";
        // begin the actual population cycle here
        gen=0;
        if(cycle_mig_sel==1) // update fitness again only if we have a migration selection life cycle
        {
            if (epistasis==2)
            {
                update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
                update_fitness_v2(nb_female_2, nb_Lblocks,dim, amp, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_1, nb_Lblocks,dim, amp, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_2, nb_Lblocks,dim, amp, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_female_1,nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_female_2, nb_Lblocks,dim, amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_1, nb_Lblocks,dim, amp, half_epistasis, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_2,nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }
        }
        while((gen<gen_Max)&((pop1_tf==0)|(pop2_tf==0)|(pop1_tfL==0)|(pop2_tfL==0)))// wait to reach gen_max or to have all markers lost or fixed
        {
            //cout << nb_male_1 << " "<< nb_male_2  <<" "<< nb_female_1<<" "<<nb_female_2<<endl;

            //cout << gen<<endl;
            gen++;


            // cout << nb_male_1-mig2to1_mal << " "<< nb_female_1-mig2to1_fem <<endl;
            // choose parents for next generation
            if (skip_mate_choice==0) // with matce choice
            {
                //    cout << "choose parents"<<endl;
                choose_parents(r, dim, parent_pop1_choosy, parent_pop1_other, n_ind, nb_male_1, nb_female_1, pop1_males, pop1_females, dist_opt_sq);// choose parent population 1
                //cout << nb_male_2 <<" "<<nb_female_2<<endl;
                choose_parents(r, dim, parent_pop2_choosy, parent_pop2_other, n_ind, nb_male_2, nb_female_2, pop2_males, pop2_females, dist_opt_sq);// choose parent population 2
                //cout <<"choose parents done\n";// choose parent population 2
            }
            else //without mate choice
            {
                //  cout << "random mating"<<endl;
                choose_parents(r, dim, parent_pop1_choosy, parent_pop1_other, n_ind, nb_male_1, nb_female_1, pop1_males, pop1_females);// choose parent population 1
                choose_parents(r, dim, parent_pop2_choosy, parent_pop2_other, n_ind, nb_male_2, nb_female_2, pop2_males, pop2_females);// choose parent population 2
            }
            //cout <<"choose parents done\n";


            // generate the new zygotes for each populations
            new_gamete(r,mutations, n_ind, mu_GW, mu_mod, parent_pop1_choosy, parent_pop1_other, pop1_males, pop1_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos, off1, off2, recomb, pos);
            new_gamete(r,mutations, n_ind, mu_GW, mu_mod, parent_pop2_choosy, parent_pop2_other, pop2_males, pop2_females, temp_2, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos, off1, off2, recomb, pos);
            //cout <<"mutation done\n";

            // update population (zygotes are replacing their parents) to give the n+1 juveniles:


            update_population(n_ind, pop1_males, pop1_females, temp_1, nb_male_1, nb_female_1);
            update_population(n_ind, pop2_males, pop2_females, temp_2, nb_male_2, nb_female_2);
            //cout <<"update population done\n";


            // check status of marker alleles. if allele is lost/fixed for the first update the variable pop* with the current generation number
            if (pop1_tf==0)
            {
                track_neutral(pop1_males, pop1_females, nb_male_1, nb_female_1, pop1_mk);
                if((pop1_mk==0)|(pop1_mk==n_ind))
                    pop1_tf=gen;
            }

            if (pop1_tfL==0)
            {
                track_linked_neutral(pop1_males, pop1_females, nb_male_1, nb_female_1, pop1_mkL);
                if(((pop1_mkL==0)|(pop1_mkL==n_ind)))
                    pop1_tfL=gen;
            }

            if (pop2_tf==0)
            {
                track_neutral(pop2_males, pop2_females, nb_male_2, nb_female_2, pop2_mk);
                if(((pop2_mk==0)|(pop2_mk==n_ind)))
                    pop2_tf=gen;
            }
            if (pop2_tfL==0)
            {
                track_linked_neutral(pop2_males, pop2_females, nb_male_2, nb_female_2, pop2_mkL);
                if(((pop2_mkL==0)|(pop2_mkL==n_ind)))
                    pop2_tfL=gen;
            }
            if ((nb_male_1==0)|(nb_male_2==0)|(nb_female_1==0)|(nb_female_2==0))
            {
                cout << "population died out"<<endl;
                break;
            }
            // update fitness of the whole population
            if (epistasis==2)
            {
                update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
                update_fitness_v2(nb_female_2, nb_Lblocks,dim, amp, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_1, nb_Lblocks,dim, amp, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_2, nb_Lblocks,dim, amp, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_female_1,nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_female_2, nb_Lblocks,dim, amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_1, nb_Lblocks,dim, amp, half_epistasis, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_2,nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }

            //cout <<pop1_mk<< " "<<pop2_mk<<endl;

        }


//***** End simulations for female *****

// save pop state
//cout << "exporting population"<<endl;
        fout << pop1_mkL << " " << pop1_tfL << " " << pop1_mk << " " << pop1_tf << " " << pop2_mkL <<" " << pop2_tfL << " " << pop2_mk << " "<< pop2_tf << endl;
    }

    //***** Main loop for male*****
    for (int rep=0; rep<nb_trials; rep++)
    {
        //cout <<rep<<endl;
        //initialize counters with 2 copy of the marker allele and titme of loss to 0. 0 is a default (and impossible to obtain value). Observing it in the exported file means the marker was never loss not fixed.
        pop1_mk=2;
        pop2_mk=2;
        pop1_tf=0;
        pop2_tf=0;
        pop1_mkL=2;
        pop2_mkL=2;
        pop1_tfL=0;
        pop2_tfL=0;

        //copy the populations from the initial templates using https://cplusplus.com/reference/algorithm/copy/
        copy(pop1_females_ini,pop1_females_ini+n_ind,pop1_females);
        copy(pop2_females_ini,pop2_females_ini+n_ind,pop2_females);
        copy(pop1_males_ini,pop1_males_ini+n_ind,pop1_males);
        copy(pop2_males_ini,pop2_males_ini+n_ind,pop2_males);

        nb_female_1=nb_female_1_ini;
        nb_female_2=nb_female_2_ini;
        nb_male_1=nb_male_1_ini;
        nb_male_2=nb_male_2_ini;

        //choose male migrants
        mig1to2_mal=int(gsl_rng_uniform(r) * nb_male_1);
        mig2to1_mal=int(gsl_rng_uniform(r) * nb_male_2);
        //migration happens here
        swap_one_individual(pop1_males, pop2_males, temp_1, mig1to2_mal, mig2to1_mal);
        //measure F1 fitness in pop 1
        if (skip_mate_choice==0)
            choose_parents_test_choosy(r, dim, parent_pop1_choosy, parent_pop1_other, n_boot, nb_male_1, nb_female_1, pop1_males, pop1_females, dist_opt_sq, mig1to2_mal, nb_F1, nb_mating);
        else
            choose_parents_test_migrant(r, dim, parent_pop1_choosy, parent_pop1_other, n_boot, nb_male_1, nb_female_1, pop1_males, pop1_females, mig1to2_mal, nb_F1, nb_mating);

        new_gamete(r, mutations, nb_F1, mu_GW, 0, parent_pop1_choosy, parent_pop1_other, pop1_males, pop1_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, 0, 0, Co_pos, off1, off2, recomb, pos);
        if (epistasis==2)
            update_fitness_v2(nb_F1, nb_Lblocks, dim, amp, opt_1, temp_1, eps,nb_dmi,pos_dmi,target_chrom);
        else
            update_fitness(nb_F1, nb_Lblocks, dim, amp, half_epistasis, opt_1, temp_1, eps, nb_dmi, pos_dmi, target_chrom);
        update_mean_metrics(nb_F1, dim, mean_abs_fit, mean_dmi, mean_fitness, temp_1); // this is only fitness of the F1 cross
        fout << "m1m\t";
        for (i =0; i<dim; i++)
        {
            fout << pop1_males[mig1to2_mal].pheno[i] <<"\t";
        }
        fout << pop1_males[mig1to2_mal].abs_fit << "\t" << pop1_males[mig1to2_mal].dmi_fit << "\t" << pop1_males[mig1to2_mal].pref <<"\t" << nb_F1 << "\t" << nb_mating << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t"<<mean_fitness<<"\t";

        //measure F1 fitness in pop 2
        if (skip_mate_choice==0)
            choose_parents_test_choosy(r, dim, parent_pop2_choosy, parent_pop2_other, n_boot, nb_male_2, nb_female_2, pop2_males, pop2_females, dist_opt_sq, mig2to1_mal, nb_F1, nb_mating);
        else
            choose_parents_test_migrant(r, dim, parent_pop2_choosy, parent_pop2_other, n_boot, nb_male_2, nb_female_2, pop2_males, pop2_females, mig2to1_mal, nb_F1, nb_mating);

        new_gamete(r, mutations, nb_F1, mu_GW, 0, parent_pop2_choosy, parent_pop2_other, pop2_males, pop2_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, 0, 0, Co_pos, off1, off2, recomb, pos);
        if (epistasis==2)
            update_fitness_v2(nb_F1, nb_Lblocks, dim, amp, opt_2, temp_1, eps,nb_dmi,pos_dmi,target_chrom);
        else
            update_fitness(nb_F1, nb_Lblocks, dim, amp, half_epistasis, opt_2, temp_1, eps, nb_dmi, pos_dmi, target_chrom);
        update_mean_metrics(nb_F1, dim, mean_abs_fit, mean_dmi, mean_fitness, temp_1); // this is only fitness of the F1 cross
        fout << "m2m\t";
        for (i =0; i<dim; i++)
        {
            fout << pop2_males[mig2to1_mal].pheno[i] <<"\t";
        }
        fout << pop2_males[mig2to1_mal].abs_fit << "\t" << pop2_males[mig2to1_mal].dmi_fit << "\t" << pop2_males[mig2to1_mal].pref <<"\t" << nb_F1 << "\t" << nb_mating << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t"<<mean_fitness<<"\t";


        // cout <<"begin loop\n";
        // begin the actual population cycle here
        gen=0;
        if(cycle_mig_sel==1) // update fitness again only if we have a migration selection life cycle
        {
            if (epistasis==2)
            {
                update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
                update_fitness_v2(nb_female_2, nb_Lblocks,dim, amp, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_1, nb_Lblocks,dim, amp, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_2, nb_Lblocks,dim, amp, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_female_1,nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_female_2, nb_Lblocks,dim, amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_1, nb_Lblocks,dim, amp, half_epistasis, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_2,nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }
        }
        while((gen<gen_Max)&((pop1_tf==0)|(pop2_tf==0)|(pop1_tfL==0)|(pop2_tfL==0)))// wait to reach gen_max or to have all markers lost or fixed
        {
            //cout << nb_male_1 << " "<< nb_male_2  <<" "<< nb_female_1<<" "<<nb_female_2<<endl;

            //cout << gen<<endl;
            gen++;


            // cout << nb_male_1-mig2to1_mal << " "<< nb_female_1-mig2to1_fem <<endl;
            // choose parents for next generation
            if (skip_mate_choice==0) // with matce choice
            {
                //    cout << "choose parents"<<endl;
                choose_parents(r, dim, parent_pop1_choosy, parent_pop1_other, n_ind, nb_male_1, nb_female_1, pop1_males, pop1_females, dist_opt_sq);// choose parent population 1
                //cout << nb_male_2 <<" "<<nb_female_2<<endl;
                choose_parents(r, dim, parent_pop2_choosy, parent_pop2_other, n_ind, nb_male_2, nb_female_2, pop2_males, pop2_females, dist_opt_sq);// choose parent population 2
                //cout <<"choose parents done\n";// choose parent population 2
            }
            else //without mate choice
            {
                //  cout << "random mating"<<endl;
                choose_parents(r, dim, parent_pop1_choosy, parent_pop1_other, n_ind, nb_male_1, nb_female_1, pop1_males, pop1_females);// choose parent population 1
                choose_parents(r, dim, parent_pop2_choosy, parent_pop2_other, n_ind, nb_male_2, nb_female_2, pop2_males, pop2_females);// choose parent population 2
            }
            //cout <<"choose parents done\n";


            // generate the new zygotes for each populations
            new_gamete(r,mutations, n_ind, mu_GW, mu_mod, parent_pop1_choosy, parent_pop1_other, pop1_males, pop1_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos, off1, off2, recomb, pos);
            new_gamete(r,mutations, n_ind, mu_GW, mu_mod, parent_pop2_choosy, parent_pop2_other, pop2_males, pop2_females, temp_2, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos, off1, off2, recomb, pos);
            //cout <<"mutation done\n";

            // update population (zygotes are replacing their parents) to give the n+1 juveniles:


            update_population(n_ind, pop1_males, pop1_females, temp_1, nb_male_1, nb_female_1);
            update_population(n_ind, pop2_males, pop2_females, temp_2, nb_male_2, nb_female_2);
            //cout <<"update population done\n";


            // check status of marker alleles. if allele is lost/fixed for the first update the variable pop* with the current generation number
            if (pop1_tf==0)
            {
                track_neutral(pop1_males, pop1_females, nb_male_1, nb_female_1, pop1_mk);
                if((pop1_mk==0)|(pop1_mk==n_ind))
                    pop1_tf=gen;
            }

            if (pop1_tfL==0)
            {
                track_linked_neutral(pop1_males, pop1_females, nb_male_1, nb_female_1, pop1_mkL);
                if(((pop1_mkL==0)|(pop1_mkL==n_ind)))
                    pop1_tfL=gen;
            }

            if (pop2_tf==0)
            {
                track_neutral(pop2_males, pop2_females, nb_male_2, nb_female_2, pop2_mk);
                if(((pop2_mk==0)|(pop2_mk==n_ind)))
                    pop2_tf=gen;
            }
            if (pop2_tfL==0)
            {
                track_linked_neutral(pop2_males, pop2_females, nb_male_2, nb_female_2, pop2_mkL);
                if(((pop2_mkL==0)|(pop2_mkL==n_ind)))
                    pop2_tfL=gen;
            }
            if ((nb_male_1==0)|(nb_male_2==0)|(nb_female_1==0)|(nb_female_2==0))
            {
                cout << "population died out"<<endl;
                break;
            }
            // update fitness of the whole population
            if (epistasis==2)
            {
                update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
                update_fitness_v2(nb_female_2, nb_Lblocks,dim, amp, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_1, nb_Lblocks,dim, amp, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_2, nb_Lblocks,dim, amp, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_female_1,nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_female_2, nb_Lblocks,dim, amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_1, nb_Lblocks,dim, amp, half_epistasis, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_2,nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
            }

            //cout <<pop1_mk<< " "<<pop2_mk<<endl;

        }


//***** End simulations for males *****

// save pop state
//cout << "exporting population"<<endl;
        fout << pop1_mkL << " " << pop1_tfL << " " << pop1_mk << " " << pop1_tf << " " << pop2_mkL <<" " << pop2_tfL << " " << pop2_mk << " "<< pop2_tf << endl;
    }
//
    fout.close();


// beginning measurement phase
//cout <<"freeing memory"<<endl;
//cout << eight_n<<endl;

    //cout <<"mutations"<<endl;
    delete [] parent_pop1_choosy;
    //cout <<"parent_pop1_choosy"<<endl;
    delete [] parent_pop1_other;
    //cout <<"parent_pop1_other"<<endl;
    delete [] parent_pop2_choosy;
    //cout <<"parent_pop2_choosy"<<endl;
    delete [] parent_pop2_other;
    //cout <<"parent_pop2_other"<<endl;
    //cout << "chrom"<<endl;
    //cout << "pop2f"<<endl;
    delete [] pop1_females;
    //cout << "pop1f"<<endl;
    delete [] pop1_males;
    //cout << "pop1m"<<endl;
    delete [] pop2_females;
    //cout << "pop2f"<<endl;
    delete [] pop2_males;
    //cout << "pop2m"<<endl;
    delete [] mutations;
    //cout << "mutations"<<endl;
    delete [] pos_dmi;
    //cout << "pos_dmi"<<endl;
    delete [] target_chrom;
    //cout << "target"<<endl;
    delete [] Co_pos;
    //cout << "co_pos"<<endl;
    delete [] temp_1;
    //cout <<"temp1"<< endl;
    delete [] temp_2;
    //cout <<"temp2"<< endl;
    //cout <<"SUCCESS";
    return 0;
}

