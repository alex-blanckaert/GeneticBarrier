#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "all_barrier.h"
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

//using gsl 2.7
//using boost 1.80

using namespace std;

// array are fixed size (compared to vector)


int main(int argc, char *argv[])
{
    //***** Definition of parameters *****
    // check that the program is given 4 or 5 arguments
    if (!((argc ==5) | (argc==6)))
    {
        cout <<"only "<< argc << " arguments given to the program instead of 5 (for new simulations) or 6 (using existing population)"<<endl;
        return 1;
    }

    //Define parameter variables.
    unsigned int seed;
    const int nb_Lblocks=5;
    int n_ind; // population size (diploid)
    int gen_Max;//Max number of generation
    int nb_loci;//Number of locus
    double mu;//Mutation rate
    double mu_mod;//Mutation rate of the modifier locus
    double mig;//Migration rate
    int dim;//Dimension of Fisher landscape
    double sig; //Variance in mutation effect
    double fit_alt;
    double epistasis; //Q parameter in Gros 2009
    double rho;
    double rho_LB;
    double rho_CL;
    double var_pref;
    double var_cost;
    int nb_dmi;
    double eps;
    double ini_pref;
    double ini_cost;
    int common_origin; // 0 is the ancestral genotype without DMI, 1 is with LA no DMI, 2 is DMI without LA, 3 is DMI + LA, 4 is importing file
    int export_metrics;
    int dmi_neutral; // 0 is no, 1 is yes
    int skip_mate_choice;
    int load_mutation_file; // 0 is no, 1 is yes
    int cycle_mig_sel; // 0 is no (ie sel then migration; 1 is yes (mig then sel)
    //const vector<double> opt_1= {1,1,1,1,1}; // Position optimum 1
    //const vector<double> opt_2= {-1,-1,-1,-1,-1}; // Position optimum 2
    int nb_F1=100;

    vector<double> opt_1; // Position optimum 1
    vector<double> opt_2; // Position optimum 2
    vector<int> pos;
    // read parameter file
    if (lireFichier(argv[1], seed, n_ind, nb_loci, dim, gen_Max, mu, mu_mod, mig, rho, rho_LB, rho_CL, sig, fit_alt, epistasis, opt_1, opt_2, var_pref, var_cost, nb_dmi, eps, dmi_neutral, ini_pref, ini_cost, skip_mate_choice, common_origin, export_metrics, cycle_mig_sel, load_mutation_file))
    {
        cout <<"wrong parameter file"<<endl;
        return 1;
    }
    // check that initial preference and initial cost are not out of bounds
    if ((ini_pref<0)|(ini_pref>1))
    {
        cout << "initial prefence must be between 0 and 1" <<endl;
        return 1;
    }
    if (ini_cost<0.5)
    {
        cout << "initial cost must be at least 0.5" <<endl;
        return 1;
    }

    if (nb_loci % 5 !=0)
    {
        cout << "number of loci must be a multiple of 5" <<endl;
        return 1;
    }
    //ignore argument 5 corresponding to a population file if common origin is not 4.
    //if ((argc==6)& (common_origin!=4))
    //{
    //    cout<< " pre-existing population file was given but ignored due to common_origin="<<common_origin<<endl;
    //}
    // initialize RNG
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    //printf ("generator type: %s\n", gsl_rng_name (r));

    // define variable that are not parameters
    if (nb_F1>n_ind)
        nb_F1=n_ind;

    int i;
    const int two_n=2*n_ind;
    //const int four_n=4*n_ind;
    //const int eight_n=8*n_ind;
    const double mu_GW=mu*nb_loci;
    //const double mu_mod=mu*10;
    const double half_epistasis=epistasis/2;
    const double half_ini_pref=ini_pref/2;
    const double half_ini_cost=ini_cost/2;
    const int nS=nb_loci*dim;
    const int two_nb_dmi=2*nb_dmi;
    double dist_opt_sq=0.0;
    int nb_loci_per_LG=nb_loci/nb_Lblocks;
    double rho_LG=rho*(nb_loci_per_LG-1);
    for (i =0; i<dim; i++)
        dist_opt_sq+=(opt_1[i]-opt_2[i])*(opt_1[i]-opt_2[i]);
    double amp=-log(fit_alt)/pow(dist_opt_sq,half_epistasis); // average effect of mutations (see Gros 2009)



    Ind * pop1_males= new Ind [n_ind];// define table  of  N individuals - note that we only need half of it - but since it fluctuates N is a better limit
    Ind * pop1_females= new Ind [n_ind];
    Ind * pop2_males= new Ind [n_ind];
    Ind * pop2_females= new Ind [n_ind];

    Ind * temp_1= new Ind [n_ind]; // placeholder to hold offspring - exactly n_ind
    Ind * temp_2= new Ind [n_ind];

    int * parent_pop1_choosy =new int [n_ind]; // define table to hold ID of the parents (N) of the next generation
    int * parent_pop1_other =new int [n_ind];
    int * parent_pop2_choosy =new int [n_ind];
    int * parent_pop2_other =new int [n_ind];
    // counter of number of males and females in each population and for the number of migrating individuals
    int nb_male_1, nb_female_1, nb_male_2, nb_female_2, gen;
    double mean_fitness, mean_pref, mean_cost,mean_dmi,mean_abs_fit;
    vector<double> mean_pheno(dim,0);
    // "mutations" will hold the effect of the 1 allele at each locus
    // on each phenotypic axis:

    double * mutations = new double [nS];
    int * pos_dmi = new int [two_nb_dmi];
    int * target_chrom= new int [two_nb_dmi];
    int * Co_pos= new int [(nb_loci_per_LG-1)];
    for (i=0; i<(nb_loci_per_LG-1); i++)
        Co_pos[i]=i;

    // if the file given in argument exist, read the file otherwise generate new mutations
    if (lireMutationFile(argv[2],r, nb_Lblocks, nb_loci_per_LG,mutations, pos_dmi, target_chrom, dim, nb_loci, nS, sig, nb_dmi, dmi_neutral,load_mutation_file))
    {
        cout <<"could not read or write mutation file "<<endl;
        return 1;
    }
    //placeholder for recombination function
    boost::dynamic_bitset<> recomb;
    boost::dynamic_bitset<> off1;
    boost::dynamic_bitset<> off2;

    off1.resize(nb_loci_per_LG);
    off2.resize(nb_loci_per_LG);
    recomb.resize(nb_loci_per_LG);

    // define output streams
    string nomFichier(argv[3]);
    ofstream fout;
    fout.open(nomFichier);
    string nomFichierM(argv[4]);
    ofstream foutM;
    foutM.open(nomFichierM);

    //***** Initialization of the simulation *****

    // Determine the number of males and females in each population, trying to be as close as possible to an equal sex ratio.
    if(n_ind%2==0)
    {
        nb_male_1=n_ind/2;
        nb_male_2=n_ind/2;
        nb_female_1=n_ind/2;
        nb_female_2=n_ind/2;
    }
    else
    {
        nb_male_1=(n_ind-1)/2;
        nb_male_2=(n_ind-1)/2;
        nb_female_1=(n_ind+1)/2;
        nb_female_2=(n_ind+1)/2;
    }
    // initialize all populations, making sure they have either the default values or the ones corresponding to the chosen origin
    // first for the 2 temporary vector
    initialize_pop(temp_1,n_ind,nb_Lblocks,nb_loci_per_LG,dim,0,0,0,opt_1,0,nb_dmi,pos_dmi,target_chrom,mutations);
    initialize_pop(temp_2,n_ind,nb_Lblocks,nb_loci_per_LG,dim,0,0,0,opt_1,0,nb_dmi,pos_dmi,target_chrom,mutations);

    // population 1 - DMI are initialized for the first locus of the pair(A)
    initialize_pop(pop1_males,n_ind,nb_Lblocks,nb_loci_per_LG,dim,half_ini_pref,half_ini_cost,common_origin,opt_1,0,nb_dmi,pos_dmi,target_chrom,mutations);
    initialize_pop(pop1_females,n_ind,nb_Lblocks,nb_loci_per_LG,dim,half_ini_pref,half_ini_cost,common_origin,opt_1,0,nb_dmi,pos_dmi,target_chrom,mutations);
    // population 2 - DMI are initialized for the first locus of the pair(B)
    initialize_pop(pop2_males,n_ind,nb_Lblocks,nb_loci_per_LG,dim,half_ini_pref,half_ini_cost,common_origin,opt_2,1,nb_dmi,pos_dmi,target_chrom,mutations);
    initialize_pop(pop2_females,n_ind,nb_Lblocks,nb_loci_per_LG,dim,half_ini_pref,half_ini_cost,common_origin,opt_2,1,nb_dmi,pos_dmi,target_chrom,mutations);
    // each of the table is updated for n_ind entries
    // load a previous population
    if (common_origin==4)
    {
        reload_population(argv[5],n_ind,nb_loci, dim, mu, mu_mod, mig, rho, rho_LB, rho_CL, sig, fit_alt, epistasis, opt_1, opt_2, var_pref, var_cost, nb_dmi, eps, nb_female_1, nb_male_1, nb_female_2, nb_male_2, pop1_females, pop1_males, pop2_females, pop2_males, nb_Lblocks, nb_loci_per_LG,cycle_mig_sel);
        amp=-log(fit_alt)/pow(dist_opt_sq,half_epistasis);
    }

    // compute min and max range for phenotype 0
    mean_pref=0;
    mean_dmi=0;
    mean_abs_fit=0;
    for (i=0; i<nS; i+=dim)
    {
        if (mutations[i]>0)
            mean_abs_fit+=mutations[i];
        else
            mean_pref+=mutations[i];
    }
    // export parameters to output file
    foutM << "seed " << seed << " pop_size " <<n_ind <<" nb_loci "<< nb_loci << " nb_dim "<< dim << " mut_rate "<< mu << " mut_rate_mod "<< mu_mod <<" mig_rate "<<mig<< " rec_rate "<<rho<< " rec_rate_between_LB "<< rho_LB << " rec_rate_for_CL "<< rho_CL << " var_mut "<<sig<< " fit_alt "<< fit_alt<< " epistasis "<<epistasis<< " var_pref "<<var_pref<< " var_cost "<<var_cost << " nb_dmi " << nb_dmi<< " dmi_e "<<eps<< " dmi_neutral "<< dmi_neutral<< " skip_mating_choice "<< skip_mate_choice <<" initial_pref "<< ini_pref <<" initial_cost "<< ini_cost << " common_origin "<< common_origin << " life_cycle_is_mig_sel " << cycle_mig_sel<< " opt_one ";

    fout  << "seed " << seed << " pop_size " <<n_ind <<" nb_loci "<< nb_loci << " nb_dim "<< dim << " mut_rate "<< mu << " mut_rate_mod " << mu_mod <<" mig_rate "<<mig<< " rec_rate "<<rho<< " rec_rate_between_LB "<< rho_LB << " rec_rate_for_CL "<< rho_CL << " var_mut "<<sig<< " fit_alt "<< fit_alt<< " epistasis "<<epistasis<< " var_pref "<<var_pref<< " var_cost "<<var_cost << " nb_dmi " << nb_dmi<< " dmi_e "<<eps<< " dmi_neutral "<< dmi_neutral<< " skip_mating_choice "<< skip_mate_choice <<" initial_pref "<< ini_pref <<" initial_cost "<< ini_cost << " common_origin "<< common_origin << " life_cycle_is_mig_sel " << cycle_mig_sel<<" opt_one ";
    for (i=0; i<dim; i++)
    {
        fout << opt_1[i]<<" ";
        foutM << opt_1[i]<<" ";
    }
    fout <<"opt_two ";
    foutM <<"opt_two ";
    for (i=0; i<dim; i++)
    {
        fout << opt_2[i]<<" ";
        foutM << opt_2[i]<<" ";
    }
    fout<<  endl;
    foutM << "min_pheno "<<mean_pref<< " max_pheno "<< mean_abs_fit <<endl;


    //***** Main loop *****for sel _mig

    //cout <<"begin loop\n";
    for (gen = 0; gen < gen_Max; gen++)
    {
        //cout << nb_male_1 << " "<< nb_male_2  <<" "<< nb_female_1<<" "<<nb_female_2<<endl;
        // Initial step migration between populations at the adult stage
        //cout <<"GEN " << gen << endl;
        if (cycle_mig_sel==1) // migration then selection
        {
            migration(r,pop1_females,pop2_females, temp_1, mig, nb_female_1,nb_female_2);
            //cout << "migration"<<endl;
            migration(r,pop1_males,pop2_males, temp_1, mig, nb_male_1,nb_male_2);
            //cout << "migration"<<endl;

            //cout << nb_male_1 << " "<< nb_male_2  <<" "<< nb_female_1<<" "<<nb_female_2<<endl;
            //cout <<"migration done\n";
            if ((nb_male_1==0)|(nb_male_2==0)|(nb_female_1==0)|(nb_female_2==0))
            {
                cout << "population died out"<<endl;
                break;
            }

            if (epistasis==2)
            {
                update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
                update_fitness_v2(nb_female_2, nb_Lblocks, dim, amp, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_1,   nb_Lblocks, dim, amp, opt_1, pop1_males,  eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_2,   nb_Lblocks, dim, amp, opt_2, pop2_males,  eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_female_1, nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_female_2, nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_1,   nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_males,  eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_2,   nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males,  eps,nb_dmi,pos_dmi,target_chrom);
            }
        }
        else if(cycle_mig_sel==0) // selection then migration
        {
            if (epistasis==2)
            {
                update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
                update_fitness_v2(nb_female_2, nb_Lblocks, dim, amp, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_1,   nb_Lblocks, dim, amp, opt_1, pop1_males,  eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness_v2(nb_male_2,   nb_Lblocks, dim, amp, opt_2, pop2_males,  eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_female_1, nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_female_2, nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_1,   nb_Lblocks, dim, amp, half_epistasis, opt_1, pop1_males,  eps,nb_dmi,pos_dmi,target_chrom);
                update_fitness(nb_male_2,   nb_Lblocks, dim, amp, half_epistasis, opt_2, pop2_males,  eps,nb_dmi,pos_dmi,target_chrom);
            }
            migration(r,pop1_females,pop2_females, temp_1, mig, nb_female_1,nb_female_2);
            //cout << "migration"<<endl;
            migration(r,pop1_males,pop2_males, temp_1, mig, nb_male_1,nb_male_2);
            //cout << "migration"<<endl;

            //cout << nb_male_1 << " "<< nb_male_2  <<" "<< nb_female_1<<" "<<nb_female_2<<endl;
            //cout <<"migration done\n";
            if ((nb_male_1==0)|(nb_male_2==0)|(nb_female_1==0)|(nb_female_2==0))
            {
                cout << "population died out"<<endl;
                break;
            }

        }

        //cout <<"update fitness done\n"

        //export metrics
        if ((gen< two_n & gen%(n_ind/10)==0)|(gen%export_metrics==0))
        {
            cout  <<gen<<"\n";
            //cout <<mig1to2_fem <<" "<< mig2to1_fem <<" "<< mig1to2_mal <<" "<< mig2to1_mal<<endl;
            update_mean_metrics(nb_female_1, dim, mean_fitness,mean_abs_fit,mean_dmi, mean_pheno, mean_pref, mean_cost, pop1_females);
            foutM << gen << "\t" <<double(nb_female_1)/n_ind << "\t" << mean_abs_fit<<"\t"<<mean_dmi << "\t" <<  mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
            for (i=0; i<dim; i++)
            {
                foutM << "\t"<< mean_pheno[i];
            }

            update_mean_metrics(nb_male_1, dim, mean_fitness, mean_abs_fit,mean_dmi, mean_pheno, mean_pref, mean_cost, pop1_males);
            foutM << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t" <<  mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
            for (i=0; i<dim; i++)
            {
                foutM << "\t"<< mean_pheno[i];
            }

            update_mean_metrics(nb_female_2, dim, mean_fitness, mean_abs_fit,mean_dmi, mean_pheno, mean_pref, mean_cost, pop2_females);
            foutM << "\t"<<double(nb_female_2)/n_ind << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t" <<  mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
            for (i=0; i<dim; i++)
            {
                foutM << "\t"<< mean_pheno[i];
            }

            update_mean_metrics(nb_male_2, dim,mean_fitness,mean_abs_fit, mean_dmi, mean_pheno, mean_pref, mean_cost, pop2_males);
            foutM << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t" << mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
            for (i=0; i<dim; i++)
            {
                foutM << "\t"<< mean_pheno[i];
            }

            // create F1s to measure hybrid load, simply by pairing females of pop1 with males of pop 2 and vice versa
            choose_parents(r, dim, parent_pop1_choosy, parent_pop2_other, nb_F1, nb_male_1, nb_female_2, pop1_males, pop2_females);// choose parent population 1
                //cout << nb_male_2 <<" "<<nb_female_2<<endl;
            choose_parents(r, dim, parent_pop2_choosy, parent_pop1_other, nb_F1, nb_male_2, nb_female_1, pop2_males, pop1_females);// choose parent population 2
            new_gamete(r,mutations, nb_F1, mu_GW, mu_mod, parent_pop1_choosy, parent_pop2_other, pop1_males, pop2_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos,off1, off2, recomb, pos);
            new_gamete(r,mutations, nb_F1, mu_GW, mu_mod, parent_pop2_choosy, parent_pop1_other, pop2_males, pop1_females, temp_2, nb_loci_per_LG, dim, nb_Lblocks,rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos, off1, off2, recomb, pos);

            if (epistasis==2)
            {
                update_fitness_v2(nb_F1, nb_Lblocks, dim, amp, opt_1, temp_1,eps,nb_dmi,pos_dmi,target_chrom ); // does not affect pref_a or pref_b
                update_fitness_v2(nb_F1, nb_Lblocks, dim, amp, opt_2, temp_2,eps,nb_dmi,pos_dmi,target_chrom);
            }
            else
            {
                update_fitness(nb_F1, nb_Lblocks, dim, amp, half_epistasis, opt_1, temp_1,eps,nb_dmi,pos_dmi,target_chrom); // does not affect pref_a or pref_b
                update_fitness(nb_F1, nb_Lblocks, dim, amp, half_epistasis, opt_2, temp_2,eps,nb_dmi,pos_dmi,target_chrom);
            }
            mean_fitness=0;
            mean_dmi=0;
            for ( i=0;i<nb_F1; i++){
                mean_fitness+=temp_1[i].dmi_fit;
                mean_dmi+=temp_2[i].dmi_fit;
            }
            mean_fitness/=nb_F1;
            mean_dmi/=nb_F1;

            foutM << "\t" << mean_fitness << "\t" << mean_dmi <<endl;
        } // end of the exporting metrics part

        // sampling parents for the next generation population 1
        //cout << nb_male_1 << " "<< nb_female_1 <<endl;
        if (skip_mate_choice==1)
        {
            choose_parents(r, dim, parent_pop1_choosy, parent_pop1_other, n_ind, nb_male_1, nb_female_1, pop1_males, pop1_females);// choose parent population 1
            //cout << nb_male_2 <<" "<<nb_female_2<<endl;
            choose_parents(r, dim, parent_pop2_choosy, parent_pop2_other, n_ind, nb_male_2, nb_female_2, pop2_males, pop2_females);// choose parent population 2
        }
        else
        {
            choose_parents(r, dim, parent_pop1_choosy, parent_pop1_other, n_ind, nb_male_1, nb_female_1, pop1_males, pop1_females, dist_opt_sq);// choose parent population 1
            //cout << nb_male_2 <<" "<<nb_female_2<<endl;
            choose_parents(r, dim, parent_pop2_choosy, parent_pop2_other, n_ind, nb_male_2, nb_female_2, pop2_males, pop2_females, dist_opt_sq);// choose parent population 2
        }
        //cout <<"choose parents done\n";

        new_gamete(r,mutations, n_ind, mu_GW, mu_mod, parent_pop1_choosy, parent_pop1_other, pop1_males, pop1_females, temp_1, nb_loci_per_LG, dim, nb_Lblocks, rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos,off1, off2, recomb, pos);
        new_gamete(r,mutations, n_ind, mu_GW, mu_mod, parent_pop2_choosy, parent_pop2_other, pop2_males, pop2_females, temp_2, nb_loci_per_LG, dim, nb_Lblocks,rho, rho_LG, rho_LB, rho_CL, var_pref, var_cost, Co_pos, off1, off2, recomb, pos);


        update_population(n_ind, pop1_males, pop1_females, temp_1, nb_male_1, nb_female_1);
        update_population(n_ind, pop2_males, pop2_females, temp_2, nb_male_2, nb_female_2);
        if ((nb_male_1==0)|(nb_male_2==0)|(nb_female_1==0)|(nb_female_2==0))
        {
            cout << "population died out gen:"<<gen <<endl;
            break;
        }

    }

    //***** End simulations *****

    // save pop state
    //cout << "exporting population"<<endl;
    //fout << "seed " << seed << " pop_size " <<n_ind <<" ... "<< endl;

    // update fitness values before exporting
    if (epistasis==2)
    {
        update_fitness_v2(nb_female_1, nb_Lblocks, dim, amp,  opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom );
        update_fitness_v2(nb_female_2, nb_Lblocks,dim,  amp,  opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness_v2(nb_male_1, nb_Lblocks,dim,  amp, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness_v2(nb_male_2, nb_Lblocks,dim,  amp, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
    }
    else
    {
        update_fitness(nb_female_1,nb_Lblocks, dim,  amp, half_epistasis, opt_1, pop1_females,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness(nb_female_2, nb_Lblocks,dim,  amp, half_epistasis, opt_2, pop2_females,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness(nb_male_1, nb_Lblocks,dim,  amp, half_epistasis, opt_1, pop1_males,eps,nb_dmi,pos_dmi,target_chrom);
        update_fitness(nb_male_2,nb_Lblocks, dim,  amp, half_epistasis, opt_2, pop2_males,eps,nb_dmi,pos_dmi,target_chrom);
    }


    // export the population

    fout <<"Female_pop1 "<< nb_female_1 <<endl;
    export_pop(fout, nb_female_1, dim, nb_Lblocks, pop1_females);
    fout <<"Male_pop1 " <<nb_male_1<< endl;
    export_pop(fout, nb_male_1, dim, nb_Lblocks, pop1_males);
    fout <<"Female_pop2 "<< nb_female_2 <<endl;
    export_pop(fout, nb_female_2, dim, nb_Lblocks, pop2_females);
    fout <<"Male_pop2 "<< nb_male_2 <<endl;
    export_pop(fout, nb_male_2, dim, nb_Lblocks, pop2_males);
    fout <<endl;

    fout.close();
//
//    //cout <<mig1to2_fem <<" "<< mig2to1_fem <<" "<< mig1to2_mal <<" "<< mig2to1_mal<<endl;
//    update_mean_metrics(nb_female_1, mean_fitness,mean_abs_fit,mean_dmi, mean_pheno, mean_pref, mean_cost, pop1_females);
//    foutM << gen << "\t" <<double(nb_female_1)/n_ind << "\t" << mean_abs_fit<<"\t"<<mean_dmi << "\t" << mean_pheno << "\t" << mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
//
//    update_mean_metrics(nb_male_1, mean_fitness, mean_abs_fit,mean_dmi, mean_pheno, mean_pref, mean_cost, pop1_males);
//    foutM << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t" << mean_pheno << "\t" << mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
//
//    update_mean_metrics(nb_female_2, mean_fitness, mean_abs_fit,mean_dmi, mean_pheno, mean_pref, mean_cost, pop2_females);
//    foutM << "\t"<<double(nb_female_2)/n_ind << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t" << mean_pheno << "\t" << mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
//
//    update_mean_metrics(nb_male_2, mean_fitness,mean_abs_fit, mean_dmi, mean_pheno, mean_pref, mean_cost, pop2_males);
//    foutM << "\t" << mean_abs_fit <<"\t"<<mean_dmi<< "\t" << mean_pheno << "\t" << mean_pref << "\t" << mean_cost<<"\t"<<mean_fitness;
//
//
//    foutM << "\tNA\tNA\tNA\tNA\t";
//
//    foutM << "\tNA\tNA\tNA\tNA\t";
//
//    foutM  << "\tNA\tNA\tNA\tNA\t";
//
//    foutM << "\tNA\tNA\tNA\tNA\t";
//
//    //mean_cost is used here for the origin
//    more_statistics(nb_female_1,gen1_females,opt_1,amp,half_epistasis,mean_fitness,mean_pheno,nb_loci,dim,mean_cost);
//    foutM << "\t" << mean_fitness<<"\t"<<mean_pheno<< "\t"<<mean_cost;
//    more_statistics(nb_male_1,gen1_males,opt_1,amp,half_epistasis,mean_fitness,mean_pheno,nb_loci,dim,mean_cost);
//    foutM << "\t" << mean_fitness<<"\t"<<mean_pheno<< "\t"<<mean_cost;
//    more_statistics(nb_female_2,gen2_females,opt_2,amp,half_epistasis,mean_fitness,mean_pheno,nb_loci,dim,mean_cost);
//    foutM << "\t" << mean_fitness<<"\t"<<mean_pheno<< "\t"<<mean_cost;
//    more_statistics(nb_male_2,gen2_males,opt_2,amp,half_epistasis,mean_fitness,mean_pheno,nb_loci,dim,mean_cost);
//    foutM << "\t" << mean_fitness<<"\t"<<mean_pheno<< "\t"<<mean_cost;


    foutM.close();
    // beginning measurement phase
    //cout <<"freeing memory"<<endl;
    //cout << eight_n<<endl;


    //cout << "pop2m"<<endl;



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

