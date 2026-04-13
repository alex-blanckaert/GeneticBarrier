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

// function to read the mutation file . mutation file must start with *  and end with *. If load_mutation_file is false it will create a new mutation file according to the parameters.
bool lireMutationFile(const char *param_file, gsl_rng* r, const int &nb_LDblocks, const int & nb_loci_per_LG, double * mutations, int * pos_dmi, int *target_chrom, const int &dim, const int& nb_loci, const int &nS, const double &sig, const int &nb_dmi, const int &dmi_neutral, const int &load_mutation_file)
{
    FILE * fichierE;
    fichierE = fopen(param_file,"r");
    int i, x;
    bool term;

    if (fichierE)
    {
        double temp;
        do
        {
            x = fgetc(fichierE);
        }
        while (!((x == '*') || (x == EOF)));
        // each parameter set must start with *
        if (x == EOF)
            term = true;
        else
        {
            term = false;
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            if (x!=dim)
                term=true;
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            if (x!=nS)
                term=true;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            if (temp!=sig)
                term=true;
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            if (x!=nb_dmi)
                term=true;
            for (i=0; i<nS; i++)
            {
                if(fscanf(fichierE,"%lf ",&mutations[i])!=1)
                    term=true;
            }
            if (nb_dmi>0)
            {
                for (i=0; i<(2*nb_dmi); i++)
                {
                    if(fscanf(fichierE,"%d ",&pos_dmi[i])!=1)
                        term=true;
                }
            }
            x = fgetc(fichierE);
            if (x != '*')
                term=true;
        }
        fclose(fichierE);
    }
    else
    {
        term=false;
        if ( load_mutation_file==1)
        {
            cout << "failed to load mutation file despite option load_mutation_file=1";
            return true;
        }
        //fclose(fichierE);
        // create new mutation file
        string nomFichier=param_file;
        ofstream fout;
        fout.open(nomFichier);
        fout << "* "<< dim<< " "<< nS << " "<<sig<<" "<<nb_dmi<<" ";
        for (i = 0; i < nS; i++)
        {
            mutations[i] = gsl_ran_gaussian(r,sig);
        }
        int * seq_nb_loci= new int [nb_loci];
        for (i=0; i<nb_loci; i++)
            seq_nb_loci[i]=i;
        if (nb_dmi>0)
        {
            gsl_ran_choose(r,pos_dmi,2*nb_dmi,seq_nb_loci,nb_loci,sizeof (int));
            gsl_ran_shuffle(r,pos_dmi,2*nb_dmi,sizeof (int));
            if (dmi_neutral==1)
            {
                for (i=0; i<(2*nb_dmi); i++)
                {
                    for (x=0; x<dim; x++)
                    {
                        mutations[pos_dmi[i]*dim+x]=0;
                    }
                }
            }
        }

        for (i = 0; i < nS; i++)
            fout << mutations[i] << " ";

        if (nb_dmi>0)
            for (i=0; i<(2*nb_dmi); i++)
                fout << pos_dmi[i]<< " ";

        fout << "*"<< endl;
        delete [] seq_nb_loci;
    }

    for (i=0; i<2*nb_dmi; i++)
    {
        //     cout << pos_dmi[i]/nb_LDblocks;
        switch(pos_dmi[i]/nb_loci_per_LG)
        {
        case 0:
            target_chrom[i]=0;
            break;
        case 1:
        {
            pos_dmi[i]-=nb_loci_per_LG;
            target_chrom[i]=1;
            break;
        }
        case 2:
        {
            pos_dmi[i]-=2*nb_loci_per_LG;
            target_chrom[i]=2;
            break;
        }
        case 3:
        {
            pos_dmi[i]-=3*nb_loci_per_LG;
            target_chrom[i]=3;
            break;
        }
        case 4:
        {
            pos_dmi[i]-=4*nb_loci_per_LG;
            target_chrom[i]=4;
            break;
        }
        }
        // cout << " "<< pos_dmi[i]<< " "<< target_chrom[i] << endl;
    }
    return term;
}

//function to load an output file and continue from its current state. will check that the old and new parameters are not conflicting (ie migration can change for example but the phenotypic dimensions cannot)
bool reload_population(const char *pop_file, const int &n_ind, const int &nb_loci, const int &dim, const double &mu, const double &mu_mod, const double &mig, const double &rho, const double &rho_LB,const double &rho_CL, const double &sig, const double &amp, const double &epistasis, const vector<double> &opt1, const vector<double> &opt2, const double &var_pref, const double &var_cost, const int &nb_dmi, const double &eps, int &nb_female_1, int &nb_male_1, int &nb_female_2, int &nb_male_2,Ind *pop1_females, Ind *pop1_males, Ind *pop2_females, Ind *pop2_males, const int &nb_LDblocks, const int & nb_loci_per_LB, const int& cycle_mig_sel)
{
    FILE * fichierE;
    fichierE = fopen(pop_file,"r");

    if (!fichierE)
    {
        cout << "Population file does not exist"<< endl;
        return 1;
    }

    unsigned int old_seed;
    int x, i,j,k;
    double temp;
    bool term;
    //char bit_set;

    x = fgetc(fichierE);
    if (x == EOF)
        term = true;
    else
    {
        // read the parameters, notifying the one that have changed and stopping the simulations oif the critical ones (number of loci, dimension, nb dmi, sigma are affected)
        term = false;
        if(fscanf(fichierE,"%*s %d ",&old_seed)!=1)
            term=true;
        if (old_seed<0)
            term=true;
        //cout <<term << " " <<seed<< " "<<fgetc(fichierE) <<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        if (x!=n_ind)
            cout << "change in population size from"<< x <<" to " << n_ind<<endl;
        //cout <<term << " "<<n_ind<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        if (x!=nb_loci)
        {
            cout << "number of loci has changed; it should not!"<<endl;
            term=true;
        }
        //cout <<term << " "<<nb_loci<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        if (x!=dim)
        {
            cout << "number of dimension has changed; it should not"<<endl;
            term=true;
        }
        //cout <<term << " "<<gen_Max<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=mu)
            cout << "mutation rate has changed from"<< temp <<" to " << mu <<endl;
        //cout <<term << " "<<mu<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=mu_mod)
            cout << "mutation rate MC has changed from"<< temp <<" to " << mu_mod <<endl;
        //cout <<term << " "<<mu<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=mig)
            cout << "migration rate has changed from"<< temp <<" to " << mig <<endl;
        //cout <<term << " "<<mig<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=rho)
            cout << "recombination rate has changed from"<< temp <<" to " <<rho <<endl;
        //cout <<term << " "<<rho<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=rho_LB)
            cout << "recombination rate between linkage group has changed from"<< temp <<" to " <<rho_LB <<endl;
        //cout <<term << " "<<rho<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=rho_CL)
            cout << "recombination rate with mate choice loci has changed from"<< temp <<" to " <<rho_CL <<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=sig)
        {
            cout << "mutation effect has changed from"<< temp <<" to " << sig <<endl;
            term=true;
        }
        //cout <<term << " "<<sig<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=amp)
            cout << "strength of selection on phenotype has changed from"<< temp <<" to " << amp <<endl;
        //cout <<term << " "<<amp<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=epistasis)
            cout << "epistasis in FL has changed from"<< temp <<" to " << epistasis  <<endl;
        //cout <<term << " "<<epistasis<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=var_pref)
            cout << "mutation size for preference has changed from"<< temp <<" to " << var_pref <<endl;
        //cout <<term << " "<<var_pref<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=var_cost)
            cout << "mutation size for cost has changed from"<< temp <<" to " << var_cost <<endl;
        //cout <<term << " "<<var_cost<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        if (x!=nb_dmi)
        {
            cout << "number of DMI has changed; it should not" <<endl;
            term=true;
        }
        //cout <<term << " "<<nb_dmi<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=eps)
            cout << "epistasis strength in DMI has changed from"<< temp <<" to " << eps <<endl;
        //cout <<term << " "<<eps<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1) //skip dmi_neutral
            term=true;
        if(fscanf(fichierE,"%*s %d ",&x)!=1) //skip skip_mate_choice
            term=true;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1) // skip initial pref
            term=true;
        //cout <<term << " "<<ini_pref<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)// skip initial cost
            term=true;
        //cout <<term << " "<<ini_cost<<endl;
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1) // skip origin
            term=true;
        //cout <<term << " "<<common_origin<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        if (x!=cycle_mig_sel)
        {
            cout << "life cycle has changed; it should not" <<endl;
            term=true;
        }
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=opt1[0])
            cout <<"optimum 1 position has changed from "<< temp <<" to " << opt1[0] <<endl;
        if (dim>1)
        {
            for (i=1; i<dim; i++)
            {
                if(fscanf(fichierE,"%lf ",&temp)!=1)
                    term=true;
                if (temp!=opt1[i])
                    cout <<"optimum 1 position has changed"<<endl;

            }
        }
        if(fscanf(fichierE,"%*s %lf ",&temp)!=1)
            term=true;
        if (temp!=opt2[0])
            cout <<"optimum 2 position has changed from "<< temp <<" to " << opt2[0] <<endl;
        if (dim>1)
        {
            for (i=1; i<dim; i++)
            {
                if(fscanf(fichierE,"%lf ",&temp)!=1)
                    term=true;
                if (temp!=opt2[i])
                    cout <<"optimum 2 position has changed"<<endl;
            }
        }
        // read the populations
        cout << "loading female population 1"<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        nb_female_1=x;
        for (i=0; i<nb_female_1; i++)
        {
            //loading the "A" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop1_females[i].neutral_a=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_females[i].pref_a=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_females[i].cost_a=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop1_females[i].chr_a[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop1_females[i].chr_a[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
            // loading the "B" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop1_females[i].neutral_b=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_females[i].pref_b=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_females[i].cost_b=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop1_females[i].chr_b[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop1_females[i].chr_b[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
        }

        cout << "loading male population 1"<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        nb_male_1=x;
        for (i=0; i<nb_male_1; i++)
        {
            //loading the "A" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop1_males[i].neutral_a=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_males[i].pref_a=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_males[i].cost_a=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop1_males[i].chr_a[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop1_males[i].chr_a[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
            // loading the "B" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop1_males[i].neutral_b=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_males[i].pref_b=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop1_males[i].cost_b=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop1_males[i].chr_b[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop1_males[i].chr_b[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
        }

        cout << "loading female population 2"<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        nb_female_2=x;
        for (i=0; i<nb_female_2; i++)
        {
            //loading the "A" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop2_females[i].neutral_a=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_females[i].pref_a=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_females[i].cost_a=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop2_females[i].chr_a[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop2_females[i].chr_a[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
            // loading the "B" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop2_females[i].neutral_b=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_females[i].pref_b=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_females[i].cost_b=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop2_females[i].chr_b[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop2_females[i].chr_b[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
        }
        cout << "loading male population 2"<<endl;
        if(fscanf(fichierE,"%*s %d ",&x)!=1)
            term=true;
        nb_male_2=x;
        for (i=0; i<nb_male_2; i++)
        {
            //loading the "A" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop2_males[i].neutral_a=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_males[i].pref_a=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_males[i].cost_a=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop2_males[i].chr_a[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop2_males[i].chr_a[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
            // loading the "B" chromosome
            if(fscanf(fichierE,"%d ",&x)!=1)
                term=true;
            pop2_males[i].neutral_b=x;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_males[i].pref_b=temp;
            if(fscanf(fichierE,"%lf ",&temp)!=1)
                term=true;
            pop2_males[i].cost_b=temp;
            for (j=0; j < nb_LDblocks; j++)
            {
                for (k=0; k< dim; k++)
                {
                    if(fscanf(fichierE,"%lf ",&temp)!=1)
                        term=true;
                    pop2_males[i].chr_b[j].pheno[k]=temp;
                }
                for (k=0; k<nb_loci_per_LB; k++)
                {
                    x = fgetc(fichierE);
                    //   cout <<x-'0'<< " ";
                    if ((x-'0') == 1)
                        pop2_males[i].chr_b[j].sel.flip((nb_loci_per_LB-k-1));
                }
            }
        }

        //cout <<term << " "<<common_origin<<endl;
        cout <<"population loaded"<<endl;
    }
    fclose(fichierE);
    return term;
}

// function to update fitness for the general case- a different function handle the case when half_epistasis =1
void update_fitness(const int &nb_ind,const int &nb_LDblocks, const int &dim, const double &amp, const double &half_epistasis,const vector<double> &opt, Ind *pop, const double &eps, const int &nb_dmi, const int *dmi_pos, const int * target_chrom)
{
    double w_max = 0;
    double d;
    int i, j, k, target, locus;
    if ((amp>0)&(eps>0)) // non-flat fitness landscape and DMIs
    {
        for (i = 0; i < nb_ind; i++)
        {
            d = 0;
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
            }
            //cout <<"calculate pheno and d_sq for ind"<<endl;
            pop[i].abs_fit = exp(-amp * pow(d,half_epistasis));// fitness of individual i:
            d=1;
            for (j=0; j<nb_dmi; j++)
            {
                locus=dmi_pos[2*j];
                target=target_chrom[2*j];
                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
                locus=dmi_pos[2*j+1];
                target=target_chrom[2*j+1];
                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];// k is counting the number of interaction between A and B

                switch(k)
                {
                case 0:// aabb
                    break;
                case 1://AbBb
                    d*=1-eps;
                    break;
                case 2://AABb or AaBB
                    d*=(1-eps)*(1-eps);
                    break;
                case 4://AABB
                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
                    break;
                }
            }
            pop[i].dmi_fit=d;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            if (w_max < (pop[i].abs_fit*pop[i].dmi_fit))
                w_max = (pop[i].abs_fit*pop[i].dmi_fit);
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
        for (i =0; i < nb_ind; i++)
        {
            pop[i].fitness = (pop[i].abs_fit*pop[i].dmi_fit)/w_max;//compute relative fitness
            //cout <<pop[i].fitness<<" ";
        }
    }
    else if ((amp>0)&(eps==0))// no DMIs
    {
        for (i = 0; i < nb_ind; i++)
        {
            d = 0;
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                //   cout <<j<<endl;
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
            }
            //cout <<"calculate pheno and d_sq for ind"<<endl;
            pop[i].abs_fit = exp(-amp * pow(d,half_epistasis));// fitness of individual i:
            pop[i].dmi_fit=1;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            if (w_max < (pop[i].abs_fit))
                w_max = (pop[i].abs_fit);
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
        for (i =0; i < nb_ind; i++)
        {
            pop[i].fitness = (pop[i].abs_fit)/w_max;
            //cout <<pop[i].fitness<<" ";
        }
    }
    else if ((amp==0)&(eps>0))// flat fitness landscape
    {
        for (i = 0; i < nb_ind; i++)
        {
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                //   cout <<j<<endl;
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
            }
            pop[i].abs_fit=1;
            //cout <<"calculate pheno and d_sq for ind"<<endl;
            d=1;
            for (j=0; j<nb_dmi; j++)
            {
                locus=dmi_pos[2*j];
                target=target_chrom[2*j];
                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
                locus=dmi_pos[2*j+1];
                target=target_chrom[2*j+1];
                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];

                switch(k)
                {
                case 0:
                    break;
                case 1:
                    d*=1-eps;
                    break;
                case 2:
                    d*=(1-eps)*(1-eps);
                    break;
                case 4:
                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
                    break;
                }
            }
            pop[i].dmi_fit=d;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            if (w_max < (pop[i].dmi_fit))
                w_max = (pop[i].dmi_fit);
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
        for (i =0; i < nb_ind; i++)
        {
            pop[i].fitness = (pop[i].dmi_fit)/w_max;
            //cout <<pop[i].fitness<<" ";
        }
    }
    else //flat fitness landscape and no DMIs
    {
        for (i = 0; i < nb_ind; i++)
        {
            d = 0;
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                //   cout <<j<<endl;
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
            }
            pop[i].abs_fit=1;
            pop[i].dmi_fit=1;
            pop[i].fitness = 1;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
    }
}

//function to choose parent with mate choice
void choose_parents(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other, const double &dist_opt_sq)
{
    int candidate_choosy, candidate_other, i, j;
    double pheno_dist=0, counter, sigma;
    for (j = 0; j < n_off; j++)
    {
        //cout << j<< endl;
        //cout <<nb_choosy << " " << nb_other<<endl;
        do
        {
            do
            {
                candidate_choosy = int(gsl_rng_uniform(r) * nb_choosy);// choose random individual
            }
            while (gsl_rng_uniform(r) > pop_choosy[candidate_choosy].fitness);   //individual is rejected on its own fitness (pheno +DMI)
            //cout << "choose choosy"<< candidate_c << endl;
            counter=0;
            do
            {
                do
                {
                    candidate_other = int(gsl_rng_uniform(r) * nb_other);
                }
                while ((gsl_rng_uniform(r) > pop_other[candidate_other].fitness));  //individual is rejected on its own fitness (pheno +DMI) - we keep an individuals if we fail the test, ie the number is below fitness
                //if (counter>0)
                //{
                //   cout<< "rejection"<<endl;
                //}
                counter++; // add +1 to the number of encounter
                pheno_dist=0.0; // compute the square of the phenotypic distance
                for (i=0; i<dim; i++)
                {
                    pheno_dist+=(pop_choosy[candidate_choosy].pheno[i]-pop_other[candidate_other].pheno[i])*(pop_choosy[candidate_choosy].pheno[i]-pop_other[candidate_other].pheno[i]);
                }
                if ((pop_choosy[candidate_choosy].pref==0)&(pheno_dist==0))
                {
                    sigma=1;// here sigma is 1 if the distance is 0
                }
                else if ((pop_choosy[candidate_choosy].pref==0)&(pheno_dist>0))
                {
                    sigma=0; // here sigma is 0 if the preference is 0 and the distance not 0;
                }
                else
                {
                    sigma=exp(pheno_dist/dist_opt_sq*log(pop_choosy[candidate_choosy].pref)); // here is the general form. sigma is 0 for really far away phenotype, because the log is negative
                }
                //cout<< "choose other "<< candidate_o<<" counter "<<counter << " cost "<<pop_choosy[candidate_c].cost<< " dist "<<pheno_dist<< " pref " << pop_choosy[candidate_c].pref<< endl;
            }
            while ((gsl_rng_uniform(r) > sigma)&(counter<=pop_choosy[candidate_choosy].cost));  //mating pair is rejected (ie the test is true when the distance is large) if the pheno distance is too large and the number of trial has not reach its limit
        }
        while (counter>pop_choosy[candidate_choosy].cost);   //the choosy individual is rejected if too many trials happened
        //cout <<"mate accepted"<<endl;
        parent_choosy[j] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, in position 2*j and 2j+1 in the pop vector.
        parent_other[j] = candidate_other;
        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
    }
    //cout <<endl;
}

//function to choose parent without mate choice
void choose_parents(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other)
{
    int candidate_choosy, candidate_other, j;
    for (j = 0; j < n_off; j++)
    {
        //cout << j<< endl;
        //cout <<nb_choosy << " " << nb_other<<endl;

        do
        {
            candidate_choosy = int(gsl_rng_uniform(r) * nb_choosy);// choose random individual
        }
        while (gsl_rng_uniform(r) > pop_choosy[candidate_choosy].fitness);   //individual is rejected on its own fitness (pheno +DMI)
        //cout << "choose choosy"<< candidate_c << endl;

        do
        {
            candidate_other = int(gsl_rng_uniform(r) * nb_other);
        }
        while ((gsl_rng_uniform(r) > pop_other[candidate_other].fitness));  //individual is rejected on its own fitness (pheno +DMI)
        //if (counter>0)
        //{
        //   cout<< "rejection"<<endl;;   //the choosy individual is rejected if too many trials happened
        //cout <<"mate accepted"<<endl;
        parent_choosy[j] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, in position 2*j and 2j+1 in the pop vector.
        parent_other[j] = candidate_other;
        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
    }
    //cout <<endl;
}

// function to update fitness if Q=2
void update_fitness_v2(const int &nb_ind,const int &nb_LDblocks, const int &dim, const double &amp, const vector<double> &opt, Ind *pop, const double &eps, const int &nb_dmi, const int *dmi_pos, const int *target_chrom)
{
    double w_max = 0;
    double d=0;
    int i, j, k, target,locus;
    if ((amp>0)&(eps>0)) //no flat landscape and DMI
    {
        for (i = 0; i < nb_ind; i++)// going through all 0 to nb_ind
        {
            d = 0;
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
            }
            //cout <<"calculate pheno and d_sq for ind"<<endl;
            pop[i].abs_fit = exp(-amp * d);// fitness of individual i:
            d=1;
            for (j=0; j<nb_dmi; j++)
            {
                locus=dmi_pos[2*j];
                target=target_chrom[2*j];
                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
                locus=dmi_pos[2*j+1];
                target=target_chrom[2*j+1];
                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];

                switch(k)
                {
                case 0:
                    break;
                case 1:
                    d*=1-eps;
                    break;
                case 2:
                    d*=(1-eps)*(1-eps);
                    break;
                case 4:
                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
                    break;
                }
            }
            pop[i].dmi_fit=d;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            if (w_max < (pop[i].abs_fit*pop[i].dmi_fit))
                w_max = (pop[i].abs_fit*pop[i].dmi_fit);
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
        for (i =0; i < nb_ind; i++)
        {
            pop[i].fitness = (pop[i].abs_fit*pop[i].dmi_fit)/w_max;
            //cout <<pop[i].fitness<<" ";
        }
    }
    else if ((amp>0)&(eps==0))// not flat landscape no DMIs
    {
        for (i = 0; i < nb_ind; i++)
        {
            d = 0;
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                //   cout <<j<<endl;
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
            }
            //cout <<"calculate pheno and d_sq for ind"<<endl;
            pop[i].abs_fit = exp(-amp * d);// fitness of individual i:
            pop[i].dmi_fit=1;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            if (w_max < (pop[i].abs_fit))
                w_max = (pop[i].abs_fit);
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
        for (i =0; i < nb_ind; i++)
        {
            pop[i].fitness = (pop[i].abs_fit)/w_max;
            //cout <<pop[i].fitness<<" ";
        }
    }
    else if ((amp==0)&(eps>0)) //flat landscape with DMIs
    {
        for (i = 0; i < nb_ind; i++)
        {
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                //   cout <<j<<endl;
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
            }
            pop[i].abs_fit=1;
            //cout <<"calculate pheno and d_sq for ind"<<endl;
            d=1;
            for (j=0; j<nb_dmi; j++)
            {
                locus=dmi_pos[2*j];
                target=target_chrom[2*j];
                // cout << "locus " << locus << " target " << target<< endl;
                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
                locus=dmi_pos[2*j+1];
                target=target_chrom[2*j+1];
                // cout << "locus " << locus << " target " << target<< endl;
                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];

                switch(k)
                {
                case 0:
                    break;
                case 1:
                    d*=1-eps;
                    break;
                case 2:
                    d*=(1-eps)*(1-eps);
                    break;
                case 4:
                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
                    break;
                }
            }
            pop[i].dmi_fit=d;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            if (w_max < (pop[i].dmi_fit))
                w_max = (pop[i].dmi_fit);
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
        for (i =0; i < nb_ind; i++)
        {
            pop[i].fitness = (pop[i].dmi_fit)/w_max;
            //cout <<pop[i].fitness<<" ";
        }
    }
    else // flat landscape and no DMIs
    {
        for (i = 0; i < nb_ind; i++)
        {
            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
            //cout <<"calculate pref for ind "<< i<<endl;
            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
            //cout <<"calculate cost for ind"<<endl;
            for (j = 0; j < dim; j++)
            {
                //   cout <<j<<endl;
                pop[i].pheno[j]=0;
                //   cout <<j<<endl;
                for (k=0; k< nb_LDblocks; k++)
                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
            }
            pop[i].abs_fit=1;
            pop[i].dmi_fit=1;
            pop[i].fitness = 1;
            //if(pop[i].dmi_fit<1)
            //{
            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
            //}
            //cout << "calculate fitness"<<endl;
            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
        }
    }
}



// function to pass down freely recombining loci, if the allele is code by an integer
void inherit_FR_loci(gsl_rng* r, int &new_value, const int &chr_1, const int &chr_2)
{
    if (chr_1==chr_2)// if identical, who cares
    {
        new_value=chr_1;
    }
    else // if different 50:50
    {
        if (gsl_rng_uniform(r)>0.5)
            new_value = chr_1;
        else
            new_value = chr_2;
    }
}

// function to pass down freely recombining loci, if the allele is code by a double
void inherit_FR_loci(gsl_rng* r, double &new_value, const double &chr_1, const double &chr_2)
{
    if (chr_1==chr_2)
    {
        new_value=chr_1;
    }
    else
    {
        if (gsl_rng_uniform(r)>0.5)
            new_value = chr_1;
        else
            new_value = chr_2;
    }
}


//function to handle gamete formation for the linkage block. strand_a remembers whether we should starts on the "a" or "b" chromosome of the parent.
// recombination between linkage block is handled at the end of the function
void gamete_linkage_block(gsl_rng* r, const double &rho_LG, const double &rho_LB, const int &dim, const int &nb_loci, const double &mu_GW,  const Ind &parent, const int & index_chrom, Chr &offspring_chrom, double * mutations, bool &strand_a,int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
{
    int index_dim, j, site, ns;
    double pheno;
    int co_nb = gsl_ran_poisson(r,rho_LG);//number of Co per LB block
    if (co_nb>(nb_loci-1)) //check for overflow
        co_nb=(nb_loci-1);
    int mut = gsl_ran_poisson(r,mu_GW); // number of mutations at selected loci on genome segment
    // no CO or same chromosome
    //if (((co_nb==0)|((parent1[parent1_ID[index_offsprg]].chr_a[index_chrom].sel==parent1[parent1_ID[index_offsprg]].chr_b[index_chrom].sel)))
     // first handle cases where we can ignore recombination because it does not happens or the two linkage blocks are identical
    if (((co_nb==0)|((parent.chr_a[index_chrom].sel==parent.chr_b[index_chrom].sel))))
    {
        //cout << "no recomb"<<endl;
        if (strand_a)// determine whether it is the "a" or "b" chromosome to be passed - note that it only matters if co_nb ==0; copy the genotype and phenotype for the relevant linkage block
        {
            offspring_chrom.sel=parent.chr_a[index_chrom].sel;
            for (index_dim = 0; index_dim < dim; index_dim++)
                offspring_chrom.pheno[index_dim] = parent.chr_a[index_chrom].pheno[index_dim];
        }
        else
        {
            offspring_chrom.sel=parent.chr_b[index_chrom].sel;
            for (index_dim = 0; index_dim < dim; index_dim++)
                offspring_chrom.pheno[index_dim] = parent.chr_b[index_chrom].pheno[index_dim];
        }

        // deals with mutations on the chromosome
        if (mut>0)
        {
            for (j = 0; j < mut; j++)
            {
                site = int(gsl_rng_uniform(r)*nb_loci);// note that this allows for double mutation at the same site (with probability 1/nb_loci²)
                //cout << "site " << site<< endl;

                ns = (nb_loci*index_chrom+site)*dim;
                // cout <<"site "<<site <<endl;
                // cout << chrom[free].sel[site]<< endl;
                if (offspring_chrom.sel.test(site))// if the current site is derived, remove the effect of the allele
                {
                    for (index_dim = 0; index_dim < dim; index_dim++)
                        offspring_chrom.pheno[index_dim] -= mutations[ns + index_dim];
                }
                else // if the current site is ancestral, add the effect of the allele
                {
                    for (index_dim = 0; index_dim < dim; index_dim++)
                    {
                        //    cout << k<< " "<<chrom[free].pheno[k]<< endl;
                        offspring_chrom.pheno[index_dim] += mutations[ns + index_dim];
                        //      cout <<"done"<<endl;
                    }
                    //    cout <<"loop done ";
                }
                //cout <<"flip "<< site<< endl;
                // cout << chrom[free].sel<<endl;
                offspring_chrom.sel.flip(site);// mutate the site
                //   cout <<"flip done"<<endl;
            }
            //cout << "new chr " << chrom[free].sel<<endl;
        }
        //    cout <<nb << " "<<free<<"\n";
        //  cout <<nb << " "<<free<<"\n";
    }
    else // deals with recombination between the 2 parental chromosomes
    {
        //cout << "recomb"<<endl;

        //cout << (*(pop[chr_ID])).nbchr<<endl;
        //cout << "parent"<< chr_ID << "-"<< chr_ID_alt << " all_chr_ent"<< (*(genotypes[chr_ID])).ID<< " "<< (*(genotypes[(chr_ID_alt)])).ID<<endl;
        //rec(r, chrom[free], (*(genotypes[chr_ID])), (*(genotypes[(chr_ID_alt)])), co_nb, nb_loci);
        //cout << parent1[parent1_ID[nb]].chr_a[i].sel<<" "<< parent1[parent1_ID[nb]].chr_b[i].sel<<endl;
        //cout <<index_offsprg <<" "<<index_chrom<< " "<<offspring[index_offsprg].chr_a[index_chrom].sel<<endl;
        if(strand_a)// if we start from chromosome "a"
            rec_v2(r, offspring_chrom, parent.chr_a[index_chrom], parent.chr_b[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec, pos);
        else
            rec_v2(r, offspring_chrom, parent.chr_b[index_chrom], parent.chr_a[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec, pos);
//cout << "new gamete done"<< endl;
        // deals with mutation
        if (mut>0)
        {
            for (j = 0; j < mut; j++)
            {
                site = int(gsl_rng_uniform(r)*nb_loci);
                offspring_chrom.sel.flip(site);// mutate the site
                //cout << "site " << site<< endl;
            }
        }
         // recompute the phenotype for the linkage blocks
        for (index_dim=0; index_dim<dim; index_dim++)
        {
            //chrom[free].pheno[i]=0;
            pheno=0;
            for (j = 0; j < nb_loci; j++)
            {
                if (offspring_chrom.sel.test(j))
                    pheno+= mutations[(index_chrom*nb_loci+j)*dim+index_dim];
                //chrom[free].pheno[i] += mutations[j*dim+i];
            }
            offspring_chrom.pheno[index_dim]=pheno;
        }
    }
    // deals with recombination between linkage blocks, strand_a is changed if their is an even number of CO and a recombination event between the blocks or an odd number of CO and no recombination events
    if ((co_nb%2!=0)^(gsl_rng_uniform(r)>rho_LB))// XOR statement
        strand_a=!strand_a;
}

//function that handles recombination between 2 linkage blocks of the parent.
void rec_v2(gsl_rng* r, Chr &res, const Chr &c1, const Chr &c2, const int &nbCo, const int &nb_loci,int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
{
    //vector<int> pos{ 0, nb_loci };
    pos.clear();
    pos.push_back(0);
    pos.push_back(nb_loci);
    int j,k;

    //cout << "1C "<<c1.sel<<endl;
    //cout << "2C "<<c2.sel<<endl;


    //res.sel.clear();
    //cout << "enter rec"<<endl;
    // vector "pos" holds the positions of cross-overs:

    //gsl_ran_shuffle(r,Co_pos,nb_loci,sizeof (int));
    j=0;
    do
    {
        k = int(gsl_rng_uniform(r)*nb_loci);
        if (find(pos.begin(), pos.end(), k) != pos.end())
        {
            pos.push_back(k);
            j++;
        }
    }
    while(j<nbCo);

    //for (j = 0; j < nbCo; j++)
    //{

    //    pos.push_back(Co_pos[j]);
    //cout<<pos[j+2]<<" ";
    //}
    //cout <<endl;

    sort(pos.begin(), pos.end());

    //cout << "nb CO "<<nbCo<<endl;
    // creates recombination mask:
    rec.reset();
    for (j = 0; j < (nbCo+1); j+=2)
    {
        //cout << pos[j] << " ";
        k=pos[j+1]-pos[j];
        rec.flip(pos[j],k);
    }
    //cout <<endl;
    //cout << "re "<<rec<<endl;

    off2 = (c2.sel & rec);// we start by c2, because we want to begin the chromosome by c1
    rec.flip();
    //cout << "fl "<<rec<<endl;
    //cout << "c1 "<<(c1.sel)<<endl;
    //cout << "te "<<(rec)<<endl;
    off1 = (c1.sel & rec);
    //cout <<"ok"<<endl;
    //cout << off1 << endl;
    //cout << off2 <<endl;
    //cout <<rec<<endl;
    //cout <<&off1<<endl;
    //cout <<&off2<<endl;
    //cout <<&rec<<endl;
    //cout << &c1.sel<<endl;
    // cout << &c2.sel<<endl;
    //cout << &res.sel<<endl;
    //cout << &pos<<endl;
    res.sel=c1.sel;
    rec=(off1 | off2);
    //cout <<rec<<endl;
    //cout << (off1 | off2)<< endl;
//   cout << "off "<<res.sel.size()<<endl;
    res.sel = rec;
    // cout <<"redundant"<<endl;
    //res.sel = (off1 | off2);


//cout << res.sel << endl;

//   for (j = 0; j < nbCo; j++)
//	cout<<pos[j+2]<<" ";


    //cout << "FC "<<res.sel<<endl;
    //  cout <<endl;
}

//function to update the population, splitting the juveniles into males and females and replacing the parents.
void update_population(const int &nb_Ind, Ind *pop_male, Ind *pop_female, Ind *temp_pop, int &nb_male, int &nb_female)
{
    nb_male=0;
    nb_female=0;
    for (int i = 0; i < nb_Ind; i++)
    {
        //cout <<i<<"\n";
        //cout <<"chr1 ori "<<(*(temp_pop[i])).neutral << " chr2 ori "<< (*(temp_pop[i+1])).neutral<<endl;
        //   cout <<"pop1 sex "<<(*(temp_1[i])).sex_chr << " "<< (*(temp_1[i+1])).sex_chr<<endl;
        //  cout <<"pop2 sex "<<(*(temp_2[i])).sex_chr << " "<< (*(temp_2[i+1])).sex_chr<<endl;
        if (temp_pop[i].sex)
        {
            //cout<<"pass\n";
            pop_male[nb_male] = temp_pop[i];
            nb_male++;
        }
        else
        {
            //cout<<"failed\n";
            pop_female[nb_female] = temp_pop[i];
            nb_female++;
        }
    }
    //cout << nb_Ind <<" " <<nb_male <<" "<< nb_female <<endl;
}
