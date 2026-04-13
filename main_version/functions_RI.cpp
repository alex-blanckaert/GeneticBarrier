#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "Ri_measure.h"
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

// this function exchange one individuals between the two populations
void swap_one_individual(Ind *pop1, Ind *pop2, Ind *temp1, const int &nb_mig_1to2, const int &nb_mig_2to1)
{
    temp1[0]=pop1[nb_mig_1to2];// copy migrant from pop1
    temp1[1]=pop2[nb_mig_2to1];// copy migrant from pop2

    //cout << "move to temp"<<endl;
    // give both migrants unique neutral markers
    temp1[0].neutral_a=2;
    temp1[0].neutral_b=2;

    temp1[1].neutral_a=3;
    temp1[1].neutral_b=3;

    temp1[0].neutral_La=2;
    temp1[0].neutral_Lb=2;

    temp1[1].neutral_La=3;
    temp1[1].neutral_Lb=3;

    //cout << "move back to pop"<<endl;
    // replace entry of the emigrating individual with the emigrating one.
    pop1[nb_mig_1to2]=temp1[1];
    pop2[nb_mig_2to1]=temp1[0];
}

// this function initialize the table of Ind object
void initialize_pop_anc(Ind * pop,const int &n_ind, const int &nb_Lblocks, const int & nb_loci_per_LG, const int &dim)
{
    int index_ind, index_chrom, index_dim;
    for (index_ind=0; index_ind<n_ind; index_ind++)
    {
        for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)
        {
            pop[index_ind].chr_a[index_chrom].sel.resize(nb_loci_per_LG);// give the right size to each linkage block
            pop[index_ind].chr_b[index_chrom].sel.resize(nb_loci_per_LG);
            pop[index_ind].chr_a[index_chrom].pheno.resize(dim);// give the right size to the phenotype vector for the corresponding linkage block
            pop[index_ind].chr_b[index_chrom].pheno.resize(dim);

        }
        pop[index_ind].pheno.resize(dim);// give the right size to the phenotype vector for the individual
        pop[index_ind].pref_a=0;// initialize preference
        pop[index_ind].pref_b=0;
        pop[index_ind].cost_a=0;// initialize cost
        pop[index_ind].cost_b=0;
    }
    //ancestral state, no LA, no DMI
    // we run it for every case, to make sure everything is initialized at 0
    for (index_ind=0; index_ind<n_ind; index_ind++)
    {
        for (index_dim=0; index_dim<dim; index_dim++)
        {
            for (index_chrom=0; index_chrom<nb_Lblocks; index_chrom++)
            {
                pop[index_ind].chr_a[index_chrom].pheno[index_dim]=0;
                pop[index_ind].chr_b[index_chrom].pheno[index_dim]=0;
            }
            pop[index_ind].pheno[index_dim]=0;
        }
    }
}

//function to read parameter file
bool lireFichier_RI(const char *param_file, unsigned int &seed, int &n_ind, int &nb_loci, int &dim, int &gen_Max, double &mu, double &mu_mod, double &mig, double &rho, double &rho_LB, double &rho_CL, double &sig, double &amp, double &epistasis, vector<double> &opt1, vector<double> &opt2, double &var_pref, double &var_cost, int &nb_dmi, double &eps, int & skip_mating, int & nb_trials, int & n_boot,int &cycle_mig_sel, int &load_mutation_file)
{

    FILE * fichierE;
    fichierE = fopen(param_file,"r");

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
        //cout <<term << " "<<rho<<endl;
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
        //cout <<term << " "<<ini_cost<<endl;
        //cout <<term << " "<<common_origin<<endl;
        if(fscanf(fichierE,"%*s %d ",&skip_mating)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %d ",&nb_trials)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %d ",&n_boot)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %d ",&cycle_mig_sel)!=1)
            term=true;
        if(fscanf(fichierE,"%*s %d ",&load_mutation_file)!=1)
            term=true;
        x = fgetc(fichierE);
        if (x != '*')
            term=true;
    }
    fclose(fichierE);
    return term;
}

// function to track the number of marker in a population for the unlinked one
void track_neutral(Ind *pop_m, Ind *pop_f, const int & nb_male, const int & nb_female, int &mk)
{
    int i;
    mk=0;
    for (i=0; i<nb_male; i++)//loop over males
    {
        if(pop_m[i].neutral_a>1)
            mk++;
        if(pop_m[i].neutral_b>1)
            mk++;
    }
    for (i=0; i<nb_female; i++)//loop over females
    {
        if(pop_f[i].neutral_a>1)
            mk++;
        if(pop_f[i].neutral_b>1)
            mk++;
    }
}

// function to track the number of marker in a population for the linked one - same as track_neutral() function
void track_linked_neutral(Ind *pop_m, Ind *pop_f, const int & nb_male, const int & nb_female, int &mk)
{
    int i;
    mk=0;
    for (i=0; i<nb_male; i++)
    {
        if(pop_m[i].neutral_La>1)
            mk++;
        if(pop_m[i].neutral_Lb>1)
            mk++;
    }
    for (i=0; i<nb_female; i++)
    {
        if(pop_f[i].neutral_La>1)
            mk++;
        if(pop_f[i].neutral_Lb>1)
            mk++;
    }
}

// function to create zygotes from parents already chosen - this function covers everything from segregation, recombination and mutation
// genetic material from parent1_ID (father) will always be store in the "a" chromosome, and the one from the mother in "b"
void new_gamete(gsl_rng* r, double * mutations, const int &n_ind,const double &mu_GW, const double &mu_mod, const int *parent1_ID, const int *parent2_ID, Ind *parent1, Ind *parent2, Ind *offspring, const int &nb_loci, const int &dim, const int &nb_LD_blocks,const double &rho, const double &rho_LG, const double &rho_LB, const double &rho_CL, const double & var_pref, const double & var_cost, int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
{

    int  index_offsprg, index_chrom;
    double mut_pref, mut_cost, has_rec;
    bool strand_a_p1=false, strand_a_p2=false;// this two variable keeps tracks on whether we are on the "a" or "b" chromosome of the current parent (p1 or p2)

    for (index_offsprg = 0; index_offsprg < n_ind; index_offsprg++) // loop over all offspring chromosome inherited from a given parent;
    {
        //cout <<"ind_nb "<<index_offsprg<<endl;

        // determine sex of the offspring
        if (gsl_rng_uniform(r)>0.5)
            offspring[index_offsprg].sex=1;
        else
            offspring[index_offsprg].sex=0;

        // inherit the freely recombining loci
        inherit_FR_loci(r,offspring[index_offsprg].cost_a,parent1[parent1_ID[index_offsprg]].cost_a,parent1[parent1_ID[index_offsprg]].cost_b);
        inherit_FR_loci(r,offspring[index_offsprg].cost_b,parent2[parent2_ID[index_offsprg]].cost_a,parent2[parent2_ID[index_offsprg]].cost_b);
        inherit_FR_loci(r,offspring[index_offsprg].neutral_a,parent1[parent1_ID[index_offsprg]].neutral_a,parent1[parent1_ID[index_offsprg]].neutral_b);
        inherit_FR_loci(r,offspring[index_offsprg].neutral_b,parent2[parent2_ID[index_offsprg]].neutral_a,parent2[parent2_ID[index_offsprg]].neutral_b);


        // decide whether the transmitted chromosome starts from the "a" or "b" parental one.
        if (gsl_rng_uniform(r)>0.5)
        {
            strand_a_p1=true;
            // cout <<"a ";
        }
        else
        {
            strand_a_p1=false;
            // cout <<"b ";
        }
        if (gsl_rng_uniform(r)>0.5)
        {
            strand_a_p2=true;
            //cout <<"c ";
        }
        else
        {
            strand_a_p2=false;
            //cout <<"d ";
        }



        if (rho_CL==0.5) // if the mate choice locus is independent, it can use the function for freely recombining loci
        {
            inherit_FR_loci(r,offspring[index_offsprg].pref_a,parent1[parent1_ID[index_offsprg]].pref_a,parent1[parent1_ID[index_offsprg]].pref_b);
            inherit_FR_loci(r,offspring[index_offsprg].pref_b,parent2[parent2_ID[index_offsprg]].pref_a,parent2[parent2_ID[index_offsprg]].pref_b);
        }
        else     // otherwise decide which preference allele is inherited
        {
            //parent 1
            has_rec=gsl_rng_uniform(r);
            if (strand_a_p1 & (has_rec>=rho_CL)) // from chromosome a and no recombination -> allele a
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_a;
            }
            else if (strand_a_p1 & (has_rec<rho_CL))     // from chromosome a and recombination -> allele b
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_b;
            }
            else if (!strand_a_p1 & (has_rec>=rho_CL))   // from chromosome b and no recombination -> allele b
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_b;
            }
            else if (!strand_a_p1 & (has_rec<rho_CL))    // from chromosome b and recombination -> allele a
            {
                offspring[index_offsprg].pref_a=parent1[parent1_ID[index_offsprg]].pref_a;
            }
            //parent 2
            has_rec=gsl_rng_uniform(r);
            if (strand_a_p2 & (has_rec>=rho_CL)) // from chromosome a and no recombination -> allele a
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_a;
            }
            else if (strand_a_p2 & (has_rec<rho_CL))    // from chromosome a and recombination -> allele b
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_b;
            }
            else if (!strand_a_p2 & (has_rec>=rho_CL))   // from chromosome b and no recombination -> allele b
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_b;
            }
            else if (!strand_a_p2 & (has_rec<rho_CL))    // from chromosome b and recombination -> allele a
            {
                offspring[index_offsprg].pref_b=parent2[parent2_ID[index_offsprg]].pref_a;
            }
        }

        //cout << "nb CO " << co_nb<<endl;
        //cout << (*(genotypes[chr_ID])).sel<< endl;
        //cout << (*(genotypes[chr_ID_alt])).sel<< endl;

        //determine whether mutation happens at the preference locus and the value of the new allele
        if(gsl_rng_uniform(r)<mu_mod)// for allele a
        {
            mut_pref=gsl_ran_gaussian(r,var_pref);
            offspring[index_offsprg].pref_a += mut_pref;
            if (offspring[index_offsprg].pref_a<0)
                offspring[index_offsprg].pref_a=0;
            if (offspring[index_offsprg].pref_a>0.5)
                offspring[index_offsprg].pref_a=0.5;
        }
        if(gsl_rng_uniform(r)<mu_mod)// for allele b
        {
            mut_pref=gsl_ran_gaussian(r,var_pref);
            offspring[index_offsprg].pref_b += mut_pref;
            if (offspring[index_offsprg].pref_b<0)
                offspring[index_offsprg].pref_b=0;
            if (offspring[index_offsprg].pref_b>0.5)
                offspring[index_offsprg].pref_b=0.5;
        }

        //same but for the cost locus
        if(gsl_rng_uniform(r)<mu_mod)//allele a
        {
            mut_cost=gsl_ran_gaussian(r,var_cost);
            offspring[index_offsprg].cost_a += mut_cost;
            if (offspring[index_offsprg].cost_a<0.5)
                offspring[index_offsprg].cost_a=0.5;
        }
        if(gsl_rng_uniform(r)<mu_mod)//allele b
        {
            mut_cost=gsl_ran_gaussian(r,var_cost);
            offspring[index_offsprg].cost_b += mut_cost;
            if (offspring[index_offsprg].cost_b<0.5)
                offspring[index_offsprg].cost_b=0.5;
        }

        // inheritance of main genomic regions - loop over the different blocks - the strand_a_p1/2 variable is used to keep track on whether we end up on the "a" or "b" chromosome and to have the information ready for when looking at the next block.
        // recombination between linkage blocks is handle by the gamete linkage block itself, and return via strand_a the position we should take at the loci of the next linkage block
        for (index_chrom=0; index_chrom<nb_LD_blocks; index_chrom++)
        {
            // in case the linkage block does NOT include the neutral linked marker
            if (index_chrom!=(nb_LD_blocks/2))
            {
                gamete_linkage_block(r, rho_LG, rho_LB, dim, nb_loci, mu_GW, parent1[parent1_ID[index_offsprg]], index_chrom, offspring[index_offsprg].chr_a[index_chrom], mutations, strand_a_p1, Co_pos, off1, off2, rec, pos); // from the first parent
                gamete_linkage_block(r, rho_LG, rho_LB, dim, nb_loci, mu_GW, parent2[parent2_ID[index_offsprg]], index_chrom, offspring[index_offsprg].chr_b[index_chrom], mutations, strand_a_p2, Co_pos, off1, off2, rec, pos);// from the second parent
            }
            else
            {
                //in case the linkage block includes the neutral marker
                gamete_linkage_block_wNM(r, rho_LG, rho_LB, dim, nb_loci, mu_GW, parent1[parent1_ID[index_offsprg]], index_chrom, offspring[index_offsprg].chr_a[index_chrom],offspring[index_offsprg].neutral_La, mutations, strand_a_p1, Co_pos, off1, off2, rec, pos);//parent 1
                //if ((parent1[parent1_ID[index_offsprg]].neutral_La+parent1[parent1_ID[index_offsprg]].neutral_Lb)>1){
                //    cout << " p1"<< parent1[parent1_ID[index_offsprg]].neutral_La << parent1[parent1_ID[index_offsprg]].neutral_Lb<< offspring[index_offsprg].neutral_La;
                //}
                gamete_linkage_block_wNM(r, rho_LG, rho_LB, dim, nb_loci, mu_GW, parent2[parent2_ID[index_offsprg]], index_chrom, offspring[index_offsprg].chr_b[index_chrom],offspring[index_offsprg].neutral_Lb, mutations, strand_a_p2, Co_pos, off1, off2, rec, pos);//parent 2

                //if ((parent2[parent2_ID[index_offsprg]].neutral_La+parent2[parent2_ID[index_offsprg]].neutral_Lb)>1){
                //    cout << " p2"<< parent2[parent2_ID[index_offsprg]].neutral_La<<parent2[parent2_ID[index_offsprg]].neutral_Lb <<offspring[index_offsprg].neutral_Lb;
                //}
            }
        }

    }
}

//function to handle recombination between two linkage block of a same individual and produce the linage block of the offspring
void rec_v2(gsl_rng* r, Chr &res, const Chr &c1, const Chr &c2, const int &nbCo, const int &nb_loci, int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos, bool & inherit_posL)
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

    // determine the position of the cross overs and store them into pos
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

    // sort the pos vector
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
    inherit_posL=true; //check whether the neutral marker is on the same chromosome than the first locus of the linkage block or on the homologous chromosome
    if (rec[nb_loci/2]==1)// if the mask in the position of the marker holds a 1, it means it is on the homologous chromosome
        inherit_posL=false;
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

//function to generate F1 offspring for a single choosy (male) individual
void choose_parents_test_choosy(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other,const double &dist_opt_sq, const int &candidate_choosy, int &nb_F1, int &nb_mating)
{
    int candidate_other, i, j;
    double pheno_dist=0, counter, sigma;
    nb_F1=0;
    nb_mating=0;
    for (j = 0; j < n_off; j++)// loop over each offspring
    {
        //cout << j<< endl;
        //cout <<nb_choosy << " " << nb_other<<endl;
        counter=0; // reset the counter for number of potential pairing - if it goes above the cost value, the individual will fail to mate
        do
        {
            do
            {
                candidate_other = int(gsl_rng_uniform(r) * nb_other);
            }
            while ((gsl_rng_uniform(r) > pop_other[candidate_other].fitness));  //individual is rejected on its own fitness (pheno +DMI)
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
            //sigma=exp(pheno_dist/dist_opt_sq*log(pop_choosy[candidate_choosy].pref));



            //cout<< "choose other "<< candidate_other<<" counter "<<counter << " cost "<<pop_choosy[candidate_choosy].cost<< " dist "<<pheno_dist<< " pref " << pop_choosy[candidate_choosy].pref<< " sigma "<< sigma <<endl;
        }
        while ((gsl_rng_uniform(r) > sigma)&(counter<=pop_choosy[candidate_choosy].cost));  //mating pair is rejected if the phenotypic distance is too large and the number of trial has not reach its limit
        //cout << " counter "<< counter <<" ";

        if (counter<=pop_choosy[candidate_choosy].cost)    //check that the choosy individual was rejected if too many trials happened
        {
            //cout <<"mate accepted"<<endl;
            parent_choosy[nb_F1] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, in position 2*j and 2j+1 in the pop vector.
            parent_other[nb_F1] = candidate_other;
            nb_F1++; // this counter may differ from j because an individual may fail to accept a mate
            nb_mating+=counter; // total number of mating to generate the F1 offspring - it is conditional on the offspring being generated
        }
        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
    }
    //cout << nb_mating << " " << nb_F1<<  endl;
}

//function to generate F1 offspring for a single non-choosy (female) individual. males are still the choosy one (using their preference for rejection), but their only choice is to mate or not. So a male is chosen and try up to "cost" to mate with females and then give up.
void choose_parents_test_other(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other,const double &dist_opt_sq, const int &candidate_other, int &nb_F1, int &nb_mating)
{
    int candidate_choosy, i, j;
    double pheno_dist=0, counter, sigma;
    nb_F1=0;
    nb_mating=0;
    for (j = 0; j < n_off; j++)
    {
        //cout << j<< endl;
        //cout <<nb_choosy << " " << nb_other<<endl;
        counter=0;// reset the counter for number of potential pairing - if it goes above the cost value, the individual will fail to mate

        do
        {
            candidate_choosy = int(gsl_rng_uniform(r) * nb_choosy);
        }
        while ((gsl_rng_uniform(r) > pop_choosy[candidate_choosy].fitness));  //individual is rejected on its own fitness (pheno +DMI)
        //if (counter>0)
        //{
        //   cout<< "rejection"<<endl;
        //}
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
        do
        {
            counter++; // add +1 to the number of encounter


            //cout<< "choose choosy "<< candidate_other<<" counter "<<counter << " cost "<<pop_choosy[candidate_choosy].cost<< " dist "<<pheno_dist<< " pref " << pop_choosy[candidate_choosy].pref<< " sigma "<< sigma <<endl;
        }
        while ((gsl_rng_uniform(r) > sigma)&(counter<=pop_choosy[candidate_choosy].cost));  //mating pair is rejected if the pheno distance is too large and the number of trial has not reach its limit
        //cout << " counter "<< counter <<" ";

        if (counter<=pop_choosy[candidate_choosy].cost)     //check that the choosy individual was rejected if too many trials happened
        {
            //cout <<"mate accepted"<<endl;
            parent_choosy[nb_F1] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, in position 2*j and 2*j+1 in the pop vector.
            parent_other[nb_F1] = candidate_other;
            nb_F1++;
            nb_mating+=counter;
        }
        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
    }
    //cout << nb_mating << " " << nb_F1<<  endl;
}

//function to generate F1 offspring in the absence of mate choice - here parent_choosy can be male or female as it does not really matter
void choose_parents_test_migrant(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other, const int &candidate_choosy, int &nb_F1, int &nb_mating)
{
    int candidate_other, j;
    for (j = 0; j < n_off; j++)
    {
        //cout << j<< endl;


        do
        {
            candidate_other = int(gsl_rng_uniform(r) * nb_other);
        }
        while ((gsl_rng_uniform(r) > pop_other[candidate_other].fitness));  //individual is rejected on its own fitness (pheno +DMI)
        //if (counter>0)
        //{
        //   cout<< "rejection"<<endl;
        //}

        //cout<< "choose other "<< candidate_o<<" counter "<<counter << " cost "<<pop_choosy[candidate_c].cost<< " dist "<<pheno_dist<< " pref " << pop_choosy[candidate_c].pref<< endl;

        parent_choosy[j] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, funf in position 2*j and 2j+1 in the pop vactor.
        parent_other[j] = candidate_other;


        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
    }
    nb_F1=n_off;
    nb_mating=n_off;
    //cout <<endl;
}

//function to compute the mean for a few metrics of a given population, fitness components, preference, phenotype and cost
void update_mean_metrics(const int & nb_ind, const int &dim, double &mean_fit_abs, double &mean_fit_dmi, double &mean_fit, double &mean_pref, double &mean_cost, vector<double> &mean_pheno, const Ind *pop)
{
    mean_fit_abs=0;
    mean_fit_dmi=0;
    mean_fit=0;
    mean_pref=0;
    mean_cost=0;
    mean_pheno.clear();
    mean_pheno.resize(dim,0);
    int i, j;
    for (i=0; i <nb_ind; i++)
    {
        mean_fit_abs+=pop[i].abs_fit;
        mean_fit_dmi+=pop[i].dmi_fit;
        mean_fit+=pop[i].abs_fit*pop[i].dmi_fit;
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
    mean_fit/=nb_ind;
    mean_pref/=nb_ind;
    mean_cost/=nb_ind;
    for (j=0; j<dim; j++)
    {
        mean_pheno[j]/=nb_ind;
    }
}

//function to compute the mean for a few metrics of a given population, fitness components
void update_mean_metrics(const int & nb_ind, const int &dim, double &mean_fit_abs, double &mean_fit_dmi, double &mean_fit, const Ind *pop)
{
    mean_fit_abs=0;
    mean_fit_dmi=0;
    mean_fit=0;
    int i;
    for (i=0; i <nb_ind; i++)
    {
        mean_fit_abs+=pop[i].abs_fit;
        mean_fit_dmi+=pop[i].dmi_fit;
        mean_fit+=pop[i].abs_fit*pop[i].dmi_fit;
        //cout << pop[i].pheno[0] <<" ";
    }
    mean_fit_abs/=nb_ind;
    mean_fit_dmi/=nb_ind;
    mean_fit/=nb_ind;
}

//function to handle gamete formation for the linkage block that has a neutral marker attached to it. strand_a remembers whether we should starts on the "a" or "b" chromosome of the parent.
// recombination between linkage block is handled at the end of the function
void gamete_linkage_block_wNM(gsl_rng* r, const double &rho_LG, const double &rho_LB,  const int &dim, const int &nb_loci, const double &mu_GW,  const Ind &parent, const int & index_chrom, Chr &offspring_chrom, int &neutral_linked, double * mutations, bool &strand_a,int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
{
    int index_dim, j, site, ns;
    double pheno;
    bool linked_marker=true;
    int co_nb = gsl_ran_poisson(r,rho_LG);//number of Co per LB block
    if (co_nb>(nb_loci-1)) //check for overflow - note that is close to impossible
        co_nb=(nb_loci-1);
    int mut = gsl_ran_poisson(r,mu_GW); // number of mutations at selected loci on genome segment - here we are not even checking for overflow givne the exepected mutation rate.

    //if (((co_nb==0)|((parent1[parent1_ID[index_offsprg]].chr_a[index_chrom].sel==parent1[parent1_ID[index_offsprg]].chr_b[index_chrom].sel)))

    // first handle cases where we can ignore recombination becaue it does not happens or the two linkage blocks are identical
    if (((co_nb==0)|((parent.chr_a[index_chrom].sel==parent.chr_b[index_chrom].sel)&(parent.neutral_La==parent.neutral_Lb))))
    {
        //cout << "no recomb"<<endl;
        if (strand_a)  // determine whether it is the "a" or "b" chromosome to be passed - note that it only matters if co_nb ==0 copy the genotype and phenotype for the relevent linkage block
        {
            offspring_chrom.sel=parent.chr_a[index_chrom].sel;
            neutral_linked=parent.neutral_La;
            for (index_dim = 0; index_dim < dim; index_dim++)
                offspring_chrom.pheno[index_dim] = parent.chr_a[index_chrom].pheno[index_dim];
        }
        else
        {
            offspring_chrom.sel=parent.chr_b[index_chrom].sel;
            neutral_linked=parent.neutral_Lb;
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

        if(strand_a) // if we start from chromosome "a"
        {
            rec_v2(r, offspring_chrom, parent.chr_a[index_chrom], parent.chr_b[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec, pos, linked_marker);
            if (linked_marker)// true means that neutral marker is inherited like the the starting position
            {
                neutral_linked=parent.neutral_La;
                //cout << "aa ";
            }
            else
            {

                neutral_linked=parent.neutral_Lb;
                //cout << "ab ";
            }
        }
        else
        {
            rec_v2(r, offspring_chrom, parent.chr_b[index_chrom], parent.chr_a[index_chrom], co_nb, nb_loci, Co_pos,off1,off2,rec,pos,linked_marker);
            if (linked_marker)
            {
                neutral_linked=parent.neutral_Lb;
                //cout << "bb ";
            }
            else
            {
                neutral_linked=parent.neutral_La;
                //cout << "ba ";
            }
        }
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
        // recompute the pnenotype for the linkage blocks
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


//
//void update_fitness_v2(const int &nb_ind,const int &nb_LDblocks, const int &dim, const double &amp, const vector<double> &opt, Ind *pop, const double &eps, const int &nb_dmi, const int *dmi_pos, const int *target_chrom)
//{
//    double w_max = 0;
//    double d=0;
//    int i, j, k, target,locus;
//    if ((amp>0)&(eps>0))
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            d = 0;
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
//            }
//            //cout <<"calculate pheno and d_sq for ind"<<endl;
//            pop[i].abs_fit = exp(-amp * d);// fitness of individual i:
//            d=1;
//            for (j=0; j<nb_dmi; j++)
//            {
//                locus=dmi_pos[2*j];
//                target=target_chrom[2*j];
//                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//                locus=dmi_pos[2*j+1];
//                target=target_chrom[2*j+1];
//                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//
//                switch(k)
//                {
//                case 0:
//                    break;
//                case 1:
//                    d*=1-eps;
//                    break;
//                case 2:
//                    d*=(1-eps)*(1-eps);
//                    break;
//                case 4:
//                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
//                    break;
//                }
//            }
//            pop[i].dmi_fit=d;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            if (w_max < (pop[i].abs_fit*pop[i].dmi_fit))
//                w_max = (pop[i].abs_fit*pop[i].dmi_fit);
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//        for (i =0; i < nb_ind; i++)
//        {
//            pop[i].fitness = (pop[i].abs_fit*pop[i].dmi_fit)/w_max;
//            //cout <<pop[i].fitness<<" ";
//        }
//    }
//    else if ((amp>0)&(eps==0))
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            d = 0;
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                //   cout <<j<<endl;
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
//            }
//            //cout <<"calculate pheno and d_sq for ind"<<endl;
//            pop[i].abs_fit = exp(-amp * d);// fitness of individual i:
//            pop[i].dmi_fit=1;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            if (w_max < (pop[i].abs_fit))
//                w_max = (pop[i].abs_fit);
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//        for (i =0; i < nb_ind; i++)
//        {
//            pop[i].fitness = (pop[i].abs_fit)/w_max;
//            //cout <<pop[i].fitness<<" ";
//        }
//    }
//    else if ((amp==0)&(eps>0))
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                //   cout <<j<<endl;
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//            }
//
//            pop[i].abs_fit=1;
//            //cout <<"calculate pheno and d_sq for ind"<<endl;
//            d=1;
//            for (j=0; j<nb_dmi; j++)
//            {
//                locus=dmi_pos[2*j];
//                target=target_chrom[2*j];
//                // cout << "locus " << locus << " target " << target<< endl;
//                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//                locus=dmi_pos[2*j+1];
//                target=target_chrom[2*j+1];
//                // cout << "locus " << locus << " target " << target<< endl;
//                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//
//                switch(k)
//                {
//                case 0:
//                    break;
//                case 1:
//                    d*=1-eps;
//                    break;
//                case 2:
//                    d*=(1-eps)*(1-eps);
//                    break;
//                case 4:
//                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
//                    break;
//                }
//            }
//            pop[i].dmi_fit=d;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            if (w_max < (pop[i].dmi_fit))
//                w_max = (pop[i].dmi_fit);
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//        for (i =0; i < nb_ind; i++)
//        {
//            pop[i].fitness = (pop[i].dmi_fit)/w_max;
//            //cout <<pop[i].fitness<<" ";
//        }
//    }
//    else
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                //   cout <<j<<endl;
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//            }
//            pop[i].abs_fit=1;
//            pop[i].dmi_fit=1;
//            pop[i].fitness = 1;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//    }
//}
//
//
//void update_fitness(const int &nb_ind,const int &nb_LDblocks, const int &dim, const double &amp, const double &half_epistasis,const vector<double> &opt, Ind *pop, const double &eps, const int &nb_dmi, const int *dmi_pos, const int * target_chrom)
//{
//    double w_max = 0;
//    double d;
//    int i, j, k, target, locus;
//    if ((amp>0)&(eps>0))
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            d = 0;
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
//            }
//            //cout <<"calculate pheno and d_sq for ind"<<endl;
//            pop[i].abs_fit = exp(-amp * pow(d,half_epistasis));// fitness of individual i:
//            d=1;
//            for (j=0; j<nb_dmi; j++)
//            {
//                locus=dmi_pos[2*j];
//                target=target_chrom[2*j];
//                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//                locus=dmi_pos[2*j+1];
//                target=target_chrom[2*j+1];
//                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//
//                switch(k)
//                {
//                case 0:
//                    break;
//                case 1:
//                    d*=1-eps;
//                    break;
//                case 2:
//                    d*=(1-eps)*(1-eps);
//                    break;
//                case 4:
//                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
//                    break;
//                }
//            }
//            pop[i].dmi_fit=d;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            if (w_max < (pop[i].abs_fit*pop[i].dmi_fit))
//                w_max = (pop[i].abs_fit*pop[i].dmi_fit);
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//        for (i =0; i < nb_ind; i++)
//        {
//            pop[i].fitness = (pop[i].abs_fit*pop[i].dmi_fit)/w_max;
//            //cout <<pop[i].fitness<<" ";
//        }
//    }
//    else if ((amp>0)&(eps==0))
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            d = 0;
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                //   cout <<j<<endl;
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//                d += (pop[i].pheno[j] - opt[j]) * (pop[i].pheno[j] - opt[j]); // "d" is the square of the distance to the optimum
//            }
//            //cout <<"calculate pheno and d_sq for ind"<<endl;
//            pop[i].abs_fit = exp(-amp * pow(d,half_epistasis));// fitness of individual i:
//            pop[i].dmi_fit=1;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            if (w_max < (pop[i].abs_fit))
//                w_max = (pop[i].abs_fit);
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//        for (i =0; i < nb_ind; i++)
//        {
//            pop[i].fitness = (pop[i].abs_fit)/w_max;
//            //cout <<pop[i].fitness<<" ";
//        }
//    }
//    else if ((amp==0)&(eps>0))
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                //   cout <<j<<endl;
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//            }
//            pop[i].abs_fit=1;
//            //cout <<"calculate pheno and d_sq for ind"<<endl;
//            d=1;
//            for (j=0; j<nb_dmi; j++)
//            {
//                locus=dmi_pos[2*j];
//                target=target_chrom[2*j];
//                k=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//                locus=dmi_pos[2*j+1];
//                target=target_chrom[2*j+1];
//                k*=pop[i].chr_a[target].sel[locus]+pop[i].chr_b[target].sel[locus];
//
//                switch(k)
//                {
//                case 0:
//                    break;
//                case 1:
//                    d*=1-eps;
//                    break;
//                case 2:
//                    d*=(1-eps)*(1-eps);
//                    break;
//                case 4:
//                    d*=(1-eps)*(1-eps)*(1-eps)*(1-eps);
//                    break;
//                }
//            }
//            pop[i].dmi_fit=d;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            if (w_max < (pop[i].dmi_fit))
//                w_max = (pop[i].dmi_fit);
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//        for (i =0; i < nb_ind; i++)
//        {
//            pop[i].fitness = (pop[i].dmi_fit)/w_max;
//            //cout <<pop[i].fitness<<" ";
//        }
//    }
//    else
//    {
//        for (i = 0; i < nb_ind; i++)
//        {
//            d = 0;
//            pop[i].pref=pop[i].pref_a+pop[i].pref_b;// calculate pref for individual - additive
//            //cout <<"calculate pref for ind "<< i<<endl;
//            pop[i].cost=pop[i].cost_a+pop[i].cost_b;// calculate cost for individual - additive
//            //cout <<"calculate cost for ind"<<endl;
//            for (j = 0; j < dim; j++)
//            {
//                //   cout <<j<<endl;
//                pop[i].pheno[j]=0;
//                //   cout <<j<<endl;
//                for (k=0; k< nb_LDblocks; k++)
//                    pop[i].pheno[j] += pop[i].chr_a[k].pheno[j]+pop[i].chr_b[k].pheno[j];
//            }
//            pop[i].abs_fit=1;
//            pop[i].dmi_fit=1;
//            pop[i].fitness = 1;
//            //if(pop[i].dmi_fit<1)
//            //{
//            //    cout <<"fit "<<pop[i].dmi_fit<<endl;
//            //}
//            //cout << "calculate fitness"<<endl;
//            //cout << gen<< " "<< opt[0]<< " "<<(*(hap[2*i])).sex_chr+ (*(hap[2*i+1])).sex_chr << " "<< pop[i].abs_fit<< " "<<pop[i].dmi_fit<< " "<< pop[i].pheno[0]<< endl;
//        }
//    }
//}
//
//
//void rec_v2(gsl_rng* r, Chr &res, const Chr &c1, const Chr &c2, const int &nbCo, const int &nb_loci, int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos)
//{
//    //vector<int> pos{ 0, nb_loci };
//    pos.clear();
//    pos.push_back(0);
//    pos.push_back(nb_loci);
//    int j,k;
//
//    //cout << "1C "<<c1.sel<<endl;
//    //cout << "2C "<<c2.sel<<endl;
//
//
//    //res.sel.clear();
//    //cout << "enter rec"<<endl;
//    // vector "pos" holds the positions of cross-overs:
//
//    //gsl_ran_shuffle(r,Co_pos,nb_loci,sizeof (int));
//    j=0;
//    do
//    {
//        k = int(gsl_rng_uniform(r)*nb_loci);
//        if (find(pos.begin(), pos.end(), k) != pos.end())
//        {
//            pos.push_back(k);
//            j++;
//        }
//    }
//    while(j<nbCo);
//
//    //for (j = 0; j < nbCo; j++)
//    //{
//
//    //    pos.push_back(Co_pos[j]);
//    //cout<<pos[j+2]<<" ";
//    //}
//    //cout <<endl;
//
//    sort(pos.begin(), pos.end());
//
//    //cout << "nb CO "<<nbCo<<endl;
//    // creates recombination mask:
//    rec.reset();
//    for (j = 0; j < (nbCo+1); j+=2)
//    {
//        //cout << pos[j] << " ";
//        k=pos[j+1]-pos[j];
//        rec.flip(pos[j],k);
//    }
//    //cout <<endl;
//    //cout << "re "<<rec<<endl;
//
//    off2 = (c2.sel & rec);// we start by c2, because we want to begin the chromosome by c1 (and the recombination vector tells which part to hide - and not to keep)
//    rec.flip();
//    //cout << "fl "<<rec<<endl;
//    //cout << "c1 "<<(c1.sel)<<endl;
//    //cout << "te "<<(rec)<<endl;
//    off1 = (c1.sel & rec);
//    //cout <<"ok"<<endl;
//    //cout << off1 << endl;
//    //cout << off2 <<endl;
//    //cout <<rec<<endl;
//    //cout <<&off1<<endl;
//    //cout <<&off2<<endl;
//    //cout <<&rec<<endl;
//    //cout << &c1.sel<<endl;
//    // cout << &c2.sel<<endl;
//    //cout << &res.sel<<endl;
//    //cout << &pos<<endl;
//    res.sel=c1.sel;
//    rec=(off1 | off2);
//    //cout <<rec<<endl;
//    //cout << (off1 | off2)<< endl;
////   cout << "off "<<res.sel.size()<<endl;
//    res.sel = rec;
//    // cout <<"redundant"<<endl;
//    //res.sel = (off1 | off2);
//
//
////cout << res.sel << endl;
//
////   for (j = 0; j < nbCo; j++)
////	cout<<pos[j+2]<<" ";
//
//
//    //cout << "FC "<<res.sel<<endl;
//    //  cout <<endl;
//}

//
//
//void choose_parents(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other, const double &dist_opt_sq)
//{
//    int candidate_choosy, candidate_other, i, j;
//    double pheno_dist=0, counter, sigma;
//    for (j = 0; j < n_off; j++)
//    {
//        //cout << j<< endl;
//        //cout <<nb_choosy << " " << nb_other<<endl;
//        do
//        {
//            do
//            {
//                candidate_choosy = int(gsl_rng_uniform(r) * nb_choosy);// choose random individual
//            }
//            while (gsl_rng_uniform(r) > pop_choosy[candidate_choosy].fitness);   //individual is rejected on its own fitness (pheno +DMI)
//            //cout << "choose choosy"<< candidate_c << endl;
//            counter=0;
//            do
//            {
//                do
//                {
//                    candidate_other = int(gsl_rng_uniform(r) * nb_other);
//                }
//                while ((gsl_rng_uniform(r) > pop_other[candidate_other].fitness));  //individual is rejected on its own fitness (pheno +DMI) - we keep an indiviuals if we fail the test, ie the number is below fitness
//                //if (counter>0)
//                //{
//                //   cout<< "rejection"<<endl;
//                //}
//                counter++; // add +1 to the number of encounter
//                pheno_dist=0.0; // compute the square of the phenotypic distance
//                for (i=0; i<dim; i++)
//                {
//                    pheno_dist+=(pop_choosy[candidate_choosy].pheno[i]-pop_other[candidate_other].pheno[i])*(pop_choosy[candidate_choosy].pheno[i]-pop_other[candidate_other].pheno[i]);
//                }
//                sigma=exp(pheno_dist/dist_opt_sq*log(pop_choosy[candidate_choosy].pref)); // here sigma is 0 for really far away phenotype, becuase the log is negative
//
//                //cout<< "choose other "<< candidate_o<<" counter "<<counter << " cost "<<pop_choosy[candidate_c].cost<< " dist "<<pheno_dist<< " pref " << pop_choosy[candidate_c].pref<< endl;
//            }
//            while ((gsl_rng_uniform(r) > sigma)&(counter<=pop_choosy[candidate_choosy].cost));  //mating pair is rejected (ie the test is true when the distance is large) if the pheno distance is too large and the number of trial has not reach its limit
//        }
//        while (counter>pop_choosy[candidate_choosy].cost);   //the choosy individual is rejected if too many trials happened
//        //cout <<"mate accepted"<<endl;
//        parent_choosy[j] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, funf in position 2*j and 2j+1 in the pop vactor.
//        parent_other[j] = candidate_other;
//        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
//    }
//    //cout <<endl;
//}
//
//
//void choose_parents(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other)
//{
//    int candidate_choosy, candidate_other, j;
//    for (j = 0; j < n_off; j++)
//    {
//        //cout << j<< endl;
//        //cout <<nb_choosy << " " << nb_other<<endl;
//
//        do
//        {
//            candidate_choosy = int(gsl_rng_uniform(r) * nb_choosy);// choose random individual
//        }
//        while (gsl_rng_uniform(r) > pop_choosy[candidate_choosy].fitness);   //individual is rejected on its own fitness (pheno +DMI)
//        //cout << "choose choosy"<< candidate_c << endl;
//
//        do
//        {
//            candidate_other = int(gsl_rng_uniform(r) * nb_other);
//        }
//        while ((gsl_rng_uniform(r) > pop_other[candidate_other].fitness));  //individual is rejected on its own fitness (pheno +DMI)
//        //if (counter>0)
//        //{
//        //   cout<< "rejection"<<endl;;   //the choosy individual is rejected if too many trials happened
//        //cout <<"mate accepted"<<endl;
//        parent_choosy[j] = candidate_choosy; //assign this individual as parent 1; the entry here corresponds to the individual j, funf in position 2*j and 2j+1 in the pop vactor.
//        parent_other[j] = candidate_other;
//        //cout<< parent_choosy[j]<< "-" <<parent_other[j] <<"; ";
//    }
//    //cout <<endl;
//}

//
//void update_population(const int &nb_Ind, Ind *pop_male, Ind *pop_female, Ind *temp_pop, int &nb_male, int &nb_female)
//{
//    nb_male=0;
//    nb_female=0;
//    for (int i = 0; i < nb_Ind; i++)
//    {
//        //cout <<i<<"\n";
//        //cout <<"chr1 ori "<<(*(temp_pop[i])).neutral << " chr2 ori "<< (*(temp_pop[i+1])).neutral<<endl;
//        //   cout <<"pop1 sex "<<(*(temp_1[i])).sex_chr << " "<< (*(temp_1[i+1])).sex_chr<<endl;
//        //  cout <<"pop2 sex "<<(*(temp_2[i])).sex_chr << " "<< (*(temp_2[i+1])).sex_chr<<endl;
//        if (temp_pop[i].sex)
//        {
//            //cout<<"pass\n";
//            pop_male[nb_male] = temp_pop[i];
//            nb_male++;
//        }
//        else
//        {
//            //cout<<"failed\n";
//            pop_female[nb_female] = temp_pop[i];
//            nb_female++;
//        }
//    }
//    //cout << nb_Ind <<" " <<nb_male <<" "<< nb_female <<endl;
//}
