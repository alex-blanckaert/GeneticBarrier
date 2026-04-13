// Header file: definitions of global variables, function prototypes

#ifndef ALL_BARRIER_H
#define ALL_BARRIER_H

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


// Global variables:

// #define fichierLecture "parametres4.txt"     // names of input
//#define fichierEcriture "resultats.txt"		// and output files

// "Chr": represents an haploid version of the genome

// TO DO: remove all pow()

// TO DO change signel chromosome to 10 linkage group, to lessen the need to recompute phenotype

// STRUCTURE

struct Chr
{
	boost::dynamic_bitset<> sel; // selected loci (chain of 0 and 1) note that they are written from right to left !!!!!
	vector<double> pheno;
};


struct Ind
{
    Chr chr_a[5];
    Chr chr_b[5];
    vector<double> pheno;
    double abs_fit;
    double dmi_fit;
    double fitness;
    bool sex;
    double pref=0.0;
    double cost=0.0;
    double cost_a=0.0;
    double pref_a=0.0;
    double cost_b=0.0;
    double pref_b=0.0;
    int neutral_a=0;
    int neutral_b=0;
    int neutral_La=0;
    int neutral_Lb=0;
};
struct compare
{
    int key;
    compare(int const &i): key(i) {}

    bool operator()(int const &i) {
        return (i == key);
    }
};

// FROM UNIQUE
bool lireFichier(const char *param_file, unsigned int &seed, int &n_ind, int &nb_loci, int &dim, int &gen_Max, double &mu, double &mu_mod, double &mig, double &rho, double &rho_LB, double &rho_CL, double &sig, double &amp, double &epistasis, vector<double> &opt1, vector<double> &opt2, double &var_pref, double &var_cost, int &nb_dmi, double &eps, int &dmi_neutral, double &ini_pref, double & ini_cost, int& skip_mate_choice, int & common_origin, int & export_metrics, int &cycle_mig_sel, int &load_mutation_file);

void add_dmi(Ind &ind, const int &dim, const int &nb_dmi, const int &nb_loci_per_LG, const int &locus, const int *dmi_pos, const int * target_chrom, const double * mutations);
void move_to_optimum(Ind &ind, const int &dim, const int &nb_Lblocks, const int & nb_loci_per_LG, const vector<double> &opt, const int * target_chrom, const double * mutations);

void initialize_pop(Ind * pop,const int &n_ind, const int &nb_Lblocks, const int & nb_loci_per_LG, const int &dim, const double &pref, const double &cost, const int &common_origin, const vector<double> &opt, const int &locus, const int &nb_dmi, const int *dmi_pos, const int * target_chrom, const double * mutations);

void migration(gsl_rng* r,Ind *pop1, Ind *pop2, Ind *temp1, const double &mig, int &nb_ind_1, int &nb_ind_2);

void new_gamete(gsl_rng* r, double * mutations, const int &n_ind,const double &mu_GW, const double &mu_mod, const int *parent1_ID, const int *parent2_ID, Ind *parent1, Ind *parent2, Ind *offspring, const int &nb_loci, const int &dim, const int &nb_LD_blocks,const double &rho, const double &rho_LG, const double &rho_LB, const double &rho_CL, const double & var_pref, const double & var_cost,int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos);

void update_population(const int &nb_Ind, Ind *pop_male, Ind *pop_female, Ind *temp_pop, int &nb_male, int &nb_female);

void export_pop(ofstream &fout, const int &nb_ind, const int &dim, const int &nb_Lblocks, Ind *pop);


//FROM SHARED
bool lireMutationFile(const char *param_file, gsl_rng* r, const int &nb_LDblocks, const int & nb_loci_per_LG, double * mutations, int * pos_dmi, int *target_chrom, const int &dim, const int& nb_loci, const int &nS, const double &sig, const int &nb_dmi, const int &dmi_neutral, const int &load_mutation_file);

bool reload_population(const char *pop_file, const int &n_ind, const int &nb_loci, const int &dim, const double &mu, const double &mu_mod, const double &mig, const double &rho, const double &rho_LB, const double &rho_CL, const double &sig, const double &amp, const double &epistasis, const vector<double> &opt1, const vector<double> &opt2, const double &var_pref, const double &var_cost, const int &nb_dmi, const double &eps, int &nb_female_1, int &nb_male_1, int &nb_female_2, int &nb_male_2,Ind *pop1_females, Ind *pop1_males, Ind *pop2_females, Ind *pop2_males, const int &nb_LDblocks, const int & nb_loci_per_LB, const int &cycle_mig_sel);

//int count_derived_allele(Ind *ind,const int nb_locus_per_LB, const int locus);


void choose_parents(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other, const double &dist_opt_sq);
void choose_parents(gsl_rng* r, const int &dim, int *parent_choosy, int *parent_other, const int &n_off, const int &nb_choosy, const int &nb_other,const Ind *pop_choosy,const Ind *pop_other);

//void rec(gsl_rng* r, Chr *res, Chr &c1, Chr &c2, const int &nbCo, const int &nb_loci);
void rec_v2(gsl_rng* r, Chr &res, const Chr &c1, const Chr &c2, const int &nbCo, const int &nb_loci,int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos);

void update_mean_metrics(const int & nb_ind, const int &dim, double &mean_fit, double &mean_fit_abs, double &mean_fit_dmi, vector <double> &mean_pheno, double &mean_pref, double &mean_cost, const Ind *pop);
//void parent_mean_metrics(const int & nb_ind, double &mean_fit_abs, double &mean_fit_dmi, double &mean_pheno, double &mean_pref, double &mean_cost, const Ind *pop, const int* parent);
//void more_statistics(const int &nb_ind, Chr **hap, const vector<double> &opt, const double &amp, const double &half_epistasis, double &mean_hom_fit, double &h_0, const double &nb_locus, const int &dim, double &mean_origin);
//void parent_stat(const int & nb_ind, const int* parent, const int & nb_parent, double &var);

void inherit_FR_loci(gsl_rng* r, int &new_value, const int &chr_1, const int &chr_2);
void inherit_FR_loci(gsl_rng* r, double &new_value, const double &chr_1, const double &chr_2);
//void inherit_FR_loci(gsl_rng* r, bool &new_value, const bool &chr_1, const bool &chr_2);

void update_fitness(const int &nb_ind,const int &nb_LDblocks, const int &dim, const double &amp, const double &half_epistasis,const vector<double> &opt, Ind *pop, const double &eps, const int &nb_dmi, const int *dmi_pos, const int * target_chrom);
void update_fitness_v2(const int &nb_ind,const int &nb_LDblocks, const int &dim, const double &amp, const vector<double> &opt, Ind *pop, const double &eps, const int &nb_dmi, const int *dmi_pos, const int * target_chrom);

void gamete_linkage_block(gsl_rng* r, const double &rho_LG, const double &rho_LB, const int &dim, const int &nb_loci, const double &mu_GW,  const Ind &parent, const int & index_chrom, Chr &offspring_chrom, double * mutations, bool &strand_a,int *Co_pos, boost::dynamic_bitset<> &off1, boost::dynamic_bitset<> &off2, boost::dynamic_bitset<> &rec, vector<int> &pos);


#endif
