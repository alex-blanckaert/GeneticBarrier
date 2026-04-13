#!/bin/bash

# DO NOT USE NUMBERS OR SIGNs IN THE TEXT; FOR OPTIMUM, further dimension should be added on the SAME line, following the "string number" pattern.
# all parameters value are defines here

#main parameters
pop_size=10000
nb_loci=500
dim=1
mu=0.000001
muMC=0.00001
mig=0.05
rho=0.001
rho_LB=0.01
rho_CL=0.1
sig=0.05
epistasis_Fisher=2
opt_one_dimI=5
opt_two_dimI=-5
var_pref=0.1 # no evolution of preference is obtained by fixing var_pref to 0; however preference still happens at the initial value
var_cost=0 # no evolution of cost is obtained by fixing var_cost to 0
nb_dmi=50
dmi_neutral=0 # 0 mutation forming DMIs have phenotypic effect; 1 mutations forming DMI dont have a phenotypic effect
ini_pref=1 # no preference is obtained by fixing ini_pref to 1
ini_cost=100 # max cost is obtained by fixing ini_cost to 1
load_mutation=1 # 1 if an exisiting file is provided, 0 if a new one must be created


#parameter directly affecting the barriers
fit_alt=0.1 # fit_alt controls indirectly the steepness of the phenotype to fitness function; it corresponds to the fitness of "best phenotype" in the alternative environment and indirectlycontrol the A parameter in the Gros (2009) equation w=exp(-A sum((z2-z1)²)^epistasis_Fisher/2)
epsilon_DMI=0.1 # epistasis of a single DMI in the double heterozygote (assuming codominance)
skip_mating_choice=0 # 0 for mate choice, 1 to turn it off
life_cycle=0 # 1 mig then sel, 0 sel then mig


#staring conditions for testing purpose
nb_gen=51
export_metrics=10 # time step for exporting metrics

#staring conditions default value
#nb_gen=500001
#export_metrics=10000 # time step for exporting metrics
start=0 # start 0 for ancestral state; 1 for diverged pop (LA, at the optimum); 2 for DMI without LA;3 for DMI and LA; 4 to load existing population



mutation_file="./mutations_500L_1D_var0d05_50DMI.txt" 

if test -f $mutation_file
  then
  echo "existing mutation file"
  else
  exit # this line stops the script if the mutation file is not found. To allow the program to generate a new one, this line must be commented. In addition, the load_mutation variable must be set to 0.
fi

name="output_default" # name of the output folder

###------no need to edit anything below"
if test -d "$name"/
  then
    echo "Directory $name/ already exists"
  else
    mkdir $name # create the folder
fi
cd $name # move into the folder

bash_file=$0
while test -f $bash_file
do
  sleep 2
  bash_file="${0%.sh}_$(date --utc +%Y%m%d_%H%M%SZ).sh"
done

cp "./.$bash_file" ./

for j in 1 2 3 4 5 6 7 8 9 10
  do

    parameter_file="parameter_$(date --utc +%Y%m%d_%H%M%SZ).txt"
    while test -f $parameter_file
    do
      sleep 2
      parameter_file="parameter_$(date --utc +%Y%m%d_%H%M%SZ).txt"
    done
    seed=$(python3 -c 'import random as R; print(R.randint(1, 2**32-1))')
    echo -e "*\nseed $seed \npop_size $pop_size \nnb_loci $nb_loci \ndim $dim \nnb_gen $nb_gen \nmu $mu \nmu_mc $muMC \nmig $mig \nrho $rho \nrho_LB $rho_LB \nrho_CL $rho_CL \nsig $sig \nfit_alt $fit_alt \nepistasis_Fisher $epistasis_Fisher \nopt_one_dimI $opt_one_dimI \nopt_two_dimI $opt_two_dimI \nvar_pref $var_pref \nvar_cost $var_cost \nnb_dmi $nb_dmi \nepsilon_DMI $epsilon_DMI \ndmi_neutral $dmi_neutral \nini_pref $ini_pref \nini_cost $ini_cost \nskip_mating_choice $skip_mating_choice \nstart $start \nexport_metrics $export_metrics \nlife_cycle $life_cycle \nmut_file $load_mutation\n*" >> $parameter_file # write the parameter file

    name_pop_file="pop_LC_${life_cycle}_LA${fit_alt}_MC$((1-$skip_mating_choice))_DMI${epsilon_DMI}_rep$j.txt" #define name output
    name_metric_file="metric_$name_pop_file" #define name output


    echo $name_pop_file
    while test -f $name_pop_file
    do
      sleep 10
      name_pop_file="pop_LA${fit_alt}_MC$((1-$skip_mating_choice))_DMI${epsilon_DMI}_rep${j}_$(date --utc +%Y%m%d_%H%M%SZ).txt"
      name_metric_file="metric_$name_pop_file"
      echo "$name_pop_file was renamed"
      echo "$name_metric_file"
    done
    ./../../main_version/multiLoci $parameter_file ./../${mutation_file} $name_pop_file $name_metric_file & # call for the program; output file has the same name than the folder plus an extra part corresponding to the recombination rate for easy identification. Note that since the parameters value are written in the file itself a mistake here, while not desirable is not critical.

#done
#cd ..
done

