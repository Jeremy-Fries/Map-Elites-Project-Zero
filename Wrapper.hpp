//
//  Wrapper.hpp
//  Mapping_Elites v1
//
//  Created by Jeremy Fries on 11/16/15.
//  Copyright © 2015 Jeremy Fries. All rights reserved.
//

class Wrapper;

#ifndef Wrapper_hpp
#define Wrapper_hpp

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "Individual.hpp"
#include "Map_space.hpp"
#include "Map_Elites.hpp"

using namespace std;

class Wrapper{
protected:
    Map_Elites ME;
    Map_Elites* pME = &ME;
    Map_space mspace;
    Map_space* pMspace =  &mspace;
    Individual individual1;
    Individual* pI = &individual1;
    
public:
    const int hidden_layer_size = 3;
// --------------------------------------------------
    void initialize_wrapper();
    void wrapper_sets_I_params(int size1, int size2, double mut_mag1, double mut_mag2, int mut_amo1, int mut_amo2);
// --------------------------------------------------
        // Loop functions combine all
    void wrapper_runs_sim(vector<double>,vector<double>);
    void fill_MAP();
    void rand_bin();
    void mutate_MAP();
    void run_single_individual();
// ----------------------------------------------------------------------------------- P0 Functions
        // P0 functions
    void build_real_set();
    void empty_soln_set();
    void fitness_P0();
    double f_1(double, double, double);
    double f_2(double, double, double);
    double real_func();
    double approx_func();
    void fill_MAP_P0();
    void new_x_pos();
    void set_x_pos(double);
    double get_x_pos();
    void display_x_pos();
    void find_phenotypes();
    void display_phenotypes();
    void graph_final();
    
    // Phenotypes: min, max,            integral output, range, calc d", # of pts d' is + or -, sum of average slopes,
    
    vector<double> best_fit;
    
private:
    int isize_1, isize_2, imutate_amount_1, imutate_amount_2;
    double imutate_mag_1,imutate_mag_2;
    int bin1, bin2;
    
// P0 variables
    vector <double> real_set;
    vector <double> soln_set;
    double real_value;
    double approx_value;
    double x_pos;
    double calc_fit_rating;
    double fit_rating;
    double phenotype_1, phenotype_2;
    
    
};
// ----------------------------------------------------------------------------------- P0 Functions
void Wrapper::build_real_set(){
    // [ a1 b1 a2 b2]
    // Creates [1 1 0 0]
    real_set.clear();
    real_set.push_back(0);
    real_set.push_back(0);
    real_set.push_back(0);
    real_set.push_back(5);
}
// ----------------------------------------------------------------------------------- P0 Functions
double Wrapper::f_1(double x, double a1, double b1){
	double f = (a1*(x*x)) + b1;
    return f;
}
double Wrapper::f_2(double x, double a2, double b2){
	double f = (a2*x)+b2;
    return f;
}
// ----------------------------------------------------------------------------------- P0 Functions
// Real Function Value - real coefficient_set's value at x_pos
double Wrapper::real_func(){
    double x = get_x_pos();                                 // P0 TODO - set x pos
    double a, b;
    a = f_1(x, real_set.at(0), real_set.at(1));
    b = f_2(x, real_set.at(2), real_set.at(3));
    real_value = (a + b);
    //cout << endl << "real value is: " << real_value << endl;
    return real_value;
}
// Approx function Value - solution_set's value at x_pos
double Wrapper::approx_func(){
    double x = get_x_pos();
    double a, b;
    a = f_1(x, soln_set.at(0), soln_set.at(1));
    b = f_2(x, soln_set.at(2), soln_set.at(3));
    approx_value = (a + b);
    //cout << endl << "Approximate Value is: " << approx_value << endl;
    return approx_value;
}
// ----------------------------------------------------------------------------------- P0 Functions
// Set x position to desired value
void Wrapper::new_x_pos(){
    double pre_x = ((double)rand() / RAND_MAX);
    double x = pre_x * 2 * 3.1416;
    x_pos = x;
}
// Set x position
void Wrapper::set_x_pos(double x){
    x_pos=x;
}
// Get x position
double Wrapper::get_x_pos(){
    return x_pos;
}
// Display x position
void Wrapper::display_x_pos(){
    cout << "\t x_pos is: " << x_pos << endl;
}
// ----------------------------------------------------------------------------------- P0 Functions
void Wrapper::fitness_P0(){
    calc_fit_rating=0;
    fit_rating=0;
    for (int q = 0; q < 20; q++){
        new_x_pos();
        calc_fit_rating -= abs(real_func() - approx_func());
    }
    fit_rating=calc_fit_rating;
}
// ----------------------------------------------------------------------------------- P0 Functions
// Phenotypes
void Wrapper::find_phenotypes(){
    phenotype_1=0;
    phenotype_2=0;
    //double x = get_x_pos();
    //double a, b;
    //a = f_1(x, soln_set.at(0), soln_set.at(1));
    //b = f_2(x, soln_set.at(2), soln_set.at(3));
    
    double r1 = ((double)rand() / RAND_MAX);
    double r2 = ((double)rand() / RAND_MAX);
    
    phenotype_1=r1*100;
    phenotype_2=r2*100;
    
    //cout << a << " , " << b << endl;
}
// --------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------
//
void Wrapper::rand_bin(){
    //rand()%(max-min)+min;
    bin1=0;
    bin2=0;
    
    int b1_min=0;
    int b2_min=0;
    int b1_max=ME.get_resolution1();
    int b2_max=ME.get_resolution2();
    
    bin1=rand()%(b1_max-b1_min)+b1_min;
    bin2=rand()%(b2_max-b2_min)+b2_min;
    
    //cout << "bin to mutate is: (" << bin1 << "," << bin2 << ")" << endl;
}
// -------------------------------------------------------------------------------------------------------- To Change Settings
void Wrapper::initialize_wrapper(){
    
    pME->set_map_params(0, 100, 0, 100, 10, 10, 100, 1000000);
    // (dim1_min, dim1_max, dim2_min, dim2_max, resolution1,2, fill generation, mutate generation)
    
    //pME->display_Map_params();
    
    wrapper_sets_I_params(4, 0, 1, 0, 2, 0); // P0 - only one genome in individual.
    // individual_size 1,2, mutate_magnitude 1,2, mutation_amount 1,2)

}
// --------------------------------------------------
void Wrapper::wrapper_sets_I_params(int size1, int size2, double mut_mag1, double mut_mag2, int mut_amo1, int mut_amo2){
    isize_1=size1;
    isize_2=size2;
    imutate_mag_1=mut_mag1;
    imutate_mag_2=mut_mag2;
    imutate_amount_1=mut_amo1;
    imutate_amount_2=mut_amo2;
    
    /// assigns values to actual individual.
    //individual1.set_individual_params(size1,size2,mut_mag1,mut_mag2,mut_amo1,mut_amo2);
    //set_individual_params(int individual_size1,int individual_size2,double mutation_magnitude1,double mutation_magnitude2,int mutation_amount1,int mutation_amount2);
}
// ----------------------------------------------------------------------------------- P0 Functions
void Wrapper::fill_MAP_P0(){
    int fill_gen = ME.get_fill_generation();
    int fill_round=0;
    build_real_set();
    
    for (int g=0; g<fill_gen; g++){
        //cout << "fill round is: " << g << endl;

        
        Individual I;
        I.set_individual_params(isize_1, isize_2, imutate_mag_1, imutate_mag_2, imutate_amount_1, imutate_amount_2);
        //I.display_individual_params();
        
        I.build_individual();
        //I.display_individual1();
        //I.display_individual2();
        
        soln_set.empty();       // only for P0
        soln_set=I.genome1;     // only for P0
        fitness_P0();
        I.set_fit_rating(fit_rating);
        //I.display_fit_rating();
        set_x_pos(10);   // at x=10, the value of a and b are the phenotypes    // only for P0
        find_phenotypes();
        I.set_phenotypes(phenotype_1,phenotype_2);
        //I.display_phenotype1();
        //I.display_phenotype2();
        
        ME.place_individual_in_map(I);
        fill_round++;
        // push_back best fit rating to vector
        
    }
    cout << "completed " << fill_round << " FILL rounds" << endl;
}
// --------------------------------------------------
void Wrapper::mutate_MAP(){                         // ROMOVE COMMENTS
    int mut_gen = ME.get_mutate_generation();
    int mutation_round=0;
    for (int g=0; g<mut_gen; g++){
        //cout << "mutation round is: " << g << endl;
        soln_set.empty();       // only for P0

        Individual I;
        I.set_individual_params(isize_1, isize_2, imutate_mag_1, imutate_mag_2, imutate_amount_1, imutate_amount_2);
        rand_bin();                             // random bin location
        ME.individual_from_map(bin1, bin2);     // copy Individual's vectors into temporary holder
        I.build_individual_1_from_another(ME.get_temp_individual1());              // pass temp vectors to new Individual, setting individual1
        //I.build_individual_1_from_another(ME.get_temp_individual2());              // pass temp vectors to new Individual, setting individual2
        I.mutate1();                            // mutate Individual vectors
        //I.mutate2();                            // mutate Individual vectors
        //I.display_individual1();                // display mutated genome
        //I.display_individual2();                // display mutated genome
        
        soln_set=I.get_individual1();     // only for P0
        fitness_P0();                           // fitness for new Individual
        I.set_fit_rating(fit_rating);
        //I.display_fit_rating();
        
        set_x_pos(10);   // at x=10, the value of a and b are the phenotypes    // only for P0
        
        find_phenotypes();
        I.set_phenotypes(phenotype_1,phenotype_2);
        //I.display_phenotype1();
        //I.display_phenotype2();
        
        ME.place_individual_in_map(I);
        mutation_round++;
        if (g== mut_gen % 25){
            cout << "25% done" << endl;
        }
        if (g== mut_gen % 50){
            cout << "50% done" << endl;
        }
        if (g== mut_gen % 75){
            cout << "75% done" << endl;
        }
    }
    cout << "completed " << mutation_round << " MUTATION rounds" << endl;
    ME.how_many_full_bins();
    ME.print_contents_of_map();
    ME.best_fit_bin_genome();
    //ME.all_parents();
}
// --------------------------------------------------
// run best Individual to graph
void Wrapper::graph_final(){
    soln_set.empty();       // only for P0
    soln_set=ME.get_best_individual1();     // when implmented error.
    
    ofstream myfile;
    myfile.open ("Graph_approx_function.txt");
    for (int xx=0; xx<100; xx++){
        set_x_pos(xx);
        double in_x=get_x_pos();
        double out_y=approx_func();
        myfile << fixed << setprecision(8) << in_x << "\t\t" << out_y << "\n";
    }
    myfile.close();
}
// --------------------------------------------------
// Run single test, modifiable weights of NN, recieve fit_rating and phenotypes.
    // TODO - User inputs all weights.
    // TODO - Able to change how weights are generated, access to different random functions.
//void Wrapper::run_single_individual(){
//    build_real_set();
//    
//    Individual I;
//    I.set_individual_params(isize_1, isize_2, imutate_mag_1, imutate_mag_2, imutate_amount_1, imutate_amount_2);
//    //I.display_individual_params();
//    
//    I.build_individual();
//    I.display_individual1();
//    I.display_individual2();
//
//    empty_soln_set();
//    soln_set=I.genome1;
//    fitness_P0();
//    I.set_fit_rating(fit_rating);
//    I.display_fit_rating();
//    set_x_pos(10);   // at x=10, the value of a and b are the phenotypes
//    //find_phenotypes();
//    I.set_phenotypes(75,5);
//    I.display_phenotype1();
//    I.display_phenotype2();
//    
//    ME.place_individual_in_map(I);
//}
    // fitness function from Sim
    // I.set_fit_rating( double );
    // get phenotypes from Sim
    // I.set_phenotypes( double , double );
    // ME.place_individual_in_map(I);

// --------------------------------------------------

                    // TODO - functions
// choose bin from map and test function
//




#endif /* Wrapper_hpp */
