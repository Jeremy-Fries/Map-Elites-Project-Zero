//
//  Map_Elites.hpp
//  Mapping_Elites v1
//
//  Created by Jeremy Fries on 10/20/15.
//  Copyright © 2015 Jeremy Fries. All rights reserved.
//

class Map_Elites;

#ifndef Map_Elites_hpp
#define Map_Elites_hpp

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>

//#include "Wrapper.hpp"

///
#include "Individual.hpp"
#include "Map_space.hpp"
//#include "Map_Elites.hpp"

///
//#include "Simulator.hpp"
//#include "State.hpp"
//#include "dynamics.h"
//#include "anglemath.h"
//#include "Craft.hpp"
//
/////
//#include "Neural Network.hpp"
//#include "Layer.hpp"
//#include "Node.hpp"


using namespace std;

class Map_Elites{
    friend class Wrapper;
    friend class Map_space;
    
protected:
    int resolution1, resolution2;
    double dim1_min, dim2_min, dim3_min,
    dim1_max, dim2_max, dim3_max;
    int num_spacing1, num_spacing2, num_spacing3;
    double spacing1,  spacing2, spacing3;
    //Individual I;
    int fill_generation, mutate_generation;
    double out_scope1, out_scope2;
    //int p1_up, p1_down, p2_up, p2_down;
    //double p1_up_pos, p1_down_pos, p2_up_pos, p2_down_pos;
    //double closest_bin_dim1, closest_bin_dim2;
    //vector <int> closest_bin_vec;
    vector<Map_space> full_bins;
    vector<double> distance_between_centerbin_phenotype, fit_ratings_in_map, best_individual1,best_individual2;
    double center_bin1, center_bin2;
    double center_dist;
    
// --------------------------------------------------
public:
            // Min of dim1
    void set_min_dim1(double);
    double get_min_dim1();
    void display_min_dim1();
            // Min of dim2
    void set_min_dim2(double);
    double get_min_dim2();
    void display_min_dim2();
            // Min of dim3
    void set_min_dim3(double);
    double get_min_dim3();
    void display_min_dim3();
            // Max of dim1
    void set_max_dim1(double);
    double get_max_dim1();
    void display_max_dim1();
            // Max of dim2
    void set_max_dim2(double);
    double get_max_dim2();
    void display_max_dim2();
            // Max of dim3
    void set_max_dim3(double);
    double get_max_dim3();
    void display_max_dim3();
// --------------------------------------------------
            // Resolution
    void set_resolution(int,int);
    int get_resolution1();
    void display_resolution1();
    int get_resolution2();
    void display_resolution2();
// --------------------------------------------------
        // Fill Generation
    void set_fill_generation(int);
    int get_fill_generation();
    void display_fill_generation();
// --------------------------------------------------
            // Mutate Generation
    void set_mutate_generation(int);
    int get_mutate_generation();
    void display_mutate_generation();
// --------------------------------------------------
    vector<double>& get_best_individual1();
// --------------------------------------------------
            // Create full bin vector
    void create_full_bin();
    void find_center_bin(int p1,int p2);
    void find_pheno_dist_to_center_bin(int p1, int p2);
// --------------------------------------------------
            // To text file
    void print_contents_of_map();

    void best_fit_bin_genome();
    
    void all_parents();
    
    void print_individual(ofstream & file);
// --------------------------------------------------    
    // how many bins are full?
    void how_many_full_bins();
// --------------------------------------------------
            // Map Parameters
    void set_map_params(double dim1_min,double dim1_max,double dim2_min,double dim2_max,int resolution1,int resolution2, int fill_generation, int mutate_generation);
    void display_Map_params();
    void set_out_scope_magnitude();
// --------------------------------------------------
            // Makes Map
    void initialize_map();
    void place_individual_in_map(Individual Passed_I);
    void individual_from_map(int p1, int p2);
    vector<double>& get_temp_individual1();
    vector<double>& get_temp_individual2();
// --------------------------------------------------
private:
    vector<Map_space> Row;
    vector< vector<Map_space> > Map;
    
    vector<double> temp_individual1, temp_individual2;
        
    
};
// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
void Map_Elites::set_min_dim1(double min1){
    dim1_min = min1;
}
double Map_Elites::get_min_dim1(){
    return dim1_min;
}
void Map_Elites::display_min_dim1(){
    cout << endl << "Min of dim1 is: " << dim1_min << endl;
}
// Min of dim2
void Map_Elites::set_min_dim2(double min2){
    dim2_min = min2;
}
double Map_Elites::get_min_dim2(){
    return dim2_min;
}
void Map_Elites::display_min_dim2(){
    cout << endl << "Min of dim2 is: " << dim2_min << endl;
}
// Min of dim3
void Map_Elites::set_min_dim3(double min3){
    dim3_min = min3;
}
double Map_Elites::get_min_dim3(){
    return dim3_min;
}
void Map_Elites::display_min_dim3(){
    cout << endl << "Min of dim3 is: " << dim3_min << endl;
}
// Max of dim1
void Map_Elites::set_max_dim1(double max1){
    dim1_max = max1;
}
double Map_Elites::get_max_dim1(){
    return dim1_max;
}
void Map_Elites::display_max_dim1(){
    cout << endl << "Max of dim1 is: " << dim1_max << endl;
}
// Max of dim2
void Map_Elites::set_max_dim2(double max2){
    dim2_max = max2;
}
double Map_Elites::get_max_dim2(){
    return dim2_max;
}
void Map_Elites::display_max_dim2(){
    cout << endl << "Max of dim2 is: " << dim2_max << endl;
}
// Max of dim3
void Map_Elites::set_max_dim3(double max3){
    dim3_max = max3;
}
double Map_Elites::get_max_dim3(){
    return dim3_max;
}
void Map_Elites::display_max_dim3(){
    cout << endl << "Max of dim3 is: " << dim3_max << endl;
}
// --------------------------------------------------
// Sets resolution of Map, will determine size of Map_space and how many.
// Resolution
void Map_Elites::set_resolution(int r1, int r2) {
    resolution1=r1;
    resolution2=r2;
    num_spacing1=r1;
    num_spacing2=r2;
    //num_spacing3=r3;
    spacing1=(dim1_max-dim1_min)/resolution1;
    spacing2=(dim2_max-dim2_min)/resolution2;
    // spacing3=(dim3_max-dim3_min)/resolution3;
}
int Map_Elites::get_resolution1(){
    return resolution1;
}
void Map_Elites::display_resolution1(){
    cout << endl << "Resolution 1 is: " << resolution1 << endl;
}
int Map_Elites::get_resolution2(){
    return resolution2;
}
void Map_Elites::display_resolution2(){
    cout << endl << "Resolution 2 is: " << resolution2 << endl;
}
// --------------------------------------------------
// Fill Generation, fills Map with individuals for desired amount.
void Map_Elites::set_fill_generation(int f){
    fill_generation=f;
}
int Map_Elites::get_fill_generation(){
    return fill_generation;
}
void Map_Elites::display_fill_generation(){
    cout << endl << "Fill generation is: " << fill_generation << endl;
}
// --------------------------------------------------
// Mutate Generation, muatates individuals in Map for desired amount.
void Map_Elites::set_mutate_generation(int m){
    mutate_generation=m;
}
int Map_Elites::get_mutate_generation(){
    return mutate_generation;
}
void Map_Elites::display_mutate_generation(){
    cout << endl << "Mutate generation is: " << mutate_generation << endl;
}
// --------------------------------------------------
// Itialization function used outside of class.
// Set Map Parameters
void Map_Elites::set_map_params(double d1_min, double d1_max, double d2_min, double d2_max, int res1, int res2, int fill_gen, int mutate_gen){
    set_min_dim1(d1_min);
    set_max_dim1(d1_max);
    set_min_dim2(d2_min);
    set_max_dim2(d2_max);
    set_resolution(res1, res2);
    set_fill_generation(fill_gen);
    set_mutate_generation(mutate_gen);
    initialize_map();
    //set_out_scope_magnitude();
}
// --------------------------------------------------
// Display Map Parameters
void Map_Elites::display_Map_params(){
    cout << endl << "------------------------------- Map Parameters" << endl;
    display_min_dim1();
    display_max_dim1();
    display_min_dim2();
    display_max_dim2();
    display_resolution1();
    display_resolution2();
    display_fill_generation();
    display_mutate_generation();
    cout << endl << "-------------------------------" << endl;
}
// --------------------------------------------------
// TODO - Make dynamic
// Check out of scope magnitude
void Map_Elites::set_out_scope_magnitude(){
    // auto set scope magnitude relative to dims min and max
}
// --------------------------------------------------
// Builds Map based on parameters, ALSO creates and sets Map_space LB and UB
void Map_Elites::initialize_map(){
    double pre_LB1=dim1_min;
    double pre_LB2=dim2_min;
    //int row_in_Map=1;
    for(int d1=0; d1<num_spacing1; d1++){
        //int num_in_row=1;
        for(int d2=0; d2<num_spacing2; d2++){
            Map_space M;
            M.previous_genome1.clear();  // clears parent vector
            M.previous_genome2.clear();  // clears parent vector
            
            M.bin1=d1;
            M.bin2=d2;

            M.set_LB1(pre_LB1);

            M.set_LB2(pre_LB2);

            double calc_UB1= pre_LB1+spacing1;
            M.set_UB1(calc_UB1);
            
            double calc_UB2= pre_LB2+spacing2;
            M.set_UB2(calc_UB2);
            
            M.build_map_space();
            Row.push_back(M);
            pre_LB2+=spacing2;
            
            //cout << endl << "number in row is: "<< num_in_row << endl;
            //num_in_row++;
        }
        Map.push_back(Row);
        pre_LB1+=spacing1;
        
        //cout << endl << "rows in Map is: "<< row_in_Map << endl;
        //row_in_Map++;
    }
    cout << endl << "Map is made" << endl;
}
// --------------------------------------------------
void Map_Elites::print_individual(ofstream & file){


}
// --------------------------------------------------
// Places individual into corresponding map_space in Map
void Map_Elites::place_individual_in_map(Individual Passed_I){
    double p1=Passed_I.get_phenotype1();
    double p2=Passed_I.get_phenotype2();
    
    // out of bounds check
    double p1_lower=dim1_min;
    double p2_lower=dim2_min;
    double p1_upper=dim1_max;
    double p2_upper=dim2_max;
    
    double pLB1=dim1_min;
    double pLB2=dim2_min;
    int row_value=0;
    int element_value=0;
    
    // to find which row phenotype belongs in
    for(int d1=0; d1<num_spacing1; d1++){                                       // TODO - table flipped wrong direction, check initilaize() as well.
        // TODO - dex 0 and dex end special case handling in bin not this one
        double pUB1= pLB1+spacing1;
        // under bounds
        if (p1<p1_lower){
            row_value=0;
            //cout << endl << "under---------------" << endl;
            break;
        }
        else if (p1>=p1_upper){
            row_value=resolution1-1;
            //cout << endl << "over----------------" << endl;
            break;
        }
        // between range of LB1 and UB1
        else if (pLB1<p1 && p1<=pUB1){
            row_value=d1;
            //cout << endl << "right on the money!" << endl;
            break;
        }
        else {
            pLB1+=spacing1;
            //cout << endl << "moving up a lvl" << endl;
        }
    }
    // to find which element in row phenotype beolongs in
    for(int d2=0; d2<num_spacing2; d2++){
        double pUB2= pLB2+spacing2;
        // under bounds
        if (p2<p2_lower){
            element_value=0;
            break;
        }
        else if (p2>=p2_upper){
            element_value=resolution2-1;
            break;
        }
        // between range of LB2 and UB2
        if (pLB2<p2 && p2<=pUB2){
            element_value=d2;
            break;
        }
        else {
            pLB2+=spacing2;
        }
    }
    Map.at(row_value).at(element_value).current_individual.push_back(Passed_I);
    
    //cout << endl << "Placed in row "<< row_value << " column " <<element_value << endl;
    
    // compare new individual in map space and erase worse
    Map.at(row_value).at(element_value).compare_new_individual();
    
    // Print fitness_rating
    // Print Individual
}
// --------------------------------------------------
// Get individual from a map_space in Map
void Map_Elites::individual_from_map(int p1, int p2){           // Changed from double to int, 1/12 - JJ
    temp_individual1.clear();
    temp_individual2.clear();
    double pLB1=dim1_min;
    double pLB2=dim2_min;
    int row_value=0;
    int element_value=0;
    // to find which row phenotype belongs in
    for(int d1=0; d1<num_spacing1; d1++){
        double pUB1= pLB1+spacing1;
        // between range of LB1 and UB1
        if (pLB1<p1 && p1<pUB1){
            row_value=d1;
        }
        else {
            pLB1+=spacing1;
        }
    }
    // to find which element in row phenotype beolongs in
    for(int d2=0; d2<num_spacing2; d2++){
        double pUB2= pLB2+spacing2;
        // between range of LB2 and UB2
        if (pLB2<p2 && p2<pUB2){
            element_value=d2;
        }
        else {
            pLB2+=spacing2;
        }
    }
    if (Map.at(row_value).at(element_value).current_individual.size()>0){
        temp_individual1 = Map.at(row_value).at(element_value).current_individual.at(0).get_individual1();
        
        temp_individual2 = Map.at(row_value).at(element_value).current_individual.at(0).get_individual2();
        //cout << endl << "individuals passed to temp" << endl;
        
    }
    else {
        //cout << endl << "No individual in: " << p1 << "," << p2 << endl;
        find_center_bin(p1, p2);
        create_full_bin();
        
        int min_index = min_element(distance_between_centerbin_phenotype.begin(), distance_between_centerbin_phenotype.end()) - distance_between_centerbin_phenotype.begin();
        
        //cout << endl << "min_index is: " << min_index << endl;
        
        temp_individual1 = full_bins.at(min_index).current_individual.at(0).get_individual1();
        temp_individual2 = full_bins.at(min_index).current_individual.at(0).get_individual2();
        
        //full_bins.at(min_index).current_individual.at(0).display_individual1();
        //full_bins.at(min_index).current_individual.at(0).display_individual2();
        
        
        //int b1 = full_bins.at(min_index).bin1;
        //int b2 = full_bins.at(min_index).bin2;
        
        //cout << endl << "closest bin is: " << b1 << "," << b2 << endl;
    }
}
// Get temp_individuals
vector<double>& Map_Elites::get_temp_individual1(){
    return temp_individual1;
}
vector<double>& Map_Elites::get_temp_individual2(){
    return temp_individual2;
}
// --------------------------------------------------
void Map_Elites::create_full_bin(){
    full_bins.clear();
    distance_between_centerbin_phenotype.clear();
    
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){ // look for enum
                full_bins.push_back(Map.at(row_value).at(element_value));       // push_back Map_space that is full
                find_pheno_dist_to_center_bin(row_value, element_value);        // calculates distance from phenotype to center of desired bin
                distance_between_centerbin_phenotype.push_back(center_dist);
                Map.at(row_value).at(element_value).bin1=row_value;    // save row and element vaue to bin
                Map.at(row_value).at(element_value).bin2=element_value;    // save row and element vaue to bin
            }
        }
    }
    //cout << endl << "Full_bin vector has  " << full_bins.size() << endl;
    //cout << endl << "Distance vector has  " << distance_between_centerbin_phenotype.size() << endl;
}
// --------------------------------------------------
// find center of bin, will be implemented in individual_from_map
void Map_Elites::find_center_bin(int p1, int p2){
    center_bin1 = 0;
    center_bin2 = 0;
    
    double center1 = spacing1/2;
    double center2 = spacing2/2;
    
    double lb1 = Map.at(p1).at(p2).get_LB1();
    double lb2 = Map.at(p1).at(p2).get_LB2();

    
    //cout << endl << "lb1 is: " << lb1 << endl;
    //cout << endl << "lb2 is: " << lb2 << endl;
    
    center_bin1 = lb1+center1;                                
    center_bin2 = lb2+center2;
    //cout << endl << "center_bin1 is: " << center_bin1 << endl;
    //cout << endl << "center_bin2 is: " << center_bin2 << endl;
}
// --------------------------------------------------
// find distance from Individual in a bin to center of desired bin
void Map_Elites::find_pheno_dist_to_center_bin(int p1, int p2){
    center_dist = 0;
    
    double d1=Map.at(p1).at(p2).current_individual.at(0).get_phenotype1();
    double d2=Map.at(p1).at(p2).current_individual.at(0).get_phenotype2();
    //cout << endl << "d1 is: " << d1 << endl;
    //cout << endl << "d2 is: " << d2 << endl;

    double d1c = center_bin1;
    double d2c = center_bin2;
    //cout << endl << "d1c is: " << d1c << endl;
    //cout << endl << "d2c is: " << d2c << endl;
    
    center_dist = sqrt((d1*d1c) + (d2*d2c));
    //cout << endl << "--------------distance is:  " << center_dist << endl;
}
// --------------------------------------------------
vector<double>& Map_Elites::get_best_individual1(){
    return best_individual1;
}
// --------------------------------------------------
void Map_Elites::how_many_full_bins(){
    int num_of_full_bins=0;
    //full_bins.clear();
    //cout << endl << "num_spacing1 is: "<< num_spacing1 << endl;
    //cout << endl << "num_spacing2 is: "<< num_spacing2 << endl;
    
    for(int row_value=0; row_value<num_spacing1; row_value++){
        //cout << endl << "row is: "<< row_value << endl;
        
        for(int element_value=0; element_value<num_spacing2; element_value++){
            
            if (Map.at(row_value).at(element_value).full_bin_check()==0){
                num_of_full_bins++;
                //cout << endl << "number of full bins: " << num_of_full_bins << endl;
            }
            else{
                //cout << endl << "not full" << endl;
            }
            //cout << endl << "element is: "<< element_value << endl;
        }
    }
    cout << endl << "Final number of full bins: " << num_of_full_bins << endl;
    //cout << endl << "Full_bin vector has  " << full_bins.size() << endl;
}
// --------------------------------------------------
// Final print contents of map,
void Map_Elites::print_contents_of_map(){
    ofstream myfile;
    myfile.open ("best_fit_ratings.txt");
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){
                myfile << Map.at(row_value).at(element_value).current_individual.at(0).get_fit_rating() << '\n';
            }
        }
    }
    myfile.close();
    cout << endl << "txt file made." << endl;
}
// --------------------------------------------------
void Map_Elites::best_fit_bin_genome(){
    full_bins.clear();
    fit_ratings_in_map.clear();
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){
                full_bins.push_back(Map.at(row_value).at(element_value));       // push_back Map_space that is full
                fit_ratings_in_map.push_back(Map.at(row_value).at(element_value).get_best_fit());
            }
        }
    }
    
    // CHECK 
    int best_fit_index = max_element(fit_ratings_in_map.begin(), fit_ratings_in_map.end()) - fit_ratings_in_map.begin();
    
    ///////
    best_individual1=  full_bins.at(best_fit_index).current_individual.at(0).get_individual1();
    //best_individual2=  full_bins.at(best_fit_index).current_individual.at(0).get_individual2();
    
    full_bins.at(best_fit_index).current_individual.at(0).display_individual1();

    cout << endl << "best fit bin is: " << best_fit_index << endl;
    cout << endl << "best fit bin has been accessed: " << full_bins.at(best_fit_index).get_counter() << " times." << endl;
    cout << endl << "best fit has " << full_bins.at(best_fit_index).previous_fit_rating.size() << " different parents." << endl;
    cout << endl << "best fit has deleted " << full_bins.at(best_fit_index).old_counter << " past parents." << endl;
    cout << endl << "best fit has deleted " << full_bins.at(best_fit_index).new_counter << " potential parents." << endl;
    
    
    ofstream myfile;
    myfile.open ("best_parent_fit_ratings.txt");
    cout << endl << "full_bins size is:" << full_bins.size() << endl;

    for(int parent=0; parent < full_bins.at(best_fit_index).previous_fit_rating.size(); parent++){
        //cout << full_bins.at(best_fit_index).previous_fit_rating.at(parent) << endl;
        
        myfile << full_bins.at(best_fit_index).previous_fit_rating.at(parent) << '\n';
    }
    
    myfile.close();
    cout << endl << "Parent txt file made." << endl;
}
// --------------------------------------------------
//// Final print contents of map,
//void Map_Elites::print_contents_of_map(){
//    ofstream myfile;
//    myfile.open ("best_fit_ratings.txt");
//    for(int row_value=0; row_value<num_spacing1; row_value++){
//        for(int element_value=0; element_value<num_spacing2; element_value++){
//            if (Map.at(row_value).at(element_value).full_bin_check()==0){
//                myfile << Map.at(row_value).at(element_value).current_individual.at(0).get_fit_rating() << '\n';
//            }
//        }
//    }
//    myfile.close();
//    cout << endl << "txt file made." << endl;
//}
// --------------------------------------------------
void Map_Elites::all_parents(){
    for(int element=0; element<full_bins.size();element++){
        cout << "parent for bin: " << element << "\taccessed count: " << full_bins.at(element).counter << endl;
        for(int parent=0; parent < full_bins.at(element).previous_fit_rating.size(); parent++){
            cout << full_bins.at(element).previous_fit_rating.at(parent) << endl;
    
        }
    }
}
// --------------------------------------------------








// --------------------------------------------------









#endif /* Map_Elites_hpp */
