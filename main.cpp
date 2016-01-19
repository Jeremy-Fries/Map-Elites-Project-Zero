//
//  main.cpp
//  Map_Elites_Testing_v1
//
//  Created by Jeremy Fries on 12/18/15.
//  Copyright © 2015 Jeremy Fries. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "Wrapper.hpp"
//#include "Individual.hpp"

// TODO - Heat map, Parent tracking, Phenotype functions

// Moment of landing
// Pheno1- x-speed
// Pheno2- orientation
// fitness- total kinetic energy
// very bad fitness for anything that doesn't land

int main() {
    srand(time(NULL));
    Wrapper W;
    
    W.initialize_wrapper();
    
    W.fill_MAP_P0();
    
    //W.run_single_individual();
    
    W.mutate_MAP();
    
    W.print_stuff();
    
    W.graph_final();
    
    return 0;
}
