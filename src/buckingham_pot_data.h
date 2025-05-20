#include "utility.h"
#include <map>
using namespace std;


map<string,vector<float>> buckingham_pot_data {
    // data used in pauli & london potentials (together called a buckingham potential)
    {"H",  {432,   4.52, 1.96 }},
    {"C",  {34000, 4.59, 15.7 }},
    {"O",  {5600,  4.59, 9.41 }},
    {"N",  {9000,  4.59, 12.7 }},
    {"F",  {2610,  4.6,  5.1  }},
    {"S",  {53526, 3.99, 70.6 }},
    {"Tl", {221248.6, 4.01164, 734.245}},
    {"In", {221248.6, 4.01164, 734.245}}
};