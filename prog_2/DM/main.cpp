#include <iostream>
#include "class_row_element.hpp"  
#include "class_row.hpp"

int main() {
    row_element<double> re1(3.0, 1);
    row_element<double> re2(2.0, 2);
    std::cout << re1;
    
    return 0;
}
