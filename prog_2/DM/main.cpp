#include <iostream>
#include "class_row_element.hpp"  
#include "class_row.hpp"

int main() {
    row<double> myrow(5.0, 2); // Create a row with initial value 5.0 at column 2
    row_element<double> re1(3.0, 1);
    myrow.insert_next_item(re1); // Insert another row element
    std::cout << "first value:" << myrow.value() << std::endl;
    return 0;
}
