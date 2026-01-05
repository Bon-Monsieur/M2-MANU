#include <iostream>
#include "class_row_element.hpp"  
#include "class_row.hpp"

int main() {
    row<double> row1(1.0, 0);
    row1.append(3.0, 2);
    row1.append(5.0, 4);
    row1.append(7.0, 6);
    
    
    row<double> row2(2.0, 1);
    row2.append(4.0, 3);
    row2.append(6.0, 4);
    row2.append(8.0, 5);
    
    row1+=row2;
    print(row1);

    return 0;
}
