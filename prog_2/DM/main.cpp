#include <iostream>
#include "class_row_element.hpp"  
#include "class_row.hpp"

int main() {
    row<double> row1(5.0, 2);
    row1.insert_first_item(3.0, 0);
    row1.insert_next_item(1.0, 1);
    row1.append(4.0, 3);
    

    print(row1);
    for (int i = 0; i <= 4; ++i) {
        std::cout << "col " << i << " = " << row1[i] << std::endl;
    }
    
    

    return 0;
}
