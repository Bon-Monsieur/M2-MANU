#include <iostream>
#include "class_row_element.hpp"  
#include "class_row.hpp"

int main() {
    row<double> row1(5.0, 0);
    row1.insert_first_item(3.0, 1);
    row1.insert_next_item(1.0, 2);
    row1.append(4.0, 3);
    
    
    for (linked_list<row_element<double>>* p = &row1; p != nullptr; p = p->p_next()) {
        std::cout << p->item() << std::endl; // utilise l'opérateur<< de row_element pour l'affichage
    }
    std::cout << "Row sum: " << row1.row_sum() << std::endl;

    return 0;
}
