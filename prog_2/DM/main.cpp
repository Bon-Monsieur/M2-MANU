#include <iostream>
#include "class_row.hpp"
#include "header/class_dynamic_vector.hpp"
#include "class_sparse_matrix.hpp"

int main(){
    sparse_matrix<double> mat(5, 3.0);

    std::cout << mat.column_number() << std::endl; 
    std::cout << mat.row_number() << std::endl;    
    std::cout << mat.order() << std::endl;         
    
    
    return 0;
}
