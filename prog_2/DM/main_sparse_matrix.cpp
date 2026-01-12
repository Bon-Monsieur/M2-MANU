#include <iostream>
#include <fstream>
#include "class_row.hpp"
#include "header/class_dynamic_vector.hpp"
#include "class_sparse_matrix.hpp"

int main(){
    std::cout << "==== Testing sparse_matrix Class ====" << std::endl << std::endl;

    // Test 1: Default constructor
    std::cout << "Test 1: Default Constructor" << std::endl;
    sparse_matrix<double> mat1(3);
    std::cout << "Created empty matrix: " << mat1.row_number() << " x " << mat1.column_number() << std::endl;
    std::cout << "Order: " << mat1.order() << std::endl << std::endl;

    // Test 2: Diagonal constructor
    std::cout << "Test 2: Diagonal Constructor" << std::endl;
    sparse_matrix<double> mat2(4, 2.5);
    std::cout << "Created 4x4 matrix with diagonal = 2.5" << std::endl;
    std::cout << "Rows: " << mat2.row_number() << ", Columns: " << mat2.column_number() << ", Order: " << mat2.order() << std::endl;
    std::cout << "Matrix 2:" << std::endl;
    print(mat2);
    std::cout << std::endl;


    // TEST 3 : MATRICE A LA MAIN
    std::cout << "Test 3: Construction matrices creuses à la main" << std::endl;
    // [[1, 2, 0],
    //  [0, 0, 3],
    //  [0, 4, 0]]
    sparse_matrix<double> MatMain1(3);
    MatMain1.item(0) = new row<double>(1.0, 0); // M(0,0) = 1
    MatMain1.item(0)->append(2.0, 1);  // M(0,1) = 2
    MatMain1.item(1) = new row<double>(3.0, 2); // M(1,0) = 3
    MatMain1.item(2)  = new row<double>(4.0, 1);  // M(1,1) = 4
    std::cout << "MatMain1:" << std::endl;
    print(MatMain1);
    std::cout << std::endl;

    // [[0,0,5],
    //  [6,0,0],
    //  [0,7,0]]
    sparse_matrix<double> MatMain2(3);
    MatMain2.item(0) = new row<double>(5.0, 2); // M(0,2) = 5
    MatMain2.item(1) = new row<double>(6.0, 0); // M(1,0) = 6
    MatMain2.item(2) = new row<double>(7.0, 1); // M(2,1) = 7
    std::cout << "MatMain2:" << std::endl;
    print(MatMain2);
    std::cout << std::endl;

    // Test 5: Scalar multiplication (operator*)
    std::cout << "Test 5: Scalar Multiplication" << std::endl;
    sparse_matrix<double> mat4 = MatMain1 * 3.0;
    std::cout << "MatMain1 * 3.0:" << std::endl;
    print(mat4);
    std::cout << std::endl;

    // Test 7: operator*= (in-place scalar multiplication)
    std::cout << "Test 7: In-place Scalar Multiplication (*=)" << std::endl;
    sparse_matrix<double> mat6 = MatMain2;
    mat6 *= 0.5;
    std::cout << "MatMain2 *= 0.5:" << std::endl;
    print(mat6);
    std::cout << std::endl;

    // Test 8: Matrix addition (operator+)
    std::cout << "Test 8: Matrix Addition" << std::endl;
    sparse_matrix<double> matSum = MatMain1 + MatMain2;
    std::cout << "MatMain1 + MatMain2:" << std::endl;
    print(matSum);
    std::cout << std::endl;

    // Test 9: operator+= (in-place addition)
    std::cout << "Test 9: In-place Addition (+=)" << std::endl;
    sparse_matrix<double> matC = MatMain1;
    matC += MatMain2;
    std::cout << "MatMain1 += MatMain2:" << std::endl;
    print(matC);
    std::cout << std::endl;

    // Test 10: Matrix subtraction (operator-)
    std::cout << "Test 10: Matrix Subtraction" << std::endl;
    sparse_matrix<double> matDiff = MatMain1 - MatMain2;
    std::cout << "MatMain1 - MatMain2:" << std::endl;
    print(matDiff);
    std::cout << std::endl;

    // Test 11: operator-= (in-place subtraction)
    std::cout << "Test 11: In-place Subtraction (-=)" << std::endl;
    sparse_matrix<double> matD = MatMain1;
    matD -= MatMain2;
    std::cout << "MatMain1 -= MatMain2:" << std::endl;
    print(matD);
    std::cout << std::endl;

    // Test 12: Matrix-Vector multiplication
    std::cout << "Test 12: Matrix-Vector Multiplication" << std::endl;
    dynamic_vector<double> vec(3, 1.0);
    dynamic_vector<double> result = MatMain1 * vec;
    std::cout << "MatMain1 * (1,1,1)" << std::endl;
    std::cout << "Result size: " << result.size() << std::endl;
    std::cout << result;
    std::cout << std::endl;

    // Test 13: Matrix-Matrix multiplication
    std::cout << "Test 13: Matrix-Matrix Multiplication" << std::endl;    
    sparse_matrix<double> matProd = MatMain1 * MatMain2;
    std::cout << "MatMain1 * MatMain2:" << std::endl;
    print(matProd);
    std::cout << std::endl;

    // Test 14: Transpose    
    std::cout << "Original matrix: MatMain1" << std::endl;
    print(MatMain1);
    std::cout << "\nTransposed matrix:" << std::endl;
    sparse_matrix<double> matTransposed = transpose(MatMain1);
    print(matTransposed);
    std::cout << std::endl;

    // Test 15: Diagonal extraction
    std::cout << "Test 15: Diagonal Extraction" << std::endl;
    sparse_matrix<double> matDiagonal = diagonal(MatMain1);
    std::cout << "Diagonal of MatMain1:" << std::endl;
    print(matDiagonal);
    std::cout << std::endl;


    // Test 16: printf function (write to file)
    std::cout << "Test 16: Printf Function (write to file)" << std::endl;
    std::ofstream outfile("sparse_matrix_output.txt");
    if (outfile.is_open()) {
        std::cout << "Writing matrix to 'sparse_matrix_output.txt'..." << std::endl;
        printf(MatMain1, outfile);
        outfile.close();
        std::cout << "File written successfully!" << std::endl;
    } else {
        std::cout << "Failed to open output file!" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "==== All Tests Completed ====" << std::endl;
    
    return 0;
}
