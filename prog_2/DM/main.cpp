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
    // Création d'une matrice 2x2 vide
    sparse_matrix<double> MatMain(2);

    // Ligne 0 : (1, 2)
    MatMain.item(0) = new row<double>(1.0, 0); // M(0,0) = 1
    MatMain.item(0)->append(2.0, 1);  // M(0,1) = 2

    // Ligne 1 : (3, 4)
    MatMain.item(1) = new row<double>(3.0, 0); // M(1,0) = 3
    MatMain.item(1)->append(4.0, 1);  // M(1,1) = 4

    // Affichage
    print(MatMain);


    // Test 4: operator() - element access
    std::cout << "Test 4: Element Access operator()" << std::endl;
    try {
        double val = mat2(0, 0);
        std::cout << "mat2(0,0) = " << val << std::endl;
    } catch (const std::runtime_error& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    // Test 5: Scalar multiplication (operator*)
    std::cout << "Test 5: Scalar Multiplication" << std::endl;
    sparse_matrix<double> mat4 = MatMain * 3.0;
    std::cout << "mat2 * 3.0:" << std::endl;
    print(mat4);
    std::cout << std::endl;

    // Test 6: Scalar multiplication (commutative)
    std::cout << "Test 6: Commutative Scalar Multiplication" << std::endl;
    sparse_matrix<double> mat5 = 2.0 * mat2;
    std::cout << "2.0 * mat2:" << std::endl;
    print(mat5);
    std::cout << std::endl;

    // Test 7: operator*= (in-place scalar multiplication)
    std::cout << "Test 7: In-place Scalar Multiplication (*=)" << std::endl;
    sparse_matrix<double> mat6 = mat2;
    mat6 *= 0.5;
    std::cout << "mat2 *= 0.5:" << std::endl;
    print(mat6);
    std::cout << std::endl;

    // Test 8: Matrix addition (operator+)
    std::cout << "Test 8: Matrix Addition" << std::endl;
    sparse_matrix<double> matA(3, 5.0);
    sparse_matrix<double> matB(3, 4.0);
    sparse_matrix<double> matSum = matA + matB;
    std::cout << "Matrix A (3x3, diag=1.0) + Matrix B (3x3, diag=2.0):" << std::endl;
    print(matSum);
    std::cout << std::endl;

    // Test 9: operator+= (in-place addition)
    std::cout << "Test 9: In-place Addition (+=)" << std::endl;
    sparse_matrix<double> matC = matA;
    matC += matB;
    std::cout << "matA += matB:" << std::endl;
    print(matC);
    std::cout << std::endl;

    // Test 10: Matrix subtraction (operator-)
    std::cout << "Test 10: Matrix Subtraction" << std::endl;
    sparse_matrix<double> matDiff = matB - matA;
    std::cout << "Matrix B - Matrix A:" << std::endl;
    print(matDiff);
    std::cout << std::endl;

    // Test 11: operator-= (in-place subtraction)
    std::cout << "Test 11: In-place Subtraction (-=)" << std::endl;
    sparse_matrix<double> matD = matB;
    matD -= matA;
    std::cout << "matB -= matA:" << std::endl;
    print(matD);
    std::cout << std::endl;

    // Test 12: Matrix-Vector multiplication
    std::cout << "Test 12: Matrix-Vector Multiplication" << std::endl;
    sparse_matrix<double> matMV(3, 2.0);
    dynamic_vector<double> vec(2, 1.0);
    dynamic_vector<double> result = MatMain * vec;
    std::cout << "Matrix (3x3, diag=2.0) * Vector (size 3, all 1.0):" << std::endl;
    std::cout << "Result size: " << result.size() << std::endl;
    std::cout << result;
    std::cout << std::endl;

    // Test 13: Matrix-Matrix multiplication
    std::cout << "Test 13: Matrix-Matrix Multiplication" << std::endl;
    sparse_matrix<double> matX(2, 6.0);
    sparse_matrix<double> matY(2, 2.0);
    
    sparse_matrix<double> matProd = MatMain * matY;
    std::cout << "Matrix MatMain (2x2, 1,2,3,4) * Matrix Y (2x2, diag=2.0):" << std::endl;
    print(matProd);
    std::cout << std::endl;

    // Test 14: Transpose    
    std::cout << "Original matrix: MatMain" << std::endl;
    print(MatMain);
    std::cout << "\nTransposed matrix:" << std::endl;
    sparse_matrix<double> matTransposed = transpose(MatMain);
    print(matTransposed);
    std::cout << std::endl;

    // Test 15: Diagonal extraction
    std::cout << "Test 15: Diagonal Extraction" << std::endl;
    sparse_matrix<double> matDiagonal = diagonal(MatMain);
    std::cout << "Diagonal of MatMain:" << std::endl;
    print(matDiagonal);
    std::cout << std::endl;

/*
    // Test 16: printf function (write to file)
    std::cout << "Test 16: Printf Function (write to file)" << std::endl;
    std::ofstream outfile("sparse_matrix_output.txt");
    if (outfile.is_open()) {
        std::cout << "Writing matrix to 'sparse_matrix_output.txt'..." << std::endl;
        printf(mat2, outfile);
        outfile.close();
        std::cout << "File written successfully!" << std::endl;
    } else {
        std::cout << "Failed to open output file!" << std::endl;
    }
    std::cout << std::endl;
*/
    std::cout << "==== All Tests Completed ====" << std::endl;
    
    return 0;
}
