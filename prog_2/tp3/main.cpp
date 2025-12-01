#include "class_point.hpp"
#include "class_node.hpp"
#include "class_element.hpp"
#include <iostream>







int main(){
    point2d pt(1.0, 2.0);
    node<point2d> n1(pt, 1, 0, 0);
    node<point2d> n2;

    std::cout << "Node n1: " << n2 << std::endl;
    n2.print();
    std::cout << "Sharing elements: " << n2.sharing_elements() << std::endl;
    n2.more_sharing_elements();
    std::cout << "After increment, sharing elements: " << n2.sharing_elements() << std::endl;



    return 0;
}