#include "class_point.hpp"
#include "class_node.hpp"
#include "class_element.hpp"
#include <iostream>







int main(){
    point2d pt1(0.0,0.0);
    point2d pt2(0.0, 1.0);
    point2d pt3(1.0, 0.0);
    point2d pt4(1.0, 1.0);
    node<point2d> n1(pt1, 1, 0, 0);
    node<point2d> n2(pt2, 2, 0, 0);
    node<point2d> n3(pt3, 3, 0, 0);
    node<point2d> n4(pt4, 4, 0, 0);

    element<point2d,3> elem(n1, n2, n3);
    std::cout << "node loc " << elem << std::endl;
    



    return 0;
}