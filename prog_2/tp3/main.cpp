#include "class_point.hpp"
#include "class_node.hpp"
#include "class_element.hpp"
#include "class_mesh.hpp"
#include <iostream>


// vertex = node<point2d>


int main(){
    vertex* a = new vertex(point2d(1.,1.));
    vertex* b = new vertex(point2d(2.,2.));
    vertex* c = new vertex(point2d(2.,0.));
    vertex* d = new vertex(point2d(3.,1.));
    vertex* e = new vertex(point2d(3.,3.));

    triangle t1(*a,*b,*c);
    triangle t2(*b,*c,*d);
    triangle t3(*b,*d,*e);

    
    int ind=1;
    t1.indexing(ind);

    std::cout << "Indices after indexing t1: (4 is expected)" << ind << std::endl;
    t1.print_vertices_indices();
    t1.reset_indices();
    t1.print_vertices_indices();

    triangulation mymesh(t1);  // maillage avec t1 comme premier élément
    mymesh.append(t2);
    mymesh.append(t3);

    mymesh.indexing();
    std::cout << "indices t1: " << std::endl;
    t1.print_vertices_indices();
    std::cout << "indices t2: " << std::endl;
    t2.print_vertices_indices();
    std::cout << "indices t3: " << std::endl;
    t3.print_vertices_indices();
    


    return 0;
}