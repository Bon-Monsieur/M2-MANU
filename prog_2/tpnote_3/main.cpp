#include "class_element.hpp"
#include "class_mesh.hpp"
#include "class_node.hpp"
#include "class_point.hpp"

#include <iostream>

int main()
{
    std::vector<vertex*> vertices;
    std::vector<triangle> triangles;

    // loads mesh data from .mesh file
    std::string meshFile = "first.mesh";
    //std::string meshFile = "ell.mesh";
    //std::string meshFile = "ex.mesh";
    bool res               = mesh_reader(meshFile, vertices, triangles);

    std::size_t n_vertices = vertices.size();
    std::size_t n_elements = triangles.size();

    // form the mesh
    triangulation ell(triangles[0]); // first item will be doubled

    // populate mesh with triangles and free copied triangles
    for (auto new_triangle : triangles)
    {
        ell.append(new_triangle);
        new_triangle.~triangle();
    }

    ell.drop_first_item(); // first item is removed

    std::cout << "index elements in mesh" << std::endl;
    int n_indexed_elements = ell.indexing_elements();
    std::cout << "number of indexed elements in mesh= " << n_indexed_elements << std::endl;
    ell.print_elements_indices();
    return 0;
}