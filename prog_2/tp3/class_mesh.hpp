#pragma once
#include "class_element.hpp"
#include "../tp2/class_linked_list.hpp"
#include "class_node.hpp"
#include "class_point.hpp"
#include <vector>
#include <string>
#include <fstream> // pour ouvrir les fichiers
#include <array>

struct coords{
    double x;
    double y;
    int index;
};
typedef std::array<int,3> itriplet;

template<class T>
class mesh : public linked_list<T>{
    protected:
        // Rien pour l'instant
    public:
        mesh(T& e):linked_list<T>(e){}   // constructor

        int indexing();   // indexing the nodes in the mesh
        void build_edges(std::vector<edge>& edges, std::size_t n_vertices);
}; 

typedef mesh<triangle> triangulation;



template<class T>
int mesh<T>::indexing(){
    for(mesh<T>*p_scan = this ; p_scan; p_scan=(mesh<T>*)p_scan->p_next())  // on reset tous les indices
        p_scan->item().reset_indices();
    int count=1;
    for(mesh<T>*p_scan = this ; p_scan; p_scan=(mesh<T>*)p_scan->p_next())  // on indexe les indices dans l'ordre ensuite
        p_scan->item().indexing(count);
    return count;
} 

bool mesh_reader_(std::string const& fname, std::vector<vertex*>& vertices , std::vector<triangle>& triangles ){
    std::ifstream fichier(fname);  // ouverture du fichier

    if (!fichier.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier !" << std::endl;
        return false;
    }

    std::string marker;

    while(marker.compare("Vertices")){
        fichier >> marker;
        assert(!fichier.eof());
    }
    int n_nodes = 0;
    fichier >> n_nodes;

    coords node;
    std::vector<point2d> points;

    for (auto in = 0; in <n_nodes; in++){
        fichier >> node.x >> node.y >> node.index;
        points.push_back(point2d(node.x, node.y));
    }   

    // Form vertices from points
    for (point2d pt : points){
        vertices.push_back( new vertex(pt));
    }
    // On assigne les indices aux vertices
    for (auto iv=0; iv<n_nodes; iv++){
        vertices[iv]->index() = iv;
    }

    while(marker.compare("Triangles")){
        fichier >> marker;
        assert(!fichier.eof());
    }

    int n_elements;
    fichier >> n_elements; 
    int n1; // useless label for the time being

    std::vector<itriplet> triplets;

    for (auto ie=0;ie<n_elements;ie++){
        itriplet tmp;
        fichier >> tmp[0] >> tmp[1] >> tmp[2];

        tmp[0]--;
        tmp[1]--;
        tmp[2]--;

        triplets.push_back(tmp);
    }

    // form triangles from triplets
    for (itriplet tri : triplets){
        triangle t( *vertices[tri[0]] , *vertices[tri[1]] , *vertices[tri[2]] );
        triangles.push_back(t);
    }

    while(marker.compare("End")){
        fichier >> marker;
        assert(!fichier.eof());
    }

    fichier.close(); // Finished reading from tje file stream
    std::cout << "mesh_reader : mesh file loaded successfully." << std::endl;

    return true;
}



template<class T>
void mesh<T>::build_edges(std::vector<edge>& edges, std::size_t n_vertices){
    int n_elements = this->length();
    std::cout << "n_vertices=" << n_vertices << std::endl;

    //data structure to store vertices-to-edges data
    list<linked_list<size_t>> hashv(n_vertices); // ve_data[i] is the linked list of edges connected to vertex i

    std::vector<std::size_t> adj_vertices(n_vertices,0); // to store adjacent vertex index when looking for existing edges
    std::size_t n_edge=0;

    // pointers to iterate through linked_list<size_t>
    linked_list<size_t>* p_cell = nullptr;
    linked_list<size_t>* p_cell_prev = nullptr;

    // pointers to vertices
    vertex* scan_min = nullptr;
    vertex* scan_max = nullptr;

    size_t ismin;
    size_t ismax;

    //looping on elements->next
    for (mesh<T>* scanner = this; scanner; scanner = (mesh<T>*)scanner->p_next()){
        for(auto kk=0; kk<3; kk++){
            size_t next = (kk+1)%3;
            if (scanner->item().vertex(kk)->index() < scanner->item().vertex(next)->index()){
                scan_min = scanner->item().vertex(kk);
                scan_max = scanner->item().vertex(next);
                ismin = scan_min->index();
                ismax = scan_max->index();
            }
            else{
                scan_min = scanner->item().vertex(next);
                scan_max = scanner->item().vertex(kk);
                ismin = scan_max->index();
                ismax = scan_min->index();
            }

            if (!hashv.item(ismin)){
                hashv.item(ismin) = new linked_list<size_t>(ismin);
            }
            p_cell_prev = hashv.item(ismin);
            p_cell = hashv.item(ismin)->p_next();

            while(p_cell && (p_cell->item() <= ismax)){
                p_cell_prev = p_cell;
                p_cell = p_cell->p_next();
            }

            if (p_cell_prev->item() < ismax){ // edge not found, we create it
                p_cell_prev->insert_next_item(ismax);

                ++adj_vertices[ismin];
                ++n_edge;

                edge.push_back( edge(*scan_min, *scan_max) );
            }
        }

    }
}


// ecrire fonction qui créer un fichier avec les coordonnées des sommets version brut-force