#pragma once
#include "header/class_list.hpp"
#include "class_row.hpp"

// Pour lire un fichier Harwell-Boeing
#include <fstream>
#include <sstream>
#include <stdexcept>


template<class T>
class sparse_matrix : public list<row<T>>{
    public:
        sparse_matrix(int Nrows=0) : list<row<T>>(Nrows){}; // constructeur par défaut
        sparse_matrix(int Nrows, T const& a);   // initialisation de la diagonale
        
        // construceur à parrtir d'un fichier comme Harwell-Boeing
        sparse_matrix(std::string const& filename);

        const T operator()(int i, int j);
        int row_number() const;
        int column_number() const;
        int order () const;


};


template<class T>
sparse_matrix<T>::sparse_matrix(int Nrows, T const& a) : list<row<T>>(Nrows)
{
    for (auto i=0; i < Nrows; ++i)
    {
        this->item(i) = new row<T>(a, i); // initialisation de la diagonale
    }
}

template<class T>
sparse_matrix<T>::sparse_matrix(std::string const& filename) : list<row<T>>(0)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier" + filename);
    }

    int Nrows, Ncols, nnz;
    file >> Nrows >> Ncols >> nnz;

    // On crée Nrows lignes vides
    for (int i = 0; i < Nrows; ++i) {
        this->append(row<T>());  // ajoute une row vide
    }

    // Lire les coefficients
    int row_idx, col_idx;
    T value;
    for (int k = 0; k < nnz; ++k) {
        file >> row_idx >> col_idx >> value;

        // Les indices en Harwell-Boeing commencent à 1 donc on les ajuste
        row_idx -= 1;
        col_idx -= 1;

        // On récupère la row correspondante et l'ajoute
        row<T>& r = this->item(row_idx);
        r.append(value, col_idx);
    }

    file.close();
}


template<class T>
const T sparse_matrix<T>::operator()(int i, int j){
    row<T>& r = this->item(i);
    if (r[j] != 0){
        return r[j];
    } else {
        throw std::runtime_error("L'élément est nul dans la matrice creuse.");
    }
}

template<class T>
int sparse_matrix<T>::row_number() const {
    return this->number();
}

template<class T>
int sparse_matrix<T>::column_number() const {
    int max_col = -1;

    for (int i = 0; i < this->number(); ++i) {

        row<T> const* p = this->item(i);
        if (!p) continue;   // ligne vide

        // Parcours récursif/itératif de la row
        while (p) {
            int col = p->column();
            if (col > max_col){
                max_col = col;
            }
            p = (row<T> const*)p->p_next();
        }
    }

    return max_col + 1; // indices commencent à 0
}



template<class T>   
int sparse_matrix<T>::order () const {
    return max(this->row_number(), this->column_number());
}