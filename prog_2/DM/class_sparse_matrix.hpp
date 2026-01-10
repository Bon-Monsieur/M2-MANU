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

        const sparse_matrix& operator+=(sparse_matrix<T> const&);
        const sparse_matrix& operator-=(sparse_matrix<T> const&);
        const sparse_matrix& operator*=(T const&);

        // Opérateurs matriciels classiques
        template<class S>
        friend const sparse_matrix<S> operator+(sparse_matrix<S> const& M1, sparse_matrix<S> const& M2);
        template<class S>
        friend const sparse_matrix<S> operator-(sparse_matrix<S> const& M1, sparse_matrix<S> const& M2);
        template<class S>
        friend const sparse_matrix<S> operator*(sparse_matrix<S> const& M, const S& t );
        template<class S>
        friend const sparse_matrix<S> operator*(const S& t , sparse_matrix<S> const& M);
        template<class S>
        friend const dynamic_vector<S> operator*(sparse_matrix<S> const& M, dynamic_vector<S> const& v);
        template<class S>
        friend const sparse_matrix<S> operator*(sparse_matrix<S> const&, sparse_matrix<S> const&);
        template<class S>
        friend const sparse_matrix<S> transpose(sparse_matrix<S> const& M);
        
        template<class S>
        friend void print(sparse_matrix<S> const& M);
        
        template<class S>
        friend void printf(sparse_matrix<S> const& M, std::ofstream& os);

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
    row<T>* r = this->item(i);
    if (!r) return T(0);

    T val = (*r)[j];
    if (val != T(0)){
        return val;
    }
    throw std::runtime_error("Element nul dans la matrice creuse");
}


template<class T>
int sparse_matrix<T>::row_number() const {
    return this->number();
}

template<class T>
int sparse_matrix<T>::column_number() const {
    int max_col = -1;

    for (int ii = 0; ii < this->number(); ++ii) {
        row<T>* r = this->item(ii);
        if (!r) continue;

        for (row<T> const* p = r; p != nullptr; p = (row<T> const*)p->p_next()) {
            if (p->column() > max_col){
                max_col = p->column();
            }
        }
    }
    return max_col + 1;
}



template<class T>   
int sparse_matrix<T>::order () const {
    return max(this->row_number(), this->column_number());
}


template<class T>
const sparse_matrix<T>& sparse_matrix<T>::operator+=(sparse_matrix<T> const& other){
    if (row_number() != other.row_number())
        throw std::runtime_error("Nombre de lignes différent");

    for (int ii = 0; ii < row_number(); ++ii){
        if (!this->item(ii))
            this->item(ii) = new row<T>();

        *(this->item(ii)) += *(other.item(ii));
    }
    return *this;
}

template<class T>
const sparse_matrix<T>& sparse_matrix<T>::operator-=(sparse_matrix<T> const& other){
    if (row_number() != other.row_number())
        throw std::runtime_error("Nombre de lignes différent");

    sparse_matrix<T> temp = other;
    temp*=-1;
    *(this) += temp;
    
    return *this;
}

template<class T>
const sparse_matrix<T>& sparse_matrix<T>::operator*=(T const& scalar){
    for (int i = 0; i < row_number(); ++i){
        if (this->item(i))
            *(this->item(i)) *= scalar;
    }
    return *this;
}



template<class S>
const sparse_matrix<S> operator+(sparse_matrix<S> const& M1, sparse_matrix<S> const& M2){
    sparse_matrix<S> temp = M1;
    temp += M2;
    return temp;
}

template<class S>
const sparse_matrix<S> operator-(sparse_matrix<S> const& M1, sparse_matrix<S> const& M2){
    sparse_matrix<S> temp = M1;
    temp -= M2;
    return temp;
}

template<class S>
const sparse_matrix<S> operator*(sparse_matrix<S> const& M, const S& t ){
    sparse_matrix<S> temp = M;
    temp *= t;
    return temp;
}

template<class S>
const sparse_matrix<S> operator*(const S& t , sparse_matrix<S> const& M){
    sparse_matrix<S> temp = M;
    temp *= t;
    return temp;
}


template<class S>
const dynamic_vector<S> operator*(sparse_matrix<S> const& M, dynamic_vector<S> const& v)
{
    if (M.column_number() != v.size())
        throw std::runtime_error("Dimensions incompatibles M x v");

    dynamic_vector<S> result(M.row_number(), S(0));

    for (int ii = 0; ii < M.row_number(); ++ii) {
        row<S>* r = M.item(ii);
        if (!r) continue;

        for (row<S> const* p = r; p != nullptr; p = (row<S> const*)p->p_next()) {
            size_t col = p->item().get_column();
            S val   = p->item().get_value();
            result(ii) += val * v[col];
        }
    }

    return result;
}


template<class S>
const sparse_matrix<S> operator*(sparse_matrix<S> const& B, sparse_matrix<S> const& A)
{
    if (B.column_number() != A.row_number()) {
        throw std::runtime_error("Dimensions incompatibles pour le produit matriciel");
    }

    sparse_matrix<S> result(B.row_number());

    for (int i = 0; i < B.row_number(); ++i) {

        row<S>* res_row = nullptr;

        row<S> const* b_row = B.item(i);
        if (!b_row) continue;

        for (row<S> const* p = b_row; p != nullptr;
             p = (row<S> const*) p->p_next()) {

            int j  = p->column();   // colonne j de B
            S coef = p->value();    // B(i,j)

            row<S>* a_row = A.item(j);
            if (!a_row) continue;

            row<S> temp = *a_row;  // copie de A_j
            temp *= coef;          // coef * A_j

            if (!res_row) {
                res_row = new row<S>(temp);  // première contribution
            } else {
                *res_row += temp;
            }
        }

        result.item(i) = res_row;
    }

    return result;
}




template<class S>
const sparse_matrix<S> transpose(sparse_matrix<S> const& M){
    sparse_matrix<S> T(M.column_number());

    for (int i = 0; i < M.row_number(); ++i){
        row<S> const* r = M.item(i);
        
        for (row<S> const* p = r; p != nullptr; p = (row<S> const*) p->p_next()){
            int j = p->column();
            S v = p->value();

            if (!T.item(j))
                T.item(j) = new row<S>(v, i);
            else
                T.item(j)->append(v, i);
        }
    }
    return T;
}


template<class S>
void print(sparse_matrix<S> const& M){
    std::cout << "Sparse Matrix (" << M.row_number() << " x " << M.column_number() << "):" << std::endl;
    
    for (int i = 0; i < M.row_number(); ++i){
        row<S> const* r = M.item(i);
        while (r){
            int col = r->column();
            S val = r->value();
            
            std::cout << "M[" << i << "][" << col << "] = " << val << std::endl;
            
            r = (row<S> const*)r->p_next();
        }
    }
}


template<class S>
void printf(sparse_matrix<S> const& M, std::ofstream& os){
    for (int i = 0; i < M.row_number(); ++i){
        row<S> const* r = M.item(i);
        while (r){
            int col = r->column();
            S val = r->value();
            
            os << i << " " << col << " " << val << std::endl;
            
            r = (row<S> const*)r->p_next();
        }
    }
}
