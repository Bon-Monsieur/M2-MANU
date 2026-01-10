#pragma once
#include "header/class_linked_list.hpp"
#include "class_row_element.hpp"
#include "header/class_dynamic_vector.hpp"  



template<typename T>
class row : public linked_list<row_element<T>> {
    public:
        // constructors
        row(T const& val = 0, int col = -1) : linked_list<row_element<T>>(row_element<T>(val, col)){};
        
        inline row_element<T> const& operator()() const { return this->item_;}
        inline T const& value() const { return this->item_.get_value();}
        inline int column() const { return this->item_.get_column();}

        void insert_next_item(T const& val , int col );
        void insert_first_item(T const& val , int col );
        void append(T const& val , int col );

        const T row_sum() const;    // somme des éléments de la ligne

        const T operator[](int ii) const; // accès à l'élément de la colonne ii (0 si inexistant)

        const row& operator*=(T const& t );
        const row& operator/=(T const& t );

        // produit scalaire entre row et dynamic_vector pour la multiplication matricielle plus tard
        const T operator*(dynamic_vector<T> const& v) const; 

        void renumber_columns(dynamic_vector<int> const& renumber);

        template<typename S>
        friend const row<S> operator*(row<S> const& r , S const& t );
        template<typename S>
        friend const row<S> operator*(S const& t , row<S> const& r );
        template<typename S>
        friend const row<S> operator/(row<S> const& r , S const& t );

        void drop_items(dynamic_vector<int> const& mask);
};



template<typename T>
void row<T>::insert_next_item(T const& val , int col ){
    row_element<T> e(val, col);
    linked_list<row_element<T>>::insert_next_item(e);
}

template<typename T>
void row<T>::insert_first_item(T const& val , int col ){
    row_element<T> e(val, col);
    linked_list<row_element<T>>::insert_first_item(e);
}

template<typename T>
void row<T>::append(T const& val , int col ){
    row_element<T> e(val, col);
    linked_list<row_element<T>>::append(e);
}


template<typename T>
const T row<T>::row_sum() const {
    T sum = this->item_.get_value();

    // si y'a un élément 
    if (this->p_next_) {
        sum += ((row*)this->p_next_)->row_sum();
    }

    return sum;
}


template<typename T>
const T row<T>::operator[](int ii) const {
    T temp = this->column();
    if (temp == ii) {
        return this->value();
    } else if (temp > ii || this->p_next_ == nullptr) {
        return 0;
    }
    else{
        return ((row*)this->p_next_)->operator[](ii);
    }
    return 0;
}


template<typename T>
const row<T>& row<T>::operator*=(T const& t ){
    this->item_.set_value( this->item_.get_value() * t );

    if (this->p_next_) {
        ((row*)this->p_next_)->operator*=(t);
    }

    return *this;
}

template<typename T>
const row<T>& row<T>::operator/=(T const& t ){
    if (t == T{}) {
        throw std::runtime_error("Division by zero in row operator/=");
    }
    else{
        this->item_.set_value( this->item_.get_value() / t );

        if (this->p_next_) {
            ((row*)this->p_next_)->operator/=(t);
        }
    }
    return *this;
}


template<typename T>
const T row<T>::operator*(dynamic_vector<T> const& v) const {
    T res = T{};

    for (row<T> const* p = this; p != nullptr; p = (row<T> const*)p->p_next()) {
        res += p->value() * v[p->column()];
    }

    return res;
}


template<typename T>
void row<T>::renumber_columns(dynamic_vector<int> const& renumber) {
    // ancienne colonne
    int old_col = this->item_.get_column();

    // nouvelle colonne
    this->item_.set_column(renumber[old_col]);

    if (this->p_next_) {
        ((row*)this->p_next_)->renumber_columns(renumber);
    }
}


template<typename S>
const row<S> operator*(row<S> const& r , S const& t ){
    row<S> res(r);
    res *= t;
    return res;
}

template<typename S>
const row<S> operator*(S const& t , row<S> const& r ) {
    return r * t;
}

template<typename S>
const row<S> operator/(row<S> const& r , S const& t ) {
    row<S> res(r);
    res /= t;
    return res;
}


template<class T>
void row<T>::drop_items(dynamic_vector<int> const& mask)
{
    // Cas de la liste suivante
    if (this->p_next()) {
        row<T>* next = (row<T>*)this->p_next();
        if (!mask[next->column()]) {
            this->drop_next_item();
            drop_items(mask); // on continue sur le même élément
        } else {
            next->drop_items(mask); // on descend récursivement
        }
    }

    // Vérification du premier élément
    if (!mask[this->column()]) {
        this->drop_first_item();
    }
}
