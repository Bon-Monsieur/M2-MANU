#pragma once
#include "header/class_linked_list.hpp"
#include "class_row_element.hpp"




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

        const T row_sum() const;

        const T operator[](int ii) const;
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