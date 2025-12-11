#pragma once
#include "header/class_linked_list.hpp"
#include "class_row_element.hpp"




template<typename T>
class row : public linked_list<row_element<T>> {
    public:
        // constructors
        row(T const& val, int col) : linked_list<row_element<T>>(row_element<T>(val, col)) {};
        
        // get first element
        inline row_element<T> const& operator()() const { return this->item();}
        inline T const& value() const { return this->item().get_value(); }

};
