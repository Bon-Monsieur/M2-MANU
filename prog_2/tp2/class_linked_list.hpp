#pragma once
#include <iostream>
#include <cassert>
#include <cmath>

template<typename T>
class linked_list {
    protected:
        T item_;
        linked_list<T>* p_next_;

    public:
        // Constructors
        linked_list();
        linked_list(T const& t, linked_list* pN=nullptr): item_(t),p_next_(pN){};
        linked_list (linked_list const& L ) : item_(L.item()), p_next_(L.p_next() ? new linked_list(*L.p_next()) : nullptr) {};
        // Destructor
        ~linked_list(){
            delete p_next_; // détruit récursivement la liste
            p_next_ = nullptr;  // assigne nullptr au pointeur après suppression
        };

        // Getters
        linked_list* p_next() const { return p_next_; }
        linked_list*& p_next() { return p_next_; }
        T const& item() const { return item_; }
        T& item() { return item_;}

        //operator
        linked_list& last();

       
        


};


template<typename T>
linked_list<T>::linked_list(){
    p_next_ = nullptr;
}

template<typename T>
linked_list<T>& linked_list<T>::last(){
    linked_list<T>* current = this;
    while (current->p_next() != nullptr) {
        current = current->p_next();
    }
    return *current;
}