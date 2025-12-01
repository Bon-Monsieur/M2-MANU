#pragma once
#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>

template<typename T>
class linked_list {
    protected:
        T item_;
        linked_list* p_next_;

    public:
        // Constructors
        linked_list();
        linked_list(T const& t, linked_list* pN=nullptr): item_(t),p_next_(pN){};
        linked_list (linked_list const& L ) : item_(L.item()), p_next_(L.p_next() ? new linked_list(*L.p_next()) : nullptr) {};
        // Destructor
        ~linked_list(){
            delete p_next_;
            p_next_ = nullptr;
            }

        // Getters
        linked_list* p_next() const { return p_next_; }
        linked_list*& p_next() { return p_next_; }
        T const& item() const { return item_; }
        T& item() { return item_;}

        //operator
        linked_list& last();
        void append(T const& t);
        size_t lenght() const;
        void insert_next_item(T &t);    // Ajoute un deuxième élément
        void insert_first_item(T &t);   // Ajoute un premier élément

        const linked_list& operator=(linked_list const& L);
        void drop_next_item();
        void drop_first_item();

        void truncate_items(double threshold); 
        
        const linked_list<T>& linked_list<T>::operator+=(linked_list<T>& L);
        linked_list<T>& linked_list<T>::order( size_t length );
};


template<typename T>
linked_list<T>::linked_list(){
    p_next_ = nullptr;
    item_ = T();
}

template<typename T>
linked_list<T>& linked_list<T>::last(){
    return p_next() ? p_next()->last() : *this; // Si p_next_ n'est pas nul, on appelle last() sur le nœud suivant sinon on renvoie l'adresse de l'élément
}

template<typename T>
size_t linked_list<T>::lenght() const{
    return p_next_ ? 1 + p_next_->lenght() : 1;
}

template<typename T>
void linked_list<T>::append(T const& t){
    last().p_next() = new linked_list<T>(t);
}


template<typename T>
void linked_list<T>::insert_next_item(T &t){
    p_next() = new linked_list<T>(t, this->p_next());
}

template<typename T>
void linked_list<T>::insert_first_item(T &t){
    insert_next_item(item_); // On rajoute le même item après le premier puis on modifie la valeur du premier
    item_ = t;
}

template<typename T>
const linked_list<T>& linked_list<T>::operator=(linked_list<T> const& L){
    if (this!=&L){
        this->item() = L.item();

        if (p_next()){  // encore de la place dans la liste de gauche ?
            if (L.p_next()){
                *p_next() = *L.p_next();
            }
            else{
                delete p_next();
                p_next() = nullptr;
            }
        }
        else{   // pas de place, on en crée une nouvelle
            p_next() = new linked_list<T>(*L.p_next());

        }
    }
    return *this;
}


template<typename T>
void linked_list<T>::drop_next_item(){  // COPIER LE CODE DE CORRECTION DU PROF POUR CETTE FONCTION
    if (p_next()){ // S'il n'y a pas qu'un seul élément dans la liste
        linked_list<T>* p_temp = p_next();
        p_next() = p_next()->p_next();
        p_temp->p_next() = nullptr;
        delete p_temp;
    }
}

template<typename T>
void linked_list<T>::drop_first_item(){
    if (p_next()==nullptr){ // S'il n'y a qu'un seul élément dans la liste
        throw std::runtime_error("La liste ne contient qu'un seul élément !");
    }
    else{
        item() = p_next()->item();
        drop_next_item();
    }
}

template<typename T>
void linked_list<T>::truncate_items(double threshold){
    if (p_next()){
        if (p_next()->item().get_value() < threshold){
            drop_next_item();
            truncate_items(threshold);
        }
        else{
            p_next()->truncate_items(threshold);
        }
    }
    if (p_next() && item().get_value() < threshold )
    {
        drop_first_item();
    }
}


template<typename T>
const linked_list<T>& linked_list<T>::operator+=(linked_list<T>& L){
    linked_list<T>* p_scan1 = this;
    linked_list<T>* p_scan2 = &L;

    if (L.item() <= item()){
        insert_first_item(L.item());
        p_scan2 = p_scan2->p_next();
    }
    for(; p_scan1->p_next( ) ; p_scan1=p_scan1->p_next()){  // Exactement comme in while
        if (p_scan2 && p_scan1->item()==p_scan2.item()){
            p_scan2 = p_scan2->p_next();
            
        } // Vérifier que les éléments ne sont pas égaux
        for(; p_scan2 && (p_scan2->item() < p_scan1->p_next()->item ( ) ) ; p_scan2 = p_scan2->p_next()){
            p_scan1->insert_next_item (p_scan2->item());
            p_scan1 = p_scan1->p_next();
        }
    }
    if (p_scan2 && p_scan1->item() == p_scan2.item()){
        p_scan2 = p_scan2->p_next();
    }
    if (p_scan2){
        p_scan1->p_next() = new linked_list<T>(*p_scan2);
    }
    return *this;
}


template<class T>
linked_list<T>& linked_list<T>::order( size_t length ){
    linked_list<T>* p_scan1 = this;
    for (auto count =0 ; count < length/2-1; ++count){
        p_scan1 = p_scan1->p_next();
    }
    linked_list<T>* p_scan2 = p_scan1->p_next();
    p_scan1->p_next() = nullptr; // important to "cut" the second half of the calling linked_list

    this->order(length/2);
    *this.operator+=(p_scan2->order(length - length/2));

    return *this;
}