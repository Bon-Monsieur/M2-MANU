#include "class_linked_list.hpp"


template<typename T>
void print(linked_list<T> const& mylist) {
    std::cout << mylist.item() << " ";
    if (mylist.p_next() != nullptr){
        print(*mylist.p_next());
    }
}




int main() {
    linked_list<double> mylist(4.);
    mylist.append(5.);
    mylist.append(6.);
    double val = 0.;
    mylist.insert_next_item(val);
    mylist.insert_first_item(val);
    print(mylist);
    std::cout << std::endl;
    //linked_list<double> mylist2(1.) ;
    //mylist2= mylist;
    //print(mylist2);

    mylist.drop_next_item();
    print(mylist);
    return 0;
}