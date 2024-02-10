#include "field.h"

#include "operations.h"

#include <fstream>

template<typename type>
boundary_condition<type>::boundary_condition(boundary<type> x_0, boundary<type> x_n, boundary<type> y_0, boundary<type> y_n) :
    x_0(x_0), x_n(x_n), y_0(y_0), y_n(y_n){}

template<typename type>
boundary_condition<type>::boundary_condition() :
    x_0({[](type, type){return type(0);}, boundary_type::Dirichlet}),
    x_n({[](type, type){return type(0);}, boundary_type::Dirichlet}),
    y_0({[](type, type){return type(0);}, boundary_type::Dirichlet}),
    y_n({[](type, type){return type(0);}, boundary_type::Dirichlet}){}

template struct boundary_condition<float>;
template struct boundary_condition<double>;

template<typename type>
bool potential<type>::isDirichlet(){
    return b_type == boundary_type::Dirichlet;
}

template<typename type>
bool potential<type>::isNeumann(){
    return b_type == boundary_type::Neumann;
}

template struct potential<float>;
template struct potential<double>;

template<typename out_type, typename type>
out_type& operator<<(out_type& out, const potential<type>& u_field){
    out << '(' << u_field.value << ',' <<u_field.derivative.x << ',' <<u_field.derivative.y << ')';
    return out;
}

template std::ofstream& operator<<(std::ofstream& out, const potential<float>& u_field);
template std::ofstream& operator<<(std::ofstream& out, const potential<double>& u_field);

template<typename type>
field<type>::field(const points<type>& grid){
    generate(grid);
}

template<typename type>
field<type>::~field(){
    destroy(data, n_x);
}

template<typename type>
void field<type>::generate(const points<type>& grid){
    n_x = grid.n_x;
    n_y = grid.n_y;
    data = ::generate<potential<type>>(grid.n_x, grid.n_y);
}

template<typename type>
potential<type>* field<type>::operator[](const size_t i){
    return data[i];
}

template<typename type>
const potential<type>* field<type>::operator[](const size_t i) const {
    return data[i];
}

template struct field<float>;
template struct field<double>;

