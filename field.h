#ifndef FIELD_H
#define FIELD_H

#include <fstream>
#include <vector>
#include <optional>

#include "point.h"

enum boundary_type{
    Non,
    Dirichlet,
    Neumann
};

template<typename type>
struct boundary{
    type (*f)(type x, type y){nullptr};
    boundary_type b_type{boundary_type::Non};
};

template<typename type>
struct boundary_condition{
    boundary<type> x_0;
    boundary<type> x_n;
    boundary<type> y_0;
    boundary<type> y_n;

    boundary_condition(boundary<type> x_0, boundary<type> x_n, boundary<type> y_0, boundary<type> y_n) :
        x_0(x_0), x_n(x_n), y_0(y_0), y_n(y_n){}

    boundary_condition<type>() :
        x_0({[](type, type){return type(0);}, boundary_type::Dirichlet}),
        x_n({[](type, type){return type(0);}, boundary_type::Dirichlet}),
        y_0({[](type, type){return type(0);}, boundary_type::Dirichlet}),
        y_n({[](type, type){return type(0);}, boundary_type::Dirichlet}){}
};

template<typename type>
struct potential{
    boundary_type b_type{Non};
    type value{type(0)};
    point<type> derivative;

    bool isDirichlet(){
        return b_type == boundary_type::Dirichlet;
    }

    bool isNeumann(){
        return b_type == boundary_type::Neumann;
    }
};

template<typename type>
std::ostream& operator<<(std::ostream& out, const potential<type>& u_field){
    out << '(' << u_field.value << ',' <<u_field.derivative.x << ',' <<u_field.derivative.y << ')';
    return out;
}

template<typename type>
using additional_conditions = std::vector<std::optional<potential<type>> (*)(type, type)>;

template<typename type>
struct field
{
    potential<type>** data{nullptr};
    size_t n_x{0};
    size_t n_y{0};

    field(const points<type>& grid){
        generate(grid);
    }

    ~field(){
        destroy(data, n_x);
    }

    void generate(const points<type>& grid){
        n_x = grid.n_x;
        n_y = grid.n_y;
        data = ::generate<potential<type>>(grid.n_x, grid.n_y);
    }

    potential<type>* operator[](const size_t i){
        return data[i];
    }

    const potential<type>* operator[](const size_t i) const {
        return data[i];
    }
};

#endif // FIELD_H
