#ifndef FIELD_H
#define FIELD_H

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

    boundary_condition(boundary<type> x_0, boundary<type> x_n, boundary<type> y_0, boundary<type> y_n);
    boundary_condition();
};

template<typename type>
struct potential{
    boundary_type b_type{Non};
    type value{type(0)};
    point<type> derivative;

    bool isDirichlet();
    bool isNeumann();
};

template<typename out_type, typename type>
out_type& operator<<(out_type& out, const potential<type>& u_field);

template<typename type>
struct field
{
    potential<type>** data{nullptr};
    size_t n_x{0};
    size_t n_y{0};

    field(const points<type>& grid);
    ~field();

    void generate(const points<type>& grid);
    potential<type>* operator[](const size_t i);
    const potential<type>* operator[](const size_t i) const;
};

#endif // FIELD_H
