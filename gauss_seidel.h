#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include <iostream>
#include <chrono>

#include "point.h"
#include "function.h"
#include "field.h"

template<typename type>
type get_u_x(type delta, const potential<type>& u_field){
    switch (u_field.b_type) {
    case boundary_type::Non:
    case boundary_type::Dirichlet:
        break;
    case boundary_type::Neumann:
        return delta * u_field.derivative_x;
    }
    return u_field.value;
}

template<typename type>
type get_u_y(type delta, const potential<type>& u_field){
    switch (u_field.b_type) {
    case boundary_type::Non:
    case boundary_type::Dirichlet:
        break;
    case boundary_type::Neumann:
        return delta * u_field.derivative_y;
    }
    return u_field.value;
}

template<typename type>
void set_x_boundary_Field(size_t i, const boundary<type>& bc, const points<type>& points, potential<type>** &u_field){
    switch (bc.b_type) {
    case boundary_type::Dirichlet:
        for(size_t j = 0; j < points.n_y; j++){
            u_field[i][j].b_type = boundary_type::Dirichlet;
            u_field[i][j].value = bc.f(points.data[i][j].x, points.data[i][j].y);
        }
        break;
    case boundary_type::Neumann:
        for(size_t j = 0; j < points.n_y; j++){
            u_field[i][j].b_type = boundary_type::Neumann;
            u_field[i][j].derivative_x = bc.f(points.data[i][j].x, points.data[i][j].y);
        }
        break;
    }
}

template<typename type>
void set_y_boundary_Field(size_t j, const boundary<type>& bc, const points<type>& points, potential<type>** &u_field){
    switch (bc.b_type) {
    case boundary_type::Dirichlet:
        for(size_t i = 0; i < points.n_x; i++){
            u_field[i][j].b_type = boundary_type::Dirichlet;
            u_field[i][j].value = bc.f(points.data[i][j].x, points.data[i][j].y);
        }
        break;
    case boundary_type::Neumann:
        for(size_t i = 0; i < points.n_x; i++){
            u_field[i][j].b_type = boundary_type::Neumann;
            u_field[i][j].derivative_y = bc.f(points.data[i][j].x, points.data[i][j].y);
        }
        break;
    }
}

template<typename type>
void adjust_field(const points<type>& points, potential<type>** &u_field){
    for(size_t i = 0; i < points.n_x; i++){
        if(type dy = points.data[i][1].y - points.data[i][0].y; u_field[i][0].b_type == boundary_type::Neumann){
            u_field[i][0].value = u_field[i][1].value - dy * u_field[i][0].derivative_y;
        }
        if(type dy = points.data[i][points.n_y - 1].y - points.data[i][points.n_y - 2].y; u_field[i][points.n_y - 1].b_type == boundary_type::Neumann){
            u_field[i][points.n_y - 1].value = u_field[i][points.n_y - 2].value + dy * u_field[i][points.n_y - 1].derivative_y;
        }

        u_field[i][0].derivative_y = (u_field[i][1].value - u_field[i][0].value) / (points.data[i][1].y - points.data[i][0].y);
        u_field[i][points.n_y - 1].derivative_y = (u_field[i][points.n_y - 1].value - u_field[i][points.n_y - 2].value) / (points.data[i][points.n_y - 1].y - points.data[i][points.n_y - 2].y);
    }
    for(size_t j = 0; j < points.n_y; j++){
        if(type dx = points.data[1][j].x - points.data[0][j].x; u_field[0][j].b_type == boundary_type::Neumann){
            u_field[0][j].value = u_field[1][j].value - dx * u_field[0][j].derivative_x;
        }
        if(type dx = points.data[points.n_x - 1][j].x - points.data[points.n_x - 2][j].x;  u_field[points.n_x - 1][j].b_type == boundary_type::Neumann){
            u_field[points.n_x - 1][j].value = u_field[points.n_x - 2][j].value + dx * u_field[points.n_x - 1][j].derivative_x;
        }

        u_field[0][j].derivative_x = (u_field[1][j].value - u_field[0][j].value) / (points.data[1][j].x - points.data[0][j].x);
        u_field[points.n_x - 1][j].derivative_x = (u_field[points.n_x - 1][j].value - u_field[points.n_x - 2][j].value) / (points.data[points.n_x - 1][j].x - points.data[points.n_x - 2][j].x);
    }

    for(size_t i = 1; i < points.n_x - 1; i++){
        for(size_t j = 1; j < points.n_y - 1; j++){
            u_field[i][j].derivative_x = (u_field[i + 1][j].value - u_field[i - 1][j].value) / (points.data[i + 1][j].x - points.data[i - 1][j].x);
            u_field[i][j].derivative_y = (u_field[i][j + 1].value - u_field[i][j - 1].value) / (points.data[i][j + 1].y - points.data[i][j - 1].y);
        }
    }
}

template<typename type>
type get_correction(size_t i, size_t j, type dx, type dy, potential<type>** &u_field){
    return  - (u_field[i - 1][j].b_type == boundary_type::Neumann ? dy * dy : 0)
            - (u_field[i + 1][j].b_type == boundary_type::Neumann ? dy * dy : 0)
            - (u_field[i][j - 1].b_type == boundary_type::Neumann ? dx * dx : 0)
            - (u_field[i][j + 1].b_type == boundary_type::Neumann ? dx * dx : 0);
}

template<typename type>
struct gauss_seidel_info{
    type eps{0};
    size_t max_it{0};
};

template<typename type>
field<type> poisson_gauss_seidel(
    const gauss_seidel_info<type>& info,
    const points<type>& points,
    const function<type>& func,
    const boundary_condition<type>& bc,
    const additional_conditions<type>& ac = {})
{
#ifdef DEBUG_INFO
    const auto start = std::chrono::high_resolution_clock::now();
#endif

    field<type> u_field(points);
    set_x_boundary_Field(0, bc.x_0, points, u_field.data);
    set_x_boundary_Field(points.n_x - 1, bc.x_n, points, u_field.data);
    set_y_boundary_Field(0, bc.y_0, points, u_field.data);
    set_y_boundary_Field(points.n_y - 1, bc.y_n, points, u_field.data);

    for(const auto& condition: ac){
        for(size_t i = 0; i < points.n_x; i++){
            for(size_t j = 0; j < points.n_y; j++){
                if(std::optional<potential<type>> seted_field = condition(points.data[i][j].x, points.data[i][j].y); seted_field.has_value()){
                    u_field.data[i][j] = seted_field.value();
                }
            }
        }
    }

    bool convergence = false;
    size_t it = 0;
    for (; !convergence && it < info.max_it; it++ ){
        type dmax = type(0);
        for (size_t p = 0; p <= points.n_x - 3 + points.n_y - 3; p++ ){
            for (size_t q = p <= points.n_x - 3 ? 0 : p - (points.n_x - 3); q <= std::min(p, std::min(points.n_x - 3, points.n_y - 3)); q++ ){
                if(size_t i = 1 + p - q, j = 1 + q; u_field.data[i][j].b_type == boundary_type::Non){
                    type dx = points.data[i + 1][j].x - points.data[i][j].x, dy = points.data[i][j + 1].y - points.data[i][j].y;
                    type temp = u_field.data[i][j].value, a = dx * dx, b = dy * dy, c = type(2) * (a + b) + get_correction(i, j, dx, dy, u_field.data);

                    u_field.data[i][j].value = (b * get_u_x(dx, u_field.data[i + 1][j]) + b * get_u_x(- dx, u_field.data[i - 1][j]) +
                                           a * get_u_y(dy, u_field.data[i][j + 1]) + a * get_u_y(- dy, u_field.data[i][j - 1]) -
                                           a * b * func.data[i][j]) / c;

                    dmax = std::max(dmax, std::abs(temp - u_field.data[i][j].value));
                }
            }
        }
        convergence = dmax < info.eps;
    }

    adjust_field(points, u_field.data);

#ifdef DEBUG_INFO
    std::cout << "poisson_gauss_seidel :"<< std::endl;
    std::cout << "   convergence = " << (convergence ? "true" : "false") << std::endl;
    std::cout << "   last_it = " << it << std::endl;
    std::cout << "   time = " <<
        std::chrono::duration<float, std::chrono::milliseconds::period>(std::chrono::high_resolution_clock::now() - start).count() << " ms" << std::endl;
#endif

    return u_field;
}

#endif
