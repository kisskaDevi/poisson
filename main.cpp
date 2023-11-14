#include <iostream>
#include <cmath>

#define DEBUG_INFO

#include "operations.h"
#include "function.h"
#include "gauss_seidel.h"

template<typename type>
void example_1(points<type>& xy_grid, function<type>& g_func, boundary_condition<type>& bc, function<type>& ex_func, additional_conditions<type>&){
    type x_0 = 0.0f, x_n =  1.0f;
    type y_0 = 0.0f, y_n =  1.0f;
    xy_grid.generate(x_0, x_n, y_0, y_n);
    g_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type, type){
                                   return type(0);
                               })));
    bc = boundary_condition<type>(
        {[](type, type y){return 100.0f - 200.0f * y;}, boundary_type::Dirichlet},
        {[](type, type y){return - 100.0f + 200.0f * y;}, boundary_type::Dirichlet},
        {[](type x, type){return 100.0f - 200.0f * x;}, boundary_type::Dirichlet},
        {[](type x, type){return - 100.0f + 200.0f * x;}, boundary_type::Dirichlet});
    ex_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                    return 100.0f * (1.0f - 2.0f * x) * (1.0f - 2.0f * y);
                                })));
}

template<typename type>
void example_2(points<type>& xy_grid, function<type>& g_func, boundary_condition<type>& bc, function<type>& ex_func, additional_conditions<type>&){
    type x_0 = -1.0f, x_n =  1.0f;
    type y_0 = -1.0f, y_n =  1.0f;
    xy_grid.generate(x_0, x_n, y_0, y_n);
    g_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                   type pi = std::atan(type(1)) * type(4);
                                   return - type(2) * pi * pi * std::sin(pi * x) * std::sin(pi * y);
                               })));
    bc = boundary_condition<type>();

    ex_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                    type pi = std::atan(type(1)) * type(4);
                                    return std::sin(pi * x) * std::sin(pi * y);
                                })));
}

template<typename type>
void example_3(points<type>& xy_grid, function<type>& g_func, boundary_condition<type>& bc, function<type>& ex_func, additional_conditions<type>&){
    type x_0 = -1.0f, x_n =  1.0f;
    type y_0 = -1.0f, y_n =  1.0f;
    xy_grid.generate(x_0, x_n, y_0, y_n);
    g_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                   type pi = std::atan(type(1)) * type(4);
                                   return - type(2) * pi * pi * std::sin(pi * x) * std::sin(pi * y);
                               })));
    bc = boundary_condition<type>(
        {[](type, type y){type pi = std::atan(type(1)) * type(4); return - pi * std::sin(pi * y);}, boundary_type::Neumann},
        {[](type, type y){type pi = std::atan(type(1)) * type(4); return - pi * std::sin(pi * y);}, boundary_type::Neumann},
        {[](type x, type){type pi = std::atan(type(1)) * type(4); return - pi * std::sin(pi * x);}, boundary_type::Neumann},
        {[](type x, type){type pi = std::atan(type(1)) * type(4); return - pi * std::sin(pi * x);}, boundary_type::Neumann});

    ex_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                    type pi = std::atan(type(1)) * type(4);
                                    return std::sin(pi * x) * std::sin(pi * y);
                                })));
}

template<typename type>
void example_4(points<type>& xy_grid, function<type>& g_func, boundary_condition<type>& bc, function<type>& ex_func, additional_conditions<type>&){
    type x_0 = -4.0f, x_n = 4.0f;
    type y_0 = -4.0f, y_n = 4.0f;
    xy_grid.generate(x_0, x_n, y_0, y_n);
    g_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                   return - type(0.25) * std::cos(type(0.5) * x * y) * (x * x + y * y);
                               })));
    bc = boundary_condition<type>(
        {[](type, type y){return type(0.5) * y * std::sin(type(2.0) * y);}, boundary_type::Neumann},
        {[](type, type y){return - type(0.5) * y * std::sin(type(2.0) * y);}, boundary_type::Neumann},
        {[](type x, type){return std::cos(type(2.0) * x);}, boundary_type::Dirichlet},
        {[](type x, type){return std::cos(type(2.0) * x);}, boundary_type::Dirichlet});

    ex_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type x, type y){
                                    return std::cos(type(0.5) * x * y);
                                })));
}

template<typename type>
void example_5(points<type>& xy_grid, function<type>& g_func, boundary_condition<type>& bc, function<type>&, additional_conditions<type>& ac){
    type x_0 = -6.0f, x_n =  6.0f;
    type y_0 = -3.0f, y_n =  3.0f;
    xy_grid.generate(x_0, x_n, y_0, y_n);
    g_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type, type){return type(0);})));
    bc = boundary_condition<type>(
        {[](type, type){return type(0);}, boundary_type::Neumann},
        {[](type, type){return type(0);}, boundary_type::Neumann},
        {[](type, type){return type(0);}, boundary_type::Neumann},
        {[](type, type){return type(0);}, boundary_type::Neumann});

    ac = additional_conditions<type>{
        [](type x, type y){ return x >= -2.0f && x <= 2.0f && y >= 0.5f && y <= 0.6f ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, 2.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();},
        [](type x, type y){ return x >= -2.0f && x <= 2.0f && y >= - 0.6f && y <= - 0.5f && y <= 0.6f ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, - 2.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();}
    };
}

template<typename type>
void example_6(points<type>& xy_grid, function<type>& g_func, boundary_condition<type>& bc, function<type>&, additional_conditions<type>& ac){
    type x_0 = -6.0f, x_n =  6.0f;
    type y_0 = -3.0f, y_n =  3.0f;
    xy_grid.generate(x_0, x_n, y_0, y_n);
    g_func.generate(xy_grid, static_cast<type(*)(type,type)>(([](type, type){return type(0.0);})));
    bc = boundary_condition<type>(
        {[](type, type){return type(0);}, boundary_type::Neumann},
        {[](type, type){return type(0);}, boundary_type::Neumann},
        {[](type, type){return type(0);}, boundary_type::Neumann},
        {[](type, type){return type(0);}, boundary_type::Neumann});

    ac = additional_conditions<type>{
        [](type x, type y){ return x >= 2.0f && x <= 3.0f && y >= - 1.5f && y <= - 1.2f && y <= - 1.2f ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, 2.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();},
        [](type x, type y){ return x >= -3.0f && x <= -2.0f && y >= - 1.5f && y <= - 1.2f ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, 2.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();},
        [](type x, type y){ return (x - 0.9) * (x - 0.9) + (y + 1.35) * (y + 1.35) < 0.1 ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, - 2.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();},
        [](type x, type y){ return (x + 0.9) * (x + 0.9) + (y + 1.35) * (y + 1.35) < 0.1 ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, - 2.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();},
        [](type x, type y){ return x * x  + (y + 1.35) * (y + 1.35) < 0.1 ?
                                        std::optional<potential<type>>({boundary_type::Dirichlet, 1.0f, {0.0f, 0.0f}}) :
                                        std::optional<potential<type>>();}
    };
}

int main()
{
    using Type = float;

    points<Type> xy_grid(101, 101);
    function<Type> g_func, ex_func;
    boundary_condition<Type> bc;
    additional_conditions<Type> ac;

    example_5(xy_grid, g_func, bc, ex_func, ac);

    field<Type> u_field = poisson_gauss_seidel(gauss_seidel_info<Type>{Type(0.000001), 50000}, xy_grid, g_func, bc, ac);

    std::filesystem::path res_path = std::filesystem::current_path() / "release/res";
    if(! std::filesystem::exists(res_path)){
        std::filesystem::create_directories(res_path);
    }
    out_if_file(res_path / "xy_grid.txt", xy_grid.data, xy_grid.n_x, xy_grid.n_y);
    out_if_file(res_path / "g_func.txt", g_func.data, g_func.n_x, g_func.n_y);
    out_if_file(res_path / "ex_func.txt", ex_func.data, ex_func.n_x, ex_func.n_y);
    out_if_file(res_path / "u_field.txt", u_field.data, u_field.n_x, u_field.n_y);

    return 0;
}
