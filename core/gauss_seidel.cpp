#include "gauss_seidel.h"

#include "cuda/cuda_gauss_seidel.h"

#include <iostream>
#include <chrono>

#define DEBUG_INFO

template<typename type>
type get_u(point<type> delta, const potential<type>& u_field){
    switch (u_field.b_type) {
    case boundary_type::Non:
    case boundary_type::Dirichlet:
        break;
    case boundary_type::Neumann:
        return dot(delta, u_field.derivative);
    }
    return u_field.value;
}

template<typename type>
void set_x_boundary(size_t i, const boundary<type>& bc, const points<type>& points, field<type>& u_field){
    switch (bc.b_type) {
    case boundary_type::Non:
        break;
    case boundary_type::Dirichlet:
        for(size_t j = 0; j < points.n_y; j++){
            u_field[i][j].b_type = boundary_type::Dirichlet;
            u_field[i][j].value = bc.f(points[i][j].x, points[i][j].y);
        }
        break;
    case boundary_type::Neumann:
        for(size_t j = 0; j < points.n_y; j++){
            u_field[i][j].b_type = boundary_type::Neumann;
            u_field[i][j].derivative.x = bc.f(points[i][j].x, points[i][j].y);
        }
        break;
    }
}

template<typename type>
void set_y_boundary(size_t j, const boundary<type>& bc, const points<type>& points, field<type>& u_field){
    switch (bc.b_type) {
    case boundary_type::Non:
        break;
    case boundary_type::Dirichlet:
        for(size_t i = 0; i < points.n_x; i++){
            u_field[i][j].b_type = boundary_type::Dirichlet;
            u_field[i][j].value = bc.f(points[i][j].x, points[i][j].y);
        }
        break;
    case boundary_type::Neumann:
        for(size_t i = 0; i < points.n_x; i++){
            u_field[i][j].b_type = boundary_type::Neumann;
            u_field[i][j].derivative.y = bc.f(points[i][j].x, points[i][j].y);
        }
        break;
    }
}

template<typename type>
void adjust_field(const points<type>& points, field<type>& u_field){
    for(size_t i = 0; i < points.n_x; i++){
        if(type dy = points[i][1].y - points[i][0].y; u_field[i][0].b_type == boundary_type::Neumann){
            u_field[i][0].value = u_field[i][1].value - dy * u_field[i][0].derivative.y;
        }
        if(type dy = points[i][points.n_y - 1].y - points[i][points.n_y - 2].y; u_field[i][points.n_y - 1].b_type == boundary_type::Neumann){
            u_field[i][points.n_y - 1].value = u_field[i][points.n_y - 2].value + dy * u_field[i][points.n_y - 1].derivative.y;
        }

        u_field[i][0].derivative.y = (u_field[i][1].value - u_field[i][0].value) / (points[i][1].y - points[i][0].y);
        u_field[i][points.n_y - 1].derivative.y = (u_field[i][points.n_y - 1].value - u_field[i][points.n_y - 2].value) / (points[i][points.n_y - 1].y - points[i][points.n_y - 2].y);
    }
    for(size_t j = 0; j < points.n_y; j++){
        if(type dx = points[1][j].x - points[0][j].x; u_field[0][j].b_type == boundary_type::Neumann){
            u_field[0][j].value = u_field[1][j].value - dx * u_field[0][j].derivative.x;
        }
        if(type dx = points[points.n_x - 1][j].x - points[points.n_x - 2][j].x;  u_field[points.n_x - 1][j].b_type == boundary_type::Neumann){
            u_field[points.n_x - 1][j].value = u_field[points.n_x - 2][j].value + dx * u_field[points.n_x - 1][j].derivative.x;
        }

        u_field[0][j].derivative.x = (u_field[1][j].value - u_field[0][j].value) / (points[1][j].x - points[0][j].x);
        u_field[points.n_x - 1][j].derivative.x = (u_field[points.n_x - 1][j].value - u_field[points.n_x - 2][j].value) / (points[points.n_x - 1][j].x - points[points.n_x - 2][j].x);
    }

    for(size_t i = 1; i < points.n_x - 1; i++){
        for(size_t j = 1; j < points.n_y - 1; j++){
            u_field[i][j].derivative.x = (u_field[i + 1][j].value - u_field[i - 1][j].value) / (points[i + 1][j].x - points[i - 1][j].x);
            u_field[i][j].derivative.y = (u_field[i][j + 1].value - u_field[i][j - 1].value) / (points[i][j + 1].y - points[i][j - 1].y);
        }
    }
}

template<typename type>
void set_additional_conditions(const additional_conditions<type>& ac, const points<type>& points, field<type>& u_field){
    for(const auto& condition: ac){
        for(size_t i = 0; i < points.n_x; i++){
            for(size_t j = 0; j < points.n_y; j++){
                if(const auto set_field = condition(points[i][j].x, points[i][j].y); set_field.has_value()){
                    u_field[i][j] = set_field.value();
                }
            }
        }
    }
}

template<typename type>
void gauss_seidel(const gauss_seidel_info<type>& info, const points<type>& points, const function<type>& func, field<type>& u_field){
    bool convergence = false;
    size_t it = 0;
    for (; !convergence && it < info.max_it; it++ ){
        type dmax = type(0);
        for (size_t p = 0; p <= points.n_x - 3 + points.n_y - 3; p++ ){
            for (size_t q = p <= points.n_x - 3 ? 0 : p - (points.n_x - 3); q <= std::min(p, std::min(points.n_x - 3, points.n_y - 3)); q++ ){
                if(size_t i = 1 + p - q, j = 1 + q; u_field[i][j].b_type == boundary_type::Non){
                    type dx = type(0.5) * (points[i + 1][j].x - points[i - 1][j].x);
                    type dy = type(0.5) * (points[i][j + 1].y - points[i][j - 1].y);
                    type a = type(1) / dx / dx, b = type(1) / dy / dy, c = type(2) * (a + b);

                    type temp = u_field[i][j].value;
                    u_field[i][j].value = (
                                              a * (get_u({ type(2) * dx, 0}, u_field[i + 1][j]) + (u_field[i + 1][j].isNeumann() ? u_field[i - 1][j].value : type(0))) +
                                              a * (get_u({-type(2) * dx, 0}, u_field[i - 1][j]) + (u_field[i - 1][j].isNeumann() ? u_field[i + 1][j].value : type(0))) +
                                              b * (get_u({0,  type(2) * dy}, u_field[i][j + 1]) + (u_field[i][j + 1].isNeumann() ? u_field[i][j - 1].value : type(0))) +
                                              b * (get_u({0, -type(2) * dy}, u_field[i][j - 1]) + (u_field[i][j - 1].isNeumann() ? u_field[i][j + 1].value : type(0))) -
                                              func[i][j]) / c;

                    dmax = std::max(dmax, std::abs(temp - u_field[i][j].value));
                }
            }
        }
        convergence = dmax < info.eps;
    }
#ifdef DEBUG_INFO
    std::cout << "   convergence = " << (convergence ? "true" : "false") << std::endl;
    std::cout << "   last_it = " << it << std::endl;
#endif
}

template<typename type>
field<type> poisson_gauss_seidel(
    const gauss_seidel_info<type>& info,
    const points<type>& points,
    const function<type>& func,
    const boundary_condition<type>& bc,
    const additional_conditions<type>& ac)
{
#ifdef DEBUG_INFO
    const auto start = std::chrono::high_resolution_clock::now();
    std::cout << "poisson_gauss_seidel :"<< std::endl;
#endif

    field<type> u_field(points);

    set_x_boundary(0, bc.x_0, points, u_field);
    set_x_boundary(points.n_x - 1, bc.x_n, points, u_field);
    set_y_boundary(0, bc.y_0, points, u_field);
    set_y_boundary(points.n_y - 1, bc.y_n, points, u_field);
    set_additional_conditions(ac, points, u_field);

    if(info.use_cuda){
        cuda::gauss_seidel(cuda::gauss_seidel_info<type>{info.eps, info.max_it}, points, func, u_field);
    } else {
        gauss_seidel(info, points, func, u_field);
    }

    adjust_field(points, u_field);

#ifdef DEBUG_INFO
    std::cout << "   time = " <<
        std::chrono::duration<float, std::chrono::milliseconds::period>(std::chrono::high_resolution_clock::now() - start).count() << " ms" << std::endl;
#endif

    return u_field;
}

template struct gauss_seidel_info<float>;
template struct gauss_seidel_info<double>;

template field<float> poisson_gauss_seidel(
    const gauss_seidel_info<float>& info,
    const points<float>& points,
    const function<float>& func,
    const boundary_condition<float>& bc,
    const additional_conditions<float>& ac);

template field<double> poisson_gauss_seidel(
    const gauss_seidel_info<double>& info,
    const points<double>& points,
    const function<double>& func,
    const boundary_condition<double>& bc,
    const additional_conditions<double>& ac);
