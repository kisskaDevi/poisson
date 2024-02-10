#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "point.h"
#include "function.h"
#include "field.h"

#include <vector>
#include <optional>

template<typename type>
struct gauss_seidel_info{
    type eps{0};
    size_t max_it{0};
};

template<typename type>
using additional_conditions = std::vector<std::optional<potential<type>> (*)(type, type)>;

template<typename type>
field<type> poisson_gauss_seidel(
    const gauss_seidel_info<type>& info,
    const points<type>& points,
    const function<type>& func,
    const boundary_condition<type>& bc,
    const additional_conditions<type>& ac = {});

#endif
