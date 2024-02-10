#ifndef FUNCTION_H
#define FUNCTION_H

#include "point.h"

template<typename type>
struct function
{
    type** data{nullptr};
    size_t n_x{0};
    size_t n_y{0};

    ~function();
    void generate(const points<type>& grid, type (*f)(type x, type y));
    type* operator[](const size_t i);
    const type* operator[](const size_t i) const;
};

#endif // FUNCTION_H
