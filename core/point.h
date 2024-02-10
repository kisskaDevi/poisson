#ifndef POINT_H
#define POINT_H

template<typename type>
struct point{
    type x{0};
    type y{0};

    point() = default;
    point(type x, type y);
    point operator+(const point& p);
    point operator*(const point<type>& b);
};

template<typename out_type, typename type>
out_type& operator<<(out_type& out, const point<type>& p);

template<typename type>
type dot(const point<type>& a, const point<type>& b);

template<typename type>
struct points{
    point<type>** data{nullptr};
    size_t n_x{0};
    size_t n_y{0};

    points(size_t n_x, size_t n_y);
    ~points();

    void generate(type x_0, type x_n, type y_0, type y_n);
    point<type>* operator[](const size_t i);
    const point<type>* operator[](const size_t i) const;
};

#endif
