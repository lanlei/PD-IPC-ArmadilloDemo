#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <unordered_map>
#include <map>
#include <array>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <type_traits>
#include <unordered_set>
#include <string>

namespace Bow {
template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;

template <typename T, int n, int m>
using Array = Eigen::Array<T, n, m, 0, n, m>;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;

template <typename DerivedV>
using Field = std::vector<DerivedV, Eigen::aligned_allocator<DerivedV>>;

template <int n, int m = 1, class Derived = Bow::Vector<double, Eigen::Dynamic>>
Bow::Field<Bow::Matrix<typename Derived::Scalar, n, m>> to_field(const Eigen::MatrixBase<Derived>& mat)
{
    const auto* p_x = reinterpret_cast<const Bow::Matrix<typename Derived::Scalar, n, m>*>(mat.derived().data());
    return Bow::Field<Bow::Matrix<typename Derived::Scalar, n, m>>(p_x, p_x + mat.size() / (n * m));
}

template <typename Scalar, int n, int m>
Bow::Vector<Scalar, Eigen::Dynamic> to_vec(const Bow::Field<Bow::Matrix<Scalar, n, m>>& x)
{
    Bow::Vector<Scalar, Eigen::Dynamic> x_vec;
    x_vec.resize(x.size() * m * n);
    memcpy(x_vec.data(), reinterpret_cast<const Scalar*>(x.data()), sizeof(Scalar) * x_vec.size());
    return x_vec;
}

template <class T, int dim>
inline static std::array<T, dim> to_std_array(const T* data)
{
    if constexpr (dim == 2) {
        return std::array<T, 2>{ data[0], data[1] };
    }
    else {
        return std::array<T, 3>{ data[0], data[1], data[2] };
    }
}

template <typename DerivedV>
typename DerivedV::Scalar dotProduct(const Bow::Field<DerivedV>& lhs, const Bow::Field<DerivedV>& rhs)
{
    typename DerivedV::Scalar result(0);
    for (size_t i = 0; i < lhs.size(); ++i) result += lhs[i].dot(rhs[i]);
    return result;
}

template <typename DerivedV>
typename DerivedV::Scalar squaredNorm(const Bow::Field<DerivedV>& data)
{
    typename DerivedV::Scalar result = 0;
    for (size_t i = 0; i < data.size(); ++i) result += data[i].squaredNorm();
    return result;
}

} // namespace Bow

template <int dim>
struct VectorHash {
    typedef Bow::Vector<int, dim> IV;
    size_t operator()(const IV& a) const
    {
        std::size_t h = 0;
        for (int d = 0; d < dim; ++d) {
            h ^= std::hash<int>{}(a(d)) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

template <typename DerivedV>
Bow::Field<DerivedV> operator-(const Bow::Field<DerivedV>& lhs, const Bow::Field<DerivedV>& rhs)
{
    Bow::Field<DerivedV> result = lhs;
    for (size_t i = 0; i < result.size(); ++i) result[i] = lhs[i] - rhs[i];
    return result;
}

template <typename DerivedV>
Bow::Field<DerivedV> operator+(const Bow::Field<DerivedV>& lhs, const Bow::Field<DerivedV>& rhs)
{
    Bow::Field<DerivedV> result = lhs;
    for (size_t i = 0; i < result.size(); ++i) result[i] = lhs[i] + rhs[i];
    return result;
}

template <typename DerivedV>
void operator+=(Bow::Field<DerivedV>& lhs, const Bow::Field<DerivedV>& rhs)
{
    for (size_t i = 0; i < lhs.size(); ++i) lhs[i] += rhs[i];
}

template <typename DerivedV>
void operator+=(Bow::Field<DerivedV>& lhs, const DerivedV& rhs)
{
    for (size_t i = 0; i < lhs.size(); ++i) lhs[i] += rhs;
}

template <typename DerivedV>
void operator-=(Bow::Field<DerivedV>& lhs, const Bow::Field<DerivedV>& rhs)
{
    for (size_t i = 0; i < lhs.size(); ++i) lhs[i] -= rhs[i];
}

template <typename DerivedV>
Bow::Field<DerivedV> operator*(const typename DerivedV::Scalar lhs, const Bow::Field<DerivedV>& rhs)
{
    Bow::Field<DerivedV> result = rhs;
    for (size_t i = 0; i < result.size(); ++i) result[i] = lhs * rhs[i];
    return result;
}

template <typename DerivedV>
Bow::Field<DerivedV> operator/(const Bow::Field<DerivedV>& lhs, const typename DerivedV::Scalar rhs)
{
    Bow::Field<DerivedV> result = lhs;
    for (size_t i = 0; i < result.size(); ++i) result[i] = lhs[i] / rhs;
    return result;
}

template <typename DerivedV>
void operator*=(Bow::Field<DerivedV>& lhs, const typename DerivedV::Scalar rhs)
{
    for (size_t i = 0; i < lhs.size(); ++i) lhs[i] *= rhs;
}

#endif