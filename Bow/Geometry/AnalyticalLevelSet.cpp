#include "AnalyticalLevelSet.h"

namespace Bow::Geometry {

template <class T, int dim>
HalfSpaceLevelSet<T, dim>::HalfSpaceLevelSet(Type type, const TV& origin, const TV& outward_normal)
    : AnalyticalLevelSet<T, dim>(type), origin(origin), outward_normal(outward_normal.normalized())
{
}

template <class T, int dim>
T HalfSpaceLevelSet<T, dim>::signed_distance(const TV& X)
{
    T result = (X - origin).dot(outward_normal);
    return result;
}

template <class T, int dim>
Vector<T, dim> HalfSpaceLevelSet<T, dim>::normal(const TV& X)
{
    return outward_normal;
}

template <class T, int dim>
AlignedBoxLevelSet<T, dim>::AlignedBoxLevelSet(Type type, const TV& min_corner, const TV& max_corner)
    : AnalyticalLevelSet<T, dim>(type)
{
    center = (min_corner + max_corner) * 0.5;
    half_edges = (max_corner - min_corner) * 0.5;
}

template <class T, int dim>
T AlignedBoxLevelSet<T, dim>::signed_distance(const TV& X)
{
    TV X_centered = X - center;
    TV d = X_centered.array().abs() - half_edges;
    T dd = d.array().maxCoeff();
    TV qq = d;
    for (int i = 0; i < dim; i++)
        if (qq(i) < (T)0)
            qq(i) = (T)0;
    // min(dd, (T)0) is to deal with inside box case
    T result = min(dd, (T)0) + qq.norm();
    return result;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template class HalfSpaceLevelSet<float, 2>;
template class HalfSpaceLevelSet<float, 3>;
#endif
#ifdef BOW_COMPILE_DOUBLE
template class HalfSpaceLevelSet<double, 2>;
template class HalfSpaceLevelSet<double, 3>;
#endif
#endif

} // namespace Bow::Geometry
