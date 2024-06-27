#include "Primitives.h"
#include <Bow/Utils/Logging.h>

namespace Bow {
namespace Geometry {
template <class T, int dim>
BOW_INLINE void cube(const Vector<int, dim> resolution, const T dx, Field<Vector<T, dim>>& points, Field<Vector<int, dim + 1>>& elements, const Vector<T, dim> center)
{
    points.clear();
    elements.clear();
    if constexpr (dim == 2) {
        int nx = resolution[0];
        int ny = resolution[1];
        points.resize((nx + 1) * (ny + 1));
        elements.clear();
        for (int ix = 0; ix <= nx; ++ix)
            for (int iy = 0; iy <= ny; ++iy) {
                int index = ix * (ny + 1) + iy;
                points[index] << (T)ix * dx, (T)iy * dx;
            }

        for (int ix = 0; ix < nx; ++ix)
            for (int iy = 0; iy < ny; ++iy) {
                int vertices[] = { ix * (ny + 1) + iy, ix * (ny + 1) + iy + 1,
                    (ix + 1) * (ny + 1) + iy, (ix + 1) * (ny + 1) + iy + 1 };
                if ((ix % 2) ^ (iy % 2)) {
                    elements.push_back(Vector<int, 3>(vertices[0], vertices[2], vertices[1]));
                    elements.push_back(Vector<int, 3>(vertices[1], vertices[2], vertices[3]));
                }
                else {
                    elements.push_back(Vector<int, 3>(vertices[0], vertices[2], vertices[3]));
                    elements.push_back(Vector<int, 3>(vertices[0], vertices[3], vertices[1]));
                }
            }
    }
    else {
        // Reference: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/48509/versions/3/previews/COMP_GEOM_TLBX/html/Divide_hypercube_6_simplices_3D.html
        int nx = resolution[0];
        int ny = resolution[1];
        int nz = resolution[2];
        points.resize((nx + 1) * (ny + 1) * (nz + 1));
        elements.clear();
        for (int ix = 0; ix <= nx; ++ix)
            for (int iy = 0; iy <= ny; ++iy)
                for (int iz = 0; iz <= nz; ++iz) {
                    int index = ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz;
                    points[index] << (T)ix * dx, (T)iy * dx, (T)iz * dx;
                }

        for (int ix = 0; ix < nx; ++ix)
            for (int iy = 0; iy < ny; ++iy)
                for (int iz = 0; iz < nz; ++iz) {
                    int vertices[] = { ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz, ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz + 1,
                        ix * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz, ix * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz + 1,
                        (ix + 1) * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz, (ix + 1) * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz + 1,
                        (ix + 1) * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz, (ix + 1) * (ny + 1) * (nz + 1) + (iy + 1) * (nz + 1) + iz + 1 };
                    elements.push_back(Vector<int, 4>(vertices[0], vertices[4], vertices[2], vertices[5]));
                    elements.push_back(Vector<int, 4>(vertices[0], vertices[2], vertices[1], vertices[5]));
                    elements.push_back(Vector<int, 4>(vertices[1], vertices[2], vertices[3], vertices[5]));
                    elements.push_back(Vector<int, 4>(vertices[2], vertices[5], vertices[4], vertices[6]));
                    elements.push_back(Vector<int, 4>(vertices[2], vertices[7], vertices[5], vertices[6]));
                    elements.push_back(Vector<int, 4>(vertices[2], vertices[3], vertices[5], vertices[7]));
                }
    }
    auto point_mat = Eigen::Map<Matrix<T, dim, Eigen::Dynamic>>(reinterpret_cast<T*>(points.data()), dim, points.size());
    point_mat.colwise() -= (point_mat.rowwise().mean() - center);
    Vector<T, dim> min_corner = point_mat.rowwise().minCoeff();
    Vector<T, dim> max_corner = point_mat.rowwise().maxCoeff();
    if constexpr (dim == 2) {
        Logging::info("Min Corner: ", min_corner(0), " ", min_corner(1));
        Logging::info("Max Corner: ", max_corner(0), " ", max_corner(1));
    }
    else {
        Logging::info("Min Corner: ", min_corner(0), " ", min_corner(1), " ", min_corner(2));
        Logging::info("Max Corner: ", max_corner(0), " ", max_corner(1), " ", max_corner(2));
    }
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template void cube(const Vector<int, 2> resolution, const float dx, Field<Vector<float, 2>>& points, Field<Vector<int, 3>>& elements, const Vector<float, 2> center);
template void cube(const Vector<int, 3> resolution, const float dx, Field<Vector<float, 3>>& points, Field<Vector<int, 4>>& elements, const Vector<float, 3> center);
#endif
#ifdef BOW_COMPILE_DOUBLE
template void cube(const Vector<int, 2> resolution, const double dx, Field<Vector<double, 2>>& points, Field<Vector<int, 3>>& elements, const Vector<double, 2> center);
template void cube(const Vector<int, 3> resolution, const double dx, Field<Vector<double, 3>>& points, Field<Vector<int, 4>>& elements, const Vector<double, 3> center);
#endif
#endif
}
} // namespace Bow::Geometry