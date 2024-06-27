#include "Query.h"
#include <algorithm>

namespace Bow::Geometry {
template <int dim>
BOW_INLINE void find_boundary(const Field<Vector<int, dim + 1>>& elements, Field<Vector<int, dim>>& boundary_elements)
{
    boundary_elements.clear();
    std::unordered_map<Vector<int, dim>, Field<Vector<int, dim>>, VectorHash<dim>> half_elements;
    Field<Vector<int, dim>> element_boundary(dim + 1);
    if constexpr (dim == 2) {
        element_boundary[0] = Vector<int, dim>(1, 0);
        element_boundary[1] = Vector<int, dim>(2, 1);
        element_boundary[2] = Vector<int, dim>(0, 2);
    }
    else {
        element_boundary[0] = Vector<int, dim>(2, 1, 0);
        element_boundary[1] = Vector<int, dim>(0, 1, 3);
        element_boundary[2] = Vector<int, dim>(1, 2, 3);
        element_boundary[3] = Vector<int, dim>(2, 0, 3);
    }

    for (const auto& indices : elements) {
        for (const auto& local_elem : element_boundary) {
            Vector<int, dim> half_elem;
            for (int i = 0; i < dim; ++i) {
                half_elem[i] = indices[local_elem[i]];
            }
            Vector<int, dim> elem = half_elem;
            std::sort(elem.data(), elem.data() + dim);
            half_elements[elem].push_back(half_elem);
        }
    }

    std::unordered_set<int> touched_indices;
    for (auto it : half_elements)
        if (it.second.size() == 1)
            boundary_elements.push_back(it.second[0]);
}

BOW_INLINE void find_boundary(const Field<Vector<int, 3>>& elements,
    Field<Vector<int, 2>>& boundary_edges,
    Field<int>& boundary_points)
{
    find_boundary(elements, boundary_edges);
    std::unordered_set<int> boundary_point_set;
    for (const auto& edge : boundary_edges) {
        boundary_point_set.insert(edge[0]);
        boundary_point_set.insert(edge[1]);
    }
    boundary_points = Field<int>(boundary_point_set.begin(), boundary_point_set.end());
}

BOW_INLINE void find_boundary(const Field<Vector<int, 4>>& elements,
    Field<Vector<int, 3>>& boundary_faces,
    Field<Vector<int, 2>>& boundary_edges,
    Field<int>& boundary_points)
{
    find_boundary(elements, boundary_faces);
    std::unordered_set<int> boundary_point_set;
    std::unordered_set<Vector<int, 2>, VectorHash<2>> boundary_edge_set;
    for (const auto& face : boundary_faces) {
        boundary_point_set.insert(face[0]);
        boundary_point_set.insert(face[1]);
        boundary_point_set.insert(face[2]);
        Vector<int, 2> edge01(face[0], face[1]);
        std::sort(edge01.data(), edge01.data() + 2);
        Vector<int, 2> edge12(face[1], face[2]);
        std::sort(edge12.data(), edge12.data() + 2);
        Vector<int, 2> edge20(face[2], face[0]);
        std::sort(edge20.data(), edge20.data() + 2);
        boundary_edge_set.insert(edge01);
        boundary_edge_set.insert(edge12);
        boundary_edge_set.insert(edge20);
    }
    boundary_points = Field<int>(boundary_point_set.begin(), boundary_point_set.end());
    boundary_edges = Field<Vector<int, 2>>(boundary_edge_set.begin(), boundary_edge_set.end());
}

#ifdef BOW_STATIC_LIBRARY
void find_boundary(const Field<Vector<int, 4>>& elements, Field<Vector<int, 3>>& boundary_elements);
void find_boundary(const Field<Vector<int, 3>>& elements, Field<Vector<int, 2>>& boundary_elements);
#endif

} // namespace Bow::Geometry