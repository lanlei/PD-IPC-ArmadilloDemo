#pragma once

#include <Bow/Types.h>
#include <Bow/Simulation/MPM/MPMGrid.h>
#include <Bow/Utils/Timer.h>
#include <MarchingCubes/LocallyConsistentIsocontour.h>
#include <algorithm>
#include <utility>
#include <queue>

namespace Bow {
namespace MPM {

template <class T, int dim>
class ParticlesLevelSet {
public:
    using SparseMask = typename MPMGrid<T, dim>::SparseMask;
    using TV = Vector<T, dim>;
    using IV = Vector<int, dim>;
    using TM = Matrix<T, dim, dim>;

    std::priority_queue<std::pair<T, std::array<int, dim>>> pq;
    MPMGrid<T, dim> phi_grid;
    T h;

    ParticlesLevelSet() = default;

    void insertHeap(IV node)
    {
        if (phi_grid[node].v_and_m(dim) > 0)
            return;
        T min_value[dim];
        for (int d = 0; d < dim; ++d) {
            min_value[d] = std::numeric_limits<T>::max();
            for (int delta = -1; delta <= 1; delta += 2) {
                IV new_node = node;
                new_node[d] += delta;
                if (phi_grid[new_node].v_and_m(dim) > 0)
                    min_value[d] = std::min(min_value[d], phi_grid[new_node].v_and_m(0));
            }
        }
        std::sort(min_value, min_value + dim);
        if (min_value[0] > std::numeric_limits<T>::max() / (T)2)
            return;
        T result = std::numeric_limits<T>::max();
        if (dim == 2) {
            T a = min_value[0], b = min_value[1];
            result = std::min(result, a + h);
            if (a + h > b) {
                T ab = ((a + b) + std::sqrt(2.0 * h * h - (a - b) * (a - b))) / 2.0;
                result = std::min(result, ab);
            }
        }
        if (dim == 3) {
            T a = min_value[0], b = min_value[1], c = min_value[2];
            result = std::min(result, a + h);
            if (a + h > b) {
                T ab = ((a + b) + std::sqrt(2.0 * h * h - (a - b) * (a - b))) / 2.0;
                result = std::min(result, ab);
                if (ab > c) {
                    T abc = 1.0 / 6.0 * (std::sqrt((-2.0 * a - 2.0 * b - 2.0 * c) * (-2.0 * a - 2.0 * b - 2.0 * c) - 12.0 * (a * a + b * b + c * c - h * h)) + 2.0 * a + 2.0 * b + 2.0 * c);
                    result = std::min(result, abc);
                }
            }
        }
        pq.push(std::make_pair(-result, to_std_array<int, dim>(node.data())));
    }

    void build(const Field<TV>& particles_X, T build_dx, T blob_radius, T erosion, T smooth_epsilon = 0)
    {
        BOW_TIMER_FLAG("build particles level set");

        h = build_dx;
        phi_grid.sortParticles(particles_X, h);
        int kernel_size = ((int)std::ceil(blob_radius / h) + 1) * 2;
        phi_grid.colored_for([&](int i) {
            auto& Xp = particles_X[i];
            // Must find base_node with linear kernel
            BSplineWeights<T, dim, 1> spline(Xp, h);
            for (int d = 0; d < dim; ++d) {
                BOW_ASSERT_INFO(spline.base_node[d] + kernel_size / 2 < phi_grid.spgrid_size, "Exceed SPGrid range");
            }
            IV region = IV::Ones() * kernel_size;
            iterateRegion(region, [&](const IV& offset) {
                IV node = spline.base_node + offset - IV::Ones() * (kernel_size / 2 - 1);
                TV pos = node.template cast<T>() * h;
                T phi = (Xp - pos).norm() - blob_radius;
                if (phi_grid[node].v_and_m(dim) > 0) {
                    phi_grid[node].v_and_m(0) = std::max(-phi, phi_grid[node].v_and_m(0));
                }
                else {
                    phi_grid[node].v_and_m(0) = -phi;
                    phi_grid[node].v_and_m(dim) = 1;
                }
            });
        });
        phi_grid.countNumNodes();
        Field<IV> ready_to_insert;
        phi_grid.iterateGridSerial([&](IV node, GridState<T, dim>& g) {
            if (g.v_and_m(0) > 0) {
                g.v_and_m(0) = 0;
                g.v_and_m(dim) = 0;
                ready_to_insert.push_back(node);
            }
        });
        for (auto& node : ready_to_insert) {
            insertHeap(node);
        }
        while (!pq.empty()) {
            auto item = pq.top();
            pq.pop();
            IV node = IV(item.second.data());
            if (phi_grid[node].v_and_m(dim) == 0) {
                phi_grid[node].v_and_m(0) = -item.first;
                phi_grid[node].v_and_m(dim) = 1;
                for (int d = 0; d < dim; ++d) {
                    for (int delta = -1; delta <= 1; delta += 2) {
                        IV new_node = node;
                        new_node[d] += delta;
                        insertHeap(new_node);
                    }
                }
            }
        }
        phi_grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.v_and_m(0) = -g.v_and_m(0) + erosion;
            if (smooth_epsilon) {
                g.v_and_m(0) = (T)-1 / (1 + std::exp(g.v_and_m(0) / smooth_epsilon)) + 0.5;
            }
        });
    }

    T queryPhi(TV& particleX)
    {
        BSplineWeights<T, dim> spline(particleX, h);
        T phi = 0;
        phi_grid.iterateKernel(spline, [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
            phi += w * g.v_and_m(0);
        });
        return phi;
    }

    TV queryNabla(TV& particleX)
    {
        BSplineWeights<T, dim> spline(particleX, h);
        TV nabla = TV::Zero();
        phi_grid.iterateKernel(spline, [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
            nabla += dw * g.v_and_m(0);
        });
        return nabla;
    }

    bool queryInside(TV& particleX, T& phi, TV& normal, T isocontour)
    {
        BSplineWeights<T, dim> spline(particleX, h);
        // firstly check inside valid blocks
        if (!phi_grid.existsKernel(spline))
            return false;
        int cnt = 0;
        // secondly check inside valid grids
        phi_grid.iterateKernel(spline, [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
            if (g.v_and_m(dim) > 0) ++cnt;
        });
        if (cnt != (1 << dim)) return false;
        phi = 0;
        normal = TV::Zero();
        phi_grid.iterateKernel(spline, [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
            phi += w * g.v_and_m(0);
            normal += dw * g.v_and_m(0);
        });
        normal.normalize();
        return phi < isocontour;
    }
};

template <class T, int dim>
class DumpSurfaceMesh {
public:
    Field<Vector<T, dim>>& m_X;
    T& dx;
    std::string filename;

    Field<Vector<T, dim>> ply_vertices;
    Field<Vector<int, 3>> ply_faces;

    DumpSurfaceMesh(Field<Vector<T, dim>>& m_X, T& dx)
        : m_X(m_X), dx(dx) {}

    void operator()()
    {
        BOW_TIMER_FLAG("dump surface mesh");
        ParticlesLevelSet<T, dim> builder;
        builder.build(m_X, dx, dx, 0.6 * dx);
        MPMGrid<T, dim>& grid = builder.phi_grid;

        ply_vertices.clear();
        ply_faces.clear();
        grid.iterateWholeGridSerial([&](Vector<int, dim> node, GridState<T, dim>& g) {
            Vector<int, dim> region = Vector<int, dim>::Ones() * 2;
            int pos_cnt = 0;
            int neg_cnt = 0;
            std::array<double, 1 << dim> scalar_field;
            std::vector<std::array<double, 3>> nodes;
            iterateRegion(region, [&](const Vector<int, dim>& offset) {
                if (grid.existsNode(node + offset) && grid[node + offset].v_and_m(dim)) {
                    T dist = grid[node + offset].v_and_m(0);
                    if constexpr (dim == 2) {
                        int b = offset[0] * 2 + offset[1];
                        scalar_field[b] = dist;
                    }
                    else {
                        int b = offset[0] * 4 + offset[1] * 2 + offset[2];
                        scalar_field[b] = dist;
                    }
                    if (dist > 0)
                        pos_cnt++;
                    else
                        neg_cnt++;
                }
            });
            if (pos_cnt + neg_cnt == (1 << dim) && pos_cnt && neg_cnt) {
                TGSL::MarchingCubes::CellIsocontour(scalar_field, nodes);
                if constexpr (dim == 2) {
                    BOW_NOT_IMPLEMENTED
                }
                else {
                    int start = (int)ply_vertices.size();
                    for (auto& i : nodes) {
                        ply_vertices.push_back((Vector<T, dim>(i.data()) + node.template cast<T>()) * dx);
                        // Logging::info("Node:", i[0], " ", i[1], " ", i[2]);
                    }
                    for (int i = 0; i < (int)nodes.size() / 3; ++i) {
                        ply_faces.push_back(Vector<int, dim>(start, start + 1, start + 2));
                        start += 3;
                    }
                }
            }
        });
        IO::write_ply(filename, ply_vertices, ply_faces);
    }
};

}
} // namespace Bow::MPM