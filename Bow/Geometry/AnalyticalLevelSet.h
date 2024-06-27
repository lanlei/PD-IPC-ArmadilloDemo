#ifndef ANALYTICAL_LEVEL_SET_H
#define ANALYTICAL_LEVEL_SET_H

#include <Bow/Types.h>
#include <Bow/Utils/Logging.h>
#include <limits>
#include <memory>

namespace Bow::Geometry {

enum Type { STICKY,
    SLIP,
    SEPARATE };

template <class T, int dim>
class AnalyticalLevelSet {
    using TV = Vector<T, dim>;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Type type;
    explicit AnalyticalLevelSet(Type type)
        : type(type) {}
    virtual T signed_distance(const TV& X) = 0;
    virtual TV normal(const TV& X) { BOW_NOT_IMPLEMENTED return TV::Unit(0); }
    virtual TV velocity(const TV& X) { return TV::Zero(); }
};

template <class T, int dim>
class HalfSpaceLevelSet : public AnalyticalLevelSet<T, dim> {
    using TV = Vector<T, dim>;
    TV origin;
    TV outward_normal;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    HalfSpaceLevelSet(Type type, const TV& origin, const TV& normal);
    T signed_distance(const TV& X) override;
    TV normal(const TV& X) override;
};

template <class T, int dim>
class AlignedBoxLevelSet : public AnalyticalLevelSet<T, dim> {
    using TV = Vector<T, dim>;
    TV center;
    TV half_edges;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    AlignedBoxLevelSet(Type type, const TV& min_corner, const TV& max_corner);
    T signed_distance(const TV& X) override;
};

template <class T, int dim>
class ScriptedLevelSet : public AnalyticalLevelSet<T, dim> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <class T, int dim>
class IndexBasedBounaryCondition {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using TV = Vector<T, dim>;
    int index;
    Type type;
    TV m_n;
    T m_v;
    IndexBasedBounaryCondition(int index, Type type, const TV& normal, const T velocity)
        : index(index), type(type), m_n(normal.normalized()), m_v(velocity) {}
    TV normal() { return m_n; }
    TV velocity() { return m_v * m_n; }
};

template <class T, int dim>
class BoundaryConditionManager {
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    // levelset based
    std::vector<std::shared_ptr<AnalyticalLevelSet<T, dim>>> level_set_objects;
    // index based
    std::map<int, std::vector<std::shared_ptr<IndexBasedBounaryCondition<T, dim>>>> node_bcs;

public:
    void clear()
    {
        level_set_objects.clear();
        node_bcs.clear();
    }
    void add(std::shared_ptr<AnalyticalLevelSet<T, dim>> ls_ptr)
    {
        level_set_objects.push_back(ls_ptr);
    }
    void add(std::shared_ptr<IndexBasedBounaryCondition<T, dim>> node_bc_ptr)
    {
        node_bcs[node_bc_ptr->index].push_back(node_bc_ptr);
    }
    void update(T time)
    {
        BOW_NOT_IMPLEMENTED
    }
    void mpm_explicit_update(const TV& X, TV& V)
    {
        for (auto& ls_ptr : level_set_objects)
            if (ls_ptr->signed_distance(X) < 0) {
                TV v = ls_ptr->velocity(X);
                TV n = ls_ptr->normal(X);
                if (ls_ptr->type == STICKY) {
                    V = v;
                }
                if (ls_ptr->type == SLIP) {
                    V -= n.dot(V - v) * n + v;
                }
                if (ls_ptr->type == SEPARATE) {
                    T dot_value = n.dot(V - v);
                    if (dot_value < 0)
                        V -= dot_value * n + v;
                }
            }
    }
    void mpm_implicit_update(const TV& X, TM& basis, int& order)
    {
        basis = Matrix<T, dim, dim>::Identity();
        order = 0;
        for (auto& ls_ptr : level_set_objects)
            if (ls_ptr->signed_distance(X) < 0) {
                TV v = ls_ptr->velocity(X);
                TV n = ls_ptr->normal(X);
                if (v.norm()) { BOW_NOT_IMPLEMENTED }
                if (ls_ptr->type == STICKY) {
                    basis = Matrix<T, dim, dim>::Identity();
                    order = dim;
                    return;
                }
                if (ls_ptr->type == SLIP) {
                    for (int i = 0; i < order; ++i)
                        n -= n.dot(basis.col(i)) * basis.col(i);
                    if (n.norm()) {
                        n.normalize();
                        basis.col(order++) = n;
                    }
                }
                if (ls_ptr->type == SEPARATE) {
                    BOW_NOT_IMPLEMENTED
                }
            }
        // no need to deal with order == 0
        if constexpr (dim == 2) {
            if (order == 1) {
                basis(0, 1) = basis(1, 0);
                basis(1, 1) = -basis(0, 0);
            }
        }
        if constexpr (dim == 3) {
            if (order == 1) {
                if (basis.col(0).dot(TV(1, 0, 0)) > 0.5)
                    basis.col(1) = basis.col(0).cross(TV(0, 1, 0));
                else
                    basis.col(1) = basis.col(0).cross(TV(1, 0, 0));
                basis.col(1).normalize();
            }
            if (order == 1 || order == 2) {
                basis.col(2) = basis.col(0).cross(basis.col(1));
                basis.col(2).normalize();
            }
        }
    }

    bool level_set_based_update(const TV& x, const T dt, TM& basis, int& order, TV& target_after_transform)
    {
        basis = Matrix<T, dim, dim>::Identity();
        order = 0;
        target_after_transform = x;
        bool valid = false;
        for (auto& ls_ptr : level_set_objects)
            if (ls_ptr->signed_distance(x) < 0) {
                valid = true;
                TV v = ls_ptr->velocity(x);
                TV n = ls_ptr->normal(x);
                if (v.norm()) { BOW_NOT_IMPLEMENTED }
                if (ls_ptr->type == STICKY) {
                    basis = Matrix<T, dim, dim>::Identity();
                    order = dim;
                    return true;
                }
                if (ls_ptr->type == SLIP) {
                    for (int i = 0; i < order; ++i)
                        n -= n.dot(basis.col(i)) * basis.col(i);
                    if (n.norm()) {
                        n.normalize();
                        basis.col(order++) = n;
                    }
                }
                if (ls_ptr->type == SEPARATE) {
                    BOW_NOT_IMPLEMENTED
                }
                target_after_transform += dt * ls_ptr->velocity(x);
            }
        // no need to deal with order == 0
        if constexpr (dim == 2) {
            if (order == 1) {
                basis(0, 1) = basis(1, 0);
                basis(1, 1) = -basis(0, 0);
            }
        }
        if constexpr (dim == 3) {
            if (order == 1) {
                if (basis.col(0).dot(TV(1, 0, 0)) > 0.5)
                    basis.col(1) = basis.col(0).cross(TV(0, 1, 0));
                else
                    basis.col(1) = basis.col(0).cross(TV(1, 0, 0));
                basis.col(1).normalize();
            }
            if (order == 1 || order == 2) {
                basis.col(2) = basis.col(0).cross(basis.col(1));
                basis.col(2).normalize();
            }
        }
        target_after_transform = basis.transpose() * target_after_transform;
        return valid;
    }

    /** Must be called after level_set_based_update */
    bool index_based_update(const int node, const TV& x, const T dt, TM& basis, int& order, TV& target_after_transform)
    {
        if (node_bcs.find(node) == node_bcs.end()) return false;
        target_after_transform = basis * target_after_transform; // transform back for further editting.
        for (auto bc : node_bcs[node]) {
            TV v = bc->velocity();
            TV n = bc->normal();
            if (v.norm()) { BOW_NOT_IMPLEMENTED }
            if (bc->type == STICKY) {
                basis = Matrix<T, dim, dim>::Identity();
                order = dim;
                return true;
            }
            if (bc->type == SLIP) {
                for (int i = 0; i < order; ++i)
                    n -= n.dot(basis.col(i)) * basis.col(i);
                if (n.norm()) {
                    n.normalize();
                    basis.col(order++) = n;
                }
            }
            if (bc->type == SEPARATE) {
                BOW_NOT_IMPLEMENTED
            }
            target_after_transform += dt * v;
        }
        // no need to deal with order == 0
        if constexpr (dim == 2) {
            if (order == 1) {
                basis(0, 1) = basis(1, 0);
                basis(1, 1) = -basis(0, 0);
            }
        }
        if constexpr (dim == 3) {
            if (order == 1) {
                if (basis.col(0).dot(TV(1, 0, 0)) > 0.5)
                    basis.col(1) = basis.col(0).cross(TV(0, 1, 0));
                else
                    basis.col(1) = basis.col(0).cross(TV(1, 0, 0));
                basis.col(1).normalize();
            }
            if (order == 1 || order == 2) {
                basis.col(2) = basis.col(0).cross(basis.col(1));
                basis.col(2).normalize();
            }
        }

        target_after_transform = basis.transpose() * target_after_transform;
        return true;
    }
};

} // namespace Bow::Geometry

#ifndef BOW_STATIC_LIBRARY
#include "AnalyticalLevelSet.cpp"
#endif

#endif