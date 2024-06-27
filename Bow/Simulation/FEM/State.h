#ifndef FEM_STATE_H
#define FEM_STATE_H
#include <Bow/Macros.h>
#include <Bow/Types.h>
#include <Bow/Physics/ConstitutiveModel.h>
namespace Bow {
namespace FEM {
template <typename T, int dim>
class State {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    // constant
    T dt = 0.02;
    TV gravity = TV::Zero();
    Field<TV> m_X; // material coordinate;
    Field<T> m_mass; // node mass
    Field<T> m_density; // element density
    Field<T> m_vol; // element volume
    Field<Vector<int, dim + 1>> m_elem; // element vertex indices
    Field<TM> m_IB;
    Field<T> m_mu;
    Field<T> m_lam;
    std::unordered_map<int, std::vector<std::pair<int, int>>> m_obj_divider;

    // state
    Field<TV> m_x;
    Field<TV> m_v;
    Field<TV> m_a;

    // control
    Field<TM> BC_basis;
    std::vector<int> BC_order;
    Field<TV> BC_target;
    Field<TV> m_f;

    // intermediate
    Field<TV> m_x_tilde;

    State() {}
    ~State() {}
    /**
     * \brief add new elements to the simulation and prepare necessary data
     */
    BOW_INLINE void append(const Field<TV>& x, const Field<Vector<int, dim + 1>>& elem, const int type, const T E, const T nu, const T density);
    /**
        \brief Serialize x, v, boundary conditions. 
               For DiffFEM, serialize should be executed immediately after update_state, 
               i.e., we serialize s^{n+1} and \Lambda^n
    */
    BOW_INLINE void save_state(std::string filename);
    BOW_INLINE void load_state(std::string filename);
};
}
} // namespace Bow::FEM

#ifndef BOW_STATIC_LIBRARY
#include "State.cpp"
#endif

#endif