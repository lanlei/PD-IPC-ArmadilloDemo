#include "State.h"
#include "Utils.h"
#include <Bow/Optimization/Newton.h>
#include <iostream>
#include <Bow/Utils/Timer.h>
#include <Bow/Physics/FixedCorotated.h>
#include <oneapi/tbb.h>
#include <fstream>

namespace Bow {
namespace FEM {
template <class T, int dim>
BOW_INLINE void State<T, dim>::append(const Field<TV>& x, const Field<Vector<int, dim + 1>>& elem, const int type, const T E, const T nu, const T density)
{
    Field<Vector<int, dim + 1>> shifted_elem = elem;
    Eigen::Map<Matrix<int, dim + 1, Eigen::Dynamic>>(&(shifted_elem[0][0]), dim + 1, shifted_elem.size()).array() += m_X.size();
    T mu, lam;
    std::tie(mu, lam) = ConstitutiveModel::lame_paramters(E, nu);

    m_density.resize(m_density.size() + elem.size(), density);
    m_X.insert(m_X.end(), x.begin(), x.end());
    m_elem.insert(m_elem.end(), shifted_elem.begin(), shifted_elem.end());
    m_mu.resize(m_mu.size() + elem.size(), mu);
    m_lam.resize(m_lam.size() + elem.size(), lam);
    m_obj_divider[type].push_back(std::make_pair(m_elem.size() - elem.size(), m_elem.size()));

    m_x.insert(m_x.end(), x.begin(), x.end());
    m_v.resize(m_v.size() + x.size(), Vector<T, dim>::Zero());
    m_a.resize(m_a.size() + x.size(), Vector<T, dim>::Zero());
    m_x_tilde.insert(m_x_tilde.end(), x.begin(), x.end());
    BC_basis.resize(m_x.size(), TM::Identity());
    BC_order.resize(m_x.size(), 0);
    BC_target.resize(m_x.size(), TV::Zero());
    m_f.resize(m_x.size(), TV::Zero());
}

template <class T, int dim>
BOW_INLINE void State<T, dim>::save_state(std::string filename)
{
    std::ofstream outstream;
    outstream.open(filename, std::ios::binary);
    outstream.write(reinterpret_cast<char*>(m_x.data()), m_x.size() * dim * sizeof(T));
    outstream.write(reinterpret_cast<char*>(m_x_tilde.data()), m_x_tilde.size() * dim * sizeof(T));
    outstream.write(reinterpret_cast<char*>(m_v.data()), m_v.size() * dim * sizeof(T));
    outstream.write(reinterpret_cast<char*>(m_a.data()), m_a.size() * dim * sizeof(T));
    outstream.write(reinterpret_cast<char*>(m_f.data()), m_f.size() * dim * sizeof(T));
    outstream.write(reinterpret_cast<char*>(BC_basis.data()), BC_basis.size() * dim * dim * sizeof(T));
    outstream.write(reinterpret_cast<char*>(BC_order.data()), BC_order.size() * sizeof(int));
    outstream.write(reinterpret_cast<char*>(BC_target.data()), BC_target.size() * dim * sizeof(T));
    outstream.close();
}

template <class T, int dim>
BOW_INLINE void State<T, dim>::load_state(std::string filename)
{
    std::ifstream instream;
    instream.open(filename, std::ios::binary);
    m_x.resize(m_X.size());
    m_x_tilde.resize(m_X.size());
    m_v.resize(m_X.size());
    m_a.resize(m_X.size());
    m_f.resize(m_X.size());
    BC_basis.resize(m_X.size());
    BC_order.resize(m_X.size());
    BC_target.resize(m_X.size());
    instream.read(reinterpret_cast<char*>(m_x.data()), m_x.size() * dim * sizeof(T));
    instream.read(reinterpret_cast<char*>(m_x_tilde.data()), m_x_tilde.size() * dim * sizeof(T));
    instream.read(reinterpret_cast<char*>(m_v.data()), m_v.size() * dim * sizeof(T));
    instream.read(reinterpret_cast<char*>(m_a.data()), m_a.size() * dim * sizeof(T));
    instream.read(reinterpret_cast<char*>(m_f.data()), m_f.size() * dim * sizeof(T));
    instream.read(reinterpret_cast<char*>(BC_basis.data()), BC_basis.size() * dim * dim * sizeof(T));
    instream.read(reinterpret_cast<char*>(BC_order.data()), BC_order.size() * sizeof(int));
    instream.read(reinterpret_cast<char*>(BC_target.data()), BC_target.size() * dim * sizeof(T));
    instream.close();
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template class State<float, 2>;
#endif
#ifdef BOW_COMPILE_3D
template class State<float, 3>;
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template class State<double, 2>;
#endif
#ifdef BOW_COMPILE_3D
template class State<double, 3>;
#endif
#endif
#endif
}
} // namespace Bow::FEM