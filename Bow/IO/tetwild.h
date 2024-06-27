#ifndef BOW_TETWILD
#define BOW_TETWILD
#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow {
namespace IO {
template <class T>
BOW_INLINE void read_mesh(const std::string filename, Field<Vector<T, 3>>& X, Field<Vector<int, 4>>& indices);
template <class T, int dim>
BOW_INLINE void write_meshing_data(const std::string filename, const Field<Vector<T, dim>>& X, const Field<Matrix<T, dim, dim>>& F);
}
} // namespace Bow::IO

#ifndef BOW_STATIC_LIBRARY
#include "tetwild.cpp"
#endif

#endif