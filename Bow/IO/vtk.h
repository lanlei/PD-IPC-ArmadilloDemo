#ifndef BOW_VTK
#define BOW_VTK
#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow {
namespace IO {
template <class T>
BOW_INLINE void write_vtk(const std::string filename, const Field<Vector<T, 3>>& xyz, const Field<Vector<int, 4>>& cells, const bool binary = true);
}
} // namespace Bow::IO

#ifndef BOW_STATIC_LIBRARY
#include "vtk.cpp"
#endif

#endif