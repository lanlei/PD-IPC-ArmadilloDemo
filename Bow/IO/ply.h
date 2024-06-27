#ifndef BOW_PLY_H
#define BOW_PLY_H
#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow {
namespace IO {
template <class T, int dim>
BOW_INLINE void read_ply(const std::string filename, Field<Vector<T, dim>>& vertices, Field<Vector<int, 3>>& faces);
template <class T, int dim>
BOW_INLINE void write_ply(const std::string filename, const Field<Vector<T, dim>>& vertices, const Field<Vector<int, 3>>& faces, const bool binary = true);
template <class T, int dim>
BOW_INLINE void write_ply(const std::string filename, const Field<Vector<T, dim>>& vertices, const bool binary = true);
}
} // namespace Bow::IO

#ifndef BOW_STATIC_LIBRARY
#include "ply.cpp"
#endif

#endif