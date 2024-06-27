#include "vtk.h"
#include <oneapi/tbb.h>
#include <sstream>
#include <fstream>

namespace Bow {
namespace IO {

namespace internal {
template <class T>
inline void swap_endian(std::vector<T>& array)
{
    // vtk binary assumes big-endian byte order
    tbb::parallel_for(size_t(0), array.size(), [&](size_t i) {
        char* data = reinterpret_cast<char*>(&array[i]);
        for (long sub_i = 0; sub_i < static_cast<long>(sizeof(T) / 2); sub_i++)
            std::swap(data[sizeof(T) - 1 - sub_i], data[sub_i]);
    });
}
} // namespace internal

template <class T>
BOW_INLINE void write_vtk(const std::string filename, const Field<Vector<T, 3>>& _xyz, const Field<Vector<int, 4>>& _cells, const bool binary)
{
    int nPoints = _xyz.size();
    int nCells = _cells.size();
    std::vector<T> xyz(3 * nPoints);
    std::vector<int> cells(5 * nCells);
    std::vector<int> cell_types(nCells);
    tbb::parallel_for(0, nPoints, [&](int i) {
        for (int d = 0; d < 3; ++d) {
            xyz[3 * i + d] = _xyz[i][d];
        }
    });
    tbb::parallel_for(0, nCells, [&](int i) {
        cells[5 * i] = 4;
        for (int d = 0; d < 4; ++d) {
            cells[5 * i + d + 1] = _cells[i][d];
        }
        cell_types[i] = 10;
    });
    std::ofstream outstream;
    if (binary) {
        internal::swap_endian(xyz);
        internal::swap_endian(cells);
        internal::swap_endian(cell_types);
        outstream.open(filename, std::ios::out | std::ios::binary);
    }
    else
        outstream.open(filename, std::ios::out);

    if (outstream.fail()) throw std::runtime_error("failed to open " + filename);
    // const std::locale & fixLoc = std::locale("C");
    // outstream_binary.imbue(fixLoc);
    const std::locale& fixLoc = std::locale("C");
    outstream.imbue(fixLoc);
    outstream << "# vtk DataFile Version 2.0\n";
    outstream << "Visulaization output file\n";
    if (binary)
        outstream << "BINARY\n";
    else
        outstream << "ASCII\n";
    outstream << "DATASET UNSTRUCTURED_GRID\n";
    outstream << "POINTS " << nPoints;
    if (std::is_same<T, double>::value)
        outstream << " double\n";
    else
        outstream << " float\n";

    if (binary)
        outstream.write(reinterpret_cast<const char*>(xyz.data()), 3 * nPoints * sizeof(T));
    else {
        outstream << xyz[0];
        for (int i = 1; i < 3 * nPoints; ++i)
            outstream << " " << xyz[i];
    }
    outstream << "\n";
    outstream << "CELLS " << nCells << " " << 5 * nCells << "\n";
    if (binary)
        outstream.write(reinterpret_cast<char*>(cells.data()), 5 * nCells * sizeof(int));
    else {
        outstream << cells[0];
        for (int i = 1; i < 5 * nCells; ++i)
            outstream << " " << cells[i];
    }
    outstream << "\n";
    outstream << "CELL_TYPES " << nCells << "\n";
    if (binary)
        outstream.write(reinterpret_cast<char*>(cell_types.data()), nCells * sizeof(int));
    else {
        outstream << cell_types[0];
        for (int i = 1; i < nCells; ++i)
            outstream << " " << cell_types[i];
    }
    outstream << "\n";
    outstream.close();
}
#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template void write_vtk(const std::string filename, const Field<Vector<float, 3>>& xyz, const Field<Vector<int, 4>>& cells, const bool binary);
#endif
#ifdef BOW_COMPILE_DOUBLE
template void write_vtk(const std::string filename, const Field<Vector<double, 3>>& xyz, const Field<Vector<int, 4>>& cells, const bool binary);
#endif
#endif
}
} // namespace Bow::IO
