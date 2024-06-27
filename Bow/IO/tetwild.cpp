#include "tetwild.h"
#include <oneapi/tbb.h>
#include <sstream>
#include <fstream>

namespace Bow {
namespace IO {

template <class T>
BOW_INLINE void read_mesh(const std::string filename, Field<Vector<T, 3>>& X, Field<Vector<int, 4>>& indices)
{
    std::ifstream in;
    in.open(filename, std::ios::in);
    Field<Vector<int, 3>>* faces = nullptr;
    // auto initial_X_size = X.size();
    // auto initial_indices_size = indices.size();

    std::string line;
    Vector<T, 3> position;
    Vector<int, 4> tet;
    Vector<int, 3> face;

    bool reading_points = false;
    bool reading_faces = false;
    bool reading_tets = false;
    size_t n_points = 0;
    size_t n_faces = 0;
    size_t n_tets = 0;

    while (std::getline(in, line)) {
        std::stringstream ss(line);
        if (line.size() == (size_t)(0)) {
            // skip empty line
        }
        else if (line[0] == '#') {
            // skip comment line
        }
        else if (line.substr(0, 3) == "End") {
            break;
        }
        else if (line.substr(0, 20) == "MeshVersionFormatted") {
            ss.ignore(128, ' ');
            int mesh_ver_formatted;
            ss >> mesh_ver_formatted;
            assert(mesh_ver_formatted == 1);
        }
        else if (line.substr(0, 9) == "Dimension") {
            ss.ignore(128, ' ');
            int dimension;
            ss >> dimension;
            assert(dimension == 3);
        }
        else if (line.substr(0, 8) == "Vertices") {
            in >> n_points;
            reading_points = true;
            reading_faces = false;
            reading_tets = false;
        }
        else if (line.substr(0, 9) == "Triangles") {
            in >> n_faces;
            reading_points = false;
            reading_faces = true && (faces != nullptr);
            reading_tets = false;
        }
        else if (line.substr(0, 10) == "Tetrahedra") {
            in >> n_tets;
            reading_points = false;
            reading_faces = false;
            reading_tets = true;
        }
        else if (reading_points) {
            for (size_t i = 0; i < 3; i++)
                ss >> position[i];
            X.emplace_back(position);
            int end_mark;
            ss >> end_mark;
            assert(end_mark == -1 || end_mark == 0);
        }
        else if (reading_faces) {
            for (size_t i = 0; i < 3; i++)
                ss >> face[i];
            face.array() -= 1;
            faces->emplace_back(face);
            int end_mark;
            ss >> end_mark;
            assert(end_mark == -1 || end_mark == 0);
        }
        else if (reading_tets) {
            for (size_t i = 0; i < 4; i++)
                ss >> tet[i];
            tet.array() -= 1;
            indices.emplace_back(tet);
            int end_mark;
            ss >> end_mark;
            assert(end_mark == -1 || end_mark == 0);
        }
    }
    in.close();
    // ZIRAN_ASSERT(n_points == X.size() - initial_X_size, "mesh read X count doesn't match.");
    // ZIRAN_ASSERT((size_t)n_tets == indices.size() - initial_indices_size, "mesh read element count doesn't match.");
}

template <class T, int dim>
BOW_INLINE void write_meshing_data(const std::string filename, const Field<Vector<T, dim>>& X, const Field<Matrix<T, dim, dim>>& F)
{
    std::ofstream outfile;
    outfile.open(filename, std::ios::binary | std::ios::out);
    int nump = X.size();
    outfile.write((const char*)&nump, sizeof(int));

    std::cout << nump << " " << filename << std::endl;
    std::cout << filename << std::endl;
    for (int i = 0; i < nump; i++) {
        auto def = F[i];
        auto Xp = X[i];
        for (int q = 0; q < dim; q++) {
            float temp = Xp(q);

            outfile.write((const char*)&temp, sizeof(float));
        }
        for (int p = 0; p < dim; p++)
            for (int q = 0; q < dim; q++) {
                float temp = def(q, p);
                outfile.write((const char*)&temp, sizeof(float));
            }
    }
    outfile.close();
}
#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template void read_mesh(const std::string filename, Field<Vector<float, 3>>& X, Field<Vector<int, 4>>& indices);
template void write_meshing_data(const std::string filename, const Field<Vector<float, 2>>& X, const Field<Matrix<float, 2, 2>>& F);
template void write_meshing_data(const std::string filename, const Field<Vector<float, 3>>& X, const Field<Matrix<float, 3, 3>>& F);
#endif
#ifdef BOW_COMPILE_DOUBLE
template void read_mesh(const std::string filename, Field<Vector<double, 3>>& X, Field<Vector<int, 4>>& indices);
template void write_meshing_data(const std::string filename, const Field<Vector<double, 2>>& X, const Field<Matrix<double, 2, 2>>& F);
template void write_meshing_data(const std::string filename, const Field<Vector<double, 3>>& X, const Field<Matrix<double, 3, 3>>& F);
#endif
#endif
}
} // namespace Bow::IO
