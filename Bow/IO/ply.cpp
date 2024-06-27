#include "ply.h"
#include "tinyply.h"
#include <oneapi/tbb.h>
#include <sstream>
#include <fstream>
#include <Bow/Utils/Logging.h>

namespace Bow {
namespace IO {

namespace internal {
inline std::vector<uint8_t> read_file_binary(const std::string& pathToFile)
{
    std::ifstream file(pathToFile, std::ios::binary);
    std::vector<uint8_t> fileBufferBytes;

    if (file.is_open()) {
        file.seekg(0, std::ios::end);
        size_t sizeBytes = file.tellg();
        file.seekg(0, std::ios::beg);
        fileBufferBytes.resize(sizeBytes);
        if (file.read((char*)fileBufferBytes.data(), sizeBytes)) return fileBufferBytes;
    }
    else
        throw std::runtime_error("could not open binary ifstream to path " + pathToFile);
    return fileBufferBytes;
}

struct memory_buffer : public std::streambuf {
    char* p_start{ nullptr };
    char* p_end{ nullptr };
    size_t size;

    memory_buffer(char const* first_elem, size_t size)
        : p_start(const_cast<char*>(first_elem)), p_end(p_start + size), size(size)
    {
        setg(p_start, p_start, p_end);
    }

    pos_type seekoff(off_type off, std::ios_base::seekdir dir, std::ios_base::openmode which) override
    {
        if (dir == std::ios_base::cur)
            gbump(static_cast<int>(off));
        else
            setg(p_start, (dir == std::ios_base::beg ? p_start : p_end) + off, p_end);
        return gptr() - p_start;
    }

    pos_type seekpos(pos_type pos, std::ios_base::openmode which) override
    {
        return seekoff(pos, std::ios_base::beg, which);
    }
};

struct memory_stream : virtual memory_buffer, public std::istream {
    memory_stream(char const* first_elem, size_t size)
        : memory_buffer(first_elem, size), std::istream(static_cast<std::streambuf*>(this)) {}
};

} // namespace internal

template <class T, int dim>
BOW_INLINE void read_ply(const std::string filepath, Field<Vector<T, dim>>& vertices_out, Field<Vector<int, 3>>& faces_out)
{
    Logging::info("Reading: ", filepath);

    std::unique_ptr<std::istream> file_stream;
    std::vector<uint8_t> byte_buffer;

    try {
        // For most files < 1gb, pre-loading the entire file upfront and wrapping it into a
        // stream is a net win for parsing speed, about 40% faster.
        byte_buffer = internal::read_file_binary(filepath);
        file_stream.reset(new internal::memory_stream((char*)byte_buffer.data(), byte_buffer.size()));

        if (!file_stream || file_stream->fail()) throw std::runtime_error("file_stream failed to open " + filepath);

        tinyply::PlyFile file;
        file.parse_header(*file_stream);
        Logging::info("\t[ply_header] Type: ", (file.is_binary_file() ? "binary" : "ascii"));

        // Because most people have their own mesh types, tinyply treats parsed data as structured/typed byte buffers.
        // See examples below on how to marry your own application-specific data structures with this one.
        std::shared_ptr<tinyply::PlyData> vertices, faces;

        // The header information can be used to programmatically extract properties on elements
        // known to exist in the header prior to reading the data. For brevity of this sample, properties
        // like vertex position are hard-coded:
        try {
            vertices = file.request_properties_from_element("vertex", { "x", "y", "z" });
        }
        catch (const std::exception& e) {
            Logging::error("tinyply exception: ", e.what());
        }

        // Providing a list size hint (the last argument) is a 2x performance improvement. If you have
        // arbitrary ply files, it is best to leave this 0.
        try {
            faces = file.request_properties_from_element("face", { "vertex_indices" }, 3);
        }
        catch (const std::exception& e) {
            Logging::error("tinyply exception: ", e.what());
        }

        file.read(*file_stream);

        if (vertices) {
            Logging::info("\tRead ", vertices->count, " total vertices ");
            vertices_out.resize(vertices->count);
            if (vertices->t == tinyply::Type::FLOAT32) {
                float* verts_data = reinterpret_cast<float*>(vertices->buffer.get());
                tbb::parallel_for(size_t(0), vertices->count, [&](size_t i) {
                    for (int d = 0; d < dim; ++d)
                        vertices_out[i](d) = verts_data[3 * i + d];
                });
            }
            else if (vertices->t == tinyply::Type::FLOAT64) {
                double* verts_data = reinterpret_cast<double*>(vertices->buffer.get());
                tbb::parallel_for(size_t(0), vertices->count, [&](size_t i) {
                    for (int d = 0; d < dim; ++d)
                        vertices_out[i](d) = verts_data[3 * i + d];
                });
            }
        }
        if (faces) {
            Logging::info("\tRead ", faces->count, " total faces (triangles) ");
            faces_out.resize(faces->count);
            if (faces->t == tinyply::Type::UINT32) {
                int* facess_data = reinterpret_cast<int*>(faces->buffer.get());
                tbb::parallel_for(size_t(0), faces->count, [&](size_t i) {
                    for (int d = 0; d < 3; ++d)
                        faces_out[i](d) = facess_data[3 * i + d];
                });
            }
            else if (faces->t == tinyply::Type::INT32) {
                int32_t* facess_data = reinterpret_cast<int32_t*>(faces->buffer.get());
                tbb::parallel_for(size_t(0), faces->count, [&](size_t i) {
                    for (int d = 0; d < 3; ++d)
                        faces_out[i](d) = facess_data[3 * i + d];
                });
            }
        }
    }
    catch (const std::exception& e) {
        Logging::error("Caught tinyply exception: ", e.what());
    }
} // namespace IO

template <class T, int dim>
BOW_INLINE void write_ply(const std::string filename, const Field<Vector<T, dim>>& _vertices, const Field<Vector<int, 3>>& faces, const bool binary)
{
    Logging::info("Writing: ", filename);
    Field<Vector<T, 3>> vertices(_vertices.size());
    if constexpr (dim == 2) {
        tbb::parallel_for((size_t)0, _vertices.size(), [&](size_t i) {
            vertices[i] = Vector<T, 3>(_vertices[i][0], _vertices[i][1], 0.0);
        });
    }
    else {
        tbb::parallel_for((size_t)0, _vertices.size(), [&](size_t i) {
            vertices[i] = _vertices[i];
        });
    }
    std::ofstream outstream_binary;
    if (binary)
        outstream_binary.open(filename, std::ios::out | std::ios::binary);
    else
        outstream_binary.open(filename, std::ios::out);
    tinyply::PlyFile file;
    if (std::is_same<T, float>::value)
        file.add_properties_to_element("vertex", { "x", "y", "z" },
            tinyply::Type::FLOAT32, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), tinyply::Type::INVALID, 0);
    else
        file.add_properties_to_element("vertex", { "x", "y", "z" },
            tinyply::Type::FLOAT64, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), tinyply::Type::INVALID, 0);
    file.add_properties_to_element("face", { "vertex_indices" },
        tinyply::Type::UINT32, faces.size(), reinterpret_cast<const uint8_t*>(faces.data()), tinyply::Type::UINT8, 3);
    if (binary)
        file.write(outstream_binary, true);
    else
        file.write(outstream_binary, false);
}
template <class T, int dim>
BOW_INLINE void write_ply(const std::string filename, const Field<Vector<T, dim>>& _vertices, const bool binary)
{
    Logging::info("Writing: ", filename);
    Field<Vector<T, 3>> vertices(_vertices.size());
    if constexpr (dim == 2) {
        tbb::parallel_for((size_t)0, _vertices.size(), [&](size_t i) {
            vertices[i] = Vector<T, 3>(_vertices[i][0], _vertices[i][1], 0.0);
        });
    }
    else {
        tbb::parallel_for((size_t)0, _vertices.size(), [&](size_t i) {
            vertices[i] = _vertices[i];
        });
    }

    std::ofstream outstream_binary;
    if (binary)
        outstream_binary.open(filename, std::ios::out | std::ios::binary);
    else
        outstream_binary.open(filename, std::ios::out);

    tinyply::PlyFile file;
    if constexpr (std::is_same<T, float>::value)
        file.add_properties_to_element("vertex", { "x", "y", "z" },
            tinyply::Type::FLOAT32, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), tinyply::Type::INVALID, 0);
    else
        file.add_properties_to_element("vertex", { "x", "y", "z" },
            tinyply::Type::FLOAT64, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), tinyply::Type::INVALID, 0);
    if (binary)
        file.write(outstream_binary, true);
    else
        file.write(outstream_binary, false);
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template void read_ply(const std::string filename, Field<Vector<float, 2>>& vertices, Field<Vector<int, 3>>& faces);
template void read_ply(const std::string filename, Field<Vector<float, 3>>& vertices, Field<Vector<int, 3>>& faces);
template void write_ply(const std::string filename, const Field<Vector<float, 2>>& vertices, const Field<Vector<int, 3>>& faces, const bool);
template void write_ply(const std::string filename, const Field<Vector<float, 3>>& vertices, const Field<Vector<int, 3>>& faces, const bool);
template void write_ply(const std::string filename, const Field<Vector<float, 2>>& vertices, const bool);
template void write_ply(const std::string filename, const Field<Vector<float, 3>>& vertices, const bool);
#endif
#ifdef BOW_COMPILE_DOUBLE
template void read_ply(const std::string filename, Field<Vector<double, 2>>& vertices, Field<Vector<int, 3>>& faces);
template void read_ply(const std::string filename, Field<Vector<double, 3>>& vertices, Field<Vector<int, 3>>& faces);
template void write_ply(const std::string filename, const Field<Vector<double, 2>>& vertices, const Field<Vector<int, 3>>& faces, const bool);
template void write_ply(const std::string filename, const Field<Vector<double, 3>>& vertices, const Field<Vector<int, 3>>& faces, const bool);
template void write_ply(const std::string filename, const Field<Vector<double, 2>>& vertices, const bool);
template void write_ply(const std::string filename, const Field<Vector<double, 3>>& vertices, const bool);
#endif
#endif
}
} // namespace Bow::IO
