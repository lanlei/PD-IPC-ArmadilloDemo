#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <Bow/Utils/Logging.h>
#include <fstream>

namespace Bow {
namespace Serialization {

class Serializer {
public:
    virtual void load(std::ifstream& fin) = 0;
    virtual void save(std::ofstream& fout) = 0;
};
static std::vector<Serializer*> serializers;
static bool prepare_phase = true;

template <class T>
class SerializerImpl : public Serializer {
    T& ref;

public:
    explicit SerializerImpl(T& ref)
        : ref(ref)
    {
        BOW_ASSERT(prepare_phase)
        serializers.push_back(this);
    }
    ~SerializerImpl()
    {
        for (size_t i = 0; i < serializers.size(); ++i)
            if (serializers[i] == this) {
                serializers.erase(serializers.begin() + i);
                break;
            }
    }
    SerializerImpl(const SerializerImpl&) = delete;
    void load(std::ifstream& fin) override
    {
        T tmp;
        tmp.emplace_back();
        int size;
        fin.read(reinterpret_cast<char*>(&size), sizeof(int));
        fin.read(reinterpret_cast<char*>(ref.data()), size * sizeof(tmp[0]));
    }
    void save(std::ofstream& fout) override
    {
        T tmp;
        tmp.emplace_back();
        int size = (int)ref.size();
        ref.resize(size);
        fout.write(reinterpret_cast<char*>(&size), sizeof(int));
        fout.write(reinterpret_cast<char*>(ref.data()), size * sizeof(tmp[0]));
    }
};

inline void serialization_load_impl(std::string filename)
{
    std::ifstream fin;
    fin.open(filename, std::ios::binary);
    BOW_ASSERT_INFO(!fin.fail(), "Restart data not exists!")
    for (auto s : Serialization::serializers) {
        s->load(fin);
    }
    prepare_phase = false;
}

inline void serialization_save_impl(std::string filename)
{
    std::ofstream fout;
    fout.open(filename, std::ios::binary);
    for (auto s : Serialization::serializers) {
        s->save(fout);
    }
    prepare_phase = false;
}
}
} // namespace Bow::Serialization

#include<boost/typeof/typeof.hpp>
#define SERIALIZATION_REGISTER(name) Bow::Serialization::SerializerImpl<BOOST_TYPEOF(name)> serializer_##name{ name };
#define SERIALIZATION_LOAD Bow::Serialization::serialization_load_impl
#define SERIALIZATION_SAVE Bow::Serialization::serialization_save_impl

#endif