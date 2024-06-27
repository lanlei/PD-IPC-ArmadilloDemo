#pragma once

#include <set>
#include <Bow/Types.h>
#include <Bow/Utils/FileSystem.h>

namespace Bow {

namespace RESULT_RECORDER {

static std::set<std::string> collection;

void record(const std::string& name, double value_a, double value_b) {
    FileSystem::create_directory("result");
    if (collection.find(name) == collection.end()) {
        collection.insert(name);
        FILE* f = fopen(("result/" + name + ".txt").c_str(), "w");
        fclose(f);
    }
    FILE* f = fopen(("result/" + name + ".txt").c_str(), "a");
    fprintf(f, "%.20f %.20f\n", value_a, value_b);
    fclose(f);
}

}

}