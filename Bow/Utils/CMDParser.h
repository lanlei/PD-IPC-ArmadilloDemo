#include <cxxopts/cxxopts.hpp>
#include <Bow/Utils/Logging.h>

namespace Bow {

namespace CMD_PARSER {

int argc;
char** argv;

void parser_initialize_impl(int argc_input, char *argv_input[]) {
    argc = argc_input;
    argv = argv_input;
}

template <class T>
T parser_get_impl(std::string name, T default_value) {
    // allocate memory and copy strings
    T value;
    try
    {
        cxxopts::Options options(argv[0]);
        options
            .allow_unrecognised_options()
            .add_options()
            (name, "", cxxopts::value<T>());
        auto result = options.parse(argc, argv);
        if (result.count(name)) {
            value = result[name].as<T>();
            Logging::info(name, " = ", value, " is read from input");
        }
        else {
            value = default_value;
            Logging::info(name, " = ", value, " set as default value");
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
    return value;
}

}

#define PARSER_INITIALIZE(argc, argv) CMD_PARSER::parser_initialize_impl(argc, argv)
#define PARSER_GET(name, default_value) CMD_PARSER::parser_get_impl(name, default_value)

}
