"""This module generates the inclusion .ipp files for the C++ ModuleLoader."""

def get_iname(name):
    """Returns the name of the interface class, given its string name. Input
    can be namespaced or not, output is similar."""
    components = name.split(":")
    components[-1] = f"I{components[-1]}"
    return ":".join(components)

def denamespace(nname):
    """Returns the last element of the name, without any of the qualifying
    namespaces."""
    return nname.split(":")[-1]

def get_pname(full_name):
    """Returns the name of the function pointer for the interface, given its
    namespaced class name."""
    return f"p_{get_iname(denamespace(full_name))}"

def get_fname(impl):
    """Returns the function name for the implementation, given its namespaced string name."""
    return f"new{denamespace(impl)}"

def headers(all_implementations, ipp_prefix, hpp_prefix):
    """Generates the moduleLoaderHeaders.ipp file."""
    with open(f"{ipp_prefix}moduleLoaderHeaders.ipp", "w", encoding="utf-8") as fil:
        for interface in all_implementations:
            i_name = get_iname(denamespace(interface["name"]))
            fil.write(f"#include \"{hpp_prefix}{i_name}.hpp\"\n")
            for impl in interface["implementations"]:
                i_name = denamespace(impl)
                fil.write(f"#include \"{hpp_prefix}{i_name}.hpp\"\n")
            # An extra line between interfaces
            fil.write("\n")

def functions(all_implementations, ipp_prefix):
    """Generates the moduleLoaderFunctions.ipp file."""
    with open(f"{ipp_prefix}moduleLoaderFunctions.ipp", "w", encoding="utf-8") as fil:
        for interface in all_implementations:
            # Define the pointer to function
            full_name = interface["name"]
            i_name = get_iname(full_name)
            p_name = get_pname(full_name)
            fil.write(f"std::unique_ptr<{i_name}> (*{p_name})();\n")
            # Define function that call the function pointer
            fil.write(
                "template<>\n"
                f"std::unique_ptr<{i_name}> ModuleLoader::getImplementation() const\n"
                "{\n"
                f"    return (*{p_name})();\n"
                "}\n"
                )
            for impl in interface["implementations"]:
                # The function that return the new unique_ptr for each implementation
                fil.write(
                    f"std::unique_ptr<{i_name}> {get_fname(impl)}()\n"
                    "{\n"
                    f"    return std::unique_ptr<{impl}>(new {impl});\n"
                    "}\n"
                    )
            # An extra line between interfaces
            fil.write("\n")

def names(all_implementations, ipp_prefix):
    """Generates the moduleLoaderNames.ipp file."""
    with open(f"{ipp_prefix}moduleLoaderNames.ipp", "w", encoding="utf-8") as fil:
        fil.write(
            "    m_availableImplementationNames = {\n"
            "        "
            )
        for interface in all_implementations:
            name = interface["name"]
            fil.write(
                "{\n"
                f"        \"{name}\", ""{\n"
                )
            for impl in interface["implementations"]:
                fil.write(f"            \"{impl}\",\n")
            fil.write(
                "            }\n"
                "        },"
                )
        fil.write("\n    };\n")

def assignments(all_implementations, ipp_prefix):
    """Generates the moduleLoaderAssignments.ipp file."""
    with open(f"{ipp_prefix}moduleLoaderAssignments.ipp", "w", encoding="utf-8") as fil:
        fil.write("        ")
        for interface in all_implementations:
            name = interface["name"]
            fil.write(
                f"if (module == \"{name}\") ""{\n"
                "            "
                )
            p_name = get_pname(name)
            for impl in interface["implementations"]:
                fil.write(
                    f"if (impl == \"{impl}\") ""{\n"
                    f"                {p_name} = &{get_fname(impl)};\n"
                    "            } else "
                    )
            fil.write(
                "{\n"
                "                throwup(module, impl);\n"
                "            }\n"
                "\n"
                "        } else "
                )
        fil.write("{ }")

def generate(all_implementations, ipp_prefix = '', hpp_prefix = ''):
    """Generates the .ipp inclusion files for ModuleLoader.cpp

    :param all_implementations: The vector of dictionaries that defines the
            interfaces and implementations thereof.
    :param ipp_prefix: A directory and file prefix that precedes the ipp file
            names to provide a path from the current working directory.
    :param hpp_prefix: A text directory and file prefix to add to the hpp file
            names to suit the locations in the build system.
    """
    headers(all_implementations, ipp_prefix, hpp_prefix)
    functions(all_implementations, ipp_prefix)
    names(all_implementations, ipp_prefix)
    assignments(all_implementations, ipp_prefix)

if __name__ == "__main__":

    import argparse
    import json

    parser =  argparse.ArgumentParser(description = "Build inclusion files for ModuleLoader.cpp")
    FILE_HELP = "JSON file containing the specification of interfaces and implementations"
    parser.add_argument("json_file", metavar = "file", nargs = '?', default = None,
                         type = argparse.FileType('r', encoding="utf-8"),
                         help = FILE_HELP)
    parser.add_argument("--ipp", dest = "ipp_prefix", default = "./",
                        help = "Path and file prefix to be added to the .ipp file names.")
    parser.add_argument("--hpp", dest = "hpp_prefix", default = "include/",
                        help = "Path to the module header file name.")
    args = parser.parse_args()

    DFILE = "modules.json"
    use_default = args.json_file is None

    with args.json_file if not use_default else open(DFILE, encoding = "utf-8") as json_file:
        alli = json.load(json_file)

        generate(alli, hpp_prefix = args.hpp_prefix, ipp_prefix = args.ipp_prefix)
