def denamespace(nname):
    """Returns the last element of the name, without any of the qualifying
    namespaces."""
    return nname.split(":")[-1]

def generator(fq_interface_name, fq_impl_names):
    """Generates the text for a module support file."""
    interface_name = denamespace(fq_interface_name)

    module_class_name = interface_name + "Module"

    # Add the include file for the module class
    print(f"#include \"include/{module_class_name}.hpp\"")
    print("")
    # List the include files from the implementation names
    for fq_impl_name in fq_impl_names:
        print(f"#include \"include/{denamespace(fq_impl_name)}.hpp\"")

    # Standard includes
    print("")
    print("#include <string>")
        
    print("")
    print("namespace Module {")

    # String constants to name the implementations
    for impl_name in map(denamespace, fq_impl_names):
        print(f"const std::string {impl_name.upper()} = \"{impl_name}\";")
    print("")
    
    # Create the functionMap from the FQ implementation and FQ module names
    print("template <>")
    print(f"Module<{fq_interface_name}>::map Module<{fq_interface_name}>::functionMap" + " = {")
    for fq_impl_name in fq_impl_names:
        print("    {" + f"{denamespace(fq_impl_name).upper()}, newImpl<{fq_interface_name}, {fq_impl_name}>" + "},")
    print("};")
    print("")
    
    # Set up the function and static pointer (FQ Module)
    print("template <>")
    print(f"Module<{fq_interface_name}>::fn Module<{fq_interface_name}>::spf = functionMap.at({denamespace(fq_impl_names[0]).upper()});")
    print("template <>")
    print(f"std::unique_ptr<{fq_interface_name}> Module<{fq_interface_name}>::staticInstance")
    print(f"= std::move(Module<{fq_interface_name}>::spf());")
    print("")

    # Module name string
    print("template <>")
    print(f"std::string Module<{fq_interface_name}>::moduleName()" + "{" + f"    return \"{denamespace(fq_interface_name)}\";" + "}")
    print("")

    # global functions (FQ module & module class names)
    print("template <>")
    print(f"{fq_interface_name}& getImplementation<{fq_interface_name}>()")
    print("{")
    print(f"    return getImplTemplate<{fq_interface_name}, {module_class_name}>();")
    print("}")

    print("template <>")
    print(f"void setImplementation<{fq_interface_name}>(const std::string& implName)")
    print("{")
    print(f"    setImplTemplate<{module_class_name}>(implName);")
    print("}")

    print("template <>")
    print(f"std::unique_ptr<{fq_interface_name}> getInstance()")
    print("{")
    print(f"    return getInstTemplate<{fq_interface_name}, {module_class_name}>();")
    print("}")

    print("} /* namespace Module */")

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description = "Write out the text for a Nextsim module class source file.")
    parser.add_argument("impl", metavar = "impls", nargs = '*', default = None, help = "Fully qualified name of the implementation classes.")
    parser.add_argument("--interface", dest = "interface", help = "Fully qualified name of the interface class.")

    args = parser.parse_args()
    
    if (len(args.impl) == 0):
        raise SystemExit

    generator(args.interface, args.impl)
