def denamespace(nname):
    """Returns the last element of the name, without any of the qualifying
    namespaces."""
    return nname.split(":")[-1]

def generator(module_class_name, fq_interface_name, fq_impl_names, help_names):
    """Generates the text for a module support file."""
    interface_name = denamespace(fq_interface_name)

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
    for fq_impl_name in fq_impl_names:
        print(f"const std::string {denamespace(fq_impl_name).upper()} = \"{fq_impl_name}\";")
    print("")
    
    # Create the functionMap from the FQ implementation and FQ module names
    print("template <>")
    print(f"Module<{fq_interface_name}>::map Module<{fq_interface_name}>::functionMap" + " = {")
    for fq_impl_name in fq_impl_names:
        print(f"    {{ {denamespace(fq_impl_name).upper()}, newImpl<{fq_interface_name}, {fq_impl_name}> }},")
    print("};")
    print("")
    
    # Set up the function and static pointer (FQ Module)
    print("template <>")
    print(f"Module<{fq_interface_name}>::fn Module<{fq_interface_name}>::spf = functionMap.at({denamespace(fq_impl_names[0]).upper()});")
    print("template <>")
    print(f"std::unique_ptr<{fq_interface_name}> Module<{fq_interface_name}>::staticInstance")
    print(f"= std::move(newImpl<{fq_interface_name}, {fq_impl_names[0]}>());")
    print("")

    # Module name string
    print("template <>")
    print(f"std::string Module<{fq_interface_name}>::moduleName(){{ return \"{fq_interface_name}\"; }}")
    print("")
    
    # global functions (FQ module & module class names)

    # Recursive help generation
    # List of uppercased implementation names
    impl_names_uc = []
    for fq_impl_name in fq_impl_names:
        impl_names_uc.append(denamespace(fq_impl_name).upper())

    print(f"template<> HelpMap& getHelpRecursive<{fq_interface_name}>(HelpMap& map, bool getAll)")
    print("{")
    print("    const std::string& pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;")
    print(f"    map[pfx].push_back({{ pfx + \".\" + Module<{fq_interface_name}>::moduleName(), ConfigType::MODULE,")
    impl_namelist_uc = ", ".join(impl_names_uc)
    print(f"        {{ {impl_namelist_uc} }}, {impl_names_uc[0]}, \"\",")
    print("        \"MODULE DESCRIPTION HERE\" });")
    for help_name in help_names:
        print(f"    {help_name}::getHelpRecursive(map, getAll);")
    print("    return map;")
    print("}")
    
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

    print(f"{module_class_name}::Constructor {module_class_name}::ctor;")
    print(f"{module_class_name}::Constructor::Constructor()")
    print("{")
    print(f"    addToConfiguredModules<{fq_interface_name}, {module_class_name}>();")
    print("}")
    print("")
    print("} /* namespace Module */")

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description = "Write out the text for a Nextsim module class source file.",
                                     epilog = "Suffix the interface or any implementation name with an asterisk (*) to include a call to getHelpRecursive().")
    parser.add_argument("impl", metavar = "impls", nargs = '*', default = None, help = "Fully qualified name of the implementation classes.")
    parser.add_argument("--interface", dest = "interface", required = True, help = "Fully qualified name of the interface class.")
    parser.add_argument("--module-prefix", dest = "modulepfx", help = "Name of the module, will be suffixed by 'Module'.")
    args = parser.parse_args()
    
    if (len(args.impl) == 0):
        raise SystemExit

    # filter all the interface and implementation names for trailing asterisks,
    # signifying that getHelpRecursive should be called on these classes.
    help_names = []
    impl_names = []
    interface_name = args.interface
    if interface_name[-1] == "*":
        interface_name = interface_name[:-1]
        help_names.append(interface_name)
    for impl_name in args.impl:
        if impl_name[-1] == "*":
            impl_names.append(impl_name[:-1])
            help_names.append(impl_names[-1])
        else:
            impl_names.append(impl_name)

    # Generate a module prefix if none is set
    if args.modulepfx is None:
        iface = denamespace(interface_name)
        if (iface[0] == "I") and iface[1].isupper():
            modulepfx = iface[1:]
        else:
            modulepfx = denamespace(interface_name)
        modulepfx += "Module"
    else:
        modulepfx = denamespace(args.modulepfx)


    generator(modulepfx, interface_name, impl_names, help_names)
