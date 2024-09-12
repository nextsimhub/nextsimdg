#
# module_builder.py
#
# 15 Nov 2023
# author: Tim Spain <timothy.spain@nersc.no>

import os

from module_strings import *

# A valid implementation section has an entry "file_prefix = "
def is_impl_section_valid(section):
    return (file_prefix_str in section)

# A file and generator header for each file
def write_file_header(stream, module_name):
    stream.write(
        f"""/*!
 * @file {module_file_prefix}.{source_suffix}
 * for {module_name}
 *
 * Generated by {os.path.basename(__file__)}
 */
"""
        )

def write_source_file(source, config, strings):
    # List over all implementation sections that contain a file_prefix entry
    valid_impl_sections = []
    for section in config:
        if section in (module_section_str, default_str,):
            continue
        else:
            if file_prefix_str in config[section]:
                valid_impl_sections.append(section)
                if is_default_str in config[section] and config[section][is_default_str]:
                    default_impl = section

    module_templ = f"Module<{strings[class_name]}>"

    # Use the provided path to the Module header file
    source.write(f"#include \"{strings[header_file_path_str]}\"\n")
    source.write("\n")
#    source.write(f"#include \"{strings[interface_prefix_str]}{strings[file_prefix_str]}.{header_suffix}\"\n")
    for section in valid_impl_sections:
        source.write(f"#include \"{os.path.join(strings[internal_header_dir], config[section][file_prefix_str])}.{header_suffix}\"\n")
    source.write("""
#include <string>

namespace Module {
""")
    impl_strings = {}
    for section in valid_impl_sections:
        impl_strings[section] = config[section][file_prefix_str].upper()
        source.write(f"const std::string {impl_strings[section]} = \"{section}\";\n")
    mapVarName = "theMap"
    source.write(f"""
template <>
const {module_templ}::map& {module_templ}::functionMap()
{{
    static const map {mapVarName} = {{
""")
    for section in valid_impl_sections:
        source.write(f"        {{ {impl_strings[section]}, newImpl<{strings[class_name]}, {section}> }},\n")
    source.write(f"""    }};
    return {mapVarName};
}}

template <>
{module_templ}::fn& {module_templ}::getGenerationFunction()
{{
    static fn ptr = functionMap().at({impl_strings[default_impl]});
    return ptr;
}}

template <> std::string {module_templ}::moduleName() {{ return \"{strings[module_class_name]}\"; }}

template <> HelpMap& getHelpRecursive<{strings[class_name]}>(HelpMap& map, bool getAll)
{{
    const std::string& pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({{ pfx + "." + {module_templ}::moduleName(), ConfigType::MODULE,
        {{ """)
    for section in valid_impl_sections:
        source.write(f"{impl_strings[section]}, ")
    source.write(f"}}, {impl_strings[default_impl]}, \"\",\n")
    source.write(f"        \"{config[module_section_str][description_str]}\" }});\n")
    for section in valid_impl_sections:
        if (has_help_str in config[section]) and (config[section][has_help_str] == true_str):
            source.write(f"    {section}::getHelpRecursive(map, getAll);\n")
    source.write(f"""
    return map;
}}
template <> {strings[class_name]}& getImplementation<{strings[class_name]}>()
{{
    return getImplTemplate<{strings[class_name]}, {strings[module_class_name]}>();
}}
template <> void setImplementation<{strings[class_name]}>(const std::string& implName)
{{
    setImplTemplate<{strings[module_class_name]}>(implName);
}}
template <> std::unique_ptr<{strings[class_name]}> getInstance()
{{
    return getInstTemplate<{strings[class_name]}, {strings[module_class_name]}>();
}}

template class {module_templ};
""")
    source.write("""
} /* namespace Module */
""")
    
'''
Main program

Parse the module.cfg file in a directory and generate the corresponding
module source file.
'''
def main():
    config = get_config_with_defaults()
    
    config_status = check_config_errors(config)
    if config_status > 0:
        return config_status
    
    # Create a dictionary of common strings
    strings = common_strings(config)
    
    source = open(module_file_prefix + "." + source_suffix, "w", encoding=file_encoding)
    
    write_file_header(source, strings[module_class_name])
    write_source_file(source, config, strings)
    source.close()
    
if __name__ == "__main__":
    main()
