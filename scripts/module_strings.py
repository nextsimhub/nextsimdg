import configparser

default_str = "DEFAULT"
module_file_str = "module.cfg"
module_section_str = "module"
module_file_prefix = "module"
header_suffix = "hpp"
source_suffix = "cpp"
name_str = "name"
file_prefix_str = "file_prefix"
interface_prefix_str = "interface_prefix"
description_str = "description"
has_help_str = "has_help"
is_default_str = "is_default"
file_encoding = "utf-8"
true_str = "true"

# Variable names for the dict of common strings
class_name = "class_name"
module_class_name = "module_class_name"

def get_config_with_defaults():
    config = configparser.ConfigParser()
    config.read(module_file_str)
    
    # Add in defaults, if optional data is absent
    if not interface_prefix_str in config[module_section_str]:
        config[module_section_str][interface_prefix_str] = ""
    
    if not description_str in config[module_section_str]:
        config[module_section_str][description_str] = ""
    
    return config

def check_config_errors(config):
    # Error reporting if the module section is invalid
    if not name_str in config[module_section_str]:
        print(f"The main '{module_section_str}' section lacks a name for the module.")
        return 1
    
    if not file_prefix_str in config[module_section_str]:
        print(f"The main '{module_section_str}' section lacks a file prefix for the source files for the module.")
        return 2
    
    return 0

def common_strings(config):
    return {
        interface_prefix_str : config[module_section_str][interface_prefix_str],
        file_prefix_str : config[module_section_str][file_prefix_str],
        class_name : config[module_section_str][name_str],
        module_class_name : config[module_section_str][file_prefix_str]+"Module"
    }
