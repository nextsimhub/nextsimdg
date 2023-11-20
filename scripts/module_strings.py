import configparser
import os.path

default_str = "DEFAULT"
module_file_str = "module.cfg"
module_section_str = "module"
module_file_prefix = "module"
header_suffix = "hpp"
source_suffix = "cpp"
name_str = "name"
file_prefix_str = "file_prefix"
interface_prefix_str = "interface_prefix"
header_dir_path_str = "header_path"
header_file_path_str = "header_file_path"
header_file_name_str = "header_file_name"
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
    
    # Default to ../include for the header file path
    if not header_dir_path_str in config[module_section_str]:
        config[module_section_str][header_dir_path_str] = "../include"
    
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
    module_name = config[module_section_str][file_prefix_str]+"Module"
    header_file_name = module_name + "." + header_suffix
    return {
        interface_prefix_str : config[module_section_str][interface_prefix_str],
        file_prefix_str : config[module_section_str][file_prefix_str],
        class_name : config[module_section_str][name_str],
        module_class_name : module_name,
        header_file_name_str : header_file_name,
        header_file_path_str : os.path.join(config[module_section_str][header_dir_path_str], module_name + "." + header_suffix)
    }
