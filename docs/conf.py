#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    subprocess.call('doxygen', shell=True)

import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

def setup(app):
    app.add_stylesheet("main_stylesheet.css")

extensions = ['breathe']
breathe_projects = { 'nextsimdg': 'xml' }
templates_path = ['_templates']
html_extra_path = ['xml']
source_suffix = '.rst'
master_doc = 'index'
project = 'nextsimdg'
copyright = '2021, Nansen Environmental and Remote Sensing Center'
author = 'Nansen Environmental and Remote Sensing Center'

exclude_patterns = []
highlight_language = 'c++'
pygments_style = 'sphinx'
todo_include_todos = False
htmlhelp_basename = 'nextsimdgdoc'
