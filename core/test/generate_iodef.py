"""Generate an XML file for XIOS from a Jinja2 template based on user input."""

import argparse
from jinja2 import Environment, FileSystemLoader
import os

# Parse user input
parser = argparse.ArgumentParser()
parser.add_argument("DGCOMP", type=int, default=6, help="Number of DG components")
parser.add_argument(
    "DGSTRESSCOMP", type=int, default=8, help="Number of DG stress components"
)
parser.add_argument("infile", type=str, help="Input filename")
parser.add_argument("outfile", type=str, help="Output filename")
parsed_args = parser.parse_args()
template_dict = vars(parsed_args)

# Set up Jinja2 environment and render template
path, infile = os.path.split(template_dict.pop("infile"))
env = Environment(loader=FileSystemLoader(path))
template = env.get_template(infile)

# Generate the XML file
outfile = template_dict.pop("outfile")
output = template.render(template_dict)
with open(outfile, "w+") as f:
    f.write(output)
