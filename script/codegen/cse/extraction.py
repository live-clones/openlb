from argparse import ArgumentParser
from source.modules_cse import extractExpressions

# Parse name of dynamics to be optimized
parser = ArgumentParser()
parser.add_argument('dynamics', help='full C++ typename of the dynamics instance')
parser.add_argument('template', help='path and filename to mako template')
parser.add_argument('input', help='path and filename to .txt file containing dynamics info')
parser.add_argument('cpp_output', help='path where the .cpp files are created')
parser.add_argument('out_output', help='path where the .out files are created')
args = parser.parse_args()

# Create .cpp files to extract expression tree
extract_filename = extractExpressions(args.dynamics, args.template, args.input, args.cpp_output, args.out_output)
