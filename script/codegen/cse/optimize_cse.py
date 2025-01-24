from source.modules_cse import executeCSE
from argparse import ArgumentParser

# Parse name of dynamics to be optimized
parser = ArgumentParser()
parser.add_argument('inputfile', help='filename containing expressions')
parser.add_argument('cse_template', help='filename of mako template to execute cse')
parser.add_argument('outputfile',   help='path and filename of output file')
args = parser.parse_args()

# Generate cse-optimized C++ code from tree and fill template
print('Generate from ' + args.inputfile)
executeCSE(args.inputfile, args.cse_template, args.outputfile)
