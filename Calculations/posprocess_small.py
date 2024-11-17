import re
import os
import argparse
#import numpy as np
import fileinput
import signal


# define command line arguments
parser = argparse.ArgumentParser(
    description='Post-Process timing/memory info and plot figures.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--prefix', metavar='prefix', default='CSL_Munich', 
                    help='A folder to look for benchmark results')
parser.add_argument('--custom_model', metavar='custom_model', default='false', 
                    help='Only plot runs with none and acegen caching')
parser.add_argument('--likwid', help='Posprocess LIKWID results', action="store_true")
args = parser.parse_args()

prefix = args.prefix if args.prefix.startswith('/') else os.path.join(os.getcwd(), args.prefix)

files = [os.path.join(prefix, k,'timings.txt') for k in os.listdir(prefix) if os.path.isfile(os.path.join(prefix, k,'timings.txt'))]


custom_model = True if args.custom_model=="true" else False



print(' Custom model: {0}'.format(custom_model))
print(' LIKWID: {0}'.format(args.likwid))

print ('Gather data from {0}'.format(prefix))
print ('found {0} files'.format(len(files)))

pattern = r'[+\-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?'

# sections in timer output:
sections = [
    'vmult (MF)',
    'vmult (Trilinos)',
    'Linear solver',
    'Assemble linear system',
    'Coarse solve level 0',
    'Setup MF: AdditionalData',
    'Setup MF: GMG setup',
    'Setup MF: MGTransferMatrixFree',
    'Setup MF: MappingQEulerian',
    'Setup MF: cache() and diagonal()',
    'Setup MF: ghost range',
    'Setup MF: interpolate_to_mg'
]

start_line = 'Total wallclock time elapsed since start'

mf2d_data = []
mf3d_data = []


for f in files:
    fin = open(f, 'r')
    ready = False

    timing = [100 for i in range(len(sections))]

    # reset CG iterations in case AMG did not have enough memory
    cg_iterations = 9999
    for line in fin:
        if 'dim   =' in line:
            dim = int(re.findall(pattern,line)[0])

        elif 'p     =' in line:
            p = int(re.findall(pattern,line)[0])

        elif 'q     =' in line:
            q = int(re.findall(pattern,line)[0])

        elif 'cells =' in line:
            cells = int(re.findall(pattern,line)[0])

        elif 'dofs  =' in line:
            dofs = int(re.findall(pattern,line)[0])

        elif 'Trilinos memory =' in line:
            tr_memory = float(re.findall(pattern,line)[0])

        elif 'MF cache memory =' in line:
            mf_memory = float(re.findall(pattern,line)[0])

        elif 'Average CG iter =' in line:
            cg_iterations = int(re.findall(pattern,line)[0])

        if '_scalar_ref' in f:
            cachings = 'scalar_ref'
        elif '_scalar' in f:
            cachings = 'scalar_def'
        elif '_tensor2' in f:
            cachings = 'tensor2'
        elif '_none' in f:
            cachings = 'none'
        elif '_acegen_cached' in f:
            cachings = 'acegen_ref'
        elif '_acegen_actual' in f:
            cachings = 'acegen_def'
        # first check tensor4_ns so that that tensor4 is not triggered
        # for this case
        elif '_tensor4_ns' in f:
            cachings = 'tensor4_ns'
        elif '_tensor4' in f:
            cachings = 'tensor4S'

        if '_amg_' in f:
            cachings = 'amg'

        if start_line in line:
            print ('dim={0} p={1} q={2} cells={3} dofs={4} tr_memory={5} mf_memory={6} cg_it={7} caching={8} file={9}'.format(dim, p, q, cells, dofs, tr_memory, mf_memory, cg_iterations, cachings, f))
            ready = True
        if not args.likwid:
            print("WITH LIKWID")
            continue

        if ready:
            # we could have sections starting from the same part
            line_split = line.split("|")
            line_sec = line_split[1].strip() if len(line_split) > 1 else ''
            for idx, s in enumerate(sections):
                if s == line_sec:
                    nums = re.findall(pattern, "".join(line.rsplit(s)))
                    # do time of a single vmult and the reset -- total time
                    n = int(nums[0]) if 'vmult (' in s else int(1)  # how many times
                    t = float(nums[1]) # total time
                    print ('  {0} {1} {2} {3}'.format(s,idx,n,t))
                    timing[idx] = t / n

    # finish processing the file, put the data
    tp = tuple((p,cachings, dofs, mf_memory, timing[0]))
    if dim == 2:
        mf2d_data.append(tp)
    elif dim == 3:
        mf3d_data.append(tp)

print("2D Data:")
headers = ["p", "cachings", "dofs", "mf_memory"] 
print("\t".join(headers))
for data in mf2d_data:
    row = list(data )
    print("\t".join(map(str, row)))


print("3D Data:")
headers = ["p", "cachings", "dofs", "mf_memory"] 
print("\t".join(headers))
for data in mf3d_data:
    row = list(data )
    print("\t".join(map(str, row)))
