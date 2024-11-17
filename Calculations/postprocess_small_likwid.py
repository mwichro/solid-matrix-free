import re
import os
import argparse
#import numpy as np
import fileinput
import signal
import glob


# define command line arguments
parser = argparse.ArgumentParser(
    description='Post-Process timing/memory info and plot figures.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--prefix', metavar='prefix', default='CSL_Munich', 
                    help='A folder to look for benchmark results')
parser.add_argument('--custom_model', metavar='custom_model', default='false', 
                    help='Only plot runs with none and acegen caching')
args = parser.parse_args()

prefix = args.prefix if args.prefix.startswith('/') else os.path.join(os.getcwd(), args.prefix)

folders =  os.listdir(prefix)
print("Folders found:")
for folder in folders :
    print(folder)



custom_model = True if args.custom_model=="true" else False



print(' Custom model: {0}'.format(custom_model))

print ('Gather data from {0}'.format(prefix))

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


for f in folders:

    files = glob.glob(os.path.join(prefix,f, '*.toutput'))
    if len(files) != 1:
        print(f"Expected exactly one .touput file in folder {f}, but found {len(files)}")
        continue
    else:
        print(f"Found exactly one .touput file in folder {f}")
    file = files[0]
    print("opening file: {0}".format(file))


    fin = open(file, 'r')
    ready = False

    if '_2d_' in f:
        dim=2
    else:
        dim=3 

    match = re.search(r'_p[0-9]q(\d+)', f)
    if match:
        q = int(match.group(1))
    else:
        q = None

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
        continue

    print ('dim={0} caching={1} q={2} file={3}'.format(dim, cachings,q,  f))
    
    quadrature_loop_data = None
    for line in fin:
        if 'Region: vmult_quadrature_loop' in line:
            quadrature_loop_data = {}
            while True:
                line = fin.readline()
                if line.startswith('Region'):
                    break
                if 'RETIRED_SSE_AVX_FLOPS_ALL' in line:
                    quadrature_loop_data['RETIRED_SSE_AVX_FLOPS_ALL'] = float(re.findall(pattern,line)[1])
                if 'Runtime (RDTSC) [s]' in line:
                    quadrature_loop_data['Runtime (RDTSC) [s]'] = float(re.findall(pattern,line)[0])
                if 'Region calls' in line:
                    quadrature_loop_data['Region calls'] = int(re.findall(pattern,line)[0])
                if 'DP [MFLOP/s]' in line:
                    quadrature_loop_data['DP MFLOPs'] = float(re.findall(pattern,line)[0])
            break

    if quadrature_loop_data:
        flops = quadrature_loop_data.get('RETIRED_SSE_AVX_FLOPS_ALL', 0)
        flops_per_qpoint = flops / ( quadrature_loop_data.get('Region calls', 1) * q**dim) 
        runtime = quadrature_loop_data.get('Runtime (RDTSC) [s]', 1)
        calls = quadrature_loop_data.get('Region calls', 1)
        flops_per_s = quadrature_loop_data.get('DP MFLOPs', 0)/1e3
        # flops_per_s = flops / runtime if runtime > 0 else 0
        print(f"Quadrature loop: GFLOPS/s = {flops_per_s}, Total FLOPS = {flops}, Number of calls = {calls} , FLOPS per qpoint = {flops_per_qpoint}")

        data_entry = (q, cachings, flops_per_qpoint, flops_per_s, runtime)
        if dim == 2:
            mf2d_data.append(data_entry)
        else:
            mf3d_data.append(data_entry)
    else:
        print("Quadrature loop data not found")



print("2D Data:")
headers = ["q", "cachings", "flops_per_qpoint", "[GFLOPS/s]", "runtime"] 
print("\t".join(headers))
for data in mf2d_data:
    row = list(data)
    print("{:<10} {:<15} {:<20.3f} {:<15.3f} {:<10.7f}".format(row[0], row[1], row[2], row[3], row[4]))


print("3D Data:")
headers = ["q", "cachings", "flops_per_qpoint", "flops_per_s", "runtime"] 
print("\t".join(headers))
for data in mf3d_data:
    row = list(data)
    print("{:<10} {:<15} {:<20.3f} {:<15.3f} {:<10.3f}".format(row[0], row[1], row[2], row[3], row[4]))
