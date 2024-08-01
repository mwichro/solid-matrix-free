import re
import os
import argparse
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np
import fileinput
import matplotlib.ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
import signal

# https://stackoverflow.com/questions/42656139/set-scientific-notation-with-fixed-exponent-and-significant-digits-for-multiple
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

# define command line arguments
parser = argparse.ArgumentParser(
    description='Post-Process timing/memory info and plot figures.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('prefix', metavar='prefix', default='CSL_Munich', nargs='?',
                    help='A folder to look for benchmark results')
args = parser.parse_args()

prefix = args.prefix if args.prefix.startswith('/') else os.path.join(os.getcwd(), args.prefix)

files = [os.path.join(prefix, k,'timings.txt') for k in os.listdir(prefix) if os.path.isfile(os.path.join(prefix, k,'timings.txt'))]

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

mf2d_data_scalar = []
mf2d_data_tensor2 = []
mf2d_data_tensor4 = []
mf2d_data_tensor4_ns = []
mb2d_data = []

mf3d_data_scalar = []
mf3d_data_tensor2 = []
mf3d_data_tensor4 = []
mf3d_data_tensor4_ns = []
mb3d_data = []


for f in files:
    fin = open(f, 'r')
    ready = False

    # timing = ['-' for i in range(2*len(sections))]
    timing = [np.nan for i in range(len(sections))]

    # reset CG iterations in case AMG did not have enough memory
    cg_iterations = np.nan
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

        if start_line in line:
            print ('dim={0} p={1} q={2} cells={3} dofs={4} tr_memory={5} mf_memory={6} cg_it={7} file={8}'.format(dim, p, q, cells, dofs, tr_memory, mf_memory, cg_iterations, f))
            ready = True

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
    tp = tuple((p, dofs, tr_memory, mf_memory, timing, cg_iterations))
    if 'MF_CG' in f:
        if '_scalar' in f:
            if dim == 2:
                mf2d_data_scalar.append(tp)
            else:
                mf3d_data_scalar.append(tp)
        elif '_tensor2' in f:
            if dim == 2:
                mf2d_data_tensor2.append(tp)
            else:
                mf3d_data_tensor2.append(tp)
        # first check tensor4_ns so that that tensor4 is not triggered
        # for this case
        elif '_tensor4_ns' in f:
            if dim == 2:
                mf2d_data_tensor4_ns.append(tp)
            else:
                mf3d_data_tensor4_ns.append(tp)
        elif '_tensor4' in f:
            if dim == 2:
                mf2d_data_tensor4.append(tp)
            else:
                mf3d_data_tensor4.append(tp)
    else:
        if dim == 2:
            mb2d_data.append(tp)
        else:
            mb3d_data.append(tp)

# now we have lists of tuples ready
# first, sort by degree:
mf2d_data_scalar.sort(key=lambda tup: tup[0])
mf2d_data_tensor2.sort(key=lambda tup: tup[0])
mf2d_data_tensor4.sort(key=lambda tup: tup[0])
mf2d_data_tensor4_ns.sort(key=lambda tup: tup[0])
mb2d_data.sort(key=lambda tup: tup[0])

mf3d_data_scalar.sort(key=lambda tup: tup[0])
mf3d_data_tensor2.sort(key=lambda tup: tup[0])
mf3d_data_tensor4.sort(key=lambda tup: tup[0])
mf3d_data_tensor4_ns.sort(key=lambda tup: tup[0])
mb3d_data.sort(key=lambda tup: tup[0])

# now get the data for printing
deg2d = [tup[0] for tup in mf2d_data_scalar]
deg3d = [tup[0] for tup in mf3d_data_scalar]

# vmult time per dof
time2d_tr    = [tup[4][1]/tup[1] for tup in mb2d_data]
time2d_sc    = [tup[4][0]/tup[1] for tup in mf2d_data_scalar]
time2d_t2    = [tup[4][0]/tup[1] for tup in mf2d_data_tensor2]
time2d_t4    = [tup[4][0]/tup[1] for tup in mf2d_data_tensor4]
time2d_t4_ns = [tup[4][0]/tup[1] for tup in mf2d_data_tensor4_ns]

time3d_tr    = [tup[4][1]/tup[1] for tup in mb3d_data]
time3d_sc    = [tup[4][0]/tup[1] for tup in mf3d_data_scalar]
time3d_t2    = [tup[4][0]/tup[1] for tup in mf3d_data_tensor2]
time3d_t4    = [tup[4][0]/tup[1] for tup in mf3d_data_tensor4]
time3d_t4_ns = [tup[4][0]/tup[1] for tup in mf3d_data_tensor4_ns]

through2d_tr    = [tup[1]/tup[4][1] for tup in mb2d_data]
through2d_sc    = [tup[1]/tup[4][0] for tup in mf2d_data_scalar]
through2d_t2    = [tup[1]/tup[4][0] for tup in mf2d_data_tensor2]
through2d_t4    = [tup[1]/tup[4][0] for tup in mf2d_data_tensor4]
through2d_t4_ns = [tup[1]/tup[4][0] for tup in mf2d_data_tensor4_ns]

through3d_tr    = [tup[1]/tup[4][1] for tup in mb3d_data]
through3d_sc    = [tup[1]/tup[4][0] for tup in mf3d_data_scalar]
through3d_t2    = [tup[1]/tup[4][0] for tup in mf3d_data_tensor2]
through3d_t4    = [tup[1]/tup[4][0] for tup in mf3d_data_tensor4]
through3d_t4_ns = [tup[1]/tup[4][0] for tup in mf3d_data_tensor4_ns]

# solver time per dof
solver2d_tr        = [tup[4][2]/tup[1] for tup in mb2d_data]
solver2d_sc        = [tup[4][2]/tup[1] for tup in mf2d_data_scalar]
solver2d_t2        = [tup[4][2]/tup[1] for tup in mf2d_data_tensor2]
solver2d_t4        = [tup[4][2]/tup[1] for tup in mf2d_data_tensor4]
solver2d_t4_ns     = [tup[4][2]/tup[1] for tup in mf2d_data_tensor4_ns]
solver2d_t4_coarse = [tup[4][4]/tup[1] for tup in mf2d_data_tensor4]

solver3d_tr        = [tup[4][2]/tup[1] for tup in mb3d_data]
solver3d_sc        = [tup[4][2]/tup[1] for tup in mf3d_data_scalar]
solver3d_t2        = [tup[4][2]/tup[1] for tup in mf3d_data_tensor2]
solver3d_t4        = [tup[4][2]/tup[1] for tup in mf3d_data_tensor4]
solver3d_t4_ns     = [tup[4][2]/tup[1] for tup in mf3d_data_tensor4_ns]
solver3d_t4_coarse = [tup[4][4]/tup[1] for tup in mf3d_data_tensor4]

assembly2d_tr = [tup[4][3]/tup[1] for tup in mb2d_data]
assembly3d_tr = [tup[4][3]/tup[1] for tup in mb3d_data]

mf_gmg_2d_t4    = [(tup[4][5]  + \
                    tup[4][6]  + \
                    tup[4][7]  + \
                    tup[4][8]  + \
                    tup[4][9]  + \
                    tup[4][10] + \
                    tup[4][11])/tup[1] for tup in mf2d_data_tensor4]

mf_gmg_3d_t4    = [(tup[4][5]  + \
                    tup[4][6]  + \
                    tup[4][7]  + \
                    tup[4][8]  + \
                    tup[4][9]  + \
                    tup[4][10] + \
                    tup[4][11])/tup[1] for tup in mf3d_data_tensor4]

mf_gmg_2d_t4_ns = [(tup[4][5]  + \
                    tup[4][6]  + \
                    tup[4][7]  + \
                    # tup[4][8]  + \  # without MappingQEulerian
                    tup[4][9]  + \
                    tup[4][10] + \
                    tup[4][11])/tup[1] for tup in mf2d_data_tensor4_ns]

mf_gmg_3d_t4_ns = [(tup[4][5]  + \
                    tup[4][6]  + \
                    tup[4][7]  + \
                    # tup[4][8]  + \ # without MappingQEulerian
                    tup[4][9]  + \
                    tup[4][10] + \
                    tup[4][11])/tup[1] for tup in mf3d_data_tensor4_ns]

# CG iterations
cg2d_tr = [tup[5] for tup in mb2d_data]
cg2d_t4 = [tup[5] for tup in mf2d_data_tensor4]

cg3d_tr = [tup[5] for tup in mb3d_data]
cg3d_t4 = [tup[5] for tup in mf3d_data_tensor4]

# Mb per dof
mem2d_tr    = [tup[2]/tup[1] for tup in mb2d_data]
mem2d_sc    = [tup[3]/tup[1] for tup in mf2d_data_scalar]
mem2d_t2    = [tup[3]/tup[1] for tup in mf2d_data_tensor2]
mem2d_t4    = [tup[3]/tup[1] for tup in mf2d_data_tensor4]
mem2d_t4_ns = [tup[3]/tup[1] for tup in mf2d_data_tensor4_ns]

mem3d_tr    = [tup[2]/tup[1] for tup in mb3d_data]
mem3d_sc    = [tup[3]/tup[1] for tup in mf3d_data_scalar]
mem3d_t2    = [tup[3]/tup[1] for tup in mf3d_data_tensor2]
mem3d_t4    = [tup[3]/tup[1] for tup in mf3d_data_tensor4]
mem3d_t4_ns = [tup[3]/tup[1] for tup in mf3d_data_tensor4_ns]

# file location
fig_prefix = os.path.join(os.getcwd(), '../doc/' + os.path.basename(os.path.normpath(prefix)) + '_')

#
#
#  MATPLOTLIB
#
#

# 6.1:
params   = {'legend.fontsize': 16,
            'font.size': 18}
# 6.2:
params2 = {'legend.fontsize': 14,
           'font.size': 16}

plt.rcParams.update(params)

#
# Plot Section 6.1
#
ratio = 0.7
timing_exp = -8 if "CSL" in args.prefix else -7

ax = plt.figure().gca()
#ax.yaxis.set_major_formatter(OOMFormatter(timing_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot(deg2d,time2d_tr, 'rs--', label='Trilinos')
plt.plot(deg2d,time2d_sc, 'bo--', label='MF scalar')
plt.plot(deg2d,time2d_t2, 'g^--', label='MF tensor2')
plt.plot(deg2d,time2d_t4, 'cv--', label='MF tensor4')
#plt.plot(deg2d,time2d_t4_ns, 'mD--', label='MF tensor4 P')
plt.xlabel('polynomial degree')
plt.ylabel('vmult wall time (s) / DoF')
leg = plt.legend(loc='best', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'timing2d.pdf', format='pdf', bbox_inches = 'tight')
# plt.show()


# clear
plt.clf()
exit()

ax = plt.figure().gca()
ax.yaxis.set_major_formatter(OOMFormatter(timing_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot(deg3d,time3d_tr, 'rs--', label='Trilinos')
plt.plot(deg3d,time3d_sc, 'bo--', label='MF scalar')
plt.plot(deg3d,time3d_t2, 'g^--', label='MF tensor2')
plt.plot(deg3d,time3d_t4, 'cv--', label='MF tensor4')
#plt.plot(deg3d,time3d_t4_ns, 'mD--', label='MF tensor4 P')
plt.xlabel('polynomial degree')
plt.ylabel('vmult wall time (s) / DoF')
leg = plt.legend(loc='best', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'timing3d.pdf', format='pdf', bbox_inches = 'tight')

# clear
thoughput_exp = 9 if "CSL" in args.prefix else 8
plt.clf()

ax = plt.figure().gca()
ax.yaxis.set_major_formatter(OOMFormatter(thoughput_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot(deg2d,through2d_tr, 'rs--', label='Trilinos')
plt.plot(deg2d,through2d_sc, 'bo--', label='MF scalar')
plt.plot(deg2d,through2d_t2, 'g^--', label='MF tensor2')
plt.plot(deg2d,through2d_t4, 'cv--', label='MF tensor4')
#plt.plot(deg2d,through2d_t4_ns, 'mD--', label='MF tensor4 P')
plt.xlabel('polynomial degree')
plt.ylabel('vmult DoF / s')
# leg = plt.legend(loc='lower right', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'throughput2d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()

ax = plt.figure().gca()
ax.yaxis.set_major_formatter(OOMFormatter(thoughput_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot(deg3d,through3d_tr, 'rs--', label='Trilinos')
plt.plot(deg3d,through3d_sc, 'bo--', label='MF scalar')
plt.plot(deg3d,through3d_t2, 'g^--', label='MF tensor2')
plt.plot(deg3d,through3d_t4, 'cv--', label='MF tensor4')
#plt.plot(deg3d,through3d_t4_ns, 'mD--', label='MF tensor4 P')
plt.xlabel('polynomial degree')
plt.ylabel('vmult DoF / s')
# leg = plt.legend(loc='lower right', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'throughput3d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()
ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.plot(deg2d,mem2d_tr, 'rs--', label='Trilinos')
plt.plot(deg2d,mem2d_sc, 'bo--', label='MF scalar')
plt.plot(deg2d,mem2d_t2, 'g^--', label='MF tensor2')
plt.plot(deg2d,mem2d_t4, 'cv--', label='MF tensor4')
#plt.plot(deg2d,mem2d_t4_ns, 'mD--', label='MF tensor4 P')
plt.xlabel('polynomial degree')
plt.ylabel('memory (Mb) / DoF')
# leg = plt.legend(loc='best', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'memory2d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()
ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.plot(deg3d,mem3d_tr, 'rs--', label='Trilinos')
plt.plot(deg3d,mem3d_sc, 'bo--', label='MF scalar')
plt.plot(deg3d,mem3d_t2, 'g^--', label='MF tensor2')
plt.plot(deg3d,mem3d_t4, 'cv--', label='MF tensor4')
#plt.plot(deg3d,mem3d_t4_ns, 'mD--', label='MF tensor4 P')
plt.xlabel('polynomial degree')
plt.ylabel('memory (Mb) / DoF')
# leg = plt.legend(loc='best', ncol=1)
y_lim = ax.get_ylim()
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'memory3d.pdf', format='pdf', bbox_inches = 'tight')

#
# Plot Section 6.2:
#

# clear
plt.clf()
ax = plt.figure().gca()
plt.yscale('log', basey=10)
plt.ylim(top=1e-2,bottom=1e-9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.rcParams.update(params2)

plt.plot(deg2d,solver2d_tr, 'rs--', label='Trilinos Solver')
# plt.plot(deg2d,solver2d_sc, 'bo--', label='MF scalar')
# plt.plot(deg2d,solver2d_t2, 'g^--', label='MF tensor2')
#plt.plot(deg2d,solver2d_t4_ns, 'mD--', label='MF Solver P')  # tensor4')
#plt.plot(deg2d,mf_gmg_2d_t4_ns,'mo--', label='MF Solver P setup')
plt.plot(deg2d,solver2d_t4, 'cv--', label='MF Solver')  # tensor4')
plt.plot(deg2d,mf_gmg_2d_t4, 'c>--', label='MF Solver setup')  # tensor4')
plt.plot(deg2d,solver2d_t4_coarse, 'g^--', label='MF Coarse Solver')
plt.xlabel('polynomial degree')
plt.ylabel('wall time (s) / DoF')
plt.plot(deg2d,assembly2d_tr, 'bp--', label='Trilinos Assembly')
leg = plt.legend(loc='lower right', ncol=1, labelspacing=0.1)
plt.savefig(fig_prefix + 'solver2d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()
ax = plt.figure().gca()
plt.yscale('log', basey=10)
plt.ylim(top=1e-2,bottom=1e-9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.rcParams.update(params2)

plt.plot(deg3d,solver3d_tr, 'rs--', label='Trilinos Solver')
# plt.plot(deg3d,solver3d_sc, 'bo--', label='MF scalar')
# plt.plot(deg3d,solver3d_t2, 'g^--', label='MF tensor2')
#plt.plot(deg3d,solver3d_t4_ns, 'mD--', label='MF Solver P')
#plt.plot(deg3d,mf_gmg_3d_t4_ns,'mo--', label='MF Solver P setup')
plt.plot(deg3d,solver3d_t4, 'cv--', label='MF Solver')
plt.plot(deg3d,mf_gmg_3d_t4, 'c>--', label='MF Solver setup')
plt.plot(deg3d,solver3d_t4_coarse, 'g^--', label='MF Coarse Solver')
plt.plot(deg3d,assembly3d_tr, 'bp--', label='Trilinos Assembly')
plt.xlabel('polynomial degree')
plt.ylabel('wall time (s) / DoF')
leg = plt.legend(loc='lower right', ncol=1, labelspacing=0.1)
plt.savefig(fig_prefix + 'solver3d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()
plt.rcParams.update(params2)

plt.plot(deg2d,cg2d_tr, 'rs--', label='Trilinos')
plt.plot(deg2d,cg2d_t4, 'cv--', label='MF')  # tensor4')
plt.xlabel('polynomial degree')
plt.ylabel('average number of CG iterations')
leg = plt.legend(loc='best', ncol=1)
plt.savefig(fig_prefix + 'cg2d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()
plt.rcParams.update(params2)

ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot(deg3d,cg3d_tr, 'rs--', label='Trilinos')
plt.plot(deg3d,cg3d_t4, 'cv--', label='MF')  # tensor4')
plt.xlabel('polynomial degree')
plt.ylabel('average number of CG iterations')
leg = plt.legend(loc='best', ncol=1)
plt.savefig(fig_prefix + 'cg3d.pdf', format='pdf', bbox_inches = 'tight')
