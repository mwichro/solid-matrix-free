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
parser.add_argument('--prefix', metavar='prefix', default='CSL_Munich', 
                    help='A folder to look for benchmark results')
parser.add_argument('--custom_model', metavar='custom_model', default='false', 
                    help='Only plot runs with none and acegen caching')
parser.add_argument('--log_scale', metavar='log_scale', default='false', 
                    help='Logarythimc scaling of y axis')
args = parser.parse_args()

prefix = args.prefix if args.prefix.startswith('/') else os.path.join(os.getcwd(), args.prefix)

files = [os.path.join(prefix, k,'timings.txt') for k in os.listdir(prefix) if os.path.isfile(os.path.join(prefix, k,'timings.txt'))]


custom_model = True if args.custom_model=="true" else False
log_scale = True if args.log_scale=="true" else False

print(' Custom model: {0}'.format(custom_model))
print(' Log scale model: {0}'.format(args.log_scale))

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

mf2d_data_scalar_ref = []
mf2d_data_scalar = []
mf2d_data_none = []
mf2d_data_acegen = []
mf2d_data_tensor2 = []
mf2d_data_tensor4 = []
mf2d_data_tensor4_ns = []
mb2d_data = []

mf3d_data_scalar_ref = []
mf3d_data_scalar = []
mf3d_data_none = []
mf3d_data_acegen = []
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
        if '_scalar_ref' in f:
            if dim == 2:
                mf2d_data_scalar_ref.append(tp)
            else:
                mf3d_data_scalar_ref.append(tp)
        elif '_scalar' in f:
            if dim == 2:
                mf2d_data_scalar.append(tp)
            else:
                mf3d_data_scalar.append(tp)
        elif '_tensor2' in f:
            if dim == 2:
                mf2d_data_tensor2.append(tp)
            else:
                mf3d_data_tensor2.append(tp)
        elif '_none' in f:
            if dim == 2:
                mf2d_data_none.append(tp)
            else:
                mf3d_data_none.append(tp)
        elif '_acegen_actual' in f:
            if dim == 2:
                mf2d_data_acegen.append(tp)
            else:
                mf3d_data_acegen.append(tp)
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
mf2d_data_none.sort(key=lambda tup: tup[0])
mf2d_data_acegen.sort(key=lambda tup: tup[0])
mf2d_data_tensor2.sort(key=lambda tup: tup[0])
mf2d_data_tensor4.sort(key=lambda tup: tup[0])
mf2d_data_tensor4_ns.sort(key=lambda tup: tup[0])
mb2d_data.sort(key=lambda tup: tup[0])

mf3d_data_scalar.sort(key=lambda tup: tup[0])
mf3d_data_scalar_ref.sort(key=lambda tup: tup[0])
mf3d_data_none.sort(key=lambda tup: tup[0])
mf3d_data_acegen.sort(key=lambda tup: tup[0])
mf3d_data_tensor2.sort(key=lambda tup: tup[0])
mf3d_data_tensor4.sort(key=lambda tup: tup[0])
mf3d_data_tensor4_ns.sort(key=lambda tup: tup[0])
mb3d_data.sort(key=lambda tup: tup[0])

# now get the data for printing
deg2d = [tup[0] for tup in mf2d_data_none]
deg3d = [tup[0] for tup in mf3d_data_none]

# vmult time per dof
time2d_tr    = [tup[4][1]/tup[1] for tup in mb2d_data]
time2d_sc    = [tup[4][0]/tup[1] for tup in mf2d_data_scalar]
time2d_nn    = [tup[4][0]/tup[1] for tup in mf2d_data_none]
time2d_ag    = [tup[4][0]/tup[1] for tup in mf2d_data_acegen]
time2d_t2    = [tup[4][0]/tup[1] for tup in mf2d_data_tensor2]
time2d_t4    = [tup[4][0]/tup[1] for tup in mf2d_data_tensor4]
time2d_t4_ns = [tup[4][0]/tup[1] for tup in mf2d_data_tensor4_ns]

time3d_tr    = [tup[4][1]/tup[1] for tup in mb3d_data]
time3d_sc    = [tup[4][0]/tup[1] for tup in mf3d_data_scalar]
time3d_nn    = [tup[4][0]/tup[1] for tup in mf3d_data_none]
time3d_ag    = [tup[4][0]/tup[1] for tup in mf3d_data_acegen]
time3d_t2    = [tup[4][0]/tup[1] for tup in mf3d_data_tensor2]
time3d_t4    = [tup[4][0]/tup[1] for tup in mf3d_data_tensor4]
time3d_t4_ns = [tup[4][0]/tup[1] for tup in mf3d_data_tensor4_ns]

through2d_tr    = [tup[1]/tup[4][1] for tup in mb2d_data]
through2d_sc    = [tup[1]/tup[4][0] for tup in mf2d_data_scalar]
through2d_nn    = [tup[1]/tup[4][0] for tup in mf2d_data_none]
through2d_ag    = [tup[1]/tup[4][0] for tup in mf2d_data_acegen]
through2d_t2    = [tup[1]/tup[4][0] for tup in mf2d_data_tensor2]
through2d_t4    = [tup[1]/tup[4][0] for tup in mf2d_data_tensor4]
through2d_t4_ns = [tup[1]/tup[4][0] for tup in mf2d_data_tensor4_ns]

through3d_tr    = [tup[1]/tup[4][1] for tup in mb3d_data]
through3d_sc    = [tup[1]/tup[4][0] for tup in mf3d_data_scalar]
through3d_scref = [tup[1]/tup[4][0] for tup in mf3d_data_scalar_ref]
through3d_nn    = [tup[1]/tup[4][0] for tup in mf3d_data_none]
through3d_ag    = [tup[1]/tup[4][0] for tup in mf3d_data_acegen]
through3d_t2    = [tup[1]/tup[4][0] for tup in mf3d_data_tensor2]
through3d_t4    = [tup[1]/tup[4][0] for tup in mf3d_data_tensor4]
through3d_t4_ns = [tup[1]/tup[4][0] for tup in mf3d_data_tensor4_ns]

# solver time per dof
solver2d_tr        = [tup[4][2]/tup[1] for tup in mb2d_data]
solver2d_sc        = [tup[4][2]/tup[1] for tup in mf2d_data_scalar]
solver2d_nn        = [tup[4][2]/tup[1] for tup in mf2d_data_none]
solver2d_ag        = [tup[4][2]/tup[1] for tup in mf2d_data_acegen]
solver2d_t2        = [tup[4][2]/tup[1] for tup in mf2d_data_tensor2]
solver2d_t4        = [tup[4][2]/tup[1] for tup in mf2d_data_tensor4]
solver2d_t4_ns     = [tup[4][2]/tup[1] for tup in mf2d_data_tensor4_ns]
solver2d_t4_coarse = [tup[4][4]/tup[1] for tup in mf2d_data_tensor4]


solver3d_tr        = [tup[4][2]/tup[1] for tup in mb3d_data]
solver3d_sc        = [tup[4][2]/tup[1] for tup in mf3d_data_scalar]
solver3d_nn        = [tup[4][2]/tup[1] for tup in mf3d_data_none]
solver3d_ag        = [tup[4][2]/tup[1] for tup in mf3d_data_acegen]
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

MB2double = 1024.0 * 1024.0 / 8.0  # Number of doubles in 1 MB

# doubles per dof
mem2d_sc    = [(tup[3]* MB2double )/tup[1] for tup in mf2d_data_scalar]
mem2d_tr    = [(tup[2]* MB2double )/tup[1] for tup in mb2d_data]
mem2d_nn    = [(tup[3]* MB2double )/tup[1] for tup in mf2d_data_none]
mem2d_ag    = [(tup[3]* MB2double )/tup[1] for tup in mf2d_data_acegen]
mem2d_t2    = [(tup[3]* MB2double )/tup[1] for tup in mf2d_data_tensor2]
mem2d_t4    = [(tup[3]* MB2double )/tup[1] for tup in mf2d_data_tensor4]
mem2d_t4_ns = [(tup[3]* MB2double )/tup[1] for tup in mf2d_data_tensor4_ns]

mem3d_tr    = [(tup[2]* MB2double)/tup[1] for tup in mb3d_data]
mem3d_sc    = [(tup[3]* MB2double)/tup[1] for tup in mf3d_data_scalar]
mem3d_nn    = [(tup[3]* MB2double)/tup[1] for tup in mf3d_data_none]
mem3d_ag    = [(tup[3]* MB2double)/tup[1] for tup in mf3d_data_acegen]
mem3d_t2    = [(tup[3]* MB2double)/tup[1] for tup in mf3d_data_tensor2]
mem3d_t4    = [(tup[3]* MB2double)/tup[1] for tup in mf3d_data_tensor4]
mem3d_t4_ns = [(tup[3]* MB2double)/tup[1] for tup in mf3d_data_tensor4_ns]

print("mem3d_ag" , mem3d_ag)
# file location
fig_prefix = os.path.join(os.getcwd(), '../doc/' + os.path.basename(os.path.normpath(prefix)) + '_')

#
#
#  MATPLOTLIB
#
#

nn_line = 'gH-.'
tr_line = 'rs--'
sc_line = 'bX:'
t4_line = 'cv:'
ag_line = 'y^-.'
t2_line = 'mD--'

tr_name = 'sparse matrix'

ne_name = 'exact log'
nP_name = 'Pade log'
ag_name = 'tensorAD'

sc_name = 'scalar curr'
sr_name = 'scalar ref.'
t4_name = 'tensor4'
t2_name = 'tensor2'
if custom_model == True:
    ne_name = 'Pade root'
    print("custom model")



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

if args.log_scale == True:
    plt.yscale('log')

plt.plot(deg2d,time2d_nn, nn_line, label='MF none')
if custom_model == False:
    plt.plot(deg2d,time2d_tr, tr_line, label=tr_name)
    plt.plot(deg2d,time2d_sc, sc_line, label=sc_name)
    plt.plot(deg2d,time2d_t4, t4_line, label=t4_name)
    plt.plot(deg2d,time2d_t2, t2_line, label=t2_name)
    # plt.plot(deg2d,time2d_t2, 'g^--', label='MF tensor2')
    # plt.plot(deg2d,time2d_t4_ns, 'mD--', label='MF tensor4 P')
else:
    plt.plot(deg2d,time2d_ag, ag_line, label= ag_name)


plt.xlabel('polynomial degree')
plt.ylabel('vmult wall time (s) / DoF')
leg = plt.legend(loc='best', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'timing2d.pdf', format='pdf', bbox_inches = 'tight')
# plt.show()


# clear
plt.clf()


ax = plt.figure().gca()
# ax.yaxis.set_major_formatter(OOMFormatter(timing_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

if args.log_scale == True:
    plt.yscale('log')

plt.plot(deg3d,time3d_nn, nn_line, label='MF none')
if custom_model == False:
    plt.plot(deg3d,time3d_tr, tr_line, label=tr_name)
    plt.plot(deg3d,time3d_sc, sc_line, label=sc_name)
    plt.plot(deg3d,time3d_t4, t4_line, label=t4_name)
    plt.plot(deg3d,time3d_t2, t2_line, label=t2_name)
    # plt.plot(deg3d,time3d_t2, 'g^--', label='MF tensor2')
    # plt.plot(deg3d,time3d_t4_ns, 'mD--', label='MF tensor4 P')
else:
    plt.plot(deg3d,time3d_ag, ag_line, label=ag_name)


plt.xlabel('polynomial degree')
plt.ylabel('vmult wall time (s) / DoF')
leg = plt.legend(loc='best', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'timing3d.pdf', format='pdf', bbox_inches = 'tight')


# comparison of impelmentations:
if custom_model == False:
    ax = plt.figure().gca()
    #ax.yaxis.set_major_formatter(OOMFormatter(timing_exp, "%1.1f"))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    if args.log_scale == True:
        plt.yscale('log')

    plt.plot(deg2d,time2d_nn, nn_line, label=ne_name)
    plt.plot(deg2d,time2d_sc, sc_line, label=sc_name)
    plt.plot(deg2d,time2d_t4, t4_line, label=t4_name)
    plt.plot(deg2d,time2d_ag, ag_line, label=ag_name)
    plt.plot(deg2d,time2d_t2, t2_line, label=t2_name)


    plt.xlabel('polynomial degree')
    plt.ylabel('vmult wall time (s) / DoF')
    leg = plt.legend(loc='best', ncol=1)
    ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
    plt.savefig(fig_prefix + 'comparison2d.pdf', format='pdf', bbox_inches = 'tight')

    ax = plt.figure().gca()
    #ax.yaxis.set_major_formatter(OOMFormatter(timing_exp, "%1.1f"))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    if args.log_scale == True:
        plt.yscale('log')

    plt.plot(deg3d,time3d_nn, nn_line, label=ne_name)
    plt.plot(deg3d,time3d_sc, sc_line, label=sc_name)
    plt.plot(deg3d,time3d_t4, t4_line, label=t4_name)
    plt.plot(deg3d,time3d_ag, ag_line, label=ag_name)
    plt.plot(deg3d,time3d_t2, t2_line, label=t2_name)


    plt.xlabel('polynomial degree')
    plt.ylabel('vmult wall time (s) / DoF')
    leg = plt.legend(loc='best', ncol=1)
    ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
    plt.savefig(fig_prefix + 'comparison3d.pdf', format='pdf', bbox_inches = 'tight')


# clear
thoughput_exp = 9 if "CSL" in args.prefix else 8
plt.clf()




ax = plt.figure().gca()
# ax.yaxis.set_major_formatter(OOMFormatter(thoughput_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

if args.log_scale == True:
    plt.yscale('log')

if custom_model == False:
    plt.plot(deg2d,through2d_nn, nn_line, label='MF none')
    plt.plot(deg2d,through2d_tr, tr_line, label= tr_name)
    plt.plot(deg2d,through2d_sc, sc_line, label= sc_name)
    plt.plot(deg2d,through2d_t4, t4_line, label= t4_name)
    plt.plot(deg2d,through2d_t2, t2_line, label= t2_name)
    plt.plot(deg2d,through2d_ag, ag_line, label= ag_name)

    # plt.plot(deg2d,through2d_t2, 'g^--', label='MF tensor2')
    # plt.plot(deg2d,through2d_t4_ns, 'mD--', label='MF tensor4 P')
else:
    plt.plot(deg2d,through2d_nn, 'gH-', label= "Halley root")
    plt.plot(deg2d,through2d_ag, 'y^-', label= "tensorAD")

plt.xlabel('polynomial degree')
plt.ylabel('vmult DoF / s')
# leg = plt.legend(loc='lower right', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'throughput2d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()

#reajustment data:

scRef_dav = [156439436.6197183, 364170491.8032787, 537218181.8181818, 56960000.0]
scalarRefSch = [4.50559284116331e8 , 10.3877703206562e8 , 13.5018642803877e8 , 17.0335570469798e8 ]

sc_curr_dav = [144249350.64935064, 347100000.00000006, 542354110.8986615, 617066666.6666666]
t4_curr_dav = [100517647.05882354, 222144000.0, 305658620.6896552, 372308379.8882681]
scalarCurrSch=  [4.386278896346e8 , 9.146905294556e8  ,12.42803877703e8  ,14.50410141685e8 ]
t4DefSch= [5.042505592841165e8, 10.39970171513795e8, 14.22967934377330e8, 17.06935123042505e8 ]

# iNH data
iNH_remputeAll= [3.3986528033456e8, 8.1230841519797e8, 11.412711168357e8, 13.568955330495e8]
iNH_TensorsDeformed = [2.8560759050389e8, 6.3747909451464e8, 9.0253829479631e8, 10.795784516927e8]
iNH_Scalar_Ref = [4.36500754147813e8, 10.108597285067873e8, 14.054298642533935e8, 16.938159879336347e8]
iNH_Scalar_Cur =  [3.8944193061840124e8, 8.226244343891402e8,11.206636500754147e8,13.318250377073907e8]

cNH_ratio_ref = [y/x for x,y in zip(scalarRefSch, scRef_dav)]
cNH_ratio_sc = [y/x for x,y in zip(scalarCurrSch, sc_curr_dav)]
cNH_ratio_t4 = [y/x for x,y in zip(t4DefSch,t4_curr_dav )]

iNH_Scalar_Ref_rescaled = [x*y for x,y in zip(iNH_Scalar_Ref, cNH_ratio_ref)]
max_iNH_Scalar_Ref_rescaled = [x * max(cNH_ratio_ref) for x in iNH_Scalar_Ref]
min_iNH_Scalar_Ref_rescaled = [x * min(cNH_ratio_ref) for x in iNH_Scalar_Ref]

iNH_Scalar_Cur_rescaled = [x*y for x,y in zip(iNH_Scalar_Cur, cNH_ratio_sc)]
max_iNH_Scalar_Curr_rescaled = [x * max(cNH_ratio_sc) for x in iNH_Scalar_Cur]
min_iNH_Scalar_Curr_rescaled = [x * min(cNH_ratio_sc) for x in iNH_Scalar_Cur]


iNH_TensorsDeformed_rescaled = [x*y for x,y in zip(iNH_TensorsDeformed, cNH_ratio_t4)]
max_iNH_TensorsDeformed_rescaled = [x * max(cNH_ratio_t4) for x in iNH_TensorsDeformed]
min_iNH_TensorsDeformed_rescaled = [x * min(cNH_ratio_t4) for x in iNH_TensorsDeformed]

iNH_remputeAll_rescaled = [x*y for x,y in zip(iNH_remputeAll, cNH_ratio_sc)]
max_iNH_remputeAll_rescaled = [x * max(cNH_ratio_sc) for x in iNH_remputeAll]
min_iNH_remputeAll_rescaled = [x * min(cNH_ratio_sc) for x in iNH_remputeAll]


ax = plt.figure().gca()
# ax.yaxis.set_major_formatter(OOMFormatter(thoughput_exp, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

if args.log_scale == True:
    plt.yscale('log')
    

# nn_line = 'gH-.'
# sc_line = 'bX:'
# t4_line = 'cv:'
# ag_line = 'y^-.'
# t2_line = 'mD--'

# ne_name = 'exact log'
# nP_name = 'Pade log'
# ag_name = 'tensorAD'


if custom_model == False:
    plt.plot(deg3d,through3d_tr, tr_line, label=tr_name)
    plt.plot(deg3d,through3d_nn, 'gH-', label= ne_name)
# if custom_model == False:
    plt.plot(deg3d,through3d_sc, sc_line, label= sc_name)
    plt.plot(deg3d,through3d_t4, t4_line, label= t4_name)
    plt.plot(deg3d,through3d_t2, t2_line, label= t2_name)
    plt.plot(deg3d,through3d_ag, ag_line, label= ag_name)
    print("through3d_scRef:", through3d_scref)
    print("through3d_sc:", through3d_sc)
    print("through3d_t4:", through3d_t4)
    # plt.plot(deg3d,through3d_t2, 'g^--', label='MF tensor2')
    # plt.plot(deg3d,through3d_t4_ns, 'mD--', label='MF tensor4 P')
else:
    plt.plot(deg3d,through3d_nn, 'gH-', label= "Halley root")
    plt.plot(deg3d,through3d_ag, 'y^-', label= "tensorAD")
    plt.plot(deg3d,max_iNH_Scalar_Curr_rescaled, "bX-.", label= "scalar curr.")
    plt.plot(deg3d,max_iNH_Scalar_Ref_rescaled, "bX:", label= "scalar ref.")
    plt.plot(deg3d,max_iNH_TensorsDeformed_rescaled, 'y^:', label= "tensor4, curr")
    # plt.fill_between(deg3d, min_iNH_TensorsDeformed_rescaled, max_iNH_TensorsDeformed_rescaled, color='yellow', alpha=0.3)
    plt.plot(deg3d,max_iNH_remputeAll_rescaled, 'gH:', label= "recompute all")
    plt.fill_between(deg3d, min_iNH_remputeAll_rescaled, max_iNH_remputeAll_rescaled, color='green', alpha=0.2)


plt.xlabel('polynomial degree')
plt.ylabel('vmult DoF / s')
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'throughput3d.pdf', format='pdf', bbox_inches = 'tight')

leg = plt.legend(loc='lower right', ncol=1)
# save the legend in a separate figure
fig_leg = figure(figsize=(10, 1))
ax_leg = fig_leg.add_subplot(111)
ax_leg.axis('off')
if custom_model == False:
    legend = ax_leg.legend(*ax.get_legend_handles_labels(), loc='center', frameon=True, ncol=3)
else:
    legend = ax_leg.legend(*ax.get_legend_handles_labels(), loc='center', frameon=True, ncol=3)
fig_leg.canvas.draw()
bbox  = legend.get_window_extent().transformed(fig_leg.dpi_scale_trans.inverted())
fig_leg.savefig(fig_prefix + 'legend-comparison.pdf', format='pdf', bbox_inches=bbox)



# clear
plt.clf()
ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# ax.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

if args.log_scale == True:
    plt.yscale('log')

plt.plot(deg2d,mem2d_nn, nn_line, label=ne_name )
plt.plot(deg2d,mem2d_ag, ag_line, label=ag_name)
if custom_model == False:
    # plt.plot(deg2d,mem2d_tr, tr_line, label='Trilinos')
    plt.plot(deg2d,mem2d_sc, sc_line, label=sc_name)
    plt.plot(deg2d,mem2d_t4, t4_line, label=t4_name)
    plt.plot(deg2d,mem2d_t2, t2_line, label=t2_name)
    # plt.plot(deg2d,mem2d_ag, ag_line, label='MF SS')
    # plt.plot(deg2d,mem2d_t4_ns, 'mD--', label='MF tensor4 P')

plt.xlabel('polynomial degree')
plt.ylabel('storage size per DoF')
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))  # Use integer format
# leg = plt.legend(loc='best', ncol=1)
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'memory2d.pdf', format='pdf', bbox_inches = 'tight')

# leg = plt.legend(loc='lower right', ncol=1)
# save the legend in a separate figure
fig_leg = figure(figsize=(10, 1))
ax_leg = fig_leg.add_subplot(111)
ax_leg.axis('off')
legend = ax_leg.legend(*ax.get_legend_handles_labels(), loc='center', frameon=True, ncol=3)
fig_leg.canvas.draw()
bbox  = legend.get_window_extent().transformed(fig_leg.dpi_scale_trans.inverted())
fig_leg.savefig(fig_prefix + 'legend-mem.pdf', format='pdf', bbox_inches=bbox)




# clear
plt.clf()
ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# ax.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f"))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

if args.log_scale == True:
    plt.yscale('log')
    
plt.plot(deg3d,mem3d_nn, nn_line, label=ne_name)
plt.plot(deg3d,mem3d_ag, ag_line, label=ag_name)
if custom_model == False:
    # plt.plot(deg3d,mem3d_tr, tr_line, label='Trilinos')
    plt.plot(deg3d,mem3d_sc, sc_line, label=sc_name)
    plt.plot(deg3d,mem3d_t4, t4_line, label=t4_name)
    plt.plot(deg3d,mem3d_t2, t2_line, label=t2_name)
    # plt.plot(deg3d,mem3d_t4_ns, 'mD--', label='MF tensor4 P')
# else:
#    plt.plot(deg3d,mem3d_ag, ag_line, label='MF SS')

plt.xlabel('polynomial degree')
plt.ylabel('storage size per DoF')
# plt.ylim(25, 125)  # Set y-axis range directly to 25-125
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))  # Use integer format

# leg = plt.legend(loc='best', ncol=1)
y_lim = ax.get_ylim()
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
plt.savefig(fig_prefix + 'memory3d.pdf', format='pdf', bbox_inches = 'tight')



#
# Plot Section 6.2:
#

# clear


if custom_model:
    exit()

plt.clf()
ax = plt.figure().gca()
plt.yscale('log', base=10)
plt.ylim(top=1e-2,bottom=1e-9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.rcParams.update(params2)

plt.plot(deg2d,assembly2d_tr, 'bp--', label='Matrix Assembly')
plt.plot(deg2d,solver2d_tr, 'rs--', label='Matrix-based solver')
# plt.plot(deg2d,solver2d_sc, 'bo--', label='MF scalar')
# plt.plot(deg2d,solver2d_t2, 'g^--', label='MF tensor2')
#plt.plot(deg2d,solver2d_t4_ns, 'mD--', label='MF Solver P')  # tensor4')
#plt.plot(deg2d,mf_gmg_2d_t4_ns,'mo--', label='MF Solver P setup')
plt.plot(deg2d,solver2d_t4, 'cv-', label='Matrix-free solver')  # tensor4')
plt.plot(deg2d,mf_gmg_2d_t4, 'm>-', label='Matrix-free setup')  # tensor4')
# plt.plot(deg2d,solver2d_t4_coarse, 'g^--', label='MF Coarse Solver')
plt.xlabel('polynomial degree')
plt.ylim(top=1e-3, bottom=1e-7)  # Set y-axis scale from 1e-7 to 1e-3
plt.ylabel('wall time (s) / DoF')
leg = plt.legend(loc='lower right', ncol=1, labelspacing=0.1)
plt.savefig(fig_prefix + 'solver2d.pdf', format='pdf', bbox_inches = 'tight')

# clear
plt.clf()
ax = plt.figure().gca()
plt.yscale('log', base=10)
plt.ylim(top=1e-2,bottom=1e-9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.rcParams.update(params2)

plt.plot(deg3d,assembly3d_tr, 'bp--', label='Matrix Assembly')
plt.plot(deg3d,solver3d_tr, 'rs--', label='Matrix-based solver')
# plt.plot(deg3d,solver3d_sc, 'bo--', label='MF scalar')
# plt.plot(deg3d,solver3d_t2, 'g^--', label='MF tensor2')
#plt.plot(deg3d,solver3d_t4_ns, 'mD--', label='MF Solver P')
#plt.plot(deg3d,mf_gmg_3d_t4_ns,'mo--', label='MF Solver P setup')
plt.plot(deg3d,solver3d_t4, 'cv-', label='Matrix-free solver')
plt.plot(deg3d,mf_gmg_3d_t4, 'm>-', label='Matrix-free setup')  
# plt.plot(deg3d,solver3d_t4_coarse, 'g^--', label='MF Coarse Solver')
plt.xlabel('polynomial degree')
plt.ylabel('wall time (s) / DoF')
plt.ylim(top=3e-3, bottom=1e-7)  # Set y-axis scale from 1e-7 to 1e-3
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
plt.plot(deg3d,cg3d_tr, tr_line, label='Trilinos')
plt.plot(deg3d,cg3d_t4, 'cv--', label='MF')  # tensor4')
plt.xlabel('polynomial degree')
plt.ylabel('average number of CG iterations')
leg = plt.legend(loc='best', ncol=1)
plt.savefig(fig_prefix + 'cg3d.pdf', format='pdf', bbox_inches = 'tight')
