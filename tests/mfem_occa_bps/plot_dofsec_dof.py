#!/usr/bin/env python
#
# -- Plot the results --
#
# Example Usage:
#   ./plot_dofsec_dof.py bp1_occa_cuda.log "BP1 - CUDA (ray)" bp3_occa_cuda.log "BP3 - CUDA (ray)"
#
# The log files are stdout redirects from running the tests in this
# directory using the go.sh script in the base directory.

from sys import stdout as out
import fileinput
import pprint
import sys
from math import ceil
import matplotlib.pyplot as plt
from os.path import splitext
from itertools import groupby
# import matplotlib as mpl
# mpl.style.use('classic')


def load_log(filename):
    '''Read all input files specified on the command line, or stdin and parse
    the content, storing it as a list of dictionaries - one dictionary for
    each test run.'''
    it = fileinput.input(filename)
    state=0
    line=''
    i=0
    mesh_p=0
    config='unknown'
    compiler='unknown'
    test='unknown'
    num_procs=0
    num_procs_node=0
    lnfmt='%05i'
    data={}
    runs=[]
    while True:
        ##
        if state%2==0:
            ##
            try:
                line=it.next()
                i=i+1
            except StopIteration:
                break
            state=state+1
            ##
        elif state==1:
            ##
            state=0
            if 'Reading configuration' in line:
                ##
                ## This is the beginning of a new file.
                ##
                # out.write('\n'+lnfmt%i+': %s'%line)
                config=line.split()[2]
                num_procs=0
                num_procs_node=0
            elif 'Setting up compiler' in line:
                # out.write(lnfmt%i+': %s'%line)
                compiler=line.split()[3]
            elif 'Reading test file: ' in line:
                # out.write(lnfmt%i+': %s'%line)
                test=line.strip().split('Reading test file: ',1)[-1]
            elif 'Running the tests using a total of' in line:
                # out.write(lnfmt%i+': %s'%line)
                num_procs=int(line.split('a total of ',1)[1].split(None,1)[0])
            elif 'tasks per node' in line:
                # out.write(lnfmt%i+': %s'%line)
                num_procs_node=int(line.split(' tasks per',1)[0].rsplit(None,1)[1])
            elif line == 'Options used:\n':
                ##
                ## This is the beginning of a new run.
                ##
                if 'cg-iteration-dps' in data:
                    runs.append(data)
                    # print
                    # pprint.pprint(data)
                data={}
                data['config']=config
                data['compiler']=compiler
                data['test']=test
                data['num-procs']=num_procs
                data['num-procs-node']=num_procs_node
                # out.write('\n'+lnfmt%i+': %s'%line)
                # out.write('*'*len(lnfmt%i)+':    --mesh-p %i\n'%mesh_p)
                state=2
            elif 'using integration rule with' in line:
                # out.write(lnfmt%i+': %s'%line)
                data['quadrature-pts']=int(line.split(' ')[5])
            elif '"DOFs/sec" in assembly' in line:
                # out.write(lnfmt%i+': %s'%line)
                data['assembly-dps']=1e6*float(line.split(' ')[3])
            elif '"DOFs/sec" in CG' in line:
                # out.write(lnfmt%i+': %s'%line)
                data['cg-iteration-dps']=1e6*float(line.split(' ')[3])
            elif 'Number of finite element unknowns' in line:
                # out.write(lnfmt%i+': %s'%line)
                data['num-unknowns']=long(line.rsplit(None,1)[1])
            elif ' elements'==line[:9]:
                # out.write(lnfmt%i+': %s'%line)
                data['num-elem']=int(line.rsplit(None,1)[1])
                ##
        elif state==3:
            ##
            if line[:5] == '   --':
                # out.write(lnfmt%i+': %s'%line)
                state=2
                opt=line[5:].strip().split(' ',1)
                if opt[0] in ('refine-serial', 'refine-parallel', 'order', 'ref-levels'):
                    data[opt[0]]=int(opt[1])
                elif opt[0] in ('mesh', 'basis-type', 'preconditioner'):
                    data[opt[0]]=opt[1]
                elif opt[0] in ('standard-version', 'hpc-version'):
                    data['std | hpc']=opt[0]
                elif opt[0] in ('assembly', 'matrix-free'):
                    data['action-type']=opt[0]
            else:
                state=1

    if 'cg-iteration-dps' in data:
        runs.append(data)
    for run in runs:
        if not 'quadrature-pts' in run:
            run['quadrature-pts']=0

    print 'Number of test runs read: %i'%len(runs)

    return runs

def make_plot(data, title=None):
    '''Make a plot from the dictionaries.'''
    fig = plt.figure()
    a = fig.gca()

    orders = list(set([d['order'] for d in data]))

    for p in orders:
        pdata = sorted(filter(lambda x: x['order'] == p, data),
                       key=lambda x: x['quadrature-pts'])

        # Loop over sets of the same number of quadrature points
        for k, g in groupby(pdata, lambda x: x['quadrature-pts']):
            glist = list(g)

            dof = [d['num-unknowns'] for d in glist]
            dps = [d['cg-iteration-dps'] for d in glist]

            line = '-'
            marker = 'o'

            if k > 0:
                qpts = int(ceil(k**(1.0/3.0)))
                if qpts == p+1:
                    linestyle = '--'
                label = "p={:d}, q={:s}".format(p, "p+"+str(qpts-p))
            else:
                label = "p={:d}".format(p)

            a.plot(dof, dps, label=label, linestyle=line, marker=marker)

    a.set_xscale('log')
    a.set_yscale('log')

    a.set_xlim((1e2, 2e7))
    a.set_ylim((1e5, 2e9))

    a.set_xlabel('Points per compute node')
    a.set_ylabel(r'[DOF $\times$ CG iterations] / [compute nodes $\times$ seconds]')

    if title is not None:
        a.set_title(title)

    a.grid(True, which='both', linestyle=':', linewidth=0.25)
    a.legend(loc='best')

    return fig

if __name__ == '__main__':
    for f, title in zip(sys.argv[1::2], sys.argv[2::2]):
        d = load_log(f)
        figfile = splitext(f)[0] + ".pdf"

        fig = make_plot(d, title)

        fig.savefig(figfile, bbox_inches='tight')
