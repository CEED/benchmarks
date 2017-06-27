# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation’s exascale computing imperative.

from sys import stdout as out
import fileinput
import pprint

#####   Read all input files specified on the command line, or stdin and parse
#####   the content, storing it as a list of dictionaries - one dictionary for
#####   each test run.

it=fileinput.input()
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
      elif line == '| NEK5000:  Open Source Spectral Element Solver            |\n':
         ##
         ## This is the beginning of a new run.
         ##
         # if 'cg-iteration-dps' in data:
         #    runs.append(data)
            # print
            # pprint.pprint(data)
         data={}
         data['config']=config
         data['compiler']=compiler
         data['test']=test
         data['num-procs']=num_procs
         data['num-procs-node']=num_procs_node
         data['mesh-order']=mesh_p
         # out.write('\n'+lnfmt%i+': %s'%line)
         # out.write('*'*len(lnfmt%i)+':    --mesh-p %i\n'%mesh_p)
      elif 'case scalar:' in line:
         # out.write(lnfmt%i+': %s'%line)
         ##  np,nx,nelt,nelgt,ndof,nppp,maxits,telaps,dofpss,titers,tppp_s
         ##   2, 3,   4,    5,   6,   7,     8,     9,    10,    11,    12
         e=line.split()
         # out.write(str(e)+'\n')
         data['order']=int(e[3])
         data['assembly-dps']=0
         data['cg-iteration-dps']=int(e[6])*int(e[8])*3/float(e[9])
         data['num-unknowns']=int(e[6])
         data['quadrature-pts']=(int(e[3])+3)**3 ##  !!! q=(p+3)^3
         data['num-elem']=int(e[5])
         data['action-type']='matrix-free'
         data['std | hpc']='hpc'
         data['case']='scalar'
         runs.append(data.copy())
      ##
      elif 'case vector:' in line:
         # out.write(lnfmt%i+': %s'%line)
         ##  np,nx,nelt,nelgt,ndof,nppp,maxitv,telapv,dofpsv,titerv,tppp_v
         ##   2, 3,   4,    5,   6,   7,     8,     9,    10,    11,    12
         e=line.split()
         # out.write(str(e)+'\n')
         data['order']=int(e[3])
         data['assembly-dps']=0
         data['cg-iteration-dps']=3*int(e[6])*int(e[8])/float(e[9])
         data['num-unknowns']=3*int(e[6])
         data['quadrature-pts']=(int(e[3])+3)**3 ##  !!! q=(p+3)^3
         data['num-elem']=int(e[5])
         data['action-type']='matrix-free'
         data['std | hpc']='hpc'
         data['case']='vector'
         runs.append(data.copy())
      ##

for run in runs:
   if not 'quadrature-pts' in run:
      run['quadrature-pts']=0
# print
# print
# pprint.pprint(runs)

print 'Number of test runs read: %i'%len(runs)
