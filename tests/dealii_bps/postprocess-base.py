# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
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
# testbed platforms, in support of the nation's exascale computing imperative.

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
      elif ' p |  q | n_elements |' in line:
         ##
         ## This is the beginning of a new run.
         ##
         data={}
         data['file']=fileinput.filename()
         data['config']=config
         data['compiler']=compiler
         data['test']=test
         data['num-procs']=num_procs
         data['num-procs-node']=num_procs_node
         data['mesh-order']=mesh_p
         data['code']="deal.II"
         # read next line and go to state 3:
         state=2
      ##
   elif state==3:
      ##
      if len(line.strip()) == 0:
         # read next line and go to state 1:
         state=0
      else:
         ##
         e=line.split()
         # p,q,n_elements,n_dofs,time/it,dofs/s/it,CG_its,time/matvec
         # 0,2,         4,     6,      8,       10,    12,         14
         # out.write(str(e)+'\n')
         data['order']=int(e[0])
         data['assembly-dps']=0
         data['cg-iteration-dps']=float(e[10])
         data['num-unknowns']=int(e[6])
         data['quadrature-pts']=int(e[2])**3
         data['num-elem']=int(e[4])
         data['action-type']='matrix-free'
         data['std | hpc']='hpc'
         data['case']='scalar' if ('/bp1.sh' in data['test']) else 'vector'
         runs.append(data.copy())
         # read next line and come back to state 3:
         state=2
      ##
   ##

# print
# print
# pprint.pprint(runs)

print 'Number of test runs read: %i'%len(runs)
