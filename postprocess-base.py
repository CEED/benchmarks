# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supportted by the Exascale Computing Project
# (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
# organizations (Office of Science and the National Nuclear Security
# Administration) responsible for the planning and preparation of a capable
# exascale ecosystem, including software, applications, hardware, advanced
# system engineering and early testbed platforms, in support of the nation's
# exascale computing imperative.

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
mesh_p=1
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
         data['file']=fileinput.filename()
         data['config']=config
         data['compiler']=compiler
         data['test']=test
         data['num-procs']=num_procs
         data['num-procs-node']=num_procs_node
         data['mesh-order']=mesh_p
         data['code']="MFEM"
         test_=test.rsplit('/',1)[-1]
         data['case']='scalar' if (('bp1' in test_) or ('bp3' in test_) or
                                   ('bp5' in test_)) else 'vector'
         # out.write('\n'+lnfmt%i+': %s'%line)
         # out.write('*'*len(lnfmt%i)+':    --mesh-p %i\n'%mesh_p)
         state=2
      elif 'setup' in line: 
         data['amg-setup'] = float(line.split()[1])
      elif 'Total iterations' in line: 
         data['iterations'] = float(line.split()[2])
         print(line.split())
      elif '"DOFs/sec" in CG' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['cg-iteration-dps']=1e6*float(line.split(' ')[3])
      elif 'Number of finite element unknowns' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['num-unknowns']=long(line.rsplit(None,1)[1])
      elif 'Number of qudrature points per element' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['quadrature-pts']=int(line.split(' = ',1)[1])
      elif 'Total number of elements:' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['num-elem']=int(line.rsplit(None,1)[1])
      ##
   elif state==3:
      ##
      if line[:5] == '   --':
         # out.write(lnfmt%i+': %s'%line)
         state=2
         opt=line[5:].strip().split(' ',1)
         if opt[0] in ('order'):
            data[opt[0]]=int(opt[1])
         elif opt[0] in ('device'):
            data['mfem-device']=opt[1]
      else:
         state=1

if 'cg-iteration-dps' in data:
   runs.append(data)
   # print
   # pprint.pprint(data)
for run in runs:
   if not 'quadrature-pts' in run:
      run['quadrature-pts']=0
# print
# print
# pprint.pprint(runs)

print 'Number of test runs read: %i'%len(runs)
