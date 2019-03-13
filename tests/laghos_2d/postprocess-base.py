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
      #out.write(lnfmt%i+'\033[37m: %s\033[m'%line)
      if 'Reading configuration' in line:
         ##
         ## This is the beginning of a new file.
         ##
         #out.write('\n'+lnfmt%i+'\033[33m: %s\033[m'%line)
         config=line.split()[2]
         num_procs=0
         num_procs_node=0
      elif 'Setting up compiler' in line:
         #out.write(lnfmt%i+'\033[33m: %s\033[m'%line)
         compiler=line.split()[3]
      elif 'Reading test file: ' in line:
         #out.write(lnfmt%i+'\033[33m: %s\033[m'%line)
         test=line.strip().split('Reading test file: ',1)[-1]
      elif 'Running the tests using a total of' in line:
         #out.write(lnfmt%i+'\033[33m: %s\033[m'%line)
         num_procs=int(line.split('a total of ',1)[1].split(None,1)[0])
      elif 'tasks per node' in line:
         #out.write(lnfmt%i+'\033[33m: %s\033[m'%line)
         num_procs_node=int(line.split(' tasks per',1)[0].rsplit(None,1)[1])
      elif line == 'Options used:\n':
         #out.write('\nBeginning of a new run:\n')
         ##
         ## This is the beginning of a new run.
         ##
         #out.write('\033[33mcg-H1-rate\033[m')
         if 'cg-H1-rate' in data:
            runs.append(data)
         #print
         #pprint.pprint(data)
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
         data['case']='scalar'
         data['action-type']='partial-assembly'
         data['quadrature-pts']=1
         data['quadrature-type']=1
         state=2
      elif 'CG (H1) rate' in line:
         #out.write(lnfmt%i+'\033[31m: %s\033[m'%line)
         data['cg-H1-rate']=1e6*float(line.split(' ')[8])
         #out.write(lnfmt%i+'\033[31m: %s\033[m'%line.split(' ')[8])
      elif 'Number of local/global kinematic' in line:
         #out.write(lnfmt%i+'\033[31m: %s\033[m'%line)
         data['kinematic-dofs']=long(line.split(' ')[7].split('/')[1])
         #print data['kinematic-dofs']
      elif 'Number of kinematic' in line:
         #out.write(lnfmt%i+'\033[31m: %s\033[m'%line)
         data['kinematic-dofs']=long(line.split(' ')[6])
         #print data['kinematic-dofs']
      elif 'Zones' in line:
         #out.write(lnfmt%i+'\033[31m: %s\033[m'%line)
         data['num-elem']=long(line.split(' ')[2])
         data['num-elem-min']=long(line.split(' ')[2])
         data['num-elem-max']=long(line.split(' ')[3])
         if (data['num-elem-min'] != data['num-elem-max']) :
            raise ValueError('num-elem-min != num-elem-max!')
      ##
   elif state==3:
      ##
      if line[:5] == '   --':
         #out.write(lnfmt%i+'\033[32m: %s\033[m'%line)
         state=2
         opt=line[5:].strip().split(' ',1)
         if opt[0] in ('refine-serial', 'refine-parallel'):
            data[opt[0]]=int(opt[1])
         elif opt[0] in ('mesh'):
            data[opt[0]]=opt[1]
         elif opt[0] in ('order-kinematic'):
            data['order']=int(opt[1])
            #out.write(lnfmt%i+'\033[32;1m: order: '+opt[1]+'\033[m\n')
         elif opt[0] in ('partial-assembly', 'full-assembly'):
            data['action-type']=opt[0]
            #out.write(lnfmt%i+'\033[32;1m: action-type: '+opt[0]+'\033[m\n')
      else:
         state=1

if 'cg-H1-rate' in data:
   runs.append(data)

#print
#pprint.pprint(data)

for run in runs:
   if not 'quadrature-pts' in run:
      run['quadrature-pts']=0
# print
#pprint.pprint(runs)

print 'Number of test runs read: %i'%len(runs)
