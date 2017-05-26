
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
         data['num-procs']=num_procs
         data['num-procs-node']=num_procs_node
         data['mesh-order']=mesh_p
         # out.write('\n'+lnfmt%i+': %s'%line)
         # out.write('*'*len(lnfmt%i)+':    --mesh-p %i\n'%mesh_p)
         state=2
      elif '"DOFs/sec" in assembly' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['assembly-dps']=1e6*float(line.split(' ')[3])
      elif '"DOFs/sec" in CG' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['cg-iteration-dps']=1e6*float(line.split(' ')[3])
      elif 'Number of finite element unknowns' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['num-unknowns']=long(line.rsplit(None,1)[1])
      elif 'MESH_P=' in line:
         # out.write('\n'+lnfmt%i+': %s'%line)
         mesh_p=int(line.partition('MESH_P=')[2].split(None,1)[0])
      elif 'Average nonzero entries per row' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['entries/row']=float(line.rsplit(' ',1)[1])
      elif 'using integration rule with' in line:
         # out.write(lnfmt%i+': %s'%line)
         data['quadrature-pts']=int(line.split('with ',1)[1].split(' ',1)[0])
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
         if opt[0] in ('refine-serial', 'refine-parallel', 'order'):
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
   # print
   # pprint.pprint(data)
for run in runs:
   if not 'quadrature-pts' in run:
      run['quadrature-pts']=0
# print
# print
# pprint.pprint(runs)

print 'Number of test runs read: %i'%len(runs)
