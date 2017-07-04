# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project
# (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
# organizations (Office of Science and the National Nuclear Security
# Administration) responsible for the planning and preparation of a capable
# exascale ecosystem, including software, applications, hardware, advanced
# system engineering and early testbed platforms, in support of the nation's
# exascale computing imperative.


#####   Load the data
execfile('postprocess-base.py')


#####   Sample plot output

from pylab import *

rcParams['font.sans-serif'].insert(0,'Noto Sans')
rcParams['font.sans-serif'].insert(1,'Open Sans')
rcParams['figure.figsize']=[10, 8] # default: 8 x 6

cm=get_cmap('Set1') # 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3'
if '_segmentdata' in cm.__dict__:
   cm_size=len(cm.__dict__['_segmentdata']['red'])
elif 'colors' in cm.__dict__:
   cm_size=len(cm.__dict__['colors'])
colors=[cm(1.*i/(cm_size-1)) for i in range(cm_size)]

# colors=['blue','green','crimson','turquoise','m','gold','cornflowerblue',
#         'darkorange']

sel_runs=runs
action_type=sel_runs[0]['action-type']
print 'Using action-type:', action_type
sel_runs=[run for run in sel_runs if run['action-type']==action_type]

for run in sel_runs:
   run['config']=run['config'].rsplit('/',1)[-1].rsplit('.sh',1)[0]
configs=[]
for run in sel_runs:
   if not (run['config'] in configs):
      configs.append(run['config'])
print 'Using configurations:', configs

for run in sel_runs:
   run['test']=run['test'].rsplit('/',1)[-1].rsplit('.sh',1)[0]
tests=[]
for run in sel_runs:
   if not (run['test'] in tests):
      tests.append(run['test'])
print 'Using tests:', tests

compilers=[]
for run in sel_runs:
   if not (run['compiler'] in compilers):
      compilers.append(run['compiler'])
print 'Using compilers:', compilers

if 'case' in sel_runs[0]:
   cases=list(set([run['case'] for run in sel_runs]))
   print 'Using cases:', cases
for run in sel_runs:
   run['vdim']=3 if ('case' in run and run['case']=='vector') else 1

codes = list(set([run['code'] for run in sel_runs]))
code  = codes[0]
sel_runs=[run for run in sel_runs if run['code']==code]

if len(configs)>1:
   key='config'
   val=configs[0]
   val2=configs[1]
elif len(tests)>1:
   key='test'
   val=tests[0]
   val2=tests[1]
elif len(compilers)>1:
   key='compiler'
   val=compilers[0]
   val2=compilers[1]
else:
   key='case'
   val='scalar'
   val2='vector'

pl_set=[(run['num-procs'],run['num-procs-node'])
        for run in sel_runs]
pl_set=sorted(set(pl_set))
print
pprint.pprint(pl_set)

for plt in pl_set:
   num_procs=plt[0]
   num_procs_node=plt[1]
   num_nodes=num_procs/num_procs_node
   pl_runs=[run for run in sel_runs
            if run['num-procs']==num_procs and
               run['num-procs-node']==num_procs_node]
   pl_runs=sorted(pl_runs)
   if len(pl_runs)==0:
      continue

   print
   print 'compute nodes: %i, number of MPI tasks = %i'%(num_nodes,num_procs)

   figure()
   i=0
   sol_p_set=sorted(set([run['order'] for run in pl_runs]))
   for sol_p in sol_p_set:
      qpts=sorted(list(set([run['quadrature-pts'] for run in pl_runs
                            if run['order']==sol_p])))
      qpts.reverse()
      print 'Order: %i, quadrature points:'%sol_p, qpts
      qpts_1d=[int(q**(1./3)+0.5) for q in qpts]

      pl_runs_1=[run for run in pl_runs if run[key]==val]
      pl_runs_2=[run for run in pl_runs if run[key]!=val]
      d1=[[run['order'],run['num-elem'],
           1.*run['num-unknowns']/num_nodes/run['vdim'],
           run['cg-iteration-dps']/num_nodes]
          for run in pl_runs_1
          if run['order']==sol_p and
             run['quadrature-pts']==qpts[0]]
      d2=[[run['order'],run['num-elem'],
           1.*run['num-unknowns']/num_nodes/run['vdim'],
           run['cg-iteration-dps']/num_nodes]
          for run in pl_runs_2
          if run['order']==sol_p and
             run['quadrature-pts']==qpts[0]]
      di=set([e[2] for e in d1]).intersection(set([e[2] for e in d2]))
      if len(di)>0:
         d=[[npts,
             max([e[3] for e in d1 if e[2]==npts]),
             max([e[3] for e in d2 if e[2]==npts])]
            for npts in di]
         d=asarray(sorted(d))
         plot(d[:,0],d[:,2]/d[:,1],'o-',color=colors[i],
              label='p=%i, q=p+%i'%(sol_p,qpts_1d[0]-sol_p))
      ##
      if len(qpts)==1:
         i=i+1
         continue
      #
      pl_runs_1=[run for run in pl_runs if run[key]==val]
      pl_runs_2=[run for run in pl_runs if run[key]!=val]
      d1=[[run['order'],run['num-elem'],
           1.*run['num-unknowns']/num_nodes/run['vdim'],
           run['cg-iteration-dps']/num_nodes]
          for run in pl_runs_1
          if run['order']==sol_p and
             run['quadrature-pts']==qpts[1]]
      d2=[[run['order'],run['num-elem'],
           1.*run['num-unknowns']/num_nodes/run['vdim'],
           run['cg-iteration-dps']/num_nodes]
          for run in pl_runs_2
          if run['order']==sol_p and
             run['quadrature-pts']==qpts[1]]
      di=set([e[2] for e in d1]).intersection(set([e[2] for e in d2]))
      if len(di)>0:
         d=[[npts,
             max([e[3] for e in d1 if e[2]==npts]),
             max([e[3] for e in d2 if e[2]==npts])]
            for npts in di]
         d=asarray(sorted(d))
         plot(d[:,0],d[:,2]/d[:,1],'s--',color=colors[i],
              label='p=%i, q=p+%i'%(sol_p,qpts_1d[1]-sol_p))
      ##
      i=i+1
   ##

   title('Config: %s %s (%i node%s, %i tasks/node), %s, %s, %s'%(
         code,str('-').join(configs),num_nodes,'' if num_nodes==1 else 's',
         num_procs_node,str('-').join(compilers),str('-').join(tests),
         'PA' if action_type=='matrix-free' else 'TA'))
   xscale('log') # subsx=[2,4,6,8]
   # yscale('log')
   # rng=arange(1e7,1.02e8,1e7)
   # yticks(rng,['%i'%int(v/1e6) for v in rng])
   # ylim(min(rng),max(rng))
   # xlim(0.5,max([run['order'] for run in pl_runs])+0.5)
   grid('on', color='gray', ls='dotted')
   grid('on', axis='both', which='minor', color='gray', ls='dotted')
   gca().set_axisbelow(True)
   xlabel('Points per compute node')
   ylabel('[%s dps] / [%s dps]'%(str(val2),str(val)))
   legend(ncol=2, loc='best')

   if 1: # write .pdf file?
      pdf_file='plot4_%s_%s_%s_%s_N%03i_pn%i.pdf'%(
               code,str('-').join(tests),str('-').join(configs),
               str('-').join(compilers),num_nodes,num_procs_node)
      print 'saving figure --> %s'%pdf_file
      savefig(pdf_file, format='pdf', bbox_inches='tight')

if 1: # show the figures?
   print '\nshowing figures ...'
   show()
