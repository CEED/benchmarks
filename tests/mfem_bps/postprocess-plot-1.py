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
rcParams['legend.fontsize']='medium' # default: 'large'

cm=get_cmap('Set1') # 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3'
if '_segmentdata' in cm.__dict__:
   cm_size=len(cm.__dict__['_segmentdata']['red'])
elif 'colors' in cm.__dict__:
   cm_size=len(cm.__dict__['colors'])
colors=[cm(1.*i/(cm_size-1)) for i in range(cm_size)]

# Down-select from all runs:
sel_runs=runs
configs=list(set([run['config'] for run in sel_runs]))
# print 'Present configurations:', configs
config=configs[0]
print 'Using configuration:', config
config_short=config.rsplit('/',1)[-1].rsplit('.sh',1)[0]
sel_runs=[run for run in sel_runs if run['config']==config]

compilers=list(set([run['compiler'] for run in sel_runs]))
compilers=compilers[0:2]
print 'Using compilers:', compilers
sel_runs=[run for run in sel_runs if run['compiler'] in compilers]

tests=list(set([run['test'] for run in sel_runs]))
# print 'Present tests:', tests
test=tests[0]
print 'Using test:', test
test_short=test.rsplit('/',1)[-1].rsplit('.sh',1)[0]
sel_runs=[run for run in sel_runs if run['test']==test]

action_type=sel_runs[0]['action-type']
print 'Using action-type:', action_type
sel_runs=[run for run in sel_runs if run['action-type']==action_type]

if 'case' in sel_runs[0]:
   cases=list(set([run['case'] for run in sel_runs]))
   case=cases[0]
   print 'Using case:', case
   sel_runs=[run for run in sel_runs if run['case']==case]

codes = list(set([run['code'] for run in sel_runs]))
code  = codes[0]
sel_runs=[run for run in sel_runs if run['code']==code]

pl_set=[run['num-procs']/run['num-procs-node'] for run in sel_runs]
pl_set=sorted(set(pl_set))
# print
# pprint.pprint(pl_set)

for plt in pl_set:
   num_nodes=plt
   pl_runs=[run for run in sel_runs
            if run['num-procs']==num_nodes*run['num-procs-node']]
   pl_runs=sorted(pl_runs)
   if len(pl_runs)==0:
      continue

   np_set=sorted(set([run['num-procs'] for run in pl_runs]))
   print
   print 'compute nodes: %i, np ='%num_nodes, np_set

   figure()
   i=0
   for np in np_set:
      rel_qpts=[int(run['quadrature-pts']**(1./3)+0.5)-run['order']
                for run in pl_runs if run['num-procs']==np]
      rel_qpts=sorted(list(set(rel_qpts)))
      rel_qpts.reverse()
      print 'q-p =', rel_qpts

      d=[[run['order'],run['num-elem'],1.*run['num-unknowns']/num_nodes,
          run['cg-iteration-dps']/num_nodes,run['compiler']]
         for run in pl_runs
         if run['num-procs']==np and
            run['quadrature-pts']==(run['order']+rel_qpts[0])**3]
      dps_idx=3
      # print
      # print 'np = %i'%np
      # pprint.pprint(sorted(d))
      d=[[order,
          max([e[dps_idx] for e in d if e[0]==order and e[4]==compilers[0]]),
          max([e[dps_idx] for e in d if e[0]==order and e[4]==compilers[-1]])]
         for order in set([e[0] for e in d])]
      d=asarray(sorted(d))
      tpn=np/num_nodes
      plot(d[:,0],d[:,1],'o-',color=colors[i],
           label='%s, %i t/n, q=p+%i'%(compilers[0],tpn,rel_qpts[0]))
      if len(compilers)==2:
         plot(d[:,0],d[:,2],'d-',color=colors[i],
              label='%s, %i t/n, q=p+%i'%(compilers[-1],tpn,rel_qpts[0]))
         fill_between(d[:,0],d[:,1],d[:,2],facecolor=colors[i],alpha=0.1)
      i=i+1
   i=0
   for np in np_set:
      rel_qpts=[int(run['quadrature-pts']**(1./3)+0.5)-run['order']
                for run in pl_runs if run['num-procs']==np]
      rel_qpts=sorted(list(set(rel_qpts)))
      rel_qpts.reverse()
      print 'q-p =', rel_qpts
      if len(rel_qpts)==1:
         i=i+1
         continue

      d=[[run['order'],run['cg-iteration-dps']/num_nodes,run['compiler']]
         for run in pl_runs
         if run['num-procs']==np and
            (run['quadrature-pts']==(run['order']+rel_qpts[1])**3)]
      d=[[order,
          max([e[1] for e in d if e[0]==order and e[2]==compilers[0]]),
          max([e[1] for e in d if e[0]==order and e[2]==compilers[-1]])]
         for order in set([e[0] for e in d])]
      d=asarray(sorted(d))
      tpn=np/num_nodes
      plot(d[:,0],d[:,1],'s--',color=colors[i],
           label='%s, %i t/n, q=p+%i'%(compilers[0],tpn,rel_qpts[1]))
      if len(compilers)==2:
         plot(d[:,0],d[:,2],'^--',color=colors[i],
              label='%s, %i t/n, q=p+%i'%(compilers[-1],tpn,rel_qpts[1]))
         fill_between(d[:,0],d[:,1],d[:,2],facecolor=colors[i],alpha=0.1)
      i=i+1

   title('Config: %s (%i node%s), %s, %s'%(
         config_short,num_nodes,'' if num_nodes==1 else 's',test_short,
         'PA' if action_type=='matrix-free' else 'TA'))
   yscale('log')
   # rng=arange(1e7,1.02e8,1e7)
   # yticks(rng,['%i'%int(v/1e6) for v in rng])
   # ylim(min(rng),max(rng))
   xlim(0.5,max([run['order'] for run in pl_runs])+0.5)
   grid('on', color='gray', ls='dotted')
   grid('on', axis='y', which='minor', color='gray', ls='dotted')
   gca().set_axisbelow(True)
   xlabel('Order, p')
   ylabel('[DOFs x CG iterations] / [compute nodes x seconds]')
   legend(ncol=2, loc='best')

   if 1: # write .pdf file?
      savefig('plot1_%s_%s_%s_N%03i.pdf'%(
              code,test_short,config_short,num_nodes),
              format='pdf', bbox_inches='tight')

if 1: # show the figures?
   show()
