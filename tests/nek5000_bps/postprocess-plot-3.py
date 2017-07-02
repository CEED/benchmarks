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

configs=list(set([run['config'].rsplit('/',1)[-1].rsplit('.sh',1)[0]
                  for run in sel_runs]))
# print 'Present configurations:', configs
config=configs[0]
print 'Using configuration:', config
config_short=config
sel_runs=[run for run in sel_runs if
          run['config'].rsplit('/',1)[-1].rsplit('.sh',1)[0]==config]

compilers=sorted(list(set([run['compiler'] for run in sel_runs])))
compilers=compilers[0:2]
print 'Using compilers:', compilers
sel_runs=[run for run in sel_runs if run['compiler'] in compilers]

tests=list(set([run['test'].rsplit('/',1)[-1].rsplit('.sh',1)[0]
                for run in sel_runs]))
# print 'Present tests:', tests
test=tests[0]
print 'Using test:', test
test_short=test
sel_runs=[run for run in sel_runs if
          run['test'].rsplit('/',1)[-1].rsplit('.sh',1)[0]==test]

if 'case' in sel_runs[0]:
   cases=list(set([run['case'] for run in sel_runs]))
   case=cases[0]
   vdim=1 if case=='scalar' else 3
   print 'Using case:', case
   sel_runs=[run for run in sel_runs if run['case']==case]

codes = list(set([run['code'] for run in sel_runs]))
code  = codes[0]
sel_runs=[run for run in sel_runs if run['code']==code]

pl_set=[(run['order'],
         run['quadrature-pts'],
         run['num-procs']/run['num-procs-node'])
        for run in sel_runs]
pl_set=sorted(set(pl_set))
print
pprint.pprint(pl_set)

for plt in pl_set:
   sol_p=plt[0]
   qpts=plt[1]
   num_nodes=plt[2]
   pl_runs=[run for run in sel_runs
            if run['order']==sol_p and
               run['quadrature-pts']==qpts and
               run['num-procs']==num_nodes*run['num-procs-node']]
   pl_runs=sorted(pl_runs)
   if len(pl_runs)==0:
      continue

   print
   print 'order: %s, qpts: %i, compute nodes: %i'%(sol_p,qpts,num_nodes)

   figure()
   i=0
   compilers_fig=list(set([run['compiler'] for run in pl_runs]))
   tpn_set=sorted(set([run['num-procs-node'] for run in pl_runs]))
   for tpn in tpn_set:
      d=[[1.*run['num-unknowns']/num_nodes/vdim,
          run['cg-iteration-dps']/num_nodes,
          run['compiler']]
         for run in pl_runs if run['num-procs-node']==tpn]
      if len(d)==0:
         i=i+1
         continue
      compilers_d=sorted(list(set([e[2] for e in d])))
      if len(compilers_d)==1:
         d=[[dpn,
             max([e[1] for e in d if e[0]==dpn and e[2]==compilers_d[0]])]
            for dpn in set([e[0] for e in d])]
      else:
         d=[[dpn,
             max([e[1] for e in d if e[0]==dpn and e[2]==compilers_d[0]]),
             max([e[1] for e in d if e[0]==dpn and e[2]==compilers_d[1]])]
            for dpn in set([e[0] for e in d])]
      d=asarray(sorted(d))
      plot(d[:,0],d[:,1],'o-',color=colors[i],
           label='%s, %i t/n'%(compilers_d[0],tpn))
      if len(compilers_d)==2:
         plot(d[:,0],d[:,2],'d-',color=colors[i],
              label='%s, %i t/n'%(compilers_d[1],tpn))
         fill_between(d[:,0],d[:,1],d[:,2],facecolor=colors[i],alpha=0.2)
      ##
      i=i+1
   ##
   y=asarray([3e5, 7e7])
   slope1=600
   slope2=6000
   plot(y/slope1,y,'k--',label='%g iter/s'%slope1)
   plot(y/slope2,y,'k-',label='%g iter/s'%slope2)

   title('Config: %s (%i node%s), p=%i, q=%i, %s, %s'%(
         config_short,num_nodes,'' if num_nodes==1 else 's', sol_p,
         int(qpts**(1./3)+0.5),test_short,
         'PA' if action_type=='matrix-free' else 'TA'))
   xscale('log') # subsx=[2,4,6,8]
   yscale('log')
   # rng=arange(1e7,1.02e8,1e7)
   # yticks(rng,['%i'%int(v/1e6) for v in rng])
   # ylim(min(rng),max(rng))
   # xlim(0.5,max([run['order'] for run in pl_runs])+0.5)
   grid('on', color='gray', ls='dotted')
   grid('on', axis='both', which='minor', color='gray', ls='dotted')
   gca().set_axisbelow(True)
   xlabel('Points per compute node')
   ylabel('[DOFs x CG iterations] / [compute nodes x seconds]')
   legend(ncol=2, loc='best')

   if 1: # write .pdf file?
      if len(compilers_fig)==1:
         pdf_file='plot3_%s_%s_%s_%s_N%03i_p%i_q%i.pdf'%(code,test_short,
                  config_short,compilers[0],
                  num_nodes,sol_p,int(qpts**(1./3)+0.5))
      else:
         pdf_file='plot3_%s_%s_%s_N%03i_p%i_q%i.pdf'%(code,test_short,
                  config_short,num_nodes,sol_p,int(qpts**(1./3)+0.5))
      print 'saving figure --> %s'%pdf_file
      savefig(pdf_file, format='pdf', bbox_inches='tight')

if 1: # show the figures?
   print '\nshowing figures ...'
   show()
