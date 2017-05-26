
#####   Load the data
execfile('postprocess-base.py')


#####   Sample plot output


from pylab import *

rcParams['font.sans-serif'].insert(0,'Noto Sans')
rcParams['font.sans-serif'].insert(1,'Open Sans')
rcParams['figure.figsize']=[10, 8] # default: 8 x 6
rcParams['legend.fontsize']='medium' # default: 'large'

cm=get_cmap('Set1') # 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3'
cm_size=len(cm.__dict__['_segmentdata']['red'])
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

action_type='matrix-free'
sel_runs=[run for run in sel_runs if run['action-type']==action_type]

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
      d=[[run['order'],run['num-elem'],1.*run['num-unknowns']/num_nodes,
          run['cg-iteration-dps']/num_nodes,run['compiler']]
         for run in pl_runs
         if run['num-procs']==np and
            run['quadrature-pts']==(run['order']+2)**3]
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
           label='%s, %i t/n, q=p+2'%(compilers[0],tpn))
      if len(compilers)==2:
         plot(d[:,0],d[:,2],'d-',color=colors[i],
              label='%s, %i t/n, q=p+2'%(compilers[-1],tpn))
         fill_between(d[:,0],d[:,1],d[:,2],facecolor=colors[i],alpha=0.1)
      i=i+1
   i=0
   for np in np_set:
      d=[[run['order'],run['cg-iteration-dps']/num_nodes,run['compiler']]
         for run in pl_runs
         if run['num-procs']==np and
            (run['quadrature-pts']==(run['order']+1)**3)]
      d=[[order,
          max([e[1] for e in d if e[0]==order and e[2]==compilers[0]]),
          max([e[1] for e in d if e[0]==order and e[2]==compilers[-1]])]
         for order in set([e[0] for e in d])]
      d=asarray(sorted(d))
      tpn=np/num_nodes
      plot(d[:,0],d[:,1],'s--',color=colors[i],
           label='%s, %i t/n, q=p+1'%(compilers[0],tpn))
      if len(compilers)==2:
         plot(d[:,0],d[:,2],'^--',color=colors[i],
              label='%s, %i t/n, q=p+1'%(compilers[-1],tpn))
         fill_between(d[:,0],d[:,1],d[:,2],facecolor=colors[i],alpha=0.1)
      i=i+1

   title('%s (%i node%s), Mass, PA'%(
         config_short,num_nodes,'' if num_nodes==1 else 's'))
   yscale('log')
   # rng=arange(1e7,1.02e8,1e7)
   # yticks(rng,['%i'%int(v/1e6) for v in rng])
   # ylim(min(rng),max(rng))
   xlim(0.5,max([run['order'] for run in pl_runs])+0.5)
   grid('on', color='gray')
   grid('on', axis='y', which='minor', color='gray')
   gca().set_axisbelow(True)
   xlabel('Order, p')
   ylabel('[DOFs x CG iterations] / [compute nodes x seconds]')
   legend(ncol=2, loc='best')

   # savefig('test_bp1_%s_N%03i.pdf'%(config_short,num_nodes),
   #         format='pdf', bbox_inches='tight')

show()
