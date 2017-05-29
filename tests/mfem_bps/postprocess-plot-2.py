
#####   Load the data
execfile('postprocess-base.py')


#####   Sample plot output


from pylab import *

rcParams['font.sans-serif'].insert(0,'Noto Sans')
rcParams['font.sans-serif'].insert(1,'Open Sans')
rcParams['figure.figsize']=[10, 8] # default: 8 x 6

cm=get_cmap('Set1') # 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3'
cm_size=len(cm.__dict__['_segmentdata']['red'])
colors=[cm(1.*i/(cm_size-1)) for i in range(cm_size)]

# colors=['blue','green','crimson','turquoise','m','gold','cornflowerblue',
#         'darkorange']

sel_runs=runs
action_type='matrix-free'
sel_runs=[run for run in sel_runs if run['action-type']==action_type]

configs=list(set([run['config'] for run in sel_runs]))
# print 'Present configurations:', configs
config=configs[0]
print 'Using configuration:', config
config_short=config.rsplit('/',1)[-1].rsplit('.sh',1)[0]
sel_runs=[run for run in sel_runs if run['config']==config]

tests=list(set([run['test'] for run in sel_runs]))
# print 'Present tests:', tests
test=tests[0]
print 'Using test:', test
test_short=test.rsplit('/',1)[-1].rsplit('.sh',1)[0]
sel_runs=[run for run in sel_runs if run['test']==test]

pl_set=[(run['compiler'],run['num-procs'],run['num-procs-node'])
        for run in sel_runs]
pl_set=sorted(set(pl_set))
print
pprint.pprint(pl_set)

for plt in pl_set:
   compiler=plt[0]
   num_procs=plt[1]
   num_procs_node=plt[2]
   num_nodes=num_procs/num_procs_node
   pl_runs=[run for run in sel_runs
            if run['compiler']==compiler and
               run['num-procs']==num_procs and
               run['num-procs-node']==num_procs_node]
   pl_runs=sorted(pl_runs)
   if len(pl_runs)==0:
      continue

   print
   print 'compiler: %s, compute nodes: %i, number of MPI tasks = %i'%(
      compiler,num_nodes,num_procs)

   figure()
   i=0
   sol_p_set=sorted(set([run['order'] for run in pl_runs]))
   for sol_p in sol_p_set:
      d=[[run['order'],run['num-elem'],1.*run['num-unknowns']/num_nodes,
          run['cg-iteration-dps']/num_nodes]
         for run in pl_runs
         if run['order']==sol_p and
            run['quadrature-pts']==(run['order']+2)**3]
      # print
      # print 'order = %i'%sol_p
      # pprint.pprint(sorted(d))
      d=[[e[2],e[3]] for e in d if e[0]==sol_p]
      # (DOFs/[sec/iter]/node)/(DOFs/node) = iter/sec
      d=[[nun,
          min([e[1] for e in d if e[0]==nun]),
          max([e[1] for e in d if e[0]==nun])]
         for nun in set([e[0] for e in d])]
      d=asarray(sorted(d))
      plot(d[:,0],d[:,2],'o-',color=colors[i],label='p=%i, q=p+2'%sol_p)
      if list(d[:,1]) != list(d[:,2]):
         plot(d[:,0],d[:,1],'o-',color=colors[i])
         fill_between(d[:,0],d[:,1],d[:,2],facecolor=colors[i],alpha=0.1)
      ##
      d=[[run['order'],run['num-elem'],1.*run['num-unknowns']/num_nodes,
          run['cg-iteration-dps']/num_nodes]
         for run in pl_runs
         if run['order']==sol_p and
            (run['quadrature-pts']==(run['order']+1)**3)]
      d=[[e[2],e[3]] for e in d if e[0]==sol_p]
      if len(d)==0:
         i=i+1
         continue
      d=[[nun,
          min([e[1] for e in d if e[0]==nun]),
          max([e[1] for e in d if e[0]==nun])]
         for nun in set([e[0] for e in d])]
      d=asarray(sorted(d))
      plot(d[:,0],d[:,2],'s--',color=colors[i],label='p=%i, q=p+1'%sol_p)
      if list(d[:,1]) != list(d[:,2]):
         plot(d[:,0],d[:,1],'s--',color=colors[i])
      ##
      i=i+1
   ##
   y=asarray([3e5, 7e7])
   slope1=600
   slope2=6000
   plot(y/slope1,y,'k--',label='%g iter/s'%slope1)
   plot(y/slope2,y,'k-',label='%g iter/s'%slope2)

   title('Config: %s (%i node%s, %i tasks/node), %s, %s, PA'%(
         config_short,num_nodes,'' if num_nodes==1 else 's',
         num_procs_node,compiler,test_short))
   xscale('log') # subsx=[2,4,6,8]
   yscale('log')
   # rng=arange(1e7,1.02e8,1e7)
   # yticks(rng,['%i'%int(v/1e6) for v in rng])
   # ylim(min(rng),max(rng))
   # xlim(0.5,max([run['order'] for run in pl_runs])+0.5)
   grid('on', color='gray')
   grid('on', axis='both', which='minor', color='gray')
   gca().set_axisbelow(True)
   xlabel('DOFs per compute node')
   ylabel('[DOFs x CG iterations] / [compute nodes x seconds]')
   legend(ncol=2, loc='best')

   # savefig('test_%s_%s_%s_N%03i_pn%i.pdf'%(
   #         test_short,config_short,compiler,num_nodes,num_procs_node),
   #         format='pdf', bbox_inches='tight')

show()