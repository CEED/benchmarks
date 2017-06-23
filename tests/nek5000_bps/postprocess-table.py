
#####   Load the data
execfile('postprocess-base.py')


#####   Sample data output

hpc={'s': 'std' , 'h': 'hpc'}
mf={'a': 'A', 'm': 'PA'}
set1=sorted(
   [(run['order'],run['quadrature-pts'],mf[run['action-type'][0]],
     hpc[run['std | hpc'][0]],run['compiler'],run['num-procs'],
     run['assembly-dps']/1e6,run['cg-iteration-dps']/1e6,
     1.*run['num-elem']*run['quadrature-pts']/run['num-unknowns'])
    for run in runs])

out.write('''\
    | quad | asse | std |  comp  |     | assembly dps | cg-iter dps | qpts per
  p |  pts | mbly | hpc |  iler  |  np |   millions   |   millions  |  unknown
----+------+------+-----+--------+-----+--------------+-------------+---------
''')
line_fmt=' %2i | %4i |  %2s  | %s | %6s | %3i | %11.6f  | %11.6f | %7.4f\n'
for run in set1:
   out.write(line_fmt%run)
