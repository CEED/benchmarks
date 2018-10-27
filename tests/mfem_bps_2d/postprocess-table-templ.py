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


#####   Load the data
execfile('postprocess-base-templ.py')


#####   Sample data output

hpc={'s': 'std' , 'h': 'hpc'}
mf={'a': 'A', 'm': 'PA'}
set1=sorted(
   [(run['order'],run['quadrature-pts'],mf[run['action-type'][0]],
     hpc[run['std | hpc'][0]],run['compiler'],run['num-procs'],
     run['num-unknowns'],
     run['assembly-dps']/1e6,run['cg-iteration-dps']/1e6,
     1.*run['num-elem']*run['quadrature-pts']/run['num-unknowns'])
    for run in runs])

out.write('''\
    | quad | asse | std |  comp  |     |    num   | assembly dps | cg-iter dps | qpts per
  p |  pts | mbly | hpc |  iler  |  np |   dofs   |   millions   |   millions  |  unknown
----+------+------+-----+--------+-----+----------+--------------+-------------+---------
''')
line_fmt=' %2i | %4i |  %2s  | %s | %6s | %3i | %8i | %11.6f  | %11.6f | %7.4f\n'
for run in set1:
   out.write(line_fmt%run)
