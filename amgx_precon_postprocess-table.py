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
execfile('amgx_precon_postprocess-base.py')

def safe_div(x,y):
    if y == 0:
        return 0
    return x / y

#####   Sample data output

set1=sorted(
   [(run['mfem-device'],
     run['order'],run['compiler'],run['num-procs'],
     run['num-unknowns'], safe_div(run['num-unknowns'],run['cost-of-precon'])/1e6 , run['iterations'],
     safe_div(run['num-unknowns'],run['time-per-cg-step'])/1e6, safe_div(run['num-unknowns'],
     run['total-cg-time'])/1e6, safe_div(run['num-unknowns'],run['amg-setup'])/1e6
     ,run['cg-iteration-dps']/1e6)
    for run in runs])

out.write('''\
    mfem   |    |  comp  |     |number of | precon dps    | number of  | cg-1-it dps |  cg-a-it dps | amg-set dps  |cg-iter dps
   device  |  p |  iler  |  np |unknowns  |  millions     | iterations | millions    |  millions    |   millions   |millions
-----------+----+--------+-----+---------+-------------+-------------+-------------+--------------+--------------+-------------
''')
line_fmt=' %9s | %2i | %6s | %3i | %6i  | %11.6f | %7i | %11.8f | %11.6f |  %11.6f |%11.6f\n'
for run in set1:
   out.write(line_fmt%run)
