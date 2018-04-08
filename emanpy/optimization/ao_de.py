#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

# ==========================================================================
## Copyright (C) 2016 Dr. Alejandro Pina Ortega
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##      http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
# ==========================================================================

"""
    Starts the optimization engine based upon differential evolution strategy
    to search for the optimum machine geometry that meets the given goals and
    constraints evaluated through analytical methods.
"""

# ==========================================================================
## Program:   ao_de.py
## Author:    ajpina
## Date:      5/1/17
## Version:   0.1.1
##
## Revision History:
##      Date     Version  Description
##  - 5/1/17:  0.1.1    Call routine aa_spm.py
##
# ==========================================================================

import sys
sys.path.append('../..')
import getopt
import json
import time
import numpy as np
import csv
from scipy.optimize import differential_evolution

import emanpy.src.aa_spm as spm
from emanpy.results.results import Result



class Usage(Exception):
    def __init__(self, msg):
        self.msg = "[Error]: %s" % ( msg )

def is_invalid_geometry(oRr, Mw, delta):
    for i in range(-1, 9, 1):
        current_magnet_edge_point = (i+1)*36*np.pi/180 + (0.5*(Mw[i]/oRr) + delta[i]*np.pi/180)
        next_magnet_edge_point =  (i+2)*36*np.pi/180 - (0.5*(Mw[i+1]/oRr) + delta[i+1]*np.pi/180)
        if (next_magnet_edge_point-current_magnet_edge_point) < 2e-3:
            return True
    return False


def cost_function(x, dir, initial_candidate ):
    w0_1 = x[0]
    w0_2 = x[1]
    w0_3 = x[2]
    w0_4 = x[3]
    w0_5 = x[4]
    w0_6 = x[5]
    w0_7 = x[6]
    w0_8 = x[7]
    w0_9 = x[8]
    w0_10 = x[9]
    w0_11 = x[10]
    w0_12 = x[11]
    SOpos_1 = x[12]
    SOpos_2 = x[13]
    SOpos_3 = x[14]
    SOpos_4 = x[15]
    SOpos_5 = x[16]
    SOpos_6 = x[17]
    SOpos_7 = x[18]
    SOpos_8 = x[19]
    SOpos_9 = x[20]
    SOpos_10 = x[21]
    SOpos_11 = x[22]
    SOpos_12 = x[23]
    Ml = x[24]
    Mw_1 = x[25]
    Mw_2 = x[26]
    Mw_3 = x[27]
    Mw_4 = x[28]
    Mw_5 = x[29]
    Mw_6 = x[30]
    Mw_7 = x[31]
    Mw_8 = x[32]
    Mw_9 = x[33]
    Mw_10 = x[34]
    delta_1 = x[35]
    delta_2 = x[36]
    delta_3 = x[37]
    delta_4 = x[38]
    delta_5 = x[39]
    delta_6 = x[40]
    delta_7 = x[41]
    delta_8 = x[42]
    delta_9 = x[43]
    delta_10 = x[44]
    current = x[45]
    iSr = 22.95e-3
    ag = 0.75e-3
    oRr = iSr - ag - Ml

    if is_invalid_geometry(oRr=oRr, Mw=[Mw_1, Mw_2, Mw_3, Mw_4, Mw_5, Mw_6, Mw_7, Mw_8, Mw_9, Mw_10],
                           delta=[delta_1, delta_2, delta_3, delta_4, delta_5, delta_6, delta_7, delta_8, delta_9, delta_10]):
        return 1e3

    geometry = initial_candidate['machine']
    geometry['stator']['slots'][0]['w0'] = w0_1
    geometry['stator']['slots'][1]['w0'] = w0_2
    geometry['stator']['slots'][2]['w0'] = w0_3
    geometry['stator']['slots'][3]['w0'] = w0_4
    geometry['stator']['slots'][4]['w0'] = w0_5
    geometry['stator']['slots'][5]['w0'] = w0_6
    geometry['stator']['slots'][6]['w0'] = w0_7
    geometry['stator']['slots'][7]['w0'] = w0_8
    geometry['stator']['slots'][8]['w0'] = w0_9
    geometry['stator']['slots'][9]['w0'] = w0_10
    geometry['stator']['slots'][10]['w0'] = w0_11
    geometry['stator']['slots'][11]['w0'] = w0_12
    geometry['stator']['slots'][0]['SOpos'] = SOpos_1
    geometry['stator']['slots'][1]['SOpos'] = SOpos_2
    geometry['stator']['slots'][2]['SOpos'] = SOpos_3
    geometry['stator']['slots'][3]['SOpos'] = SOpos_4
    geometry['stator']['slots'][4]['SOpos'] = SOpos_5
    geometry['stator']['slots'][5]['SOpos'] = SOpos_6
    geometry['stator']['slots'][6]['SOpos'] = SOpos_7
    geometry['stator']['slots'][7]['SOpos'] = SOpos_8
    geometry['stator']['slots'][8]['SOpos'] = SOpos_9
    geometry['stator']['slots'][9]['SOpos'] = SOpos_10
    geometry['stator']['slots'][10]['SOpos'] = SOpos_11
    geometry['stator']['slots'][11]['SOpos'] = SOpos_12
    geometry['rotor']['magnets']['dimension'][0]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][1]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][2]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][3]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][4]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][5]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][6]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][7]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][8]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][9]['Ml'] = Ml
    geometry['rotor']['magnets']['dimension'][0]['Mw'] = Mw_1
    geometry['rotor']['magnets']['dimension'][1]['Mw'] = Mw_2
    geometry['rotor']['magnets']['dimension'][2]['Mw'] = Mw_3
    geometry['rotor']['magnets']['dimension'][3]['Mw'] = Mw_4
    geometry['rotor']['magnets']['dimension'][4]['Mw'] = Mw_5
    geometry['rotor']['magnets']['dimension'][5]['Mw'] = Mw_6
    geometry['rotor']['magnets']['dimension'][6]['Mw'] = Mw_7
    geometry['rotor']['magnets']['dimension'][7]['Mw'] = Mw_8
    geometry['rotor']['magnets']['dimension'][8]['Mw'] = Mw_9
    geometry['rotor']['magnets']['dimension'][9]['Mw'] = Mw_10
    geometry['rotor']['magnets']['dimension'][0]['delta'] = delta_1
    geometry['rotor']['magnets']['dimension'][1]['delta'] = delta_2
    geometry['rotor']['magnets']['dimension'][2]['delta'] = delta_3
    geometry['rotor']['magnets']['dimension'][3]['delta'] = delta_4
    geometry['rotor']['magnets']['dimension'][4]['delta'] = delta_5
    geometry['rotor']['magnets']['dimension'][5]['delta'] = delta_6
    geometry['rotor']['magnets']['dimension'][6]['delta'] = delta_7
    geometry['rotor']['magnets']['dimension'][7]['delta'] = delta_8
    geometry['rotor']['magnets']['dimension'][8]['delta'] = delta_9
    geometry['rotor']['magnets']['dimension'][9]['delta'] = delta_10
    geometry['rotor']['magnets']['dimension'][0]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][1]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][2]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][3]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][4]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][5]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][6]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][7]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][8]['iMr'] = oRr
    geometry['rotor']['magnets']['dimension'][9]['iMr'] = oRr
    geometry['rotor']['oRr'] = oRr
    initial_candidate['analysis']['ripple']['current'] = current

    tmp_filename = "%s/candidate_12s10p.haf" % (dir)

    with open(tmp_filename, 'w') as heman_outfile:
        json.dump(initial_candidate, heman_outfile)

    results = spm.main(["aa_spm.py", "-d", dir, "-f", "candidate_12s10p.haf", "-r", "1", "-l", "0"])

    p1 = results.pressure_radial_ol[1]
    p2 = results.pressure_radial_ol[2]
    cost = p1*100 + p2*80
    t_avg = np.average(results.torque_ripple_y)
    if t_avg < 6.5:
        cost = cost + 500
    else:
        cost = cost + 1.0/t_avg

    results_filename = "%s/optimization_results_12s10p.csv" % (dir)

    with open(results_filename,'ab') as f:
        fWriter = csv.writer(f,delimiter=',',quoting=csv.QUOTE_NONE)
        data = [w0_1,w0_2,w0_3,w0_4,w0_5,w0_6,w0_7,w0_8,w0_9,w0_10,w0_11,w0_12,
                SOpos_1,SOpos_2,SOpos_3,SOpos_4,SOpos_5,SOpos_6,SOpos_7,SOpos_8,SOpos_9,SOpos_10,SOpos_11,SOpos_12,
                Ml,Mw_1,Mw_2,Mw_3,Mw_4,Mw_5,Mw_6,Mw_7,Mw_8,Mw_9,Mw_10,
                delta_1,delta_2,delta_3,delta_4,delta_5,delta_6,delta_7,delta_8,delta_9,delta_10,current,oRr,
                p1,p2,t_avg,cost]
        fWriter.writerow(data)

    return cost


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hd:f:o:", ["help","dir","file-machine","optim-file"])
        except getopt.error, msg:
             raise Usage(msg)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print 'ao_de.py -d [dir_name] -f [machine_file_name] -o [optim_file_name]'
                sys.exit()
            elif opt in ("-d", "--dir"):
                dir = arg
            elif opt in ("-f", "--file-machine"):
                afile = arg
            elif opt in ("-o", "--optim-file"):
                ofile = arg

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

    afilename = "%s/%s" % (dir, afile)
    ofilename = "%s/%s" % (dir, ofile)


    start1 = time.clock()

    with open(afilename) as heman_infile:
        candidate = json.load(heman_infile)

    b = [ (0.1e-3, 3e-3),   # Width Slot Opening 1
          (0.1e-3, 3e-3),   # Width Slot Opening 2
          (0.1e-3, 3e-3),   # Width Slot Opening 3
          (0.1e-3, 3e-3),   # Width Slot Opening 4
          (0.1e-3, 3e-3),   # Width Slot Opening 5
          (0.1e-3, 3e-3),   # Width Slot Opening 6
          (0.1e-3, 3e-3),   # Width Slot Opening 7
          (0.1e-3, 3e-3),   # Width Slot Opening 8
          (0.1e-3, 3e-3),   # Width Slot Opening 9
          (0.1e-3, 3e-3),   # Width Slot Opening 10
          (0.1e-3, 3e-3),   # Width Slot Opening 11
          (0.1e-3, 3e-3),   # Width Slot Opening 12
          (4, 26),          # Pos Slot Opening 1
          (34, 56),         # Pos Slot Opening 2
          (64, 86),         # Pos Slot Opening 3
          (94, 116),        # Pos Slot Opening 4
          (124, 146),       # Pos Slot Opening 5
          (154, 176),       # Pos Slot Opening 6
          (184, 206),       # Pos Slot Opening 7
          (214, 236),       # Pos Slot Opening 8
          (244, 266),       # Pos Slot Opening 9
          (274, 296),       # Pos Slot Opening 10
          (304, 326),       # Pos Slot Opening 11
          (334, 356),       # Pos Slot Opening 12
          (2.7e-3, 3.5e-3), # Magnet Length
          (10e-3, 12e-3),   # Width Magnet 1
          (10e-3, 12e-3),   # Width Magnet 2
          (10e-3, 12e-3),   # Width Magnet 3
          (10e-3, 12e-3),   # Width Magnet 4
          (10e-3, 12e-3),   # Width Magnet 5
          (10e-3, 12e-3),   # Width Magnet 6
          (10e-3, 12e-3),   # Width Magnet 7
          (10e-3, 12e-3),   # Width Magnet 8
          (10e-3, 12e-3),   # Width Magnet 9
          (10e-3, 12e-3),   # Width Magnet 10
          (-5.0, 5.0),      # Pos Magnet 1
          (-5.0, 5.0),      # Pos Magnet 2
          (-5.0, 5.0),      # Pos Magnet 3
          (-5.0, 5.0),      # Pos Magnet 4
          (-5.0, 5.0),      # Pos Magnet 5
          (-5.0, 5.0),      # Pos Magnet 6
          (-5.0, 5.0),      # Pos Magnet 7
          (-5.0, 5.0),      # Pos Magnet 8
          (-5.0, 5.0),      # Pos Magnet 9
          (-5.0, 5.0),      # Pos Magnet 10
        (100, 140) ]        # Stator Current

    optim_res = differential_evolution(cost_function, bounds=b, args=(dir, candidate),
                                       strategy='best1bin', maxiter=80000, popsize=1000, mutation=(0.5, 0.6),
                                       recombination=0.9, init='random')

    print optim_res.x, optim_res.fun


    finish = time.clock()
    print 'Total optimization time is ', finish - start1, 's'



if __name__ == '__main__':
    sys.exit(main())