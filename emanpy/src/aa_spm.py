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
    Creates a surface permanent magnet machines (SPM)
    and computes its performance analytically.
"""

# ==========================================================================
## Program:   aa_spm.py
## Author:    ajpina
## Date:      12/23/15
## Version:   0.1.1
##
## Revision History:
##      Date     Version  Author    Description
##  - 12/23/15:  0.1.1              Calculate Flux Density, Cogging and Ripple
##
# ==========================================================================

import getopt
import json
import sys
sys.path.append('../..')
import time
import logging

from emanpys.magnets import *
from emanpys.materials import *
from emanpys.rotors import *
from emanpys.slots import *
from emanpys.stators import *
from emanpys.results import *
from emanpys.windings import *

from emanpys import spm as machine
from emanpys.src.constants import *


#from database import db_connector as db


class Usage(Exception):
    def __init__(self, msg):
        self.msg = "[Error]: %s" % ( msg )


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hd:f:r:l:", ["help","dir","file","result-id","log"])
        except getopt.error, msg:
             raise Usage(msg)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print 'aa_spm.py -d [dir_name] -f [file_name] -r [id] -l [level]'
                sys.exit()
            elif opt in ("-d", "--dir"):
                dir = arg
            elif opt in ("-f", "--file"):
                file = arg
            elif opt in ("-r", "--result-id"):
                rid = arg
            elif opt in ("-l", "--log"):
                loglevel = int(arg)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

    filename = "%s/%s" % (dir, file)

    start1 = time.clock()
    with open(filename) as heman_file:
        heman_settings = json.load(heman_file)

    logfile = "%s/%s_%s_%s.log" % (dir, heman_settings['UID'],
                                        heman_settings['SID'],
                                        heman_settings['AID'])

    if loglevel >= LOG_ALL:
        loglevel = LOG_ALL
        logging.basicConfig(filename=logfile, level=logging.DEBUG,
                            format='%(asctime)s - [%(name)s] %(levelname)s: %(message)s')
    elif loglevel == LOG_INFO:
        logging.basicConfig(filename=logfile, level=logging.INFO,
                            format='%(asctime)s - [%(name)s] %(levelname)s: %(message)s')
    elif loglevel == LOG_WARN:
        logging.basicConfig(filename=logfile, level=logging.WARNING,
                            format='%(asctime)s - [%(name)s] %(levelname)s: %(message)s')
    elif loglevel == LOG_ERROR:
        logging.basicConfig(filename=logfile, level=logging.ERROR,
                            format='%(asctime)s - [%(name)s] %(levelname)s: %(message)s')
    elif loglevel == LOG_CRITICAL:
        logging.basicConfig(filename=logfile, level=logging.CRITICAL,
                            format='%(asctime)s - [%(name)s] %(levelname)s: %(message)s')


    statorObj = heman_settings['machine']['stator']
    rotorObj = heman_settings['machine']['rotor']
    analysisObj = heman_settings['analysis']
    slts = []
    for i,slot in enumerate(statorObj['slots']):
        slts.insert(i, Type0(h0=slot['h0'], h1=slot['h1'], h2=slot['h2'], h3=slot['h3'],
                              w0=slot['w0'], w1=slot['w1'], w2=slot['w2'],
                              SOp=slot['SOpos'], Sp=slot['Spos']))

    pms = []
    for i,magnet in enumerate(rotorObj['magnets']['dimension']):
        pms.insert(i, Magnet(Ml=magnet['Ml'], Mw=magnet['Mw'], Mps=magnet['Mps'], iMr=magnet['iMr'],
                                P=2*(rotorObj['pp']), magType=rotorObj['magnets']['magType'],
                                delta=magnet['delta']))

    connMatrix = np.zeros((statorObj['winding']['NoPhases']*statorObj['winding']['Layers'],
                           statorObj['Ns']), dtype=int)
    layers = ['CM1', 'CM2', 'CM3', 'CM4']
    phases = ['A', 'B', 'C', 'D', 'E', 'F']
    row = 0
    for la in range(0,statorObj['winding']['Layers']):
        for ph in range(0,statorObj['winding']['NoPhases']):
            for ns in range(0,statorObj['Ns']):
                connMatrix[row,ns] = statorObj['winding']['CM'][layers[la]][ns][phases[ph]]
            row +=1

    wm = Material(rho=statorObj['winding']['material']['resistivity'])
    wd = Concentrated(phases=statorObj['winding']['NoPhases'],
                      layers=statorObj['winding']['Layers'], conn=statorObj['winding']['Type'],
                      C=connMatrix, mat=wm, series=statorObj['winding']['Cseries'],
                      parallel=statorObj['winding']['Cparallel'], coilTurns=statorObj['winding']['Cturns'],
                      wih=statorObj['winding']['wih'], condDiam=statorObj['winding']['condDiam'],
                      coilPitch=statorObj['winding']['Cpitch'], Ns=statorObj['Ns'], slot=slts[0],
                      iSr=statorObj['iSr'], Sl=statorObj['Sl'])

    st = StandardOuterStator(Ns=statorObj['Ns'], slots=slts, iSr=statorObj['iSr'],
                             oSr=statorObj['oSr'], Sl=statorObj['Sl'], winding=wd)

    rt = InnerRotorSPM(pp=rotorObj['pp'], magnets=pms, iRr=rotorObj['iRr'], oRr=rotorObj['oRr'],
                       Rl=rotorObj['Rl'], init_pos=rotorObj['init_pos'])

    analysis = Study(analysisObj)

    a = machine.InnerRotor(st, rt)
    start2 = time.clock()
    res = a.solve(analysis)
    finish = time.clock()
    log_msg = "[AA_SPM] Solving time is %fsec" % (finish - start2)
    logging.info(log_msg)
    log_msg = "[AA_SPM] Total time is %fsec" % (finish - start1)
    logging.info(log_msg)


    #db_config = db.read_db_config(filename='./database/heman_config.ini')

    sql1 = """
        UPDATE Results SET Resistance_Stator=%s, Ke=%s, status=%s
        WHERE id=%s
    """
    value1 = (res.stator_phase_resistance,
              0.003,
              1,
              rid)



    #db.update_register(db_config, sql1, value1)

    sql2 = """
        INSERT INTO Points_Graph (Results, Data_Set, number, x, y)
        VALUES (%s, %s, %s, %s, %s)
    """

    value2 = []

    zipped_cogging = zip(res.cogging_torque_x, res.cogging_torque_y)
    i=1
    for x,y in zipped_cogging:
        value2.append((rid,1,i,"{:.6f}".format(x),"{:.6f}".format(y)))
        i+=1

    zipped_ripple = zip(res.torque_ripple_x, res.torque_ripple_y)
    i=1
    for x,y in zipped_ripple:
        value2.append((rid,2,i,"{:.6f}".format(x),"{:.6f}".format(y)))
        i+=1

    zipped_static_torque = zip(res.static_torque_x, res.static_torque_y)
    i=1
    for x,y in zipped_static_torque:
        value2.append((rid,3,i,"{:.6f}".format(x),"{:.6f}".format(y)))
        i+=1

    zipped_nl_Bg_r = zip(res.nl_Bg_theta, res.nl_Bg_r)
    i=1
    for x,y in zipped_nl_Bg_r:
        value2.append((rid,4,i,"{:.6f}".format(x),"{:.6f}".format(y)))
        i+=1

    zipped_nl_Bg_t = zip(res.nl_Bg_theta, res.nl_Bg_t)
    i=1
    for x,y in zipped_nl_Bg_t:
        value2.append((rid,4,i,"{:.6f}".format(x),"{:.6f}".format(y)))
        i+=1


    #db.multiple_insert(db_config, sql2, value2)

    sql3 = """
        UPDATE Results SET status=%s
        WHERE id=%s
    """
    value3 = (2,rid)

    #db.update_register(db_config, sql3, value3)

    log_msg = "[AA_SPM] Magnet Flux = %f" % (res.magnet_flux)
    logging.info(log_msg)
    log_msg = "[AA_SPM] Lmq = %f" % (res.Lmq)
    logging.info(log_msg)
    log_msg = "[AA_SPM] Self Inductance (AirGap) = %f" % (res.self_inductance_ag)
    logging.info(log_msg)
    log_msg = "[AA_SPM] Mutual Inductance (AirGap) = %f" % (res.mutual_inductance_ag)
    logging.info(log_msg)
    log_msg = "[AA_SPM] End winding Leakage = %f" % (res.end_winding_leakage)
    logging.info(log_msg)

    if loglevel == LOG_ALL:

        import matplotlib.pyplot as plt
        plt.figure(1)
        # plt.ion()
        plt.title('Air Gap Flux Density')
        plt.subplot(211)
        plt.plot(res.nl_Bg_theta, res.nl_Bg_r, label='No Load')
        plt.plot(res.nl_Bg_theta, res.ol_Bg_r, label='On Load')
        plt.legend()
        plt.subplot(212)
        plt.plot(res.nl_Bg_theta, res.nl_Bg_t, label='No Load')
        plt.plot(res.nl_Bg_theta, res.ol_Bg_t, label='On Load')


        plt.figure(2)
        plt.plot(res.cogging_torque_x, res.cogging_torque_y, '-ro')
        plt.title('Cogging Torque')

        plt.figure(3)
        plt.plot(res.torque_ripple_x, res.torque_ripple_y, '-ro')
        plt.title('Torque Ripple')

        plt.figure(4)
        plt.plot(res.static_torque_x, res.static_torque_y, '-ro')
        plt.title('Static Torque')

        plt.figure(5)
        plt.title('Winding Function - Turns Density')
        plt.subplot(311)
        plt.plot(res.wf[0,:], label='Winding Function')
        plt.plot(res.td[0, :], label='Turns Density')
        plt.legend()
        plt.subplot(312)
        plt.plot(res.wf[1,:])
        plt.plot(res.td[1, :])
        plt.subplot(313)
        plt.plot(res.wf[2, :])
        plt.plot(res.td[2, :])


        plt.figure(6)
        plt.title('No Load Flux Linkage')
        plt.plot(res.nl_flux_linkage_x, res.nl_flux_linkage_y[0, :], '-o', color='blue', label='a')
        plt.plot(res.nl_flux_linkage_x, res.nl_flux_linkage_y[1, :], '-o', color='green', label='b')
        plt.plot(res.nl_flux_linkage_x, res.nl_flux_linkage_y[2, :], '-o', color='orange', label='c')
        plt.legend()

        plt.figure(7)
        plt.title('On Load Flux Linkage')
        plt.plot(res.ol_flux_linkage_x, res.ol_flux_linkage_y[0, :]-res.nl_flux_linkage_y[0, :], '-o', color='blue', label='a')
        plt.plot(res.ol_flux_linkage_x, res.ol_flux_linkage_y[1, :]-res.nl_flux_linkage_y[1, :], '-o', color='green', label='b')
        plt.plot(res.ol_flux_linkage_x, res.ol_flux_linkage_y[2, :]-res.nl_flux_linkage_y[2, :], '-o', color='orange', label='c')
        plt.legend()

        plt.figure(8)
        x=np.linspace(0,50,51)

        bar1 = res.wh
        bar2 = res.kw_v
        plt.subplot(211)
        plt.title('Winding Harmonics')
        plt.bar(x,bar1[0:51],width=0.5)
        #plt.bar(x+0.25, bar2[0:51], width=0.5)
        plt.subplot(212)
        plt.title('Winding Factors')
        plt.bar(x+1, bar2[0:51])

        plt.figure(9)
        plt.title('Phase Current')
        plt.plot(res.phase_current_x, res.phase_current_y[0, :], '-o', color='blue', label='a')
        plt.plot(res.phase_current_x, res.phase_current_y[1, :], '-o', color='green', label='b')
        plt.plot(res.phase_current_x, res.phase_current_y[2, :], '-o', color = 'orange', label='c')
        plt.legend()

        plt.figure(10)
        plt.title('Back EMF')
        plt.plot(res.bemf_x, res.bemf_y[0, :], '-o', color='blue', label='a')
        plt.plot(res.bemf_x, res.bemf_y[1, :], '-o', color='green', label='b')
        plt.plot(res.bemf_x, res.bemf_y[2, :], '-o', color='orange', label='c')
        plt.legend()
        #np.savetxt('bemf.csv', res.bemf_y.transpose(), delimiter=',')
        #np.savetxt('bemf2.csv', res.bemf_x.transpose(), delimiter=',')

        plt.figure(11)
        plt.title('Radial Pressure Harmonics')
        x = np.arange(0, 5 * GCD(len(slts),len(pms)) + 1, 1)
        plt.bar(x - 0.2, res.pressure_radial_nl[x], width=0.4, label='No Load')
        plt.bar(x + 0.2, res.pressure_radial_ol[x], width=0.4, label='On Load')
        plt.legend()


        plt.show()

    logging.shutdown()
    return res



if __name__ == '__main__':
    sys.exit(main())
