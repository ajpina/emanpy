#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

# ==========================================================================
# Copyright (C) 2016 Dr. Alejandro Pina Ortega
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==========================================================================

"""
    Analysis of permanent magnet synchronous motors with surface-mounted magnets
    in axial flux topology, two outer rotors with coreless stator in between.
"""

# ==========================================================================
# Program:   spm_or_af_subdomain.py
# Author:    ajpina
# Date:      12/23/16
# Version:   0.2
#
# Revision History:
#      Date     Version    Author      Description
#  - 12/23/16:  0.1.1      ajpina      Calculate Air Gap Flux Density
#  - 12/19/16:  0.1.1      ajpina      Calculate Cogging, Ripple
#  - 10/01/17:  0.1.1      ajpina      Calculate Inductances
#  - 05/01/23:  0.2        ajpina      New SD model for axial flux
# ==========================================================================

import scipy.interpolate as si
import scipy.linalg as sl
import logging
from scipy.fftpack import fft
from emanpy.src.constants import *

from emanpy.results import Result
from .magnetization import get_fs_coeff_linear
from emanpy.analysis import Analysis

class SPMOuterRotorAxialFluxSubDomain(Analysis):

    __n_harms = 60
    __pp = 3
    __Ns = 9
    __oRr = 15.2e-3
    __Rm = 18.2e-3
    __Rs = 19.2e-3
    __Rt = 20e-3
    __Rsb = 30e-3
    __Aso = np.ones(9) * 0.1
    __As = np.ones(9) * 0.3
    __alpha_so = np.linspace(0, 360, 9, endpoint=False) * DEG2RAD
    __alpha_s = np.linspace(0, 360, 9, endpoint=False) * DEG2RAD
    __alpha_p_v = np.zeros(6)
    __delta_v = np.zeros(6)
    __M = np.ones(6) * 1.2 / MU0
    __mur = np.ones(6) * 1.1
    __magnetisation = 'Parallel'
    __init_pos = 0
    __C = np.zeros((6, 9))
    __Nph = 3
    __Ct = 20
    __Cp = 3
    __Cs = 1
    __WiH = 1
    __gcd = 3
    __lcm = 18
    __Sl = 30e-3
    __Rl = 30e-3
    __Lx = 20e-3
    __y1 = 2e-3
    __y2 = 3e-3
    __y3 = 4e-3
    __y4 = 5e-3
    __y5 = 6e-3
    __y6 = 7e-3
    __coil_cs = 1e-6
    __r = 20e-3
    __rfactor = 16.0

    def __init__(self, analysis_settings, spm):
        self.type = "Surface Permanent Magnet Machine - Inner Rotor - Subdomain"
        self.__r = spm.slice*(spm.stator.outer_radius - spm.stator.inner_radius) + spm.stator.inner_radius
        self.__pp = spm.rotor.pp
        self.__Ns = spm.stator.slots_number
        self.__oRr = spm.rotor.outer_radius
        self.__Rm = spm.rotor.magnets[0].length + spm.rotor.outer_radius
        self.__init_pos = spm.rotor.rotor_position * DEG2RAD
        self.__Rs = spm.stator.inner_radius
        self.__Aso = np.zeros(4*self.__Ns)
        self.__As = np.zeros(4*self.__Ns)
        self.__alpha_so = np.zeros(4*self.__Ns)
        self.__alpha_s = np.zeros(4*self.__Ns)
        self.__gcd = GCD(self.__Ns, 2 * self.__pp)
        self.__lcm = LCM(self.__Ns, 2 * self.__pp)
        self.__Sl = spm.stator.outer_radius - spm.stator.inner_radius
        self.__Rl = spm.rotor.outer_radius - spm.rotor.inner_radius 
        self.__Lx = 2.0 * np.pi * self.__r
        self.__y1 = spm.rotor.pos_y + spm.rotor.length_y
        self.__y2 = self.__y1 + spm.rotor.magnets[0].length 
        self.__y3 = spm.stator.pos_y
        self.__y4 = self.__y3 + spm.stator.length_y
        self.__y5 = self.__y4 + (self.__y3 - self.__y2)
        self.__y6 = self.__y5 + spm.rotor.magnets[0].length
        self.__results = Result()

        self.__set_slot_geometry(spm)
        self.__set_magnet_geometry(spm)
        self.__set_winding_definition(spm)
        self.__configure(analysis_settings)

        log_msg = "[SPM] Magnet Widths: " + str(self.__alpha_p_v)
        logging.debug(log_msg)
        log_msg = "[SPM] Magnet Shift: " + str(self.__delta_v * RAD2DEG)
        logging.debug(log_msg)
        log_msg = "[SPM] Slot Opening Positions: " + str(self.__alpha_so * RAD2DEG)
        logging.debug(log_msg)
        log_msg = "[SPM] Slot Positions: " + str(self.__alpha_s * RAD2DEG)
        logging.debug(log_msg)



    def __set_slot_geometry(self, spm):
        slot = spm.stator.slots[0]
        if slot.get_type() == 'Type0':
            self.__Rt = spm.stator.inner_radius + (slot.h0 + slot.h1) / 2.0
            self.__Rsb = self.__Rt + slot.h2 + slot.h3
            for i, sl in enumerate(spm.stator.slots):
                self.__Aso[i] = (sl.w0)
                self.__As[i] = (((sl.w1 + sl.w2) / 2.0))
                self.__alpha_so[i] = sl.so_position
                self.__alpha_s[i] = sl.s_position
        else:
            self.__Rt = spm.stator.inner_radius + (slot.h0 + slot.h1) / 2.0
            self.__Rsb = self.__Rt + slot.h2 + slot.h3
            for i, sl in enumerate(spm.stator.slots):
                self.__Aso[i] = (sl.w0)
                self.__As[i] = (((sl.w1 + sl.w2) / 2.0))
                self.__alpha_so[i] = sl.so_position
                self.__alpha_s[i] = sl.s_position
        

    def __set_magnet_geometry(self, spm):
        magnet = spm.rotor.magnets[0]
        if magnet.get_type() == 'Rectangular':
            self.__alpha_p_v = np.zeros(2 * self.__pp)
            self.__delta_v = np.zeros(2 * self.__pp)
            self.__M = np.zeros(2 * self.__pp)
            self.__mur = np.zeros(2 * self.__pp)
            for i, mn in enumerate(spm.rotor.magnets):
                self.__alpha_p_v[i] = 0.5 * mn.width
                self.__delta_v[i] = mn.deviation
                self.__M[i] = mn.material.Br / MU0
                self.__mur[i] = mn.material.mur

            self.__magnetisation = mn.magnetisation
        else:
            self.__alpha_p_v = np.zeros(2 * self.__pp)
            self.__delta_v = np.zeros(2 * self.__pp)
            self.__M = np.zeros(2 * self.__pp)
            self.__mur = np.zeros(2 * self.__pp)
            for i, mn in enumerate(spm.rotor.magnets):
                self.__alpha_p_v[i] = 0.5 * mn.width
                self.__delta_v[i] = mn.deviation
                self.__M[i] = mn.material.Br / MU0
                self.__mur[i] = mn.material.mur

            self.__magnetisation = mn.magnetisation

    def __set_winding_definition(self, spm):
        winding = spm.stator.winding
        if winding.get_type() == 'concentrated_coreless':
            self.__Nph = winding.phases
            self.__Ct = winding.turns_coil
            self.__Cp = winding.coil_parallel
            self.__Cs = winding.coil_series
            self.__C = winding.conn_matrix
            self.__WiH = winding.wires_in_hand
        else:
            self.__Nph = winding.phases
            self.__Ct = winding.turns_coil
            self.__Cp = winding.coil_parallel
            self.__Cs = winding.coil_series
            self.__C = winding.conn_matrix
            self.__WiH = winding.wires_in_hand

        # Assuming slot 1 is the average coil side width
        self.__coil_cs = winding.slot_fill_factor * (self.__y4 - self.__y3)*self.__Lx*self.__As[1] / self.__Ns
      
        Cpositions = self.__alpha_s
        #self.__results.stator_coil_resistance, self.__results.stator_phase_resistance = winding.get_resistances()
        self.__results.td = winding.turns_density(m=3, Ns=self.__Ns, Spos=Cpositions )
        self.__results.wf = winding.winding_function(m=3, pp=self.__pp, Ns=self.__Ns, Spos=Cpositions )
        self.__results.wh = winding.winding_harmonics(m=3, pp=self.__pp, Ns=self.__Ns, Spos=Cpositions )
        self.__results.kw_v = winding.winding_factors(m=3, pp=self.__pp, Ns=self.__Ns, Spos=Cpositions )

    def __configure(self, analysis_settings):
        self.noload_bemf = analysis_settings['noload'].get('bemf', False)
        self.noload_speed = analysis_settings['noload'].get('speed', 1000)
        self.cogging = analysis_settings['noload'].get('cogging', False)
        self.noload_pressure = analysis_settings['noload'].get('pressure', False)
        self.rippe = analysis_settings['load'].get('ripple', False)
        self.load_speed = analysis_settings['load'].get('speed', 1000)
        self.load_current = analysis_settings['load'].get('current', 100)
        self.load_voltage = analysis_settings['load'].get('voltage', None)
        self.load_gamma = analysis_settings['load'].get('gamma', 0)
        self.load_losses = analysis_settings['load'].get('losses', False)
        self.load_pressure = analysis_settings['load'].get('pressure', False)
        self.inductance_max_current = analysis_settings['inductance'].get('max_current', 10)
        self.inductance_steps = analysis_settings['inductance'].get('steps', 5)
        self.torque_speed_voltage_dc = analysis_settings['torque_speed'].get('voltage_dc', 12.0)
        self.torque_speed_m = analysis_settings['torque_speed'].get('m_index', 1.0)
        self.torque_speed_steps = analysis_settings['torque_speed'].get('speed_steps', 5)
        self.winding_function = analysis_settings['winding'].get('winding_function', False)
        self.winding_harmonics = analysis_settings['winding'].get('winding_harmonics', False)
        self.winding_factors = analysis_settings['winding'].get('winding_factors', False)

    def __assembly_A(self):
        muJ = 1
        muR = self.__mur[0]
        
        Nharms = np.arange(1, self.__n_harms)
        N = len(Nharms)
        Q = 4*self.__Ns   # Account for 2 coilsides + 2 coil separation
        tao_p = 2 * np.pi / (2*self.__pp)
        r = self.__r * self.__rfactor
        B = Nharms * np.pi / tao_p
        B_nn = np.diag(B)
        N_nn = np.diag(Nharms)
        G1_nn = np.diag(np.exp((B/r) * self.__y1))
        Inv_G1_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y1), dtype=float))
        G2_nn = np.diag(np.exp((B/r) * self.__y2))
        Inv_G2_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y2), dtype=float)) 
        G3_nn = np.diag(np.exp((B/r) * self.__y3))
        Inv_G3_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y3), dtype=float))
        G4_nn = np.diag(np.exp((Nharms/r) * self.__y3))
        Inv_G4_nn = np.diag(np.reciprocal(np.exp((Nharms/r) * self.__y3), dtype=float))
        In = np.eye(N, dtype=int)
        NpB = Nharms + B
        Inv_NpB = np.diag(np.reciprocal(NpB, dtype=float))
        alpha_nn = (-2j/(tao_p)) * np.dot(Inv_NpB, np.diag(np.exp(1j*tao_p*NpB)) - In)
        G5_nn = np.diag(np.exp((Nharms/r) * self.__y4))
        Inv_G5_nn = np.diag(np.reciprocal(np.exp((Nharms/r) * self.__y4), dtype=float))
        G6_nn = np.diag(np.exp((B/r) * self.__y4))
        Inv_G6_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y4), dtype=float))
        G7_nn = np.diag(np.exp((B/r) * self.__y5))
        Inv_G7_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y5), dtype=float))
        G8_nn = np.diag(np.exp((B/r) * self.__y6))
        Inv_G8_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y6), dtype=float))
        Zero_nn = np.zeros((N, N))
        


        Eq2_Xm = np.dot(G2_nn, Inv_G1_nn**2) - Inv_G2_nn
        Eq2_Wg = -muR * G2_nn
        Eq2_Xg = muR * Inv_G2_nn
        Eq2_Wj = Zero_nn
        Eq2_Xj = Zero_nn
        Eq2_Wg2 = Zero_nn
        Eq2_Xg2 = Zero_nn
        Eq2_Xm2 = Zero_nn


        Eq3_Xm = np.dot(G2_nn, Inv_G1_nn**2) + Inv_G2_nn
        Eq3_Wg = -G2_nn
        Eq3_Xg = -Inv_G2_nn
        Eq3_Wj = Zero_nn
        Eq3_Xj = Zero_nn
        Eq3_Wg2 = Zero_nn
        Eq3_Xg2 = Zero_nn
        Eq3_Xm2 = Zero_nn

        Eq4_Xm = Zero_nn
        Eq4_Wg = (np.pi/tao_p) * G3_nn
        Eq4_Xg = -(np.pi/tao_p) * Inv_G3_nn
        Eq4_Wj = -(1.0/muJ) * np.dot(alpha_nn, G4_nn)
        Eq4_Xj = (1.0/muJ) * np.dot(alpha_nn, Inv_G4_nn)
        Eq4_Wg2 = Zero_nn
        Eq4_Xg2 = Zero_nn
        Eq4_Xm2 = Zero_nn

        Eq5_Xm = Zero_nn
        Eq5_Wg = Zero_nn
        Eq5_Xg = Zero_nn
        Eq5_Wj = -(1.0/muJ) * np.dot(alpha_nn, G5_nn)
        Eq5_Xj = (1.0/muJ) * np.dot(alpha_nn, Inv_G5_nn)
        Eq5_Wg2 = (np.pi/tao_p) * G6_nn
        Eq5_Xg2 = -(np.pi/tao_p) * Inv_G6_nn
        Eq5_Xm2 = Zero_nn
 
        Eq6_Xm = Zero_nn
        Eq6_Wg = 1j*np.dot(B_nn, G3_nn)
        Eq6_Xg = 1j*np.dot(B_nn, Inv_G3_nn)
        Eq6_Wj = -1j*np.dot(np.dot(N_nn,G4_nn), alpha_nn)
        Eq6_Xj = -1j*np.dot(np.dot(N_nn, Inv_G4_nn), alpha_nn)
        Eq6_Wg2 = Zero_nn
        Eq6_Xg2 = Zero_nn
        Eq6_Xm2 = Zero_nn

        Eq7_Xm = Zero_nn
        Eq7_Wg = Zero_nn
        Eq7_Xg = Zero_nn
        Eq7_Wj = -1j*np.dot(np.dot(N_nn,G5_nn), alpha_nn)
        Eq7_Xj = -1j*np.dot(np.dot(N_nn, Inv_G5_nn), alpha_nn)
        Eq7_Wg2 = 1j*np.dot(B_nn, G6_nn)
        Eq7_Xg2 = 1j*np.dot(B_nn, Inv_G6_nn)
        Eq7_Xm2 = Zero_nn
        

        Eq8_Xm2 = np.dot(G7_nn, Inv_G8_nn**2) - Inv_G7_nn
        Eq8_Wg2 = -muR * G7_nn
        Eq8_Xg2 = muR * Inv_G7_nn
        Eq8_Wj = Zero_nn
        Eq8_Xj = Zero_nn
        Eq8_Wg = Zero_nn
        Eq8_Xg = Zero_nn
        Eq8_Xm = Zero_nn


        Eq9_Xm2 = np.dot(G7_nn, Inv_G8_nn**2) + Inv_G7_nn
        Eq9_Wg2 = -G7_nn
        Eq9_Xg2 = -Inv_G7_nn
        Eq9_Wj = Zero_nn
        Eq9_Xj = Zero_nn
        Eq9_Wg = Zero_nn
        Eq9_Xg = Zero_nn
        Eq9_Xm = Zero_nn


        Eq2 = np.concatenate(
            (Eq2_Xm, Eq2_Wg, Eq2_Xg, Eq2_Wj, Eq2_Xj, Eq2_Wg2, Eq2_Xg2, Eq2_Xm2), 
            axis=1)
        Eq3 = np.concatenate(
            (Eq3_Xm, Eq3_Wg, Eq3_Xg, Eq3_Wj, Eq3_Xj, Eq3_Wg2, Eq3_Xg2, Eq3_Xm2),  
            axis=1)
        Eq4 = np.concatenate(
            (Eq4_Xm, Eq4_Wg, Eq4_Xg, Eq4_Wj, Eq4_Xj, Eq4_Wg2, Eq4_Xg2, Eq4_Xm2), 
            axis=1)
        Eq5 = np.concatenate(
            (Eq5_Xm, Eq5_Wg, Eq5_Xg, Eq5_Wj, Eq5_Xj, Eq5_Wg2, Eq5_Xg2, Eq5_Xm2),  
            axis=1)
        Eq6 = np.concatenate(
            (Eq6_Xm, Eq6_Wg, Eq6_Xg, Eq6_Wj, Eq6_Xj, Eq6_Wg2, Eq6_Xg2, Eq6_Xm2), 
            axis=1)
        Eq7 = np.concatenate(
            (Eq7_Xm, Eq7_Wg, Eq7_Xg, Eq7_Wj, Eq7_Xj, Eq7_Wg2, Eq7_Xg2, Eq7_Xm2), 
            axis=1)
        Eq8 = np.concatenate(
            (Eq8_Xm, Eq8_Wg, Eq8_Xg, Eq8_Wj, Eq8_Xj, Eq8_Wg2, Eq8_Xg2, Eq8_Xm2), 
            axis=1)
        Eq9 = np.concatenate(
            (Eq9_Xm, Eq9_Wg, Eq9_Xg, Eq9_Wj, Eq9_Xj, Eq9_Wg2, Eq9_Xg2, Eq9_Xm2), 
            axis=1)

        self._A_ = np.concatenate((Eq2, Eq3, Eq4, Eq5, Eq6, Eq7, Eq8, Eq9), axis=0)


    def __assembly_b_nl(self, pos = [0.0]):
        Nharms = np.arange(1, self.__n_harms)
        N = len(Nharms)

        r = self.__r * self.__rfactor

        Eq2 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq3 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq4 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq5 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq6 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq7 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq8 = np.zeros((N, len(pos)), dtype=np.complex)
        Eq9 = np.zeros((N, len(pos)), dtype=np.complex)

        tao_p = 2 * np.pi / (2*self.__pp)
        B = Nharms * np.pi / tao_p
        Inv_G1_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y1), dtype=float))
        G2_nn = np.diag(np.exp((B/r) * self.__y2))
        G7_nn = np.diag(np.exp((B/r) * self.__y5))
        Inv_G8_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y6), dtype=float))

        for p in range(0, len(pos)):
            self._Mrn, self._Mtn = get_fs_coeff_linear(self.__magnetisation, self.__n_harms, self.__pp, pos[p], self.__M,
                         self.__alpha_p_v, self.__delta_v, self.__Lx, 0)

            # Mrycn = 2.0 * np.real(self._Mrn)
            # Mrysn = -2.0 * np.imag(self._Mrn)
            # Mrxcn = 2.0 * np.real(self._Mtn)
            # Mrxsn = -2.0 * np.imag(self._Mtn)

            for i in range(0, N):

                
                Eq2[i,p] = (1.0/B[i]) * r * (MU0 - G2_nn[i,i]*Inv_G1_nn[i,i]) * self._Mtn[i+1]
                Eq3[i,p] = -r * MU0 * (1.0/B[i]) * (1j*self._Mrn[i+1] + G2_nn[i,i]*Inv_G1_nn[i,i]*self._Mtn[i+1])
                Eq8[i,p] = (1.0/B[i]) * r * (MU0 - G7_nn[i,i]*Inv_G8_nn[i,i]) * self._Mtn[i+1]
                Eq9[i,p] = -r * MU0 * (1.0/B[i]) * (1j*self._Mrn[i+1] + G7_nn[i,i]*Inv_G8_nn[i,i]*self._Mtn[i+1])
                

        self._b_ = np.concatenate((Eq2, Eq3, Eq4, Eq5, Eq6, Eq7, Eq8, Eq9), axis=0)


    def __assembly_b_ol(self, pos = [0.0], samePos = True, I = [0.0, 0.0, 0.0]):
        if(samePos==True):
            b = np.tile(self._b_, 2)
            self._b_ = b
        else:
            Nharms = np.arange(1, self.__n_harms)

            N = len(Nharms)

            r = self.__r * self.__rfactor
            tao_p = 2 * np.pi / (2*self.__pp)

            B = Nharms * np.pi / tao_p
            Inv_G1_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y1), dtype=float))
            G2_nn = np.diag(np.exp((B/r) * self.__y2))
            G7_nn = np.diag(np.exp((B/r) * self.__y5))
            Inv_G8_nn = np.diag(np.reciprocal(np.exp((B/r) * self.__y6), dtype=float))

            Eq2 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq3 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq4 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq5 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq6 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq7 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq8 = np.zeros((N, len(pos)), dtype=np.complex)
            Eq9 = np.zeros((N, len(pos)), dtype=np.complex)


            for p in range(0, len(pos)):
                self._Mrn, self._Mtn = get_fs_coeff_linear(self.__magnetisation, self.__n_harms, self.__pp, pos[p], self.__M,
                                                    self.__alpha_p_v, self.__delta_v, self.__Lx, 0)

                # Mrycn = 2.0 * np.real(self._Mrn)
                # Mrysn = -2.0 * np.imag(self._Mrn)
                # Mrxcn = 2.0 * np.real(self._Mtn)
                # Mrxsn = -2.0 * np.imag(self._Mtn)

                three_phase_td = np.dot(I[p,:],self.__results.wf[:,0:-1])
                N_half = (np.rint(len(three_phase_td)/2.0)).astype(int)
                td_fft = fft(three_phase_td)*np.exp(1j*np.pi/2)/len(three_phase_td)

                for i in range(0, N):
                    
                    Eq2[i,p] = (1.0/B[i]) * r * (MU0 - G2_nn[i,i]*Inv_G1_nn[i,i]) * self._Mtn[i+1]
                    Eq3[i,p] = -r * MU0 * (1.0/B[i]) * (1j*self._Mrn[i+1] + G2_nn[i,i]*Inv_G1_nn[i,i]*self._Mtn[i+1])

                    n = Nharms[i]
                    b = n*np.pi / tao_p
                    npb = n+b
                    alpha = -2j*(np.exp(1j*tao_p*npb) - 1)/(tao_p*npb)
                    
                    # TO-DO: harmonics must be less than 180 due to Turns Density size
                    Jn = 2.0*td_fft[n]/self.__coil_cs
                    Eq6[i,p] = -(r**2) * MU0 * self.__mur[0] * alpha * Jn / n
                    Eq7[i,p] = -(r**2) * MU0 * self.__mur[0] * alpha * Jn / n
                    Eq8[i,p] = (1.0/B[i]) * r * (MU0 - G7_nn[i,i]*Inv_G8_nn[i,i]) * self._Mtn[i+1]
                    Eq9[i,p] = -r * MU0 * (1.0/B[i]) * (1j*self._Mrn[i+1] + G7_nn[i,i]*Inv_G8_nn[i,i]*self._Mtn[i+1])

            b = np.concatenate((Eq2, Eq3, Eq4, Eq5, Eq6, Eq7, Eq8, Eq9), axis=0)
            self._b_ = np.concatenate((self._b_, b), axis=1)


    def get_self_and_mutual_inductance(self):
        """Get self and mutual inductance by means of winding function

        Compute self and mutual inductance according to Paul et al. 'Drive response modeling of dual
        wound surface permanent magnet machines' in IEMDC (2017)

        Args:
            posNL:      Rotor positions to be computed under no load

        Returns:
            Lself:      Self Inductance
            Lmutual:    Mutual Inductance
        """
        tao_s = 2 * PI * self.__Rs / self.__Ns
        mur = self.__mur
        Ml = self.__Rm - self.__oRr
        g = self.__Rs - self.__Rm
        gc = g + Ml / mur
        w0 = self.__Aso[0] * self.__Rs
        kcs = tao_s / (tao_s - (2 * w0 / PI) * (np.arctan2(w0, 2 * gc) -
                                          (gc / w0) * np.log(1 + (w0 / (2 * gc))**2)))
        ge = gc * kcs

        #Ps = self._stator._winding.slot_permeances()

        #print (gc, kcs, ge, Ps)
        Lself = 0.0
        Lmutual = 0.0

        return Lself, Lmutual


    def get_ag_flux_density(self, posNL = [0.0], posOL = [0.0], samePos = True, current = [0.0, 0.0, 0.0],
                            psi = np.linspace(0, 2*PI, 360)):
        """Get the flux density in the air gap by sub-domain method

        Compute air gap flux density according to Pina Ortega et al. 'Analytical Model for
        predicting effects of manufacturing variations on cogging torque in surface-mounted
        permanent magnet motors' in IEEE Transactions on Industry Application (2016), and 'Analytical
        prediction of torque ripple in surface-mounted permanent magnet motors due to manufacturing
        variations' in IEEE Transactions on Energy Conversion (2016)

        Args:
            posNL:      Rotor positions to be computed under no load
            posOL:      Rotor positions to be computed under load
            samePos:    No load and load positions are the same
            current:    Current vector
            psi:        Spatial position in the air gap

        Returns:
            Bg_r:        Radial component for each position evaluated
            Bg_t:        Tangential component for each position evaluated
        """
        self.__assembly_A()
        self.__assembly_b_nl(pos=posNL)
        self.__assembly_b_ol(pos=posOL,samePos=samePos, I=current)
        unk = sl.solve(self._A_, self._b_)
        Nharms = np.arange(1, self.__n_harms)
        N = len(Nharms)
        Wg = unk[1 * N:2 * N, : ]
        Xg = unk[2 * N:3 * N, : ]
        Wg2 = unk[5 * N:6 * N, : ]
        Xg2 = unk[6 * N:7 * N, : ]
        g = self.__y3 - self.__y2
        g2 = self.__y5 - self.__y4
        # Flux Density Close to stator for better prediction of flux linkage and
        # Radial forces
        y = self.__y3 - g / 2.0
        y2 = self.__y5 - g2 / 2.0
        Bg_r = np.zeros((len(posNL)+len(posOL),len(psi)))
        Bg_t = np.zeros((len(posNL)+len(posOL),len(psi)))
        Bg2_r = np.zeros((len(posNL)+len(posOL),len(psi)))
        Bg2_t = np.zeros((len(posNL)+len(posOL),len(psi)))
        prd = 2*PI/(self.__Lx*self.__pp)
        tao_p = 2 * np.pi / (2*self.__pp)
        r = self.__r * self.__rfactor
        for p in range(0, len(posNL)+len(posOL)):
            for i in range(0, N):
                n = Nharms[i]
                b = n * np.pi / tao_p
                e_pkny = np.exp((b/r)*y)
                e_mkny = np.exp(-(b/r)*y)
                term1 = (b/r) * (Wg[i,p] * e_pkny - Xg[i,p] * e_mkny)
                term2 = (1j*b/r) * (Wg[i,p] * e_pkny + Xg[i,p] * e_mkny)
                
                e_pkny2 = np.exp((b/r)*y2)
                e_mkny2 = np.exp(-(b/r)*y2)
                term3 = (b/r) * (Wg2[i,p] * e_pkny2 - Xg2[i,p] * e_mkny2)
                term4 = (1j*b/r) * (Wg2[i,p] * e_pkny2 + Xg2[i,p] * e_mkny2)
                Bg_t[p,:] = np.real(Bg_t[p,:] + 2.0 * term1 * np.exp(1j*b*psi*prd))
                Bg_r[p,:] = np.real(Bg_r[p,:] + 2.0 * term2 * np.exp(1j*b*psi*prd))
                Bg2_t[p,:] = np.real(Bg2_t[p,:] + 2.0 * term3 * np.exp(1j*b*psi*prd))
                Bg2_r[p,:] = np.real(Bg2_r[p,:] + 2.0 * term4 * np.exp(1j*b*psi*prd))

        return Bg_r, Bg_t, Bg2_r, Bg2_t



    def __get_settings_cogging(self):
        final_pos = 360 / self.__lcm
        steps = 7
        posVector = np.linspace(0, final_pos, steps, endpoint=False) * DEG2RAD  + self.__init_pos
        return posVector


    def __get_settings_ripple(self):
        final_pos = 360 / self.__lcm
        steps = 7
        posVector = np.linspace(0, final_pos, steps, endpoint=False) * DEG2RAD + self.__init_pos
        thetaElec = self.__pp * np.linspace(0, final_pos, steps, endpoint=False) * DEG2RAD
        Is = self.load_current / self.__Cp 
        IVector = np.transpose([Is * np.sin(thetaElec + self.load_gamma*PI/180.0), Is * np.sin(thetaElec - 2.0*PI/3.0 + self.load_gamma*PI/180.0),
                                Is * np.sin(thetaElec + 2.0*PI/3.0 + self.load_gamma*PI/180.0)])
        return posVector, IVector

    def __get_settings_static_torque(self,  pos=np.array([]), I=np.array([])):
        final_pos = 360 / self.__pp
        steps = 9
        posVector = np.linspace(0, final_pos, steps) * DEG2RAD + self.__init_pos
        Is = self.load_current / self.__Cp 
        IVector = np.transpose([Is * np.ones(steps), -0.5 * Is * np.ones(steps),
                                -0.5 * Is * np.ones(steps)])
        return posVector, IVector

    def __get_settings_inductance_calc(self):
        posVector = np.array([self.__init_pos])
        IVector = np.array([self.load_current, 0.0, 0.0]) / self.__Cp 
        return posVector, IVector

    def __get_settings_torque_speed(self):
        posVector = np.array([self.__init_pos, self.__init_pos])
        Is = self.torque_speed_current / self.__Cp 
        IVector = np.array([[0.0, Is * np.sin(-2.0*PI/3.0), Is * np.sin(-2.0*PI/3.0)],
            [Is * np.sin(PI/2.0), Is * np.sin(-PI/6.0), Is * np.sin(7.0*PI/6.0)]])
        return posVector, IVector

    def __get_settings_bemf(self, pos=[]):
        init_pos = pos[-1]
        steps = 3
        final_pos = 2*PI / (3 * self.__pp) + self.__init_pos        # Max span 120 electric degrees
        step_size = np.abs(init_pos - final_pos) / steps
        posVector = np.linspace(init_pos + step_size, final_pos - step_size, steps, endpoint=False)
        return posVector, self.noload_speed

    def solve(self):
        if self.cogging:
            pos1 = self.__get_settings_cogging()
        if self.noload_bemf:
            pos2, wr_nl = self.__get_settings_bemf(pos=pos1)

        pos3, I1 = self.__get_settings_ripple()
        pos4, I2 = self.__get_settings_static_torque()
        pos5, I3 = self.__get_settings_inductance_calc()
        #pos6, I4 = self.__get_settings_torque_speed()
        

        posNL = np.hstack((pos1, pos2))

        posOL = np.hstack((pos3, pos4, pos5))
        I = np.vstack((I1, I2, I3))


        Bg_r, Bg_t, Bg2_r, Bg2_t = self.get_ag_flux_density(posNL=posNL, posOL=posOL, samePos=False,
                                              current=I,psi=np.linspace(0,self.__Lx,360))
        
        r = self.__r
        Lstk = (self.__Sl + self.__Rl) / 2.0
        torque = (Lstk * (r**2) / MU0) * np.trapz(Bg_r * Bg_t, axis=1, dx=PI/180)
        torque2 = (Lstk * (r**2) / MU0) * np.trapz(Bg2_r * Bg2_t, axis=1, dx=PI/180)
        pos_nl_cs = (pos1-self.__init_pos)*RAD2DEG
        cogging_cs = si.CubicSpline(pos_nl_cs,
                                    torque[0:len(pos1)],
                                    extrapolate='periodic')
        cogging_acs = np.linspace(0,360/self.__pp,180)
        cogging_cs2 = si.CubicSpline(pos_nl_cs,
                                    torque2[0:len(pos1)],
                                    extrapolate='periodic')
        
        pos_ol_cs = (pos3 - self.__init_pos) * RAD2DEG
        ripple_cs = si.CubicSpline(pos_ol_cs,
                                    torque[len(posNL):len(posNL)+len(pos3)],
                                    extrapolate='periodic')
        ripple_acs = np.linspace(0, 360 / self.__pp, 180)
        ripple_cs2 = si.CubicSpline(pos_ol_cs,
                                    torque2[len(posNL):len(posNL)+len(pos3)],
                                    extrapolate='periodic')

        pos_st_cs = (pos4 - self.__init_pos) * RAD2DEG
        staticTq_cs = si.CubicSpline(pos_st_cs,
                                   torque[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4)])
        staticTq_acs = np.linspace(0, 360 / self.__pp, 45)
        staticTq_cs2 = si.CubicSpline(pos_st_cs,
                                   torque2[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4)])


        self.__results.cogging_torque_x = cogging_acs
        self.__results.cogging_torque2_y = cogging_cs2(cogging_acs)
        self.__results.cogging_torque_y = cogging_cs(cogging_acs)
        self.__results.torque_ripple_x = ripple_acs
        self.__results.torque_ripple2_y = ripple_cs2(ripple_acs)
        self.__results.torque_ripple_y = ripple_cs(ripple_acs)
        self.__results.static_torque_x = staticTq_acs
        self.__results.static_torque2_y = staticTq_cs2(staticTq_acs)
        self.__results.static_torque_y = staticTq_cs(staticTq_acs)
        self.__results.nl_Bg_theta = np.linspace(0, self.__Lx, 360)
        self.__results.nl_Bg_r = Bg_r[0,:]
        self.__results.nl_Bg_t = Bg_t[0,:]
        self.__results.ol_Bg_r = Bg_r[len(posNL), :]
        self.__results.ol_Bg_t = Bg_t[len(posNL), :]
        self.__results.nl_Bg2_r = Bg2_r[0,:]
        self.__results.nl_Bg2_t = Bg2_t[0,:]
        self.__results.ol_Bg2_r = Bg2_r[len(posNL), :]
        self.__results.ol_Bg2_t = Bg2_t[len(posNL), :]


        nrows = len(posNL) + len(pos3) + len(pos4)

        fl = np.zeros((nrows,3))
        fl2 = np.zeros((nrows,3))

        for i in range(0, nrows):
            fl[i] = (Lstk * r / self.__Cp) * np.trapz(Bg_r[i,:] * self.__results.wf,  dx=PI/180)
            fl2[i] = (Lstk * r / self.__Cp) * np.trapz(Bg2_r[i,:] * self.__results.wf,  dx=PI/180)
 
        fl_a = np.concatenate((fl[0:len(posNL), 0], fl[0:len(posNL), 2], fl[0:len(posNL), 1], fl[0:len(posNL), 0],
                               fl[len(posNL):len(posNL)+len(pos3), 0], fl[len(posNL):len(posNL)+len(pos3), 2],
                               fl[len(posNL):len(posNL)+len(pos3), 1], fl[len(posNL):len(posNL)+len(pos3), 0],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0] ),
                              axis=0)
        fl_b = np.concatenate((fl[0:len(posNL), 1], fl[0:len(posNL), 0], fl[0:len(posNL), 2], fl[0:len(posNL), 1],
                               fl[len(posNL):len(posNL) + len(pos3), 1], fl[len(posNL):len(posNL) + len(pos3), 0],
                               fl[len(posNL):len(posNL) + len(pos3), 2], fl[len(posNL):len(posNL) + len(pos3), 1],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1]),
                              axis=0)
        fl_c = np.concatenate((fl[0:len(posNL), 2], fl[0:len(posNL), 1], fl[0:len(posNL), 0], fl[0:len(posNL), 2],
                               fl[len(posNL):len(posNL)+len(pos3), 2], fl[len(posNL):len(posNL)+len(pos3), 1],
                               fl[len(posNL):len(posNL)+len(pos3), 0], fl[len(posNL):len(posNL)+len(pos3), 2],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2] ),
                              axis=0)
        fl2_a = np.concatenate((fl2[0:len(posNL), 0], fl2[0:len(posNL), 2], fl2[0:len(posNL), 1], fl2[0:len(posNL), 0],
                               fl2[len(posNL):len(posNL)+len(pos3), 0], fl2[len(posNL):len(posNL)+len(pos3), 2],
                               fl2[len(posNL):len(posNL)+len(pos3), 1], fl2[len(posNL):len(posNL)+len(pos3), 0],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0] ),
                              axis=0)
        fl2_b = np.concatenate((fl2[0:len(posNL), 1], fl2[0:len(posNL), 0], fl2[0:len(posNL), 2], fl2[0:len(posNL), 1],
                               fl2[len(posNL):len(posNL) + len(pos3), 1], fl2[len(posNL):len(posNL) + len(pos3), 0],
                               fl2[len(posNL):len(posNL) + len(pos3), 2], fl2[len(posNL):len(posNL) + len(pos3), 1],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1]),
                              axis=0)
        fl2_c = np.concatenate((fl2[0:len(posNL), 2], fl2[0:len(posNL), 1], fl2[0:len(posNL), 0], fl2[0:len(posNL), 2],
                               fl2[len(posNL):len(posNL)+len(pos3), 2], fl2[len(posNL):len(posNL)+len(pos3), 1],
                               fl2[len(posNL):len(posNL)+len(pos3), 0], fl2[len(posNL):len(posNL)+len(pos3), 2],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 1],
                               fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 0], fl[len(posNL)+len(pos3):len(posNL)+len(pos3)+len(pos4), 2] ),
                              axis=0)

        pos_nl_bemf_cs = (posNL - self.__init_pos) * RAD2DEG

        pos_nl_abc = np.concatenate( (self.__pp*pos_nl_bemf_cs,
                                      self.__pp*pos_nl_bemf_cs + 120,
                                      self.__pp*pos_nl_bemf_cs + 240,
                                      self.__pp*pos_nl_bemf_cs + 360 ), axis=0 )

        fl_a_nl_cs = si.CubicSpline(pos_nl_abc,
                                    fl_a[0:4*len(posNL)])
        fl_b_nl_cs = si.CubicSpline(pos_nl_abc,
                                    fl_b[0:4*len(posNL)])
        fl_c_nl_cs = si.CubicSpline(pos_nl_abc,
                                    fl_c[0:4*len(posNL)])
        fl2_a_nl_cs = si.CubicSpline(pos_nl_abc,
                                    fl2_a[0:4*len(posNL)])
        fl2_b_nl_cs = si.CubicSpline(pos_nl_abc,
                                    fl2_b[0:4*len(posNL)])
        fl2_c_nl_cs = si.CubicSpline(pos_nl_abc,
                                    fl2_c[0:4*len(posNL)])

        self.__results.nl_flux_linkage_y = np.array( [fl_a_nl_cs(THETA_e_DEG),
                                         fl_b_nl_cs(THETA_e_DEG),
                                         fl_c_nl_cs(THETA_e_DEG)] )
        self.__results.nl_flux_linkage2_y = np.array( [fl2_a_nl_cs(THETA_e_DEG),
                                         fl2_b_nl_cs(THETA_e_DEG),
                                         fl2_c_nl_cs(THETA_e_DEG)] )
        self.__results.nl_flux_linkage_x = THETA_e_DEG


        # Extract fundamental of flux linkage due to magnet
        fft_magnet_flux = fft(self.__results.nl_flux_linkage_y[0,:])/len(self.__results.nl_flux_linkage_y[0,:])
        self.__results.magnet_flux = 2.0 * np.abs(fft_magnet_flux[1])
        fft_magnet_flux2 = fft(self.__results.nl_flux_linkage2_y[0,:])/len(self.__results.nl_flux_linkage2_y[0,:])
        self.__results.magnet_flux2 = 2.0 * np.abs(fft_magnet_flux2[1])

        bemf_a = np.zeros( EL_STEPS )
        bemf_b = np.zeros(EL_STEPS)
        bemf_c = np.zeros(EL_STEPS)
        bemf2_a = np.zeros( EL_STEPS )
        bemf2_b = np.zeros(EL_STEPS)
        bemf2_c = np.zeros(EL_STEPS)
        for i in range(1,30,2):
            bemf_a = ( bemf_a + (1.0/EL_STEPS) * ( 1j*i*wr_nl*self.__pp*(PI/30) ) *
                       ( fft_magnet_flux[i]*np.exp( 1j*i*THETA_e_RAD ) -
                        np.conjugate(fft_magnet_flux[i])*np.exp( -1j*i*THETA_e_RAD ) ) )
            bemf_b = ( bemf_b + (1.0/EL_STEPS) * ( 1j*i*wr_nl*self.__pp*(PI/30) ) *
                       ( fft_magnet_flux[i]*np.exp( 1j*i*(THETA_e_RAD - PI_2by3) ) -
                        np.conjugate(fft_magnet_flux[i])*np.exp( -1j*i*(THETA_e_RAD - PI_2by3) ) ) )
            bemf_c = (bemf_c + (1.0 / EL_STEPS) * (1j * i * wr_nl * self.__pp * (PI / 30)) *
                      (fft_magnet_flux[i] * np.exp(1j * i * (THETA_e_RAD + PI_2by3)) -
                       np.conjugate(fft_magnet_flux[i]) * np.exp(-1j * i * (THETA_e_RAD + PI_2by3))))
            bemf2_a = ( bemf2_a + (1.0/EL_STEPS) * ( 1j*i*wr_nl*self.__pp*(PI/30) ) *
                       ( fft_magnet_flux2[i]*np.exp( 1j*i*THETA_e_RAD ) -
                        np.conjugate(fft_magnet_flux2[i])*np.exp( -1j*i*THETA_e_RAD ) ) )
            bemf2_b = ( bemf2_b + (1.0/EL_STEPS) * ( 1j*i*wr_nl*self.__pp*(PI/30) ) *
                       ( fft_magnet_flux2[i]*np.exp( 1j*i*(THETA_e_RAD - PI_2by3) ) -
                        np.conjugate(fft_magnet_flux2[i])*np.exp( -1j*i*(THETA_e_RAD - PI_2by3) ) ) )
            bemf2_c = (bemf2_c + (1.0 / EL_STEPS) * (1j * i * wr_nl * self.__pp * (PI / 30)) *
                      (fft_magnet_flux2[i] * np.exp(1j * i * (THETA_e_RAD + PI_2by3)) -
                       np.conjugate(fft_magnet_flux2[i]) * np.exp(-1j * i * (THETA_e_RAD + PI_2by3))))

        self.__results.bemf_y = np.array( [ np.real( bemf_a ), np.real( bemf_b ), np.real( bemf_c )] )
        self.__results.bemf2_y = np.array( [ np.real( bemf2_a ), np.real( bemf2_b ), np.real( bemf2_c )] )
        self.__results.bemf_x = THETA_e_DEG

        pos_ol_abc = np.concatenate( (self.__pp*pos_ol_cs,
                                      self.__pp*pos_ol_cs + 120,
                                      self.__pp*pos_ol_cs + 240,
                                      self.__pp*pos_ol_cs + 360 ), axis=0 )

        fl_a_ol_cs = si.CubicSpline(pos_ol_abc,
                                    fl_a[4*len(posNL):4*(len(posNL)+len(pos3))])
        fl_b_ol_cs = si.CubicSpline(pos_ol_abc,
                                    fl_b[4*len(posNL):4*(len(posNL)+len(pos3))])
        fl_c_ol_cs = si.CubicSpline(pos_ol_abc,
                                    fl_c[4*len(posNL):4*(len(posNL)+len(pos3))])
        fl2_a_ol_cs = si.CubicSpline(pos_ol_abc,
                                    fl2_a[4*len(posNL):4*(len(posNL)+len(pos3))])
        fl2_b_ol_cs = si.CubicSpline(pos_ol_abc,
                                    fl2_b[4*len(posNL):4*(len(posNL)+len(pos3))])
        fl2_c_ol_cs = si.CubicSpline(pos_ol_abc,
                                    fl2_c[4*len(posNL):4*(len(posNL)+len(pos3))])

        ePos4 = (pos4 - self.__init_pos)*self.__pp

        fl_a_st_cs = (fl_a[4*(len(posNL)+len(pos3)):4*(len(posNL)+len(pos3))+len(pos4)]+
                        fl2_a[4*(len(posNL)+len(pos3)):4*(len(posNL)+len(pos3))+len(pos4)])
        fl_b_st_cs = (fl_b[4*(len(posNL)+len(pos3)):4*(len(posNL)+len(pos3))+len(pos4)]+
                        fl2_b[4*(len(posNL)+len(pos3)):4*(len(posNL)+len(pos3))+len(pos4)])
        fl_c_st_cs = (fl_c[4*(len(posNL)+len(pos3)):4*(len(posNL)+len(pos3))+len(pos4)]+
                        fl2_c[4*(len(posNL)+len(pos3)):4*(len(posNL)+len(pos3))+len(pos4)])



        k_qd_pos4 = (2.0/3.0)*np.array([[np.sin(ePos4), np.sin(ePos4-2.0*PI/3.0), np.sin(ePos4+2.0*PI/3.0)],
                                        [-np.cos(ePos4), -np.cos(ePos4-2.0*PI/3.0), -np.cos(ePos4+2.0*PI/3.0)]])
        iq_st = self.__Cp * np.diag(np.dot( np.atleast_2d(np.transpose(k_qd_pos4[0,:,:])), np.atleast_2d(np.transpose(I2))))
        id_st = self.__Cp *np.diag(np.dot( np.atleast_2d(np.transpose(k_qd_pos4[1,:,:])), np.atleast_2d(np.transpose(I2))))
        fq_st = np.diag(np.dot( np.atleast_2d(np.transpose(k_qd_pos4[0,:,:])), np.vstack((fl_a_st_cs , fl_b_st_cs, fl_c_st_cs))))
        fd_st = np.diag(np.dot( np.atleast_2d(np.transpose(k_qd_pos4[1,:,:])), np.vstack((fl_a_st_cs , fl_b_st_cs, fl_c_st_cs))))
        
        pos0 = (np.where((np.rint(ePos4*RAD2DEG)).astype(int)==90))[0][0]
        pos90 = (np.where((np.rint(ePos4*RAD2DEG)).astype(int)==0))[0][0]
        iq_st_0 = iq_st[pos0]
        iq_st_90 = iq_st[pos90]
        id_st_0 = id_st[pos0]
        id_st_90 = id_st[pos90]
        fq_st_0 = fq_st[pos0]
        fq_st_90 = fq_st[pos90]
        fd_st_0 = fd_st[pos0]
        fd_st_90 = fd_st[pos90]
        Lqd = fq_st_0 / iq_st_0
        fm = self.__results.magnet_flux + self.__results.magnet_flux2 

        rpm_max = (30.0/(self.__pp*PI))*self.torque_speed_m*self.torque_speed_voltage_dc/(1.731*fd_st_90)
        rpm_base = (30.0/(self.__pp*PI))*(self.torque_speed_m*self.torque_speed_voltage_dc/1.731)*(
            1.0/np.sqrt(fm**2 + Lqd**2 * iq_st_0**2)
        )
        self.__results.torque_speed_rpm = []
        self.__results.torque_speed_tq = []
        self.__results.torque_speed_iq = []
        self.__results.torque_speed_id = []
        for i in np.linspace(0.0 , 1.0, self.torque_speed_steps, endpoint=False):
            iq_i = (1.0-i) * self.load_current
            id_i = -i * self.load_current
            gamma_i = np.arctan2(-id_i, iq_i) * RAD2DEG + 90
            rpm_i =  (30.0/(self.__pp*PI))*(self.torque_speed_m*self.torque_speed_voltage_dc/1.731)*(
                1.0/np.sqrt(fm**2 + Lqd**2 * self.load_current**2 + 2*fm*Lqd*id_i)
            )
            Tq_i = -staticTq_cs(gamma_i/self.__pp)+staticTq_cs2(gamma_i/self.__pp)
            if i==0:
                self.__results.torque_speed_rpm.append(0)
                self.__results.torque_speed_tq.append(Tq_i)
                self.__results.torque_speed_iq.append(iq_i)
                self.__results.torque_speed_id.append(id_i)
            self.__results.torque_speed_rpm.append(rpm_i)
            self.__results.torque_speed_tq.append(Tq_i)
            self.__results.torque_speed_iq.append(iq_i)
            self.__results.torque_speed_id.append(id_i)
        self.__results.torque_speed_rpm.append(rpm_max)
        self.__results.torque_speed_tq.append(0)
        self.__results.torque_speed_iq.append(0)
        self.__results.torque_speed_id.append(id_st_90)

        self.__results.ol_flux_linkage_y = np.array( [fl_a_ol_cs(THETA_e_DEG),
                                         fl_b_ol_cs(THETA_e_DEG),
                                         fl_c_ol_cs(THETA_e_DEG)] )
        self.__results.ol_flux_linkage2_y = np.array( [fl2_a_ol_cs(THETA_e_DEG),
                                         fl2_b_ol_cs(THETA_e_DEG),
                                         fl2_c_ol_cs(THETA_e_DEG)] )
        self.__results.ol_flux_linkage_x = THETA_e_DEG

        i_a = np.concatenate((I1[:, 0], I1[:, 2], I1[:, 1], I1[:, 0]))
        i_b = np.concatenate((I1[:, 1], I1[:, 0], I1[:, 2], I1[:, 1]))
        i_c = np.concatenate((I1[:, 2], I1[:, 1], I1[:, 0], I1[:, 2]))
        i_a_cs = si.CubicSpline(pos_ol_abc, i_a)
        i_b_cs = si.CubicSpline(pos_ol_abc, i_b)
        i_c_cs = si.CubicSpline(pos_ol_abc, i_c)

        self.__results.phase_current_y = self.__Cp * np.array( [  i_a_cs(THETA_e_DEG),
                                                    i_b_cs(THETA_e_DEG),
                                                   i_c_cs(THETA_e_DEG)])
        self.__results.phase_current_x = THETA_e_DEG

        fl_only_phA_current = (Lstk * r / self.__Cp) * np.trapz(Bg_r[-1, :] * self.__results.wf, dx=PI / 180)

        self.__results.self_inductance_ag = (np.abs((fl_only_phA_current[0] - self.__results.nl_flux_linkage_y[0,0]) /
                                      (self.load_current)))
        self.__results.mutual_inductance_ag = -(np.abs((fl_only_phA_current[1] - self.__results.nl_flux_linkage_y[1,0]) /
                                         (self.load_current)))


        iq = ( K_QD[0, 0, :] * self.__results.phase_current_y[0, :] +
                         K_QD[0, 1, :] * self.__results.phase_current_y[1, :] +
                         K_QD[0, 2, :] * self.__results.phase_current_y[2, :])
        id = ( K_QD[1, 0, :] * self.__results.phase_current_y[0, :] +
                        K_QD[1, 1, :] * self.__results.phase_current_y[1, :] +
                        K_QD[1, 2, :] * self.__results.phase_current_y[2, :])

        fl_abc_noMagnet = - self.__results.nl_flux_linkage_y + self.__results.ol_flux_linkage_y

        fq = np.average(K_QD[0, 0, :] * fl_abc_noMagnet[0, :] +
                        K_QD[0, 1, :] * fl_abc_noMagnet[1, :] +
                        K_QD[0, 2, :] * fl_abc_noMagnet[2, :])
        fd = np.average(K_QD[1, 0, :] * fl_abc_noMagnet[0, :] +
                        K_QD[1, 1, :] * fl_abc_noMagnet[1, :] +
                        K_QD[1, 2, :] * fl_abc_noMagnet[2, :])

        fq_ol = (K_QD[0, 0, :] * self.__results.ol_flux_linkage_y[0, :] +
                        K_QD[0, 1, :] * self.__results.ol_flux_linkage_y[1, :] +
                        K_QD[0, 2, :] * self.__results.ol_flux_linkage_y[2, :])
        fd_ol = (K_QD[1, 0, :] * self.__results.ol_flux_linkage_y[0, :] +
                        K_QD[1, 1, :] * self.__results.ol_flux_linkage_y[1, :] +
                        K_QD[1, 2, :] * self.__results.ol_flux_linkage_y[2, :])
        
        fq_nl = (K_QD[0, 0, :] * self.__results.nl_flux_linkage_y[0, :] +
                        K_QD[0, 1, :] * self.__results.nl_flux_linkage_y[1, :] +
                        K_QD[0, 2, :] * self.__results.nl_flux_linkage_y[2, :])
        fd_nl = (K_QD[1, 0, :] * self.__results.nl_flux_linkage_y[0, :] +
                        K_QD[1, 1, :] * self.__results.nl_flux_linkage_y[1, :] +
                        K_QD[1, 2, :] * self.__results.nl_flux_linkage_y[2, :])

        self.__results.ol_flux_linkage_y = np.concatenate((self.__results.ol_flux_linkage_y,np.atleast_2d(fq_ol),np.atleast_2d(fd_ol)),axis=0)
        self.__results.nl_flux_linkage_y = np.concatenate((self.__results.nl_flux_linkage_y,np.atleast_2d(fq_nl),np.atleast_2d(fd_nl)),axis=0)
        self.__results.phase_current_y = np.concatenate((self.__results.phase_current_y,np.atleast_2d(iq),np.atleast_2d(id)),axis=0)

        # Only Flux in q-axis
        self.__results.Lmq = fq / ( np.average(iq) )
        self.__results.Lmd = self.__results.Lmq

        self.__results.self_inductance = ( self.__results.self_inductance_ag +
                                           self.__results.self_inductance_slot_leakage +
                                           self.__results.self_inductance_end_winding_leakage )
        self.__results.mutual_inductance = ( self.__results.mutual_inductance_ag +
                                             self.__results.mutual_inductance_slot_leakage +
                                             self.__results.mutual_inductance_end_winding_leakage )

        self.__results.Ld = self.__results.self_inductance - self.__results.mutual_inductance
        self.__results.Lq = self.__results.Ld
 
        Pr = (1 / (2 * MU0)) * (Bg_r * Bg_r-
                                Bg_t * Bg_t)

        # Conversion to MPa
        Pr_fft_nl = np.average( fft(Pr[0:len(pos1),:]) / 1e6 , axis=0 )
        Pr_fft_ol = np.average( fft(Pr[len(posNL)+1:len(posNL)+1+len(pos3),:]) / 1e6 , axis=0 )

        self.__results.pressure_radial_nl = (2.0 / len(Pr_fft_nl)) * np.abs(Pr_fft_nl[0: 5 * self.__gcd + 1])
        self.__results.pressure_radial_ol = (2.0 / len(Pr_fft_ol)) * np.abs(Pr_fft_ol[0: 5 * self.__gcd + 1])


        self.get_self_and_mutual_inductance()


        return True


    def get_results(self):
        return self.__results



