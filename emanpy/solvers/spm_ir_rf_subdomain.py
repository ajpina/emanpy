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
    and inner rotor structure.
"""

# ==========================================================================
# Program:   spm_ir_rf_subdomain.py
# Author:    ajpina
# Date:      12/23/16
# Version:   0.1.1
#
# Revision History:
#      Date     Version    Author      Description
#  - 12/23/16:  0.1.1      ajpina      Calculate Air Gap Flux Density
#  - 12/19/16:  0.1.1      ajpina      Calculate Cogging, Ripple
#  - 10/01/17:  0.1.1      ajpina      Calculate Inductances
#
# ==========================================================================

import scipy.interpolate as si
import scipy.linalg as sl
import logging
from scipy.fftpack import fft
from emanpy.src.constants import *

from emanpy.results import Result
from .magnetization import get_fs_coeff
from emanpy.analysis import Analysis


class SPMInnerRotorRadialFluxSubDomain(Analysis):

    __n_harms = 180
    __m_harms = 20
    __l_harms = 20
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
    __gcd = 3
    __lcm = 18
    __Sl = 30e-3
    __Rl = 30e-3

    def __init__(self, analysis_settings, spm):
        self.type = "Surface Permanent Magnet Machine - Inner Rotor - Subdomain"
        self.__pp = spm.rotor.pp
        self.__Ns = spm.stator.slots_number
        self.__oRr = spm.rotor.outer_radius
        self.__Rm = spm.rotor.magnets[0].length + spm.rotor.outer_radius
        self.__init_pos = spm.rotor.rotor_position * DEG2RAD
        self.__Rs = spm.stator.inner_radius
        self.__Aso = np.zeros(self.__Ns)
        self.__As = np.zeros(self.__Ns)
        self.__alpha_so = np.zeros(self.__Ns)
        self.__alpha_s = np.zeros(self.__Ns)
        self.__gcd = GCD(self.__Ns, 2 * self.__pp)
        self.__lcm = LCM(self.__Ns, 2 * self.__pp)
        self.__Sl = spm.stator.stack_length
        self.__Rl = spm.rotor.stack_length
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
                self.__Aso[i] = (sl.w0 / spm.stator.inner_radius)
                self.__As[i] = (((sl.w1 + sl.w2) / 2.0) / ((self.__Rt + self.__Rsb) / 2.0))
                self.__alpha_so[i] = sl.so_position * DEG2RAD
                self.__alpha_s[i] = sl.s_position * DEG2RAD
        else:
            self.__Rt = spm.stator.inner_radius + (slot.h0 + slot.h1) / 2.0
            self.__Rsb = self.__Rt + slot.h2 + slot.h3
            for i, sl in enumerate(spm.stator.slots):
                self.__Aso[i] = (sl.w0 / spm.stator.inner_radius)
                self.__As[i] = (((sl.w1 + sl.w2) / 2.0) / ((self.__Rt + self.__Rsb) / 2.0))
                self.__alpha_so[i] = sl.so_position * DEG2RAD
                self.__alpha_s[i] = sl.s_position * DEG2RAD

    def __set_magnet_geometry(self, spm):
        magnet = spm.rotor.magnets[0]
        if magnet.get_type() == 'Arc':
            self.__alpha_p_v = np.zeros(2 * self.__pp)
            self.__delta_v = np.zeros(2 * self.__pp)
            self.__M = np.zeros(2 * self.__pp)
            self.__mur = np.zeros(2 * self.__pp)
            for i, mn in enumerate(spm.rotor.magnets):
                # IT'S MULTIPLIED BY POLE NUMBER
                self.__alpha_p_v[i] = 2 * self.__pp * mn.mean_arc_angle
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
                # IT'S MULTIPLIED BY POLE NUMBER
                self.__alpha_p_v[i] = 2 * self.__pp * mn.mean_arc_angle
                self.__delta_v[i] = mn.deviation
                self.__M[i] = mn.material.Br / MU0
                self.__mur[i] = mn.material.mur

            self.__magnetisation = mn.magnetisation

    def __set_winding_definition(self, spm):
        winding = spm.stator.winding
        if winding.get_type() == 'concentrated':
            self.__Nph = winding.phases
            self.__Ct = winding.turns_coil
            self.__Cp = winding.coil_parallel
            self.__Cs = winding.coil_series
            self.__C = winding.conn_matrix
        else:
            self.__Nph = winding.phases
            self.__Ct = winding.turns_coil
            self.__Cp = winding.coil_parallel
            self.__Cs = winding.coil_series
            self.__C = winding.conn_matrix

        #self.__results.stator_coil_resistance, self.__results.stator_phase_resistance = winding.get_resistances()
        self.__results.td = winding.turns_density(m=3, Ns=self.__Ns)
        self.__results.wf = winding.winding_function(m=3, pp=self.__pp, Ns=self.__Ns)
        self.__results.wh = winding.winding_harmonics(m=3, pp=self.__pp, Ns=self.__Ns)
        self.__results.kw_v = winding.winding_factors(m=3, pp=self.__pp, Ns=self.__Ns)

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
        self.winding_function = analysis_settings['winding'].get('winding_function', False)
        self.winding_harmonics = analysis_settings['winding'].get('winding_harmonics', False)
        self.winding_factors = analysis_settings['winding'].get('winding_factors', False)

    def __assembly_A(self):
        Nharms = np.arange(2, self.__n_harms+1)
        N = len(Nharms)
        Mharms = np.arange(1, self.__m_harms)
        M = len(Mharms)
        Lharms = np.arange(1, self.__l_harms)
        L = len(Lharms)
        MNs = M*self.__Ns
        LNs = L*self.__Ns

        In = np.eye(N, dtype=int)
        Imns = np.eye(MNs, dtype=int)
        Ilns = np.eye(LNs, dtype=int)
        NN = np.zeros((N, N))
        Zero_nn = np.zeros((N, N))
        Zero_mn = np.zeros((MNs, N))
        Zero_nm = np.zeros((N, MNs))
        Zero_ln = np.zeros((LNs, N))
        Zero_nl = np.zeros((N, LNs))
        Zero_ml = np.zeros((MNs, LNs))

        rowVectorL = np.ones((1, L))
        colVectorM = np.ones((M, 1))

        RM = np.diag((self.__oRr / self.__Rm) ** (Nharms - 1))
        MG = np.diag((self.__Rm / self.__Rs) ** (Nharms - 1))
        GS = np.zeros((MNs, MNs))
        F = np.zeros((MNs, MNs))
        SOS = np.zeros((LNs, LNs))
        E = np.zeros((LNs, LNs))
        eta = np.zeros((N, MNs))
        xi = np.zeros((N, MNs))
        sigma = np.zeros((MNs, N))
        tau = np.zeros((MNs, N))
        phi = np.zeros((MNs, LNs))
        V = np.zeros((MNs, LNs))

        for i in range(0,N):
            n = Nharms[i] - 1
            NN[i][i] = n

            for iNs in range(0, self.__Ns):
                Fm = Mharms * PI / self.__Aso[iNs]
                aux_term1 = ((Fm**2.0) - n**2.0)**(-1)
                aux_term2 = ( -1**Mharms )*np.sin( n*self.__alpha_so[iNs] + n*self.__Aso[iNs] / 2.0)
                aux_term3 = np.sin(n*self.__alpha_so[iNs] - n*self.__Aso[iNs] / 2.0) * np.ones(M)
                eta[i,iNs*M:(iNs+1)*M] = (-1.0/PI)*n*aux_term1*(aux_term2 - aux_term3)

                aux_term4 = ( -1**Mharms )*np.cos( n*self.__alpha_so[iNs] + n*self.__Aso[iNs] / 2.0)
                aux_term5 = np.cos(n * self.__alpha_so[iNs] - n * self.__Aso[iNs] / 2.0) * np.ones(M)
                xi[i,iNs*M:(iNs+1)*M] = (1.0 / PI) * n * aux_term1 * (aux_term4 - aux_term5)

                sigma[iNs*M:(iNs+1)*M,i] = (2.0*PI / self.__Aso[iNs])*eta[i,iNs*M:(iNs+1)*M]
                tau[iNs*M:(iNs+1)*M,i] = (2.0*PI / self.__Aso[iNs])*xi[i,iNs*M:(iNs+1)*M]

        for iNs in range(0, self.__Ns):
            Fm = Mharms * PI / self.__Aso[iNs]
            GS[iNs*M:(iNs+1)*M,iNs*M:(iNs+1)*M] = np.diag((self.__Rs / self.__Rt)**Fm)
            F[iNs*M:(iNs+1)*M,iNs*M:(iNs+1)*M] = np.diag(Fm)

            El = Lharms * PI / self.__As[iNs]
            SOS[iNs*L:(iNs+1)*L,iNs*L:(iNs+1)*L] = np.diag((self.__Rt / self.__Rsb) ** El)
            E[iNs*L:(iNs+1)*L,iNs*L:(iNs+1)*L] = np.diag(El)

            aux_term6 = (np.dot(np.transpose([Fm**2]),rowVectorL) - np.dot(colVectorM,[El**2]))**(-1)
            aux_term7 = np.sin(El*(self.__alpha_so[iNs] - self.__alpha_s[iNs] +
                                          (self.__As[iNs] - self.__Aso[iNs]) / 2.0))
            aux_term8 =  np.dot(colVectorM, [aux_term7])
            aux_term9 = np.sin(El * (self.__alpha_so[iNs] - self.__alpha_s[iNs] +
                                     (self.__As[iNs] + self.__Aso[iNs]) / 2.0))
            aux_term10 = np.dot(np.transpose([-1 ** Mharms]), [aux_term9])

            phi[iNs*M:(iNs+1)*M,iNs*L:(iNs+1)*L] = ( (2.0 / self.__As[iNs])*
                                                     np.dot(colVectorM, [El]) * aux_term6 *
                                                     (aux_term8 - aux_term10) )
            V[iNs*M:(iNs+1)*M,iNs*L:(iNs+1)*L] = ( (self.__As[iNs] / self.__Aso[iNs]) *
                                                   phi[iNs*M:(iNs+1)*M,iNs*L:(iNs+1)*L] )

        Eq1_Wm = RM
        Eq1_Xm = -In

        Eq2_Ym = -RM
        Eq2_Zm = In

        Eq7_Wm = In
        Eq7_Xm = -RM
        Eq7_Wg = -self.__mur[0] * MG
        Eq7_Xg = self.__mur[0] * In

        Eq8_Wm = -In
        Eq8_Xm = -RM
        Eq8_Wg = MG
        Eq8_Xg = In

        Eq9_Ym = -In
        Eq9_Zm = RM
        Eq9_Yg = self.__mur[0] * MG
        Eq9_Zg = -self.__mur[0] * In

        Eq10_Ym = -In
        Eq10_Zm = -RM
        Eq10_Yg = MG
        Eq10_Zg = In

        Eq11_Wg = -NN
        Eq11_Xg = np.dot(NN, MG)
        Eq11_Yso = np.dot(np.dot(eta, F), GS)
        Eq11_Zso = -np.dot(eta, F)

        Eq12_Yg = -NN
        Eq12_Zg = np.dot(NN, MG)
        Eq12_Yso = np.dot(np.dot(xi, F), GS)
        Eq12_Zso = -np.dot(xi, F)

        Eq13_Wg = sigma
        Eq13_Xg = np.dot(sigma, MG)
        Eq13_Yg = tau
        Eq13_Zg = np.dot(tau, MG)
        Eq13_Yso = -GS
        Eq13_Zso = -Imns

        Eq14_Yso = Imns
        Eq14_Zso = GS
        Eq14_Zs = -np.dot(V, (SOS**2 + Ilns))

        Eq15_Yso = np.dot( phi.transpose(), F)
        Eq15_Zso = -np.dot( phi.transpose(), np.dot(F, GS))
        Eq15_Zs = -np.dot(E, (SOS**2 - Ilns))

        Eq1 = np.concatenate(
            (Eq1_Wm, Eq1_Xm, Zero_nn, Zero_nn, Zero_nn, Zero_nn, Zero_nn, Zero_nn, Zero_nm, Zero_nm, Zero_nl),
            axis=1)
        Eq2 = np.concatenate(
            (Zero_nn, Zero_nn, Eq2_Ym, Eq2_Zm, Zero_nn, Zero_nn, Zero_nn, Zero_nn, Zero_nm, Zero_nm, Zero_nl),
            axis=1)
        Eq7 = np.concatenate(
            (Eq7_Wm, Eq7_Xm, Zero_nn, Zero_nn, Eq7_Wg, Eq7_Xg, Zero_nn, Zero_nn, Zero_nm, Zero_nm, Zero_nl),
            axis=1)
        Eq8 = np.concatenate(
            (Eq8_Wm, Eq8_Xm, Zero_nn, Zero_nn, Eq8_Wg, Eq8_Xg, Zero_nn, Zero_nn, Zero_nm, Zero_nm, Zero_nl),
            axis=1)
        Eq9 = np.concatenate(
            (Zero_nn, Zero_nn, Eq9_Ym, Eq9_Zm, Zero_nn, Zero_nn, Eq9_Yg, Eq9_Zg, Zero_nm, Zero_nm, Zero_nl),
            axis=1)
        Eq10 = np.concatenate(
            (Zero_nn, Zero_nn, Eq10_Ym, Eq10_Zm, Zero_nn, Zero_nn, Eq10_Yg, Eq10_Zg, Zero_nm, Zero_nm, Zero_nl),
            axis=1)
        Eq11 = np.concatenate(
            (Zero_nn, Zero_nn, Zero_nn, Zero_nn, Eq11_Wg, Eq11_Xg, Zero_nn, Zero_nn, Eq11_Yso, Eq11_Zso, Zero_nl),
            axis=1)
        Eq12 = np.concatenate(
            (Zero_nn, Zero_nn, Zero_nn, Zero_nn, Zero_nn, Zero_nn, Eq12_Yg, Eq12_Zg, Eq12_Yso, Eq12_Zso, Zero_nl),
            axis=1)
        Eq13 = np.concatenate(
            (Zero_mn, Zero_mn, Zero_mn, Zero_mn, Eq13_Wg, Eq13_Xg, Eq13_Yg, Eq13_Zg, Eq13_Yso, Eq13_Zso,Zero_ml),
            axis=1)
        Eq14 = np.concatenate(
            (Zero_mn, Zero_mn, Zero_mn, Zero_mn, Zero_mn, Zero_mn, Zero_mn, Zero_mn, Eq14_Yso, Eq14_Zso, Eq14_Zs),
            axis=1)
        Eq15 = np.concatenate(
            (Zero_ln, Zero_ln, Zero_ln, Zero_ln, Zero_ln, Zero_ln, Zero_ln, Zero_ln, Eq15_Yso, Eq15_Zso, Eq15_Zs),
            axis=1)

        self._A_ = np.concatenate((Eq1, Eq2, Eq7, Eq8, Eq9, Eq10, Eq11, Eq12, Eq13, Eq14, Eq15), axis=0)


    def __assembly_b_nl(self, pos = [0.0]):
        Nharms = np.arange(2, self.__n_harms + 1)
        N = len(Nharms)
        Mharms = np.arange(1, self.__m_harms)
        M = len(Mharms)
        Lharms = np.arange(1, self.__l_harms)
        L = len(Lharms)
        MNs = M * self.__Ns
        LNs = L * self.__Ns

        Zero_n = np.zeros((N, len(pos)))
        Zero_m = np.zeros((MNs, len(pos)))
        Zero_l = np.zeros((LNs, len(pos)))
        Eq1 = np.zeros((N, len(pos)))
        Eq2 = np.zeros((N, len(pos)))
        Eq7 = np.zeros((N, len(pos)))
        Eq8 = np.zeros((N, len(pos)))
        Eq9 = np.zeros((N, len(pos)))
        Eq10 = np.zeros((N, len(pos)))

        for p in range(0, len(pos)):
            self._Mrn, self._Mtn = get_fs_coeff(self.__magnetisation, self.__n_harms, self.__pp, pos[p], self.__M,
                         self.__alpha_p_v, self.__delta_v, 0)

            Mrcn = 2.0 * np.real(self._Mrn)
            Mrsn = -2.0 * np.imag(self._Mrn)
            Mtcn = 2.0 * np.real(self._Mtn)
            Mtsn = -2.0 * np.imag(self._Mtn)

            for i in range(0, N):
                n = Nharms[i]-1
                if (np.abs(n)==1):
                    factor_Rr = np.log(self.__oRr) / 2.0
                    factor_Rm = np.log(self.__Rm) / 2.0
                else:
                    factor_Rr = factor_Rm = 1.0 / (n**2 - 1.0)

                Eq1[i,p] = MU0 * self.__oRr * factor_Rr * (Mrsn[i+1] - n * Mtcn[i+1])
                Eq2[i,p] = MU0 * self.__oRr * factor_Rr * (Mrcn[i+1] + n * Mtsn[i+1])
                Eq7[i,p] = MU0 * self.__Rm * factor_Rm * (Mrsn[i+1] - n * Mtcn[i+1])
                Eq8[i,p] = MU0 * self.__Rm * factor_Rm * (Mtcn[i+1] - n * Mrsn[i+1])
                Eq9[i,p] = MU0 * self.__Rm * factor_Rm * (Mrcn[i+1] + n * Mtsn[i+1])
                Eq10[i,p] = MU0 * self.__Rm * factor_Rm * (Mtsn[i+1] + n * Mrcn[i+1])

        self._b_ = np.concatenate((Eq1, Eq2, Eq7, Eq8, Eq9, Eq10, Zero_n, Zero_n, Zero_m,
                                      Zero_m, Zero_l), axis=0)


    def __assembly_b_ol(self, pos = [0.0], samePos = True, I = [0.0, 0.0, 0.0]):
        if(samePos==True):
            b = np.tile(self._b_, 2)
            self._b_ = b
        else:
            Nharms = np.arange(2, self.__n_harms + 1)
            N = len(Nharms)
            Mharms = np.arange(1, self.__m_harms)
            M = len(Mharms)
            Lharms = np.arange(1, self.__l_harms)
            L = len(Lharms)
            MNs = M * self.__Ns
            LNs = L * self.__Ns

            Zero_n = np.zeros((N, len(pos)))
            Zero_m = np.zeros((MNs, len(pos)))
            Zero_l = np.zeros((LNs, len(pos)))
            rowVectorL = np.ones((1, L))
            colVectorM = np.ones((M, 1))
            Ilns = np.eye(LNs)
            Eq1 = np.zeros((N, len(pos)))
            Eq2 = np.zeros((N, len(pos)))
            Eq7 = np.zeros((N, len(pos)))
            Eq8 = np.zeros((N, len(pos)))
            Eq9 = np.zeros((N, len(pos)))
            Eq10 = np.zeros((N, len(pos)))
            Eq11 = np.zeros((N, len(pos)))
            Eq12 = np.zeros((N, len(pos)))
            Eq14 = np.zeros((MNs, len(pos)))
            Eq15 = np.zeros((LNs, len(pos)))

            eta_o = np.zeros((N, self.__Ns))
            xi_o = np.zeros((N, self.__Ns))
            SOS = np.zeros((LNs, LNs))
            E = np.zeros((LNs, LNs))
            phi = np.zeros((MNs, LNs))
            V = np.zeros((MNs, LNs))
            V_o = np.zeros((LNs, 1))
            Jn = np.zeros((LNs, 1))

            C = self.__C

            for p in range(0, len(pos)):
                self._Mrn, self._Mtn = get_fs_coeff(self.__magnetisation, self.__n_harms, self.__pp, pos[p], self.__M,
                                                    self.__alpha_p_v, self.__delta_v, 0)

                Mrcn = 2.0 * np.real(self._Mrn)
                Mrsn = -2.0 * np.imag(self._Mrn)
                Mtcn = 2.0 * np.real(self._Mtn)
                Mtsn = -2.0 * np.imag(self._Mtn)

                for i in range(0, N):
                    n = Nharms[i] - 1
                    if (np.abs(n) == 1):
                        factor_Rr = np.log(self.__oRr) / 2.0
                        factor_Rm = np.log(self.__Rm) / 2.0
                    else:
                        factor_Rr = factor_Rm = 1.0 / (n ** 2 - 1.0)

                    Eq1[i, p] = MU0 * self.__oRr * factor_Rr * (Mrsn[i + 1] - n * Mtcn[i + 1])
                    Eq2[i, p] = MU0 * self.__oRr * factor_Rr * (Mrcn[i + 1] + n * Mtsn[i + 1])
                    Eq7[i, p] = MU0 * self.__Rm * factor_Rm * (Mrsn[i + 1] - n * Mtcn[i + 1])
                    Eq8[i, p] = MU0 * self.__Rm * factor_Rm * (Mtcn[i + 1] - n * Mrsn[i + 1])
                    Eq9[i, p] = MU0 * self.__Rm * factor_Rm * (Mrcn[i + 1] + n * Mtsn[i + 1])
                    Eq10[i, p] = MU0 * self.__Rm * factor_Rm * (Mtsn[i + 1] + n * Mrcn[i + 1])

                    eta_o[i,:] = (2.0 / (n * PI)) * (np.cos(n * self.__alpha_so) *
                                                        np.sin(n * self.__Aso / 2.0))
                    xi_o[i, :] = (2.0 / (n * PI)) * (np.sin(n * self.__alpha_so) *
                                                         np.sin(n * self.__Aso / 2.0))

                ###### This is valid only for concentrated winding and 3 phases
                I1_p_I2 = np.dot(I[p,:], C[0:self.__Nph, 0:self.__Ns] + C[self.__Nph:, 0:self.__Ns])
                I1_m_I2 = np.dot(I[p,:], C[0:self.__Nph, 0:self.__Ns] - C[self.__Nph:, 0:self.__Ns])
                Jo = I1_p_I2 / self.__Aso

                for iNs in range(0, self.__Ns):
                    Fm = Mharms * PI / self.__Aso[iNs]
                    El = Lharms * PI / self.__As[iNs]
                    SOS[iNs * L:(iNs + 1) * L, iNs * L:(iNs + 1) * L] = np.diag((self.__Rt / self.__Rsb) ** El)
                    E[iNs * L:(iNs + 1) * L, iNs * L:(iNs + 1) * L] = np.diag(El)
                    aux_term6 = (np.dot(np.transpose([Fm ** 2]), rowVectorL) - np.dot(colVectorM, [El ** 2])) ** (-1)
                    aux_term7 = np.sin(El * (self.__alpha_so[iNs] - self.__alpha_s[iNs] +
                                             (self.__As[iNs] - self.__Aso[iNs]) / 2.0))
                    aux_term8 = np.dot(colVectorM, [aux_term7])
                    aux_term9 = np.sin(El * (self.__alpha_so[iNs] - self.__alpha_s[iNs] +
                                             (self.__As[iNs] + self.__Aso[iNs]) / 2.0))
                    aux_term10 = np.dot(np.transpose([-1 ** Mharms]), [aux_term9])

                    phi[iNs * M:(iNs + 1) * M, iNs * L:(iNs + 1) * L] = ((2.0 / self.__As[iNs]) *
                                                                         np.dot(colVectorM, [El]) * aux_term6 *
                                                                         (aux_term8 - aux_term10))
                    V[iNs * M:(iNs + 1) * M, iNs * L:(iNs + 1) * L] = ((self.__As[iNs] / self.__Aso[iNs]) *
                                                                       phi[iNs * M:(iNs + 1) * M,
                                                                       iNs * L:(iNs + 1) * L])
                    aux_term11 = 4.0 * ((Lharms * PI) ** (-1))
                    aux_term12 = np.sin(Lharms * PI * self.__Aso[iNs] / (self.__As[iNs] * 2.0))
                    aux_term13 = np.cos(Lharms * PI / 2.0 + Lharms * PI *
                                        (self.__alpha_so[iNs] - self.__alpha_s[iNs]) / 2.0)
                    phi_o = aux_term11 * aux_term12 * aux_term13
                    V_o[iNs * L:(iNs + 1) * L,:] = np.transpose([phi_o]) * Jo[iNs]
                    aux_term14 = np.sin(Lharms * PI / 2.0) / (Lharms * PI * self.__As[iNs])
                    Jn[iNs * L:(iNs + 1) * L, :] = (8.0 * I1_m_I2[iNs] /
                                                    (self.__Rsb ** 2 - self.__Rt ** 2 ))*\
                                                   np.transpose([aux_term14])

                Eq11[:,p] = -MU0 * np.dot(eta_o, Jo)
                Eq12[:,p] = -MU0 * np.dot(xi_o, Jo)
                tmp1 = np.diag((np.diag(E ** 2 - 4.0 * Ilns)) ** (-1))
                tmp2 = Ilns * (self.__Rt ** 2) - (2.0 * (self.__Rsb ** 2) *
                                                 np.dot(np.diag(E.diagonal() ** (-1)), SOS))
                Eq14[:,p] = MU0 * V.dot(tmp1).dot(tmp2).dot(Jn).flatten()
                tmp3 = Ilns * (self.__Rt ** 2) - (self.__Rsb ** 2) * SOS
                Eq15[:,p] = (-MU0 * V_o + 2.0 * MU0 * tmp1.dot(tmp3).dot(Jn)).flatten()



            b = np.concatenate((Eq1, Eq2, Eq7, Eq8, Eq9, Eq10, Eq11, Eq12, Zero_m,
                                Eq14, Eq15), axis=0)
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
        Nharms = np.arange(2, self.__n_harms + 1)
        N = len(Nharms)
        Wg = unk[4 * N:5 * N, : ]
        Xg = unk[5 * N:6 * N, : ]
        Yg = unk[6 * N:7 * N, : ]
        Zg = unk[7 * N:8 * N, : ]
        g = self.__Rs - self.__Rm
        # Flux Density Close to stator for better prediction of flux linkage and
        # Radial forces
        r = self.__Rs - g / 2.0
        Bg_r = np.zeros((len(posNL)+len(posOL),len(psi)))
        Bg_t = np.zeros((len(posNL)+len(posOL),len(psi)))
        for p in range(0, len(posNL)+len(posOL)):
            for i in range(0, N):
                n = Nharms[i] - 1
                term1 = -(n / r) * ((Wg[i,p] * (r / self.__Rs) ** n) + (Xg[i,p] * (r / self.__Rm) ** (-n)))
                term2 = (n / r) * ((Yg[i,p] * (r / self.__Rs) ** n) + (Zg[i,p] * (r / self.__Rm) ** (-n)))
                term3 = -(n / r) * ((Wg[i,p] * (r / self.__Rs) ** n) - (Xg[i,p] * (r / self.__Rm) ** (-n)))
                term4 = -(n / r) * ((Yg[i,p] * (r / self.__Rs) ** n) - (Zg[i,p] * (r / self.__Rm) ** (-n)))
                Bg_r[p,:] = Bg_r[p,:] + term1 * np.sin(n * psi) + term2 * np.cos(n * psi)
                Bg_t[p,:] = Bg_t[p,:] + term3 * np.cos(n * psi) + term4 * np.sin(n * psi)

        return Bg_r, Bg_t



    def __get_settings_cogging(self):
        final_pos = 360 / self.__lcm
        steps = 7
        posVector = np.linspace(0, final_pos, steps) * DEG2RAD  + self.__init_pos
        return posVector


    def __get_settings_ripple(self):
        final_pos = 360 / self.__lcm
        steps = 7
        posVector = np.linspace(0, final_pos, steps) * DEG2RAD + self.__init_pos
        thetaElec = self.__pp * np.linspace(0, final_pos, steps) * DEG2RAD
        Is = self.load_current / self.__Cp
        IVector = np.transpose([Is * np.sin(thetaElec), Is * np.sin(thetaElec - 2.0*PI/3.0),
                                Is * np.sin(thetaElec + 2.0*PI/3.0)])
        return posVector, IVector

    def __get_settings_static_torque(self,  pos=np.array([]), I=np.array([])):
        final_pos = 360 / self.__pp
        steps = 7
        posVector = np.linspace(0, final_pos, steps) * DEG2RAD + self.__init_pos
        Is = self.load_current / self.__Cp
        IVector = np.transpose([Is * np.ones(steps), -0.5 * Is * np.ones(steps),
                                -0.5 * Is * np.ones(steps)])
        return posVector, IVector

    def __get_settings_inductance_calc(self):
        posVector = np.array([self.__init_pos])
        IVector = np.array([self.load_current, 0.0, 0.0]) / self.__Cp
        return posVector, IVector

    def __get_settings_bemf(self, pos=[]):
        init_pos = pos[-1]
        steps = 3
        final_pos = 2*PI / (3 * self.__pp) + self.__init_pos        # Max span 120 electric degrees
        step_size = np.abs(init_pos - final_pos) / steps
        posVector = np.linspace(init_pos + step_size, final_pos - step_size, steps)
        return posVector, self.noload_speed

    def solve(self):
        if self.cogging:
            pos1 = self.__get_settings_cogging()
        if self.noload_bemf:
            pos2, wr_nl = self.__get_settings_bemf(pos=pos1)

        pos3, I1 = self.__get_settings_ripple()
        pos4, I2 = self.__get_settings_static_torque()
        pos5, I3 = self.__get_settings_inductance_calc()


        posNL = np.hstack((pos1, pos2))

        posOL = np.hstack((pos3, pos4, pos5))
        I = np.vstack((I1, I2, I3))


        Bg_r, Bg_t = self.get_ag_flux_density(posNL=posNL, posOL=posOL, samePos=False,
                                              current=I)
        g = self.__Rs - self.__Rm
        r = (self.__Rs - g / 7.0)
        Lstk = (self.__Sl + self.__Rl) / 2.0
        torque = (Lstk * (r**2) / MU0) * np.trapz(Bg_r * Bg_t, axis=1, dx=PI/180)
        pos_nl_cs = (pos1-self.__init_pos)*RAD2DEG
        cogging_cs = si.CubicSpline(pos_nl_cs,
                                    torque[0:len(pos1)],
                                    extrapolate='periodic')
        cogging_acs = np.linspace(0,360/self.__pp,180)

        pos_ol_cs = (pos3 - self.__init_pos) * RAD2DEG
        ripple_cs = si.CubicSpline(pos_ol_cs,
                                    torque[len(posNL):len(posNL)+len(pos3)],
                                    extrapolate='periodic')
        ripple_acs = np.linspace(0, 360 / self.__pp, 180)

        pos_st_cs = (pos4 - self.__init_pos) * RAD2DEG
        staticTq_cs = si.CubicSpline(pos_st_cs,
                                   torque[len(posNL)+len(pos3):len(posNL)+len(posOL)-len(pos5)])
        staticTq_acs = np.linspace(0, 360 / self.__pp, 45)


        self.__results.cogging_torque_x = cogging_acs
        self.__results.cogging_torque_y = cogging_cs(cogging_acs)
        self.__results.torque_ripple_x = ripple_acs
        self.__results.torque_ripple_y = ripple_cs(ripple_acs)
        self.__results.static_torque_x = staticTq_acs
        self.__results.static_torque_y = staticTq_cs(staticTq_acs)
        self.__results.nl_Bg_theta = np.linspace(0, 360, 360)
        self.__results.nl_Bg_r = Bg_r[0,:]
        self.__results.nl_Bg_t = Bg_t[0,:]
        self.__results.ol_Bg_r = Bg_r[len(posNL), :]
        self.__results.ol_Bg_t = Bg_t[len(posNL), :]


        nrows = len(posNL) + len(pos3)

        fl = np.zeros((nrows,3))

        for i in range(0, nrows):
            fl[i] = (Lstk * r / self.__Cp) * np.trapz(Bg_r[i,:] * self.__results.wf,  dx=PI/180)

        fl_a = np.concatenate((fl[0:len(posNL), 0], fl[0:len(posNL), 2], fl[0:len(posNL), 1], fl[0:len(posNL), 0],
                               fl[len(posNL):len(posNL)+len(pos3), 0], fl[len(posNL):len(posNL)+len(pos3), 2],
                               fl[len(posNL):len(posNL)+len(pos3), 1], fl[len(posNL):len(posNL)+len(pos3), 0] ),
                              axis=0)
        fl_b = np.concatenate((fl[0:len(posNL), 1], fl[0:len(posNL), 0], fl[0:len(posNL), 2], fl[0:len(posNL), 1],
                               fl[len(posNL):len(posNL) + len(pos3), 1], fl[len(posNL):len(posNL) + len(pos3), 0],
                               fl[len(posNL):len(posNL) + len(pos3), 2], fl[len(posNL):len(posNL) + len(pos3), 1]),
                              axis=0)
        fl_c = np.concatenate((fl[0:len(posNL), 2], fl[0:len(posNL), 1], fl[0:len(posNL), 0], fl[0:len(posNL), 2],
                               fl[len(posNL):len(posNL)+len(pos3), 2], fl[len(posNL):len(posNL)+len(pos3), 1],
                               fl[len(posNL):len(posNL)+len(pos3), 0], fl[len(posNL):len(posNL)+len(pos3), 2] ),
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

        self.__results.nl_flux_linkage_y = np.array( [fl_a_nl_cs(THETA_e_DEG),
                                         fl_b_nl_cs(THETA_e_DEG),
                                         fl_c_nl_cs(THETA_e_DEG)] )
        self.__results.nl_flux_linkage_x = THETA_e_DEG


        # Extract fundamental of flux linkage due to magnet
        fft_magnet_flux = fft(self.__results.nl_flux_linkage_y[0,:])
        self.__results.magnet_flux = (2.0 / len(self.__results.nl_flux_linkage_y[0,:])) * np.abs(fft_magnet_flux[1])

        bemf_a = np.zeros( EL_STEPS )
        bemf_b = np.zeros(EL_STEPS)
        bemf_c = np.zeros(EL_STEPS)
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

        self.__results.bemf_y = np.array( [ np.real( bemf_a ), np.real( bemf_b ), np.real( bemf_c )] )
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

        self.__results.ol_flux_linkage_y = np.array( [fl_a_ol_cs(THETA_e_DEG),
                                         fl_b_ol_cs(THETA_e_DEG),
                                         fl_c_ol_cs(THETA_e_DEG)] )
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


        iq = np.average( K_QD[0, 0, :] * self.__results.phase_current_y[0, :] +
                         K_QD[0, 1, :] * self.__results.phase_current_y[1, :] +
                         K_QD[0, 2, :] * self.__results.phase_current_y[2, :])
        id = np.average( K_QD[1, 0, :] * self.__results.phase_current_y[0, :] +
                        K_QD[1, 1, :] * self.__results.phase_current_y[1, :] +
                        K_QD[1, 2, :] * self.__results.phase_current_y[2, :])

        fl_abc_noMagnet = - self.__results.nl_flux_linkage_y + self.__results.ol_flux_linkage_y

        fq = np.average(K_QD[0, 0, :] * fl_abc_noMagnet[0, :] +
                        K_QD[0, 1, :] * fl_abc_noMagnet[1, :] +
                        K_QD[0, 2, :] * fl_abc_noMagnet[2, :])
        fd = np.average(K_QD[1, 0, :] * fl_abc_noMagnet[0, :] +
                        K_QD[1, 1, :] * fl_abc_noMagnet[1, :] +
                        K_QD[1, 2, :] * fl_abc_noMagnet[2, :])

        # Only Flux in q-axis
        self.__results.Lmq = fq / ( iq )
        self.__results.Lmd = self.__results.Lmq

        self.__results.self_inductance = ( self.__results.self_inductance_ag +
                                           self.__results.self_inductance_slot_leakage +
                                           self.__results.self_inductance_end_winding_leakage )
        self.__results.mutual_inductance = ( self.__results.mutual_inductance_ag +
                                             self.__results.mutual_inductance_slot_leakage +
                                             self.__results.mutual_inductance_end_winding_leakage )

        self.__results.Ld = self.__results.self_inductance - self.__results.mutual_inductance
        self.__results.Lq = self.__results.Ld
        #Sb = self.__oSr - (self.__iSr + (self._stator._slots[0]).get_slot_total_height())
        #res.end_winding_leakage = self._stator._winding.end_winding_permeance(Sb)


        # Calculate average pressure on different rotor positions
        #Pr = (1 / (2*MU0)) * (Bg_r[0:len(posNL)+1:len(posNL),:] * Bg_r[0:len(posNL)+1:len(posNL),:] -
        #                      Bg_t[0:len(posNL)+1:len(posNL),:] * Bg_t[0:len(posNL)+1:len(posNL),:])
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



