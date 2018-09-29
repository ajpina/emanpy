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
    Creates a Magnetization Vector (Radial, Halback or Parallel) 
    and computes their complex Fourier series coefficients.
"""

# ==========================================================================
## Program:   magnetization.py
## Author:    ajpina
## Date:      12/23/16
## Version:   0.1.1
##
## Revision History:
##      Date     Version  Author    Description
##  - 12/23/16:  0.1.1              Calculate Radial, Halbach Magnetization
##  - 12/19/16:  0.1.1              Calculate Parallel Magnetization
##
# ==========================================================================

import numpy as np

def get_fs_coeff(type_mzt,n_harm,pp,rotor_pos,m_v,alpha_p_v,delta_v,dp):
    """Get Fourier series coefficients of Magnetization Vector

    Compute Fourier series coefficients of Magnetization Vector, 
    the implementation is done trhough complex Fourier Series

    Args:
        type_mzt:   Type of magnetization:  (0) Radial
                                            (1) Halbach
                                            (2) Parallel
        n_harm:     Number of harmonics to computer
        pp:         Rotor pole pairs
        rotor_pos:  Rotor position
        m_v:        Magnitude magnetization per magnet piece
        alpha_p_v:  Arc length ratio per magnet piece
        delta_v:    Tangential shift per magnet piece
        dp:         Display plot [true,false]

    Returns:
        Mrn:        Fourier coefficients of radial magnetization
        Mtn:        Fourier coefficients of tangential magnetization
    """


    Mrn = np.zeros(n_harm, dtype=np.complex)
    Mtn = np.zeros(n_harm, dtype=np.complex)
    theta_r = -rotor_pos
    RCo, RCn, TCo, TCn = 0, 0, 0, 0
    for v in range(0,2*pp,1):
        polarity = ((-1)**(v+2))*m_v[v]
        theta_init =(np.pi/(2*pp))*(2*v+1-alpha_p_v[v])+delta_v[v]
        theta_end =(np.pi/(2*pp))*(2*v+1+alpha_p_v[v])+delta_v[v]
        psi_init = (np.pi/(2*pp))*(2*v+1)+delta_v[v]

        if type_mzt == "Radial":
            RCo = (1/2/np.pi)*polarity*2*alpha_p_v[v]
            TCo = 0
        elif type_mzt == "Halbach":
            RCo = polarity*(1j/(4*np.pi*pp))*\
                (   np.exp(1j*pp*psi_init)*(np.exp(-1j*pp*theta_end)-\
                    np.exp(-1j*pp*theta_init))-np.exp(-1j*pp*psi_init)*\
                    (np.exp(1j*pp*theta_end)-np.exp(1j*pp*theta_init))  )

            TCo = polarity*(1/(4*np.pi*pp))*\
                (   np.exp(1j*pp*psi_init)*(np.exp(-1j*pp*theta_end)-\
                    np.exp(-1j*pp*theta_init))+np.exp(-1j*pp*psi_init)*\
                    (np.exp(1j*pp*theta_end)-np.exp(1j*pp*theta_init))  )
        elif type_mzt == "Parallel":
            RCo = polarity*(1j/(4*np.pi))*\
                (   np.exp(1j*psi_init)*(np.exp(-1j*theta_end)-\
                    np.exp(-1j*theta_init))-np.exp(-1j*psi_init)*\
                    (np.exp(1j*theta_end)-np.exp(1j*theta_init))    )
            TCo = polarity*(1/(4*np.pi))*\
                (   np.exp(1j*psi_init)*(np.exp(-1j*theta_end)-\
                    np.exp(-1j*theta_init))+np.exp(-1j*psi_init)*\
                    (np.exp(1j*theta_end)-np.exp(1j*theta_init))    )
        else:
            RCo = 0
            TCo = 0

        Mrn[0] = Mrn[0] + RCo
        Mtn[0] = Mtn[0] + TCo

        for n in range(1,n_harm,1):
            if type_mzt == "Radial":
                RCn = (1j/(2*np.pi*n))*polarity*\
                    (np.exp(-1j*n*theta_end)-np.exp(-1j*n*theta_init))*\
                    np.exp(1j*n*(theta_r))
                TCn = 0
            elif type_mzt == "Halbach":
                if np.abs(n) != pp:
                    RCn = polarity*(1j/(4*np.pi))*\
                        (   np.exp(1j*pp*psi_init)*(np.exp(-1j*(pp+n)*theta_end)-\
                            np.exp(-1j*(pp+n)*theta_init))/(pp+n)-np.exp(-1j*pp*psi_init)*\
                            (np.exp(1j*(pp-n)*theta_end)-np.exp(1j*(pp-n)*theta_init))/(pp-n)    )*\
                        np.exp(1j*n*(theta_r))
                    TCn = polarity*(1/(4*np.pi))*\
                        (   np.exp(1j*pp*psi_init)*(np.exp(-1j*(pp+n)*theta_end)-\
                            np.exp(-1j*(pp+n)*theta_init))/(pp+n)+np.exp(-1j*pp*psi_init)*\
                            (np.exp(1j*(pp-n)*theta_end)-np.exp(1j*(pp-n)*theta_init))/(pp-n)   )*\
                        np.exp(1j*n*(theta_r))
                elif n == pp:
                    RCn = polarity*(1/(4*np.pi))*\
                        (   np.exp(-1j*pp*psi_init)*(theta_end-theta_init)+1j*np.exp(1j*pp*psi_init)*\
                            (np.exp(-1j*2*pp*theta_end)-np.exp(-1j*2*pp*theta_init))/(2*pp)     )*\
                        np.exp(1j*n*(theta_r))
                    TCn = polarity*(1j/(4*np.pi))*\
                        (   np.exp(-1j*pp*psi_init)*(theta_end-theta_init)-1j*np.exp(1j*pp*psi_init)*\
                            (np.exp(-1j*2*pp*theta_end)-np.exp(-1j*2*pp*theta_init))/(2*pp)     )*\
                        np.exp(1j*n*(theta_r))
                elif n == -pp:
                    RCn = polarity*(1/(4*np.pi))*\
                        (   np.exp(1j*pp*psi_init)*(theta_end-theta_init)-1j*np.exp(-1j*pp*psi_init)*\
                            (np.exp(1j*2*pp*theta_end)-np.exp(1j*2*pp*theta_init))/(2*pp)       )*\
                        np.exp(1j*n*(theta_r))
                    TCn = -polarity*(1j/(4*np.pi))*\
                        (   np.exp(1j*pp*psi_init)*(theta_end-theta_init)+1j*np.exp(-1j*pp*psi_init)*\
                            (np.exp(1j*2*pp*theta_end)-np.exp(1j*2*pp*theta_init))/(2*pp)       )*\
                        np.exp(1j*n*(theta_r))
            elif type_mzt == "Parallel":
                if np.abs(n) != 1: 
                    RCn = polarity*(1j/(4*np.pi))*\
                        (   np.exp(1j*psi_init)*(np.exp(-1j*(1+n)*theta_end)-\
                            np.exp(-1j*(1+n)*theta_init))/(1+n)-np.exp(-1j*psi_init)*\
                            (np.exp(1j*(1-n)*theta_end)-np.exp(1j*(1-n)*theta_init))/(1-n)      )*\
                        np.exp(1j*n*(theta_r))
                    TCn = polarity*(1/(4*np.pi))*\
                        (   np.exp(1j*psi_init)*(np.exp(-1j*(1+n)*theta_end)-\
                            np.exp(-1j*(1+n)*theta_init))/(1+n)+np.exp(-1j*psi_init)*\
                            (np.exp(1j*(1-n)*theta_end)-np.exp(1j*(1-n)*theta_init))/(1-n)      )*\
                        np.exp(1j*n*(theta_r))
                elif n == 1: 
                    RCn = polarity*(1/(4*np.pi))*\
                        (   np.exp(-1j*psi_init)*(theta_end-theta_init)+1j*np.exp(1j*psi_init)*\
                            (np.exp(-1j*2*theta_end)-np.exp(-1j*2*theta_init))/(2)  )*\
                        np.exp(1j*n*(theta_r))
                    TCn = polarity*(1j/(4*np.pi))*\
                        (   np.exp(-1j*psi_init)*(theta_end-theta_init)-1j*np.exp(1j*psi_init)*\
                            (np.exp(-1j*2*theta_end)-np.exp(-1j*2*theta_init))/(2)  )*\
                        np.exp(1j*n*(theta_r))
                elif n == -1: 
                    RCn = polarity*(1/(4*np.pi))*\
                        (   np.exp(1j*psi_init)*(theta_end-theta_init)-1j*np.exp(-1j*psi_init)*\
                            (np.exp(1j*2*theta_end)-np.exp(1j*2*theta_init))/(2)    )*\
                        np.exp(1j*n*(theta_r))
                    TCn = -polarity*(1j/(4*np.pi))*\
                        (   np.exp(1j*psi_init)*(theta_end-theta_init)+1j*np.exp(-1j*psi_init)*\
                            (np.exp(1j*2*theta_end)-np.exp(1j*2*theta_init))/(2)    )*\
                        np.exp(1j*n*(theta_r))
            else:
                RCn = 0
                TCn = 0

            Mrn[n] = Mrn[n] + RCn
            Mtn[n] = Mtn[n] + TCn

            
    if dp == 1:
        theta = np.linspace(0,2*np.pi,360)
        Mr = Mrn[0]*np.ones(theta.size)
        Mt = Mtn[0]*np.ones(theta.size)
        for n in range(1,n_harm,1):
            Mr = Mr + Mrn[n]*np.exp(1j*n*theta) + np.conjugate(Mrn[n])*np.exp(-1j*n*theta)
            Mt = Mt + Mtn[n]*np.exp(1j*n*theta) + np.conjugate(Mtn[n])*np.exp(-1j*n*theta)

        import matplotlib.pyplot as plt
        plt.figure(1)
        #plt.ion()
        plt.subplot(221)
        plt.plot(theta,Mr.real)
        plt.subplot(222)
        plt.plot(theta,Mt.real)
        plt.subplot(223)
        plt.plot(theta,Mr.imag)
        plt.subplot(224)
        plt.plot(theta,Mt.imag)
        plt.show()

    return Mrn, Mtn








