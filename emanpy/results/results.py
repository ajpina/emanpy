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

__author__ = 'ajpina'

class Result:

    cogging_torque_x = []
    cogging_torque_y = []
    torque_ripple_x = []
    torque_ripple_y = []
    static_torque_x = []
    static_torque_y = []
    nl_Bg_r = []
    nl_Bg_t = []
    ol_Bg_r = []
    ol_Bg_t = []
    nl_Bg_theta = []
    stator_phase_resistance = 0.0
    stator_coil_resistance = 0.0
    wf = []
    td = []
    wh = []
    kw_v = []
    nl_flux_linkage_x = []
    nl_flux_linkage_y = []
    ol_flux_linkage_x = []
    ol_flux_linkage_y = []
    bemf_y = []
    bemf_x = []
    phase_current_x = []
    phase_current_y = []
    self_inductance = 0.0
    mutual_inductance = 0.0
    self_inductance_ag = 0.0
    mutual_inductance_ag = 0.0
    self_inductance_slot_leakage = 0.0
    self_inductance_end_winding_leakage = 0.0
    mutual_inductance_slot_leakage = 0.0
    mutual_inductance_end_winding_leakage = 0.0
    end_winding_leakage = 0.0
    Lmd = 0.0
    Lmq = 0.0
    magnet_flux = 0.0
    pressure_radial_nl = []
    pressure_radial_ol = []

    def __init__(self):
            pass

