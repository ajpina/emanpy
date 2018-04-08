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

import numpy as np
import fractions


LOG_CRITICAL = 0
LOG_ERROR = 1
LOG_WARN = 2
LOG_INFO = 3
LOG_ALL = 4


MU0 = np.pi*4e-7
DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi
PI = np.pi
EL_STEPS = 90
THETA_e_DEG = np.linspace(0, 360, EL_STEPS)
THETA_e_RAD = THETA_e_DEG * DEG2RAD
PI_2by3 = 2.0 * np.pi / 3.0
K_QD = (2.0 / 3.0) * np.array(([[np.sin(THETA_e_DEG * DEG2RAD),
                                np.sin((THETA_e_DEG - 120) * DEG2RAD),
                                np.sin((THETA_e_DEG + 120) * DEG2RAD)],
                               [np.cos(THETA_e_DEG * DEG2RAD),
                                np.cos((THETA_e_DEG - 120) * DEG2RAD),
                                np.cos((THETA_e_DEG + 120) * DEG2RAD)]]))

def LCM(a,b):
    return abs(a * b) / fractions.gcd(a, b) if a and b else 0

def GCD(a,b):
    return fractions.gcd(a, b)

