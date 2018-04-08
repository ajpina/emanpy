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

import abc


class Analysis(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, analysis_settings):
        self._solver = analysis_settings['solver']

    def set_machine(self, RotatingMachine):
        self._rotating_machine = RotatingMachine

    def solve(self):
        if self._RotatingMachine._type == 'SPM':
            if self._solver == 'subdomain':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                spm = SPMInnerRotorRadialFluxSubDomain()
                success = spm.solve()
                if success:
                    self._results = spm.get_results()
                    return True
                else:
                    return False
            elif self._solver == 'reluctance_network':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                spm = SPMInnerRotorRadialFluxSubDomain()
                success = spm.solve()
                if success:
                    self._results = spm.get_results()
                    return True
                else:
                    return False
            else:
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                spm = SPMInnerRotorRadialFluxSubDomain()
                success = spm.solve()
                if success:
                    self._results = spm.get_results()
                    return True
                else:
                    return False
        else:
            if self._solver == 'subdomain':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                spm = SPMInnerRotorRadialFluxSubDomain()
                success = spm.solve()
                if success:
                    self._results = spm.get_results()
                    return True
                else:
                    return False
            elif self._solver == 'reluctance_network':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                spm = SPMInnerRotorRadialFluxSubDomain()
                success = spm.solve()
                if success:
                    self._results = spm.get_results()
                    return True
                else:
                    return False
            else:
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                spm = SPMInnerRotorRadialFluxSubDomain()
                success = spm.solve()
                if success:
                    self._results = spm.get_results()
                    return True
                else:
                    return False

    def get_results(self):
        return self._results
