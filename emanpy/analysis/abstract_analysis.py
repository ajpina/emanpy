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



class Analysis(object):

    def __init__(self, analysis_settings, rotating_machine):
        self.solver_str = analysis_settings['solver']
        if rotating_machine.get_type() == 'SPM':
            if self.solver_str == 'subdomain':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                self.solver_instance = SPMInnerRotorRadialFluxSubDomain(analysis_settings, rotating_machine)
            elif self.solver_str == 'reluctance_network':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                self.solver_instance = SPMInnerRotorRadialFluxSubDomain(analysis_settings, rotating_machine)
            else:
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                self.solver_instance = SPMInnerRotorRadialFluxSubDomain(analysis_settings, rotating_machine)
        else:
            if self.solver_str == 'subdomain':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                self.solver_instance = SPMInnerRotorRadialFluxSubDomain(analysis_settings, rotating_machine)
            elif self.solver_str == 'reluctance_network':
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                self.solver_instance = SPMInnerRotorRadialFluxSubDomain(analysis_settings, rotating_machine)
            else:
                from emanpy.solvers import SPMInnerRotorRadialFluxSubDomain
                self.solver_instance = SPMInnerRotorRadialFluxSubDomain(analysis_settings, rotating_machine)

    def solve(self):
        return self.solver_instance.solve()

    def get_results(self):
        return self.solver_instance.get_results()




