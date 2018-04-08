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

class Study:


    def __init__(self, settings = []):
        if not settings:
            pass
        else:
            self._noLoadSpeed = settings['noLoad']['speed']
            self._rippleSpeed = settings['ripple']['speed']
            self._rippleCurrent = settings['ripple']['current']
            self._rippleVoltage = settings['ripple']['voltage']
            self._rippleGamma = settings['ripple']['gamma']