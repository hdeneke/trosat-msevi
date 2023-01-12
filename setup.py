#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of the TROPOS Satellite Suite (TROOSAT) developed 
# by the satellite group at TROPOS.
#
# TROSAT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Authors/Copyright(2012-2020):
# -Hartwig Deneke (deneke@tropos.de)

import glob
from setuptools import setup
import versioneer

setup(
    name         = 'trosat-msevi',
    version      = versioneer.get_version(),
    cmdclass     = versioneer.get_cmdclass(),
    author       = 'Hartwig Deneke', 
    author_email = 'deneke@tropos.de',
    scripts      = glob.glob('bin/*'),
    #scripts     =  ['bin/ms15_std2hres.py'],
    packages     = [ 'trosat', 'trosat.msevi', 'trosat.msevi.io' ],
    package_dir  = { '': 'src'},
    package_data = { 'trosat.msevi': ['share/*','lib/*']}
)
