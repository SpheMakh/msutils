#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 SKA South Africa
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import subprocess
import os
from setuptools import setup

pkg='MSUtils'

build_root=os.path.dirname(__file__)

def readme():
    with open(os.path.join(build_root, 'README.md')) as f:
        return f.read()

def requirements():
    with open(os.path.join(build_root, 'requirements.txt')) as f:
        return [pname.strip() for pname in f.readlines()]

def src_pkg_dirs(pkg_name):
    mbdir = os.path.join(build_root, pkg_name)
    # Ignore
    pkg_dirs = []
    l = len(mbdir) + len(os.sep)
    exclude = ['docs', '.git', '.svn', 'CMakeFiles']
    for root, dirs, files in os.walk(mbdir, topdown=True):
        # Prune out everything we're not interested in
        # from os.walk's next yield.
        dirs[:] = [d for d in dirs if d not in exclude]

        for d in dirs:
            # Take everything after that ('src...') and
            # append a '/*.*' to it
            pkg_dirs.append(os.path.join(root[l:], d, '*.*'))
    return pkg_dirs

setup(name='msutils',
      version="0.9.5",
      description='Tools for playing with Measurement sets.',
      long_description=readme(),
      url='https://github.com/SpheMakh/msutils',
      classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
#        "License :: GNU GPL v2",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Astronomy"],
      author='Sphesihle Makhathini',
      author_email='sphemakh@gmail.com',
      license='GNU GPL v2',
      packages=[pkg],
      install_requires=requirements(),
      package_data={pkg: src_pkg_dirs(pkg)},
      include_package_data=True,
)
