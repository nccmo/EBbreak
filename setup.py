#!/usr/bin/env python

from distutils.core import setup

setup(name='bp_utils',
      version='0.1.0',
      description='Python tools for various analysis for breakpoints from cancer genome sequencing data.',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/bp_utils',
      package_dir = {'': 'lib'},
      packages=['bp_utils'],
      scripts=['bp_utils'],
      license='GPL-3'
     )

