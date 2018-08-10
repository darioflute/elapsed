#!/usr/bin/env python
"""Setup script for installing elapsed."""

from setuptools import setup 

config = {
    'name': 'elapsed',
    'version': '0.0.1',
    'description': 'GUI',
    'author': 'Dario Fadda',
    'author email': 'darioflute@gmail.com',
    'url': 'https://github.com/darioflute/elapsed.git',
    'download_url': 'https://github.com/darioflute/elapsed',
    'license': 'GPL3',
    'packages': ['elapsed'],
    'scripts': ['bin/elapsed'],
    'include_package_data': True,
    'package_data': {'elapsed': ['icons/*.png', 'icons/*.gif', 'yellow.stylesheet', 'copyright.txt']}
}

setup(**config)
