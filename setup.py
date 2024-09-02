#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

setup(
    entry_points={
        'console_scripts': [
            'acore=acore.cli:main',
        ],
    },
    include_package_data=True,
    keywords='acore',
    packages=find_packages(include=['acore', 'acore.*']),
    url='https://github.com/Multiomics-Analytics-Group/acore',
    zip_safe=False,
)
