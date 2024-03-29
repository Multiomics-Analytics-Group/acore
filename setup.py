#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Alberto Santos Delgado",
    author_email='albsad@dtu.dk',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    description="A Python package with statistical functions to analyse multimodal molecular data",
    entry_points={
        'console_scripts': [
            'acore=acore.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='acore',
    name='acore',
    packages=find_packages(include=['acore', 'acore.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/Multiomics-Analytics-Group/acore',
    version='0.1.0',
    zip_safe=False,
)
