#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['click>=7.0',
                'numpy == 1.23.2',
                'pandas == 2.0.2',
                'scipy == 1.10.1',
                'networkx == 3.1',
                'biopython == 1.81',
                'combat == 0.3.3',
                'gseapy == 1.0.4',
                'kmapper == 2.0.1',
                'lifelines == 0.27.7',
                'pingouin == 0.5.3',
                'python-louvain == 0.16',
                'PyWGCNA == 1.16.8',
                'snfpy == 0.2.2',
                'umap-learn == 0.5.3',
                # not in requirements_dev:
                'statsmodels',
                'combat',
                ]

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
