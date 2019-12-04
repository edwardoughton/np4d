#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup file for np4d

"""
from glob import glob
from os.path import basename, splitext

from setuptools import find_packages
from setuptools import setup


def readme():
    """Read README contents
    """
    with open('README.md', encoding='utf8') as f:
        return f.read()


setup(
    name='np4d',
    use_scm_version=True,
    license='MIT License',
    description='Network Planning in Four Dimensions',
    long_description=readme(),
    long_description_content_type='text/markdown',
    author='Edward J. Oughton',
    author_email='edward.oughton@ouce.ox.ac.uk',
    url='https://github.com/edwardoughton/np4d',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Communications :: Telephony',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    setup_requires=[
        'setuptools_scm'
    ],
    install_requires=[
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
        'fiona',
        'shapely',
        'numpy',
        'rtree',
        'pytest',
        'imageio',
        'contextily',
        'pygifsicle',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            # eg: 'cdcam = cdcam.cli:main',
        ]
    },
)
