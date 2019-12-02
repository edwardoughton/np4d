Advanced Network Planning in Four Dimensions (anp4d)
=========================================================

Description
-----------

The **Advanced Network Planning in Four Dimensions** (`anp4d`) can quantify the demand
and supply for digital connectivity in Four Dimensions (x, y, z, coordinates + time).

Setup and configuration
-----------------------

All code for **anp4d** is written in
Python (Python>=3.5) and has a number of dependencies.
See `requirements.txt` for a full list.

Using conda
-----------

The recommended installation method is to use [conda](http://conda.pydata.org/miniconda.html),
which handles packages and virtual environments, along with the `conda-forge` channel which
has a host of pre-built libraries and packages.

Create a conda environment called `anp4d`:

    conda create --name anp4d python=3.7

Activate it (run each time you switch projects)::

    activate pysim5g

First, install required packages including `fiona`, `shapely`, `numpy`, and `rtree`:

    conda install fiona shapely numpy rtree

For development purposes, run this command once per machine:

    python setup.py develop

To install pysim5g permanently:

    python setup.py install

To generate results run:

    python scripts/run.py

Contributors
------------

- Ed Oughton (University of Oxford)
- Tom Russell (University of Oxford)
