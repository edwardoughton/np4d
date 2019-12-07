Network Planning in Four Dimensions (np4d)
==========================================

Description
-----------

The **Network Planning in Four Dimensions** (`np4d`) can quantify the demand
and supply for digital connectivity in Four Dimensions (x, y, z, coordinates + time).

We demonstrate the capability of this software using vehicle traffic flow estimates,
and focusing on the demand and supply of data for connected vehicles.

Capacity margin estimation: The City of Oxford, UK
--------------------------------------------------

![Example](/movie_capacity_margin.gif)


Setup and configuration
-----------------------

All code for **np4d** is written in Python and has a number of dependencies.
See `requirements.txt` for a full list.

Using conda
-----------

The recommended installation method is to use [conda](http://conda.pydata.org/miniconda.html),
which handles packages and virtual environments, along with the `conda-forge` channel which
has a host of pre-built libraries and packages.

Create a conda python environment called `np4d` with `gdal`:

    conda create --name np4d python gdal

Activate it (run each time you switch projects)::

    activate np4d

First, install required packages including `fiona`, `shapely`, `numpy`, `rtree`:

    conda install fiona shapely numpy rtree pytest

It helps to deactive and then reactivate the env using:

    deactivate

And then:

    conda activate np4d

Then for development purposes, run this command once per machine:

    python setup.py develop

Or, to install np4d permanently:

    python setup.py install

To generate results run:

    python scripts/run.py

Contributors
------------

- Ed Oughton (University of Oxford)
- Tom Russell (University of Oxford)
