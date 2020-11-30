Installation
------------

Prerequisites
~~~~~~~~~~~~~
You will need AnnData and ScanPy before proceeding with this package. This in turn requires that you install Python.
Follow the instructions for ScanPy_ and then return here to continue


You also need a Java compiler and virtual machine to be installed

.. _ScanPy: https://scanpy.readthedocs.io/en/stable/installation.html


Installing SCherlock
~~~~~~~~~~~~~~~~~~~~

The package consists of two parts; command-line tools for pre-processing, and the Python package to interpret
the generated content. Some routines in the Python package are useful without any of the pre-processors.

First pull the code:

.. code-block:: sh

    $ git clone git@github.com:henriksson-lab/isocounter.git


Then in this directory, run:

.. code-block:: sh

    $ make installpy


This will both generate the command-line tools (in the directory "build") and a python package,
which it will also install.
