SCherlock – Detailed Single-Cell Analysis in Python
===================================================

.. image:: _static/scherlock_overview.svg

SCherlock is a set of utilities to make the most out of your single-cell data.
It is supposed to be used in conjuntion with `anndata <https://anndata.readthedocs.io>`__ and
`scanpy <https://scanpy.readthedocs.io>`__.

The main features of SCherlock are all related to detailed qualitative and quantitative analysis
of reads - are there are changes in isoforms, promoters or enhancers? It supports a new file format,
CellPile, which enables quick pile-up generation for each clustering you perform. It can also
take output from Isocounter, a more versatile feature counter.


.. include:: _links.rst
.. include:: _key_contributors.rst

.. role:: small
.. role:: smaller

* Get started by browsing :doc:`tutorials <tutorials>`,

.. _GitHub: https://github.com/henriksson-lab/isocounter


Latest additions
----------------

.. include:: release-latest.rst

.. put references first so all references are resolved
.. NO! there is a particular meaning to this sequence
.. toctree::
   :maxdepth: 1
   :hidden:

   tutorials
   usage-principles
   installation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
