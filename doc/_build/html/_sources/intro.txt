.. _introduction:

opacplot2
#########

Python package for manipulating Equation of State (EoS) and Opacity data.

``Opacplot2`` comes with an EoS Table conversion tool named ``opac-convert``.
It also comes with an EoS Table comparison tool named ``opac-error``.
Both can be found in the :ref:`clts`

Dependencies
************

``opacplot2``'s dependencies include:

* numpy
* six
* pytables
* matplotlib
* scipy
* periodictable
* hedp (https://github.com/luli/hedp)

They can be installed as follows:

``````shell
pip install numpy six pytables matplotlib scipy periodictable
pip install git+https://github.com/luli/hedp
``````

Installation
************

This module requires Python 2.7 or 3.5. The latest version can be installed with::

   pip install git+https://github.com/flash-center/opacplot2

If you have the Propaceos Python reader, in order to include it in the
installation, you must install ``opacplot2`` as follows::

   git clone https://github.com/flash-center/opacplot2
   cp /path/to/opg_propaceos.py opacplot2/opacplot2/
   cd opacplot2
   python setup.py install