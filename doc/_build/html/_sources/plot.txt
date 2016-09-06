Plotting Capabilities
#####################

The plotting capabilities of ``opacplot2`` rely on the ``matplotlib`` python
module. It can be installed with ``pip``. It can also be installed with ``conda``
if you are using Anaconda Python. For information on how to install ``matplotlib``, see
their `website <http://matplotlib.org/users/installing.html>`__.

Configuring ``matplotlib`` to display inline
********************************************

In order to get plots to display inline with the python interpreter,
it is recommended to use ``matplotlib``
in conjunction with `Project Jupyter's QtConsole <http://jupyter.org/qtconsole/stable/>`_.

Steps:

#. Install ``jupyter`` using your preferred Python package manager.
#. Start the console with ``jupyter qtconsole``.
#. Set the mode to suit inline ``matplotlib`` plots using::

       %matplotlib inline

At any point during your Python session on ``qtconsole``, you can view what your
current plots look like by typing the name of your figure.
    

Troubleshooting ``matplotlib`` on OSX
*************************************

In order to display your plots, ``matplotlib`` uses a variety of `backends <http://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`_.
If you have your plots displaying inline on the ``qtconsole``, it is using a backend
specific to the ``qtconsole`` that is separate from other standalone backends.
Thus, for those of you using ``qtconsole``, you can ignore this section.
The recommended backend to use for OSX is ``macosx``. This can be set in ``matplotlib``'s
cofiguration file ``~/.matplotlib/matplotlibrc`` with the line::

    backend: macosx

It can also be set in an interactive session **before** you import ``matplotlib.pyplot``
with::

    matplotlib.rcParams['backend'] = 'macosx'

For more information on the ``matplotlibrc`` file, see ``matplotlib``'s `website <http://matplotlib.org/users/customizing.html#the-matplotlibrc-file>`__.


For Double Implementation Errors
================================

If you are receiving double implementation errors in Python, it is probably due
to the backend chosen for ``matplotlib``. This error is confirmed for OSX users
using the ``Tk`` backends. To fix this, switch your backend to ``macosx`` as shown
above.

For Framework Errors
====================

From ``matplotlib``'s website:

    ''On OSX, two different types of Python Builds exist: a regular build and a
    framework build. In order to interact correctly with OSX through some GUI
    frameworks you need a framework build of Python. At the time of writing the
    macosx, WX and WXAgg backends require a framework build to function correctly.
    Unfortunately virtualenv creates a non framework build even if created from a
    framework build of Python. Conda environments are framework builds. From
    Matplotlib 1.5 onwards the macosx backend checks that a framework build is
    available and fails if a non framework build is found. WX has a similar
    check build in.''

    -- `Matplotlib's documentation <http://matplotlib.org/faq/virtualenv_faq.html?highlight=jupyter#osx>`_

