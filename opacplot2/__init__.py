"""Package for manipulating Equation of State (EoS) and Opacity data.
"""

__version__ = '1.0.0'

# eos_ prefix
from .eos_plotter   import *

# opg_ prefix
from .opg_hdf5      import *
from .opg_ionmix    import *
from .opg_multi     import *
# The Propaceos module is not distributed. More details can be found in
# opg_propaceos-note.py.
#from .opg_propaceos import *
from .opg_qeos      import *
from .opg_sesame    import *
from .opg_tabop     import *

# plotting
from .opac_plotter  import *
from .histogram     import *

# Constants.
from .constants     import *

# Utilities.
from . import utils
