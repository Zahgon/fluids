# type: ignore
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2017, 2018, 2019, 2020, 2021, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


__all__ = ["u", "wraps_numpydoc"]

import functools
import inspect
import re
import sys
import types
from collections.abc import Iterable
from copy import copy
from inspect import cleandoc

import numpy as np

import fluids
import fluids.vectorized

ndarray = np.ndarray
try:
    from pint import _DEFAULT_REGISTRY as u
    from pint import DimensionalityError

except ImportError: # pragma: no cover
    raise ImportError("The unit handling in fluids requires the installation "
                      "of the package pint, available on pypi or from "
                      "https://github.com/hgrecco/pint")


"""See fluids.units.rst for documentation for this module.
"""

try:
    doc_stripped = sys.flags.optimize == 2
except:
    doc_stripped = False
# is_critical_flow is broken

def get_docstring(f):
    """Returns the docstring of a function, working in -OO mode also.
    """
    pass


def func_args(func):
    """Basic function which returns a tuple of arguments of a function or
    method.
    """
    pass

u.autoconvert_offset_to_baseunit = True


expr = re.compile("Parameters *\n *-+\n +")
expr2 = re.compile("Returns *\n *-+\n +")
match_sections = re.compile("\n *[A-Za-z ]+ *\n *-+")
match_section_names = re.compile("\n *[A-Za-z]+ *\n *-+")
variable = re.compile("[a-zA-Z_0-9]* : ")
variable_type = re.compile("[a-zA-Z_0-9]* : .*")
match_units = re.compile(r"\[[a-zA-Z0-9().\/*^\- ]*\]")


parse_numpydoc_variables_units_cache = {}



def make_dimensionless_units(unit_str):
    pass

def parse_numpydoc_variables_units(func, replace=None):
    pass

def parse_numpydoc_variables_units_docstring(text):
    pass

def check_args_order(func):
    """Reads a numpydoc function and compares the Parameters and Other
    Parameters with the input arguments of the actual function signature. Raises
    an exception if not correctly defined.

    getargspec is used for Python 2.7 compatibility and is deprecated in Python
    3.

    >>> check_args_order(fluids.core.Reynolds)
    """
    pass

def check_module_docstring_parameters(module, bad_names={"__getattr__", "all_submodules"}):
    """Reads all functions in a module and compares their Parameters and Other
    Parameters sections from their numpydoc docstrings with the actual function
    signatures. Raises an exception if any functions have mismatched definitions.

    Parameters
    ----------
    module : module
        The Python module whose functions should be checked
    bad_names : set[str], optional
        Set of function names to skip during checking [-]

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If any function's signature does not match its docstring parameters

    Examples
    --------
    >>> import fluids
    >>> check_module_docstring_parameters(fluids)
    """
    pass




pint_expression_cache = {}





in_vars_cache = {}
in_units_cache = {}
out_vars_cache = {}
out_units_cache = {}




class UnitAwareClass:
    wrapped = None
    ureg = u
    strict = True
    property_units = {} # for properties and attributes only
    method_units = {}

    def __repr__(self):
        """Called only on the class instance, not any instance - ever.
        https://stackoverflow.com/questions/10376604/overriding-special-methods-on-an-instance
        """
        return self.wrapped.__repr__()

    def __add__(self, other):
        new_obj = self.wrapped.__add__(other.wrapped)
        new_instance = copy(self)
        new_instance.wrapped = new_obj
        return new_instance

    def __sub__(self, other):
        new_obj = self.wrapped.__sub__(other.wrapped)
        new_instance = copy(self)
        new_instance.wrapped = new_obj
        return new_instance

    def __init__(self, *args, **kwargs):
        args_base, kwargs_base =  self.input_units_to_dimensionless("__init__", *args, **kwargs)
        self.wrapped = self.wrapped(*args_base, **kwargs_base)




    def __getattr__(self, name):
        instance = True
        if name in self.class_methods or name in self.static_methods:
            instance = False
        try:
            value = getattr(self.wrapped, name)
        except Exception as e:
            raise AttributeError(f"Failed to get property {name!s} with error {e!s}")
        if value is not None:
            if name in self.property_units:
                if type(value) == dict:
                    d = {}
                    unit = self.property_units[name]
                    for key, val in value.items():
                        d[key] = val*unit
                    return d
                try:
                    return value*self.property_units[name]
                except:
                    # Not everything is going to work. The most common case here
                    # is returning a list, some of the values being None and so
                    # it cannot be wrapped.
                    return value
            else:
                if hasattr(value, "__call__"):

#                    if not instance:
#                        @functools.wraps(value)
#                        # Special case where self needs to be passed in specifically
#                        def call_func_with_inputs_to_SI(*args, **kwargs):
#                            args_base, kwargs_base = self.input_units_to_dimensionless(self, name, *args, **kwargs)
#                            result = value(*args_base, **kwargs_base)
#                            if name == '__init__':
#                                return result
#                            _, _, _, out_vars, out_units = self.method_units[name]
#                            if not out_units:
#                                return
#                            return convert_output(result, out_units, out_vars, self.ureg)
#
#                    else:

                    return call_func_with_inputs_to_SI
                raise AttributeError("Error: Property does not yet have units attached")
        else:
            return value

    _another_getattr = classmethod(__getattr__)






def match_parse_units(doc, i=-1):
    pass

def convert_input(val, unit, ureg, strict=True):
    pass

def parse_expression_cached(unit, ureg):
    pass

def convert_output(result, out_units, out_vars, ureg):
    # Attempt to handle multiple return values
    # Must be able to convert all values to a pint expression
    pass

def wraps_numpydoc(ureg, strict=True):
    pass

def clean_parsed_info(parsed_info):
    pass

def wrap_numpydoc_obj(obj_to_wrap):
    pass

def kwargs_to_args(args, kwargs, signature):
    """Accepts an *args and **kwargs and a signature
    like ['rho', 'mu', 'nu'] which is an ordered list of
    all accepted arguments.

    Returns a list containing all the arguments, sorted, and
    left as None if not specified
    """
    pass
def A_multiple_hole_cylinder(Do, L, holes):
    pass

def V_multiple_hole_cylinder(Do, L, holes):
    pass

def variable_output_wrapper(func, wrapped_basic_func, output_signatures, input_length):
    pass



__pint_wrapped_functions = {}

for name in dir(fluids):
    if "RectangularOffsetStripFinExchanger" in name:
        continue
    if "ParticleSizeDistribution" in name:
        continue
    if name in ("__getattr__", "__test__"):
        continue
    obj = getattr(fluids, name)
    if isinstance(obj, types.FunctionType):
        obj = wraps_numpydoc(u)(obj)
    elif type(obj) == type:
        obj = wrap_numpydoc_obj(obj)
    elif type(obj) is types.ModuleType:
        # Functions accessed with the namespace like friction.friction_factor
        # would call the original function - leads to user confusion if they are exposed
        continue
    elif isinstance(obj, str):
        continue
    if name == "__all__":
        continue
    __all__.append(name)
    __pint_wrapped_functions.update({name: obj})

globals().update(__pint_wrapped_functions)
__all__.extend([
    "UnitAwareClass",
    "check_args_order",
    "convert_input",
    "convert_output",
    "match_parse_units",
    "parse_numpydoc_variables_units",
    "wrap_numpydoc_obj",
    "wraps_numpydoc",
])




variable_output_unit_funcs = {
    # True: arg should be present; False: arg should be None
    "nu_mu_converter": ({(True, False, True): [u.Pa*u.s],
                        (True, True, False): [u.m**2/u.s],
                        }, 3),
    "differential_pressure_meter_solver": ({(True, True, True, True, False, True, True, True): [u.m],
                                            (True, True, True, True, True, False, True, True): [u.Pa],
                                            (True, True, True, True, True, True, False, True): [u.Pa],
                                            (True, True, True, True, True, True, True, False): [u.kg/u.s],
                                            }, 8),
    "isothermal_gas": ({(True, True, False, True, True, True, True): [u.Pa],
                        (True, True, True, False, True, True, True): [u.Pa],
                        (True, True, True, True, False, True, True): [u.m],
                        (True, True, True, True, True, False, True): [u.m],
                        (True, True, True, True, True, True, False): [u.kg/u.s],
                        }, 7)
}

simple_compressible_variable_output = ({(True, True, False, True, True, True, True): [u.m],
                                        (True, True, True, False, True, True, True): [u.m],
                                        (True, True, True, True, False, True, True): [u.Pa],
                                        (True, True, True, True, True, False, True): [u.Pa],
                                        (True, True, True, True, True, True, False): [u.m**3/u.s],
                                        }, 7)
for f in ["Panhandle_A", "Panhandle_B", "Weymouth", "Spitzglass_high", "Spitzglass_low", "Oliphant", "Fritzsche"]:
    variable_output_unit_funcs[f] = simple_compressible_variable_output

IGT_Muller_variable_output = ({(True, True, True, False, True, True, True, True): [u.m],
                               (True, True, True, True, False, True, True, True): [u.m],
                               (True, True, True, True, True, False, True, True): [u.Pa],
                               (True, True, True, True, True, True, False, True): [u.Pa],
                               (True, True, True, True, True, True, True, False): [u.m**3/u.s],
                               }, 8)

for f in ["Muller", "IGT"]:
    variable_output_unit_funcs[f] = IGT_Muller_variable_output


for name, val in variable_output_unit_funcs.items():
    globals()[name] = variable_output_wrapper(getattr(fluids, name),
            __pint_wrapped_functions[name], val[0], val[1])
