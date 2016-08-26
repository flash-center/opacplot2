import numpy as np
import sys

def _repr_arr_error(obj_self, idx, var, msg, dens_key, temps_key):
    """
    Print an error message when array failed to pass consistency checks.
    
    Parameters
    ----------
    obj_self : object
        The EoS/opacity table object that is being tested.
    idx : array
        Indices of points that failed the test.
    var : str
        Variable being tested.
    msg : str
        Error message.
    dens_key : str
        Density key specific to the EoS/opacity table object.
    temps_key : str
        Temperature key specific to the EoS/opacity table object.
    """
    out = ['-'*80]
    out.append(' Test failed for {0}/{1} points: {2}'.format(idx.sum(),idx.size, msg))
    out.append('-'*80)
    arr2str = np.array2string

    if idx.size >= obj_self[temps_key][:].size * obj_self[dens_key][:].size:
        out.append(' == density mask ==')
        out.append(arr2str(np.nonzero(idx)[0]))
        out.append(arr2str(obj_self[dens_key][np.nonzero(idx)[0]]))

        out.append(' == temperature mask ==')
        out.append(arr2str(np.nonzero(idx)[1]))
        out.append(arr2str(obj_self[temps_key][np.nonzero(idx)[1]]))

    out.append(' ==     var     ==')
    if var.ndim == idx.ndim:
        out.append(arr2str(var[idx]))
    elif var.ndim == 3 and idx.ndim == 2:
        # Probably something to do with the "ionfrac sums to 1" test.
        out.append(arr2str(np.sum(var, axis=-1)[idx]))


    out.append('-'*80)
    return '\n'.join(out)

def testsuite(test_keys, var, cond=None, mode='short'):
    """
    Test suite decorator to check consistency of an opacity table

    Parameter
    ---------
    test_keys : list
        ['dens_key', 'temp_key'] where 'dens_key' and 'temp_key' are the
        class's specific density & temperature keys.
    cond : None, bool, or lambda expression
        Lambda function takes obj_self and returns a mask in the DT domain.
        In case of None, the existence of the var parameters is taken as 
        condition.
    var :  str, lambda func
        Variable to print if something goes wrong.
    mode : str 
        Type of the test: 'short' or 'full'.
    """
    def func_outerwrap(func):
        def func_wrapper(obj_self, req_mode):
            try:
                if type(var) is str:
                    var_expr = lambda x: x[var]
                else:
                    var_expr = var # If var is lambda expression.
                if cond is None:
                    # Throws KeyError if var not in obj_self.
                    # cond_expr true if obj_self[var] is not all false/zero.
                    cond_expr = np.any(var_expr(obj_self)) 
                elif type(cond) is bool:
                    cond_expr = cond 
                else:
                    # For lambda expression, check whether cond(obj_self)
                    # is not all false/zero.
                    cond_expr = np.any(cond(obj_self))
            except KeyError:
                cond_expr = False  # var is not in obj_self.
            except:
                raise # @testsuite() not called correctly.
            
            # Only run if cond_expr==True.
            # If run_testsuite(mode='full'), run all tests.
            # If run_testsuite(mode='short'), only run tests such that
            # @testsuite(mode='short').
            if cond_expr and \
                    (req_mode == 'full' or req_mode == mode == 'short'):
                mask =  ~func(obj_self, mode=req_mode)
                # If test has failed.
                if len(np.nonzero(mask)[0]):
                    sys.stdout.write('F\n')
                    print()
                    print(_repr_arr_error(obj_self, mask, 
                                          var_expr(obj_self), 
                                          func.__doc__, test_keys[0], 
                                          test_keys[1]))
                # If test passes.
                else:
                    sys.stdout.write('.')
            # If we skip this test.
            else:
                sys.stdout.write('S')
            sys.stdout.flush()
            return None
        return func_wrapper
    return func_outerwrap
    
    
    
    
    
    
    
    