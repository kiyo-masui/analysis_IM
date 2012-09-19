"""This parser is my system for reading input files for large programs.  The idea
is that the only argument for the program should always be the input file and
all the parameters are read from that file.  The input file will have plain
python syntax.  I've found this to have the best flexibility, avoiding the
need to have many versions of the same code.  However, because any python
statements are executed when the input file is read, this system is not
appropriate if security is an issue (it's an arbitrary code exceution security
hole).

This is purposely written as a set of functions rather than a class.  I can
think of no reason that you would want the parser to stick around after being
called and the output dictionaries are pretty self contained.

Revision History:
  KM August '10 - Wrote initial code (fileparser and dictparser).
                - Later converted fileparser to just parse, which is an
                  interface for both files and dicts.
  KM Oct. '10   - Added write_params
  KM Mar. '11   - Changed checking argument to feedback and type_check.
"""

import custom_exceptions as ce

def parse(ini_data, params, return_undeclared=False, prefix='',
          feedback=2, type_check=False, checking=-1):
    """Parses a python file or dictionary to get parameters.
    
    This function accepts a filename and a dictionary of keys and pre typed
    values. It returns a dictionary of the same keys with values read from
    file.  It optionally performs type checking.

    Parameters
    ----------
        ini_data: a string, containing a python file name, or a dictionary.  The
            file must contain a script (not function) that defines parameter
            values in the local namespace.  Alternately, if ini is a
            dictionary, then parameters are read from the dictionary.
            Variables must have names and
            types corresponding to the params dictionary argument.
        params: a dictionary of keys and corresponding to variable names to be
            read from file and values corresponding to defaults if the
            corresponding variable is not found in the file.
        return_undeclared: Bool default False.  Whether to return a second
            dictionary of with variables found in the parameter file but not in
            the in params argument.
        prefix: String default ''.  A prefix added to parameter names (defined
            in the keys of params) when read from the input file or dictionary.
            The prefix is not added to the returned output dictionary.
        feedback: integer 1 to 10, default 2.  Desired feedback level,
            controling what to pring to the standard out.
        type_check: Boolian default False. Whethar to raise an exception if the
            recived value for a parameter is a different type than the default
            value.
        checking (deprecated, use feedback and typecheck):
            Perform various checks:
            1's digit: perform type checking on the values in the file and in
                passed params:
                    0 not at all
                    2 print warning (default)
                    3 (or greater) raise an Exception
            10s digit: parameter feedback:
                    0 none
                    2 print message when parameters remain default value
                        (default)
                    3 print all parameters and whether they've defaulted

    Returns
    -------
        out_params: A dictionary with the same keys as argument params but 
            with values read from file.
        undeclared: Optional. A dictionary that holds any key found in the 
            file but not in params. Returned if return_undeclared=True.
    """
    
    # Deal with deprecated checking variable.
    if checking != -1 and feedback == 2 and type_check == False:
        old_typecheck = checking%10
        parcheck = (checking - old_typecheck)//10
        if old_typecheck >= 3 :
            type_check = True
        if old_typecheck >= 2 and parcheck > 2 :
            feedback = 2
        else :
            feedback = parcheck
        
    if isinstance(ini_data, str) :
        if feedback > 0 :
            print 'Reading parameters from file: '+ ini_data
        # Convert local variables defined in python script to dictionary.
        # This is in a separate function to avoid namespace issues.
        dict_to_parse = _execute_parameter_file(ini_data)
    elif isinstance(ini_data, dict) :
        if feedback > 0 :
            print 'Reading parameters from dictionary.'
        dict_to_parse = ini_data
    elif ini_data is None :
        if feedback > 0 :
            print 'No input, all parameters defaulted.'
        if return_undeclared :
            return dict(params), {}
        else :
            return dict(params)
    else :
        raise TypeError("Argument ini must be a dictionary, file name, "
                        "or None (to accept defaults).")
    
    return parse_dict(dict_to_parse, params, return_undeclared, prefix,
                      feedback, type_check)

def parse_dict(dict_to_parse, params, return_undeclared=False, prefix='',
               feedback=2, type_check=False):
    """Same as parse_ini.parse except parameters read from only dictionary.
    
    This function is intended for internal use.  All of it's functionality is
    availble from the parse function.

    This function accepts an input dictionary and a dictionary of keys 
    and pre typed
    values. It returns a dictionary of the same keys with values read from
    the input dictionary.  See the docstring for parse for more
    information, the only difference is the first argument must be a
    dictionary.

    Arguments:
        dict_to_parse: A dictionary containing keys and values to be read as
            parameters.  Entries should have keys and
            types corresponding to the pars dictionary argument (depending on
            level of checking requested).
      """
    
    # Same keys as params but for checking but contains only a flag to indicate
    # if parameter retained it's default value.
    defaulted_params = {}
    for key in params.iterkeys():
        defaulted_params[key] = True
    # Make dictionaries for outputs
    undeclared = {} # For keys found in dict_to_parse and not in params
    out_params = dict(params)

    # Loop over both input dictionaries and look for matching keys
    for inkey, invalue in dict_to_parse.iteritems():
        found_match_flag = False
        for key, value in params.iteritems():
            # Check for matching keys. Note stripping.
            if prefix + key.strip() == inkey.strip():
                if type(value) != type(invalue):
                    if type_check:
                        raise ce.FileParameterTypeError(
                            "Tried to assign an input "
                            "parameter to the value of the wrong type " 
                            "and asked for strict type checking. "
                            "Parameter name: " + key)
                    elif feedback > 1:
                        print ("Warning: Assigned an input "
                            "parameter to the value of the wrong type. "
                            "Parameter name: " + key)
                out_params[key] = invalue
                found_match_flag = True
                defaulted_params[key]=False
                # There shouldn't be another matching key so:
                break
        if not found_match_flag :
            # Value found in dict_to_parse was not found in params
            undeclared[inkey]=invalue
    # Check if parameters have remained a default value and print information
    # about the parameters that were set. Depending on feedback level.
    if feedback > 1 :
        print "Parameters set."
        for key, value in out_params.iteritems():
            if defaulted_params[key] :
                print "parameter: "+key+" defaulted to value: "+str(value)
            elif feedback > 2 :
                print "parameter: "+key+" obtained value: "+str(value)

    if return_undeclared :
        return out_params, undeclared
    else :
        return out_params

def _execute_parameter_file(this_parameter_file_name):
    """
    Executes python script in named file and returns dictionary of variables
    declared in that file.
    """
    
    # Only a few locally defined variables and all have a long name to avoid
    # namespace conflicts.

    # Execute the filename which presumably holds a python script. This will
    # bring the parameters defined there into the local scope.
    try:
        exec(open(this_parameter_file_name).read())
    except Exception as E:
        print locals()
        msg = ("Execution of parameter file " + this_parameter_file_name +
               " caused an error.  The error message was: " + repr(E))
        raise ce.ParameterFileError(msg)
    # Store the local scope as a dictionary.
    out = locals()
    # Delete all entries of out that correspond to variables defined in this
    # function (i.e. not in the read file).
    del out['this_parameter_file_name']
    # Return the dictionary of parameters read from file.
    return out

def write_params(params, file_name, prefix='', mode='w') :
    """Write a parameter dictionary to file.

    Given a dictionary of parameters, such as one of the ones read br the parse
    function, this program writes a compatible ini file.
    
    This should work if the parameters are built in types, but no promises for
    other types. Basically if the out put of 'print param' looks like it could
    go on the rhs of the assignment operator, you are in good shape.

    arguments:
        params : dictionary of parameter names and values to be written to
            file.
        file_name: sting. File_name to write to.
        prefix : prefix for teh parameter names when written to file.
        mode: 'a' or 'w'.  Whether to open the file in write or append mode.
    """
    
    if not (mode == 'w' or mode == 'a') :
        raise ValueError("Params can be written with mode either 'w' or 'a'.")
    file = open(file_name, mode)
    for par_name, value in params.iteritems() :
        line_str = prefix + par_name + ' = '
        try :
            line_str = line_str + repr(value)
        except SyntaxError :
            try :
                line_str = line_str + repr(value)
            except SyntaxError :
                line_str = line_str + "'not representable'"
        line_str = line_str + '\n'
        file.write(line_str)
    file.close()
        
