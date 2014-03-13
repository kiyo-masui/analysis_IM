"""Unit tests for parse_ini.py"""

import unittest
import os

import parse_ini
import custom_exceptions as ce

test_ini_file_name = 'testfile_parse_ini.ini'
# This is a dictionary of parameters found in the above name file, copied
# by hand.
test_ini_dict = {'a_string' : 'string',
                 'a_list' : [1,2,3],
                 'a_mixed_tuple' : ('A', 1, 2*3.14159),
                 'pi' : 3.14159,
                 'a' : 1,
                 'b' : 2,
                 'd' : 3
                 }


class TestInternal(unittest.TestCase) :
    
    def test__execute_parameter_file(self) :
        file_ini_params = parse_ini._execute_parameter_file(
                              test_ini_file_name)
        for key, value in test_ini_dict.iteritems() :
            self.assertTrue(file_ini_params.has_key(key))
            if file_ini_params.has_key(key) :
                self.assertEqual(file_ini_params[key], value)

class TestReads(unittest.TestCase) :
    
    def setUp(self) :
        self.template_dict = dict(test_ini_dict)
        for key in self.template_dict.iterkeys() :
            self.template_dict[key] = 0
    
    def test_parse_dict_reads(self) :
        self.out_dict = parse_ini.parse_dict(test_ini_dict, self.template_dict,
                                             return_undeclared=False, 
                                             feedback=0)

    def test_parse_readsfile(self) :
        self.out_dict = parse_ini.parse(test_ini_file_name, self.template_dict,
                                        return_undeclared=False, checking=0)

    def test_parse_readsdict(self) :
        self.out_dict = parse_ini.parse(test_ini_dict, self.template_dict,
                                        return_undeclared=False, checking=0)

    def test_defaults_None_ini(self) :
        self.out_dict = parse_ini.parse(None, test_ini_dict, 
                                        return_undeclared=False, checking=0)

    def tearDown(self) :
        for key in self.template_dict.iterkeys() :
            self.assertEqual(self.out_dict[key], test_ini_dict[key])
            self.assertEqual(self.template_dict[key], 0)

class TestDefaults(unittest.TestCase) :
    
    def setUp(self) :
        self.template_dict = dict(test_ini_dict)
        self.extra_params = {'q' : 10,
                             'u' : 11,
                             'another' : ['a', 1, 3*3.14159]
                             }
        self.template_dict.update(self.extra_params)

    def test_parse_dict_defaults(self) :
        self.out_dict = parse_ini.parse(test_ini_dict, self.template_dict,
                                        checking=01)

    def test_parse_file_defaults(self) :
        self.out_dict = parse_ini.parse(test_ini_file_name, self.template_dict,
                                        checking=01)

    def tearDown(self) :
        for key, value in self.template_dict.iteritems() :
            self.assertEqual(self.out_dict[key], value)

class TestReturnsUndeclared(unittest.TestCase) :
    
    def setUp(self) :
        self.template_dict = {'a':1}
        self.template_backup = dict(self.template_dict)
    
    def test_parse_dict_undeclared(self) :
        self.out_dict, self.undeclared = parse_ini.parse(test_ini_dict,
                                             self.template_dict,
                                             return_undeclared=True,
                                             checking=01)
    
    def test_parse_file_undeclared(self) :
        self.out_dict, self.undeclared = parse_ini.parse(test_ini_file_name,
                                             self.template_dict,
                                             return_undeclared=True,
                                             checking=01)

    def tearDown(self) :
        for key, value in test_ini_dict.iteritems() :
            if self.template_dict.has_key(key) :
                self.assertTrue(self.template_backup.has_key(key))
            else :
                self.assertTrue(self.undeclared.has_key(key))
                if self.undeclared.has_key(key) :
                    self.assertEqual(self.undeclared[key], value)

class TestWriteParams(unittest.TestCase) :
    
    def test_circle(self) :
        changed_ini = dict(test_ini_dict)
        changed_ini['a_string'] = 'a different string'
        changed_ini['a'] = 5
        changed_ini['a_list'] = [5,6,7]
        copy_changed = dict(changed_ini)
        parse_ini.write_params(changed_ini, 'temp.ini')
        read_ini = parse_ini.parse('temp.ini', test_ini_dict, checking = 01)
        os.remove('temp.ini')
        for key, value in copy_changed.iteritems() :
            self.assertTrue(read_ini.has_key(key))
            self.assertEqual(value, read_ini[key])
            self.assertEqual(value, changed_ini[key])

class TestExceptions(unittest.TestCase) :
    
    def test_parse_typechecking(self) :
        template_dict = dict(test_ini_dict)
        for swapped in [1, 0, None, 'a', [1,2,3], ('a',3.14), 3.14159] :
            for key, value in test_ini_dict.iteritems() :
                ini_dict = dict(test_ini_dict)
                ini_dict[key] = swapped
                if not type(swapped) is type(value) :
                    self.assertRaises(ce.FileParameterTypeError,
                                      parse_ini.parse, ini_dict,
                                      template_dict, checking=03)

class TestPrefix(unittest.TestCase) :
    

    def setUp(self) :
        self.prefixed_ini_dict = {'tt_a_string' : 'string',
                             'tt_a_list' : [1,2,3],
                             'tt_a_mixed_tuple' : ('A', 1, 2*3.14159),
                             'tt_pi' : 3.14159,
                             'tt_a' : 1,
                             'b' : 2,
                             'tt_d' : 3
                             }
        self.params_init = {'a_string' : 'blah',
                       'pi' : 3.0,
                       'a' : 5,
                       'b' : 10
                        }

    def test_prefix(self) :
        out_dict, undeclared = parse_ini.parse(self.prefixed_ini_dict,
                                             self.params_init, prefix='tt_',
                                             return_undeclared=True,
                                             checking=01)
        self.assertEqual(out_dict['a_string'], 'string')
        self.assertEqual(out_dict['a'], 1)
        self.assertEqual(out_dict['b'], 10)
        self.assertAlmostEqual(out_dict['pi'], 3.14159)
        self.assertEqual(undeclared['b'], 2)
        self.assertEqual(undeclared['tt_d'], 3)
        
    def test_circle(self) :
        changed_ini = dict(self.params_init)
        changed_ini['a'] = 15
        parse_ini.write_params(changed_ini, 'temp.ini', prefix='tt_')
        read_ini = parse_ini.parse('temp.ini', self.params_init, prefix='tt_',
                                   checking = 01)
        self.assertEqual(read_ini['a'], 15)
        self.assertEqual(read_ini['b'], 10)
        os.remove('temp.ini')


if __name__ == '__main__' :
    unittest.main()
