import unittest
import sys
import os

baseDir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(baseDir)
from common_functions import (
    sanity_dict_clean,
    load_configfile
)

class TestSanityDictClean(unittest.TestCase):
    """ test set for cf.sanity_dict_clean helper function """

    def setUp(self):
        self.dictionary = {
            "maindir": "/path/to/maindir",
            "workflow": "workflow_name",
            "other_key": 1
        }

    def test_output_function(self):
        """ tests if function cleans dictionary """
        self.assertEqual({"other_key": 1}, sanity_dict_clean(self.dictionary))

    def test_list_argument(self):
        """ tests if TypeError is raised if not dict is given as argument """
        self.assertRaises(TypeError, sanity_dict_clean, list(self.dictionary))
        self.assertRaises(TypeError, sanity_dict_clean, list(self.dictionary)[0])

class TestLoadConfigfile(unittest.TestCase):
    """ test set for cf.load_configfile helper function """

    def setUp(self):
        self.key_names = [
            'snakemake_executable',
            'tmpDir',
            'pipeline',
            'outdir',
            'configFile',
            'maxJobs',
            'indir',
            'genome',
            'ext',
            'reads',
            'aligner',
            'trim',
            'trimmer',
            'trimmerOptions',
            'dedup',
            'mapq',
            'properPairs',
            'mateOrientation',
            'insertSizeMax',
            'alignerOpts',
            'verbose'
        ]
        self.config = load_configfile("config/defaults.yaml", False)

    def test_type_return(self):
        """ tests whether function returns a dictionary """
        self.assertIsInstance(self.config, dict)

    def test_total_key_number(self):
        """ tests whether keys of dictionary match expected number """
        self.assertEqual(len(self.config), len(self.key_names))

    def test_presence_single_keys(self):
        """ tests whether all keys are present one by one """
        for k in self.key_names:
            self.assertTrue(k in self.key_names)

    def test_config_not_present(self):
        """ tests whether a FileNotFound exception is raised if wrong path is provided """ 
        self.assertRaises(FileNotFoundError, load_configfile, "config/defaults.yam", False)


if __name__ == '__main__':
    unittest.main()

