#!/usr/bin/env python

# general functions for workflow #############################################
##############################################################################
import os
import yaml


def load_configfile(configFiles, verbose, info='Config'):
    with open(configFiles, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = config

    if verbose:
        print("\n--- " + info + " ---------------------------------------------------------------------")
        print("config file: {}".format(configFiles))
        for k, v in sorted(config.items()):
            print("{}: {}".format(k, v))
        print("-" * 80, "\n")
    return config

def setDefaults():
    """
    Set a number of variables used in the wrappers and the defaults
    """
    # workflowDir is located where this current file is
    workflowDir = os.path.dirname(__file__)

    # defaults
    defaults = load_configfile(os.path.join(workflowDir, "config", "defaults.yaml"), False)

    return workflowDir, defaults
