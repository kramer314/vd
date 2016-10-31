# Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

# Name of default JSON configuration file; set to None to force manual
# configuration
default_config_fname = "./config_json.default"

# Flag to specify custom build configuration
AddOption("--config-json", dest="custom_json_fname", default=default_config_fname)

# Standard library imports
import os
import json


# Helper functions
def abspath(path):
    return os.path.abspath(path) + "/"


# Parse build configuration from JSON
custom_json_fname = GetOption("custom_json_fname")
if custom_json_fname is None:
  raise ValueError("JSON configuration file specification required.")

custom_json_f = open(custom_json_fname, "r")
custom_json = json.load(custom_json_f)
custom_json_f.close()

# Merge custom and default configurations
default_json = dict()
if default_config_fname is not None and default_config_fname != custom_json_fname:
  default_json_f = open(default_json_fname, "r")
  default_json = json.load(default_json_f)
  default_json_f.close()

config_json = default_json.copy()
config_json.update(custom_json)


# Extract source and build directories
source_dir = abspath(config_json["src_path"])
build_dir = abspath(config_json["build_path"])


# Extract library names and paths
lib_names = []
lib_paths = []
lib_path_root = config_json["lib_path_root"]
libs_dict = config_json["libs"]

for lib in libs_dict:
  lib_names.append(lib["name"])

  lib_path = lib["path"]

  if (lib["use_root"]):
    lib_path = lib_path_root + lib_path

  lib_paths.append(abspath(lib_path))


# Extract compiler flags
flags = ""

debug = bool(config_json["debug"])

flags_dict = config_json["flags"]

for flag_cat, flag_params in flags_dict.items():
    flag_list = flag_params["list"]
    if (debug and bool(flag_params["debug"])) or (not debug and bool(flag_params["prod"])):
        flags += " ".join(flag for flag in flag_list) + " "


# Begin build
env = DefaultEnvironment(ENV = os.environ, TOOLS = ['default', "gfortran"])

env.Replace(F90FLAGS = flags)
env.Replace(LINKFLAGS = flags)
env.Replace(FORTRANMODDIRPREFIX = "-J ")
env.Replace(FORTRANMODDIR = build_dir)
env.Replace(F90PATH = [lib_paths, build_dir])

Export("env")
Export("lib_names")
Export("lib_paths")

SConscript(source_dir+"SConscript", variant_dir=build_dir, duplicate=1)


# For whatever reason, we can't use duplicate=0 and have *.mod files in the
# build directory. But, if we duplicate the source tree into the build
# directory SCons doesn't automatically clean the source files, so we have to
# manually define the entire build directory as a cleaning target.
Clean(".", build_dir)
