#!/usr/bin/env python
"""ReFrame base module for Asimov-CCS tests"""

import os

import reframe as rfm
import reframe.utility.sanity as sn


#  class BuildCCS(rfm.CompileOnlyRegressionTest, SetupDependencies):
class BuildCCS(rfm.CompileOnlyRegressionTest):
  """
  Compiles CCS, making it available as a reframe fixture:
  https://reframe-hpc.readthedocs.io/en/stable/tutorial.html#test-fixtures
  """

  build_system = "Make"
  sourcedir = "asimov-ccs"

  local = True
  build_locally = False

  # loads ccs/CMP-2026-06 module file
  valid_prog_environs = ["PrgEnv-gnu", "PrgEnv-cray"]

  # Reframe only supports https git cloning using the sourcedir route which doesn't work with
  # private repos, prebuild_cmds allows to use the user's credentials
  prebuild_cmds = ["git clone git@git.ecdf.ed.ac.uk:asimov/asimov-ccs.git", "cd asimov-ccs"]

  @run_before("compile")
  def build_options(self):
    """Setup environment for build"""
    # General options
    self.build_system.max_concurrency = 8
    self.build_system.options = ["ccs_app"]
    self.build_system.srcdir = "asimov-ccs"

  @sanity_function
  def sanity_executable_exists(self):
    """Confirms that the executable was built"""
    build_dir = f"{self.stagedir}/asimov-ccs"
    return sn.path_exists(os.path.join(build_dir, "ccs_app"))


@rfm.simple_test
class RunCCS(rfm.RunOnlyRegressionTest):
  """Run CCS"""

  valid_systems = ["archer2:compute", "cirrus-ex:compute"]
  valid_prog_environs = ["PrgEnv-cray", "PrgEnv-gnu"]
  executable = "ccs_app"
  stream_binary = fixture(BuildCCS, scope="environment")
  executable_opts = [
          "--ccs_m 32",
          "--ccs_case TaylorGreenVortex",
          ]


  sourcesdir = "inputs"
  keep_files = ["residulas.log", "timers.csv"]

  strict_check = True
  tags = {"full_application"}


  @run_after("setup")
  def set_executable(self):
    """sets up executable"""
    # stage_dir is the base folder for the build job;
    # sourcedir is the folder where the executable was compiled
    # self.executable is defined in this class, above
    self.executable = os.path.join(
            self.stream_binary.stagedir,
            self.stream_binary.sourcedir,
            self.executable,
            )


  @sanity_function
  def assert_finished(self):
    """Sanity check that simulation finished successfully"""
    file_len = len(self.keep_files[0])
    return sn.assert_ge(file_len, 10)
