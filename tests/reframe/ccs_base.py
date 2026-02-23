#!/usr/bin/env python
"""ReFrame base module for Asimov-CCS tests"""

import os

import reframe as rfm
import reframe.utility.sanity as sn


class BuildCCS(rfm.CompileOnlyRegressionTest):
  """
  Compiles CCS, making it available as a reframe fixture:
  https://reframe-hpc.readthedocs.io/en/stable/tutorial.html#test-fixtures
  """

  build_system = "Make"
  modules = []
  sourcesdir = "https://git.ecdf.ed.ac.uk/asimov/asimov-ccs.git"
  #  sourcepath = "src"
  local = True
  build_locally = False

  valid_prog_environs = ["PrgEnv-gnu", "PrgEnv-cray"]

  @run_before("compile")
  def prepare_build(self):
    """Setup environment for build"""
    self.build_system.max_concurrency = 8
    self.build_system.build_dir = f"{self.stagedir}/ccs_build"
    self.build_system.options = ["all"]

    ccs_deps = "/work/e609/e609/shared/ccs-deps/cray"
    self.env_vars["PETSC_DIR"] = ccs_deps + "/petsc"
    self.env_vars["FYAMLC_DIR"] = ccs_deps + "/fyaml-c"
    self.env_vars["adios2_DIR"] = ccs_deps + "/adios2"
    self.env_vars["PARHIP_DIR"] = ccs_deps + "/parhip"
    self.env_vars["PARMETIS_DIR"] = ccs_deps + "/parmetis"
    self.env_vars["RCMF90_DIR"] = ccs_deps + "/rcmf90"
    self.env_vars["PATH"] = ccs_deps + "/.../makedepf90/bin" + ":" + os.getenv("PATH")
    if self.current_environ == "PrgEnv-Cray":
      self.env_vars["CMP"] = "cray"
    elif self.current_environ == "PrgEnv-gnu":
      self.env_vars["CMP"] = "gnu"

  @sanity_function
  def sanity_executable_exists(self):
    """Confirms that the executable was built"""
    build_dir = f"{self.stagedir}/ccs_build"
    return sn.path_exists(os.path.join(build_dir, "ccs_app"))


@rfm.simple_test
class RunCCS(rfm.RunOnlyRegressionTest):
  """Run CCS"""

  valid_systems = ["archer2:compute"]
  valid_prog_environs = ["PrgEnv-cray"]
  executable = "ccs_app"
  stream_binary = fixture(BuildCCS, scope="environment")

  keep_files = ["residulas.log", "timers.csv"]

  strict_check = True
  tags = {"full_application"}

  @run_before("compile")
  def prepare_build(self):
    """Setup environment for build"""
    self.build_system.build_dir = f"{self.stagedir}/ccs_build"
    self.build_system.options = ["all"]

    ccs_deps = "/work/e609/e609/shared/ccs-deps/cray"
    self.env_vars["PETSC_DIR"] = ccs_deps + "/petsc"
    self.env_vars["FYAMLC_DIR"] = ccs_deps + "/fyaml-c"
    self.env_vars["adios2_DIR"] = ccs_deps + "/adios2"
    self.env_vars["PARHIP_DIR"] = ccs_deps + "/parhip"
    self.env_vars["PARMETIS_DIR"] = ccs_deps + "/parmetis"
    self.env_vars["RCMF90_DIR"] = ccs_deps + "/rcmf90"
    self.env_vars["PATH"] = ccs_deps + "/.../makedepf90/bin" + ":" + os.getenv("PATH")
    if self.current_environ == "PrgEnv-Cray":
      self.env_vars["CMP"] = "cray"
    elif self.current_environ == "PrgEnv-gnu":
      self.env_vars["CMP"] = "gnu"

    @run_after("setup")
    def set_executable(self):
        """sets up executable"""
        self.executable = os.path.join(self.stream_binary.build_system.builddir, "ccs_app")


  @sanity_function
  def assert_finished(self):
    """Sanity check that simulation finished successfully"""
    file_len = len(self.keep_files[0])
    return sn.assert_ge(file_len, 10)
