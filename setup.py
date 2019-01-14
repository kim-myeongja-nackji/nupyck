from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess

class build(build_py):
    def run(self):
        subprocess.check_call("make")
        build_py.run(self)

setup(
  name = 'nupyck',
  version = "1.0.0",
  packages=find_packages(),
  package_data = {'nupyck':['nupack.so', 'parameters/*']},
  include_package_data = True,
  install_requires = ["numpy>=1.15.4,<2"],
  cmdclass = {"build_py": build}
)
