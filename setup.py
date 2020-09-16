from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from reportfunk import __version__, _program

setup(name='reportfunk',
      version=__version__,
      packages=find_packages(),
      scripts=["reportfunk/funks/funkIO.py",
            "reportfunk/funks/funkR.py",
            "reportfunk/funks/custom_logger.py",
            "reportfunk/funks/funkT.py"],
      install_requires=[
            "biopython>=1.70",
            "matplotlib>=3.2.1"
        ],
      description='snipit',
      url='https://github.com/cov-ert/reportfunk',
      author='Aine OToole, Verity Hill',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = reportfunk.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
