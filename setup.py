from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from reportfunk import __version__, _program

setup(name='reportfunk',
      version=__version__,
      packages=find_packages(),
      scripts=["reportfunk/funks/baltic.py",
            "reportfunk/funks/custom_logger.py",
            "reportfunk/funks/log_handler_handle.py",
            "reportfunk/funks/time_functions.py",
            "reportfunk/funks/prep_data_functions.py",
            "reportfunk/funks/report_functions.py",
            "reportfunk/funks/parsing_functions.py",
            "reportfunk/funks/io_functions.py",
            "reportfunk/funks/class_definitions.py",
            "reportfunk/funks/tree_functions.py",
            "reportfunk/funks/table_functions.py"],
      install_requires=[
            "biopython>=1.70",
            "matplotlib>=3.2.1",
            "epiweeks>=2.1.1"
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
