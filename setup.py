import os
from setuptools import setup, find_packages
from src.version import __version__

cur_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(cur_path, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
      name='SVision-pro',
      version=__version__,

      description='SV/CSV callers',
      long_description=long_description,

      url='https://github.com/songbowang125/',

      author='Songbo Wang',
      author_email='songbowang125@163.com',

      license='GPLv3',
      classifiers=[
      'Operating System :: POSIX :: Linux',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 3.7'
      ],

      keywords=['SVision-pro', 'Structural variants(SV)', 'Complex structural variants', 'De no SV/CSV', 'Somatic SV/CSV', 'Single molecular sequencing'],

      packages = ['src', 'src/hash_realign', 'src/network', 'src/network/utils', 'src/pre_process'],
      data_files = [("", ["LICENSE"])],

      zip_safe=False,

      scripts=['SVision-pro'],
      )