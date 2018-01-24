from distutils.core import setup
from setuptools import find_packages


setup(name = 'pyreference',
      packages = find_packages(),
      version = '0.1',
      description = 'Library for working with reference genomes',
      author = 'David Lawrence',
      author_email = 'davmlaw@gmail.com',
      url = 'https://bitbucket.org/sacgf/pyreference',
      keywords = ['genomics', 'gtf', 'gff', 'genome', 'genes'],
      classifiers = [],
      install_requires=[
        'biopython',
        'configargparse',
        'deprecated',
        'HTSeq',
        'lazy',
        'pyfasta',
      ],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*',
      scripts=['bin/pyreference_gtf_to_json.py'],)

