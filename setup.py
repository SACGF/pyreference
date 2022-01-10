from distutils.core import setup
from setuptools import find_packages

setup(name='pyreference',
      packages=find_packages(),
      version='0.6.1',
      description='Library for working with reference genomes and gene GTF/GFFs',
      author='David Lawrence',
      author_email='davmlaw@gmail.com',
      url='https://github.com/SACGF/pyreference',
      keywords=['genomics', 'gtf', 'gff', 'genome', 'genes'],
      classifiers=[],
      install_requires=[
          'biopython',
          'configargparse',
          'deprecation',
          'HTSeq',
          'lazy',
          'pysam',
      ],
      python_requires='>=2.7, >=3.5',
      scripts=['bin/pyreference_gff_to_json.py',
               'bin/pyreference_biotype.py'])
