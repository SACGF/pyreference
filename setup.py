from distutils.core import setup
from setuptools import find_packages
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def _get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(name='pyreference',
      packages=find_packages(),
      version=_get_version("pyreference/__init__.py"),
      description='Library for working with reference genomes and gene GTF/GFFs',
      long_description_content_type="text/markdown",
      long_description=open("README.md").read(),
      author='David Lawrence',
      author_email='davmlaw@gmail.com',
      url='https://github.com/SACGF/pyreference',
      keywords=['genomics', 'gtf', 'gff', 'genome', 'genes'],
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
      ],
      install_requires=[
          'numpy',
          'biopython',
          'bioutils',
          'configargparse',
          'deprecation',
          'HTSeq',
          'lazy',
          'pysam',
          'pandas',
          'seaborn',
      ],
      python_requires='>=2.7, >=3.5',
      scripts=['bin/pyreference_biotype.py'])
