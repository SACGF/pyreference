from distutils.core import setup

setup(name = 'pyreference',
      packages = ['pyreference'], # this must be the same as the name above
      version = '0.1',
      description = 'Library for working with reference genomes',
      author = 'David Lawrence',
      author_email = 'davmlaw@gmail.com',
      url = 'https://bitbucket.org/sacgf/pyreference',
      keywords = ['genomics', 'gtf', 'gff', 'genome', 'genes'],
      classifiers = [],
      install_requires=[
          'configargparse',
      ])

