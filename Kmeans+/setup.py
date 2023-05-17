from setuptools import setup, Extension
setup(name='mykmeanssp',
      version='1.0',
      description='kmeanspp for hw2',
      ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])])
