from setuptools import setup

setup(name='bringreduce',
      version='2017.00',
      description='bRING reduction code',
      url='https://github.com/BringLeiden/bringreduce.git',
      author='Talens',
      author_email='talens@strw.leidenuniv.nl',
      packages=['bringreduce'],
      scripts=['scripts/bringreduce.py'],
      zip_safe=False)