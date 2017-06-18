from setuptools import setup

setup(name='lsreduce',
      version='2017.06',
      description='La Silla reduction code',
      url='https://gitlab.strw.leidenuniv.nl/talens/lsreduce.git',
      author='Talens',
      author_email='talens@strw.leidenuniv.nl',
      packages=['lsreduce'],
      scripts=['scripts/lsreduce.py'],
      zip_safe=False)