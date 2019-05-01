from setuptools import setup


RELEASE = False
VERSION_MAJOR = 3
VERSION_MINOR = 0
VERSION_PATCH = 0
VERSION = '.'.join(map(str, [VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH]))
if not RELEASE:
    VERSION += '-dev'


setup(name='normalizer',
      version=VERSION,
      packages=['normalizer'],
      description='HGVS variant description normalizer',
      entry_points={
          'console_scripts': ['normalizer=normalizer.cli:main']}
      )
