from setuptools_scm import get_version
version = get_version()

if '+' in version:
	__version__ = version.split('+')[0]
else:
	__version__ = version
