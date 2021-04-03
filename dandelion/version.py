from setuptools_scm import get_version
try:
	version = get_version()
except:
	import os
	os.chdir('..')
	version = get_version()

if '+' in version:
	__version__ = version.split('+')[0]
else:
	__version__ = version
