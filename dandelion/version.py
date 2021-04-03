from .logging import version

if '+' in version:
	__version__ = version.split('+')[0]
else:
	__version__ = version
