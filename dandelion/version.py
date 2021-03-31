exec(open('dandelion/logging/_scmtag.py').read())

if '+' in version:
	__version__ = version.split('+')[0]
else:
	__version__ = version
