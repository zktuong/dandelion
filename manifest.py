#!/usr/bin/env python
'''
Generates the version number file as manifest
'''

import json
exec(open('dandelion/version.py').read())

def main():
    data = {}
    data['version'] = []
    data['version'].append({'version': __version__})
    with open('manifest.json', 'w') as outfile:
        json.dump(data, outfile)


if __name__ == "__main__":    
    main()
