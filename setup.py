from setuptools import setup, find_packages
import pkg_resources

with open("README.md", "r") as readme_file:
    readme = readme_file.read()
exec(open('dandelion/version.py').read())

with open("requirements_dev.txt", "rt", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh.readlines()]

setup(
    name="sc-dandelion",
    version=__version__,
    author="zktuong",
    author_email="kt16@sanger.ac.uk",
    description="sc-BCR analysis tool",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/zktuong/dandelion/",
    packages=find_packages(),
    setup_requires=['setuptools_scm'],
    install_requires=requirements,
    package_data={'dandelion': ['bin/tigger-genotype.R']},
    data_files=[('bin', ['bin/tigger-genotype.R'])],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
    ],
    zip_safe=False
)

