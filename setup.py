import sys
from setuptools import setup

if sys.version_info[:2] < (3, 5):
    print("Python >= 3.5 is required.")
    sys.exit(-1)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pyfmd',
#    version='0.1',
    description='PyFMD provides an object-oriented interface for interacting with the core part of FMD.',
    license='GNU GPL v3',
    url='https://github.com/a-amouei/fmd',
    keywords="fmd molecular-dynamics physics",
    platforms=['Unix'],
    long_description=long_description,
    long_description_content_type='text/markdown',
#    author='',
#    author_email='',
    packages=['pyfmd'],  # 'src'
    install_requires=['numpy', 'ase', 'periodictable']  #
)