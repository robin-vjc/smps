from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='smps',
    version='1.0',
    description='SMPS files parsing utilities.',
    # url='http://github.com/storborg/funniest',
    author='Robin Vujanic',
    author_email='vjc.robin@gmail.com',
    install_requires = required,
    packages=find_packages(),
)
