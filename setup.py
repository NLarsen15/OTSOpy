import sys, os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install


setup(
    name='OTSO',
    version='1.0.0',
    author='Nicholas Larsen',
    author_email='nlarsen1505@gmail.com',
    description='Geomagnetic Cutoff Computation Tool',
    #long_description=open('README.md').read(),
    #long_description_content_type='text/markdown',
    #url='https://github.com/yourusername/your-repo',
    packages=find_packages(),
    #ext_modules=ext_modules,
    include_package_data=True,
    entry_points={
            'console_scripts': [
                'OTSO.setup=OTSO:setup',
            ],
        },

    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12, <3.13',
    install_requires=[
        'psutil==7.0.0',     # Common dependency
    ],
    extras_require={
        ':python_version>="3.10"': [
            'numpy>=2.0.0',
            'pandas>=2.2.2',
            'requests==2.32.3',
        ],
    },
    )
