import sys
from setuptools import setup, find_packages

# Determine the Python version
python_version = f"{sys.version_info.major}{sys.version_info.minor}"

# Define package data based on the Python version
package_data = {}

if sys.platform == "linux":
    package_data = {
        'LinuxOTSOModule': [f'OTSO/Parameters/functions/MiddleMan.cpython-{python_version}-x86_64-linux-gnu.so'],
        'StationList': ['OTSO/Parameters/functions/StationList.csv'],
        'TSY04Parameters': ['OTSO/Parameters/functions/Parameters.par'],
    }
elif sys.platform == "win32":
    package_data = {
        'WindowsOTSOModule': [f'OTSO/Parameters/functions/MiddleMan.cp{python_version}-win_amd64.pyd'],
        'StationList': ['OTSO/Parameters/functions/StationList.csv'],
        'TSY04Parameters': ['OTSO/Parameters/functions/Parameters.par'],
    }

setup(
    name='OTSO',
    version='0.1',
    author='Nicholas Larsen',
    author_email='nlarsen1505@gmail.com',
    description='Geomagnetic Cutoff Computation Tool',
    #long_description=open('README.md').read(),
    #long_description_content_type='text/markdown',
    #url='https://github.com/yourusername/your-repo',
    packages=find_packages(),
    #ext_modules=ext_modules,
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10, <=3.12.9',
    install_requires=[
        'numpy== 2.2.3',
        'psutil==7.0.0',
        'pandas==2.2.3',
        'requests==2.32.3',
        # Add any dependencies here
    ],
    package_data=package_data,
)
