import sys, os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

class CustomInstallCommand(install):
    """Custom command to modify rpath after installation"""

    def run(self):
        # Determine the Python version
        python_version = f"{sys.version_info.major}{sys.version_info.minor}"

        package_data = {}

        if sys.platform == "linux":  # Linux or WSL environments
            # Adjust for Python 3.7 for Linux (use 'm' in the filename)
            if python_version == '37':
                package_data = {
                    'OTSO': [
                        f'Parameters/functions/MiddleMan.cpython-{python_version}m-x86_64-linux-gnu.so',  # For Python 3.7 on Linux/WSL
                        'Parameters/functions/StationList.csv',
                        'Parameters/functions/Parameters.par',
                    ],
                }
            else:
                package_data = {
                    'OTSO': [
                        f'Parameters/functions/MiddleMan.cpython-{python_version}-x86_64-linux-gnu.so',  # For other versions on Linux/WSL
                        'Parameters/functions/StationList.csv',
                        'Parameters/functions/Parameters.par',
                    ],
                }

        elif sys.platform == "win32":
            # Adjust for Windows platform (using cp37 for Python 3.7, cp38 for Python 3.8, etc.)
            if python_version == '37':
                package_data = {
                    'OTSO': [
                        f'Parameters/functions/MiddleMan.cp37-win_amd64.pyd',  # For Python 3.7 on Windows
                        'Parameters/functions/StationList.csv',
                        'Parameters/functions/Parameters.par',
                    ],
                }
            else:
                package_data = {
                    'OTSO': [
                        f'Parameters/functions/MiddleMan.cp{python_version}-win_amd64.pyd',  # For other versions on Windows
                        'Parameters/functions/StationList.csv',
                        'Parameters/functions/Parameters.par',
                    ],
                }

        elif sys.platform == "darwin":
            package_data = {
                'OTSO': [
                    f'Parameters/functions/MiddleMan.cpython-{python_version}-darwin.so',
                    'Parameters/functions/StationList.csv',
                    'Parameters/functions/Parameters.par',
                ],
            }

        # Assign the package_data to the dist object
        self.distribution.package_data = package_data

        # Run the standard install process
        super().run()

        # Modify rpath after the package is installed (for macOS)
        if sys.platform == "darwin":
            # Path where the package is installed (can be customized)
            install_dir = self.install_lib
            so_file = f"{install_dir}/OTSO/Parameters/functions/MiddleMan.cpython-312-darwin.so"

            # Use install_name_tool to add the rpath (@loader_path will allow finding dependencies)
            subprocess.check_call(['install_name_tool', '-add_rpath', '@loader_path', so_file])
            print(f"Updated rpath for {so_file} to include @loader_path.")

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
    python_requires='>=3.7, <=3.12.9',
    install_requires=[
        'psutil==7.0.0',     # Common dependency
    ],
    extras_require={
        ':python_version>="3.10"': [
            'numpy>=2.0.0',
            'pandas>=2.2.2',
            'requests==2.32.3',
        ],
        ':python_version=="3.9"': [
            'numpy<=1.25.1',
            'pandas<=1.4.4',
            'requests==2.32.3',
        ],
        ':python_version=="3.8"': [
            'numpy<=1.25.1',
            'pandas<=1.4.4',
            'requests==2.32.3',
        ],
        ':python_version<="3.7"': [
            'numpy<=1.21.6',
            'pandas<=1.3.5',
            'requests<=2.31.0'
        ],
    },
    cmdclass={
        'install': CustomInstallCommand     # Custom install command to modify rpath
    },
    )
