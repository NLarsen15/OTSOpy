import sys, os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install



# class PostInstallCommand(install):
#     """Custom installation to adjust rpath of .so files."""
#     def run(self):
#         install.run(self)  # Run the standard install process
#         so_files = []  # Collect paths to your .so files
#         for root, _, files in os.walk(os.path.join(self.install_lib, 'OTSO')):
#             for file in files:
#                 if file.endswith('darwin.so'):
#                     so_files.append(os.path.join(root, file))
#         for so_file in so_files:
#             # Update rpath using install_name_tool
#             subprocess.run(['install_name_tool', '-add_rpath', '@loader_path', so_file], check=True)

setup(
    name='OTSO',
    version='0.1.1',
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
    python_requires='>=3.7, <3.13',
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
    #cmdclass={
    #    'install': PostInstallCommand     # Custom install command to modify rpath
    #},
    )
