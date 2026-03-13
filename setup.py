from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='OTSO',
    version='1.2.3',
    author='Nicholas Larsen',
    author_email='nlarsen1505@gmail.com',
    description='Geomagnetic Cutoff Computation Tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/NLarsen15/OTSOpy',
    packages=find_packages(),
    #ext_modules=ext_modules,
    include_package_data=True,
    entry_points={
            'console_scripts': [
                'OTSO.clean=OTSO.otso_cli:clean',
                'OTSO.addstation=OTSO.otso_cli:addstation',
                'OTSO.removestation=OTSO.otso_cli:removestation',
                'OTSO.liststations=OTSO.otso_cli:liststations',
                'OTSO.IGRFupdate=OTSO.otso_cli:IGRFupdate',
                'OTSO.serverdownload=OTSO.otso_cli:serverdownload',
            ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10, <3.15',
    install_requires=[
        'python-dateutil',   # Common dependency across all versions
        'six',              # Common dependency across all versions
        'tzdata',           # Common dependency across all versions
        'psutil',
        'tqdm',
        'requests'
    ],
    extras_require={
        # Python 3.10 requirements
        ':python_version=="3.10"': [
            'meson',
            'numpy>=2.2.6, <2.3.0',
            'packaging',
            'pandas',
            'pytz',
        ],
        # Python 3.11 requirements
        ':python_version=="3.11"': [
            'meson',
            'numpy>=2.4.2, <2.5.0',
            'packaging',
            'pandas',
        ],
        # Python 3.12 requirements
        ':python_version=="3.12"': [
            'meson',
            'numpy>=2.4.2, <2.5.0',
            'packaging',
            'pandas',
            'wheel>=0.46.3',
        ],
        # Python 3.13 requirements (same as 3.12 for now)
        ':python_version=="3.13"': [
            'meson',
            'numpy>=2.4.2, <2.5.0',
            'packaging',
            'pandas',
            'wheel>=0.46.3',
        ],
        # Python 3.14 requirements
        ':python_version=="3.14"': [
            'meson',
            'numpy>=2.4.2, <2.5.0',
            'pandas',
        ],
    },
    )
