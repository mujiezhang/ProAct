from setuptools import setup, find_packages

setup(
    name='proact',
    version='1.0.0',
    author='Mujie Zhang',
    author_email='zhangmujie@sjtu.edu.cn',
    description='ProAct: Provirus Activity Detector',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/mujiezhang/ProAct',
    license='GPL-3.0-or-later',
    packages=find_packages(include=['proact*']),
    package_data={
        'proact.extract_gtdb_mg': ['hmm/**/*'],
    },
    include_package_data=True,
    python_requires='>=3.6',
    install_requires=[
        'pandas',
        'biopython',
        'pysam',
    ],
    entry_points={
        'console_scripts': [
            'proact=proact.pipeline:main',
            'proact-generate-depth=proact.generate_depth:main',
            'proact-extract-gtdb-mg=proact.extract_gtdb_mg.extract_gtdb_mg:main',
            'proact-get-depth=proact.get_MG_and_phage_depth:main',
            'proact-calculate-ptoh=proact.calculate_PtoH:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
