#!/usr/bin/env python3

from setuptools import setup, find_packages

# Include README as long description
with open("README.md", 'r') as f:
    long_description = f.read()

# Version management
from src import __version__ as VERSION

setup(
    name='crisprware',
    version=VERSION,
    author='Your Name',
    author_email='ericmalekos@gmail.com',
    url='https://github.com/ericmalekos/crisprware',
    license="Your License",
    keywords="genome editing, CRISPR, bioinformatics",
    description="Tools for CRISPR-based genome editing analysis.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),  # Automatically find and include all packages
    entry_points={
        'console_scripts': [
            'preprocess_annotation = src.preprocess_annotation:main',
            'generate_guides = src.generate_guides:main',
            'index_genome = src.index_genome:main',
            'score_guides = src.score_guides:main',
            'rank_guides = src.rank_guides:main'
        ],
    },
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Your License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[
        'pybedtools>=0.9',
        'rs3>=0.0.15',
        'biopython>=1.8',
        'lightgbm==3.3.5',
        'pybigwig',
        'pandas',
        'numpy'
    ],
    scripts=['scripts/allele_specific_guides.py', 'scripts/bigwig_to_signalwindow.py', 'scripts/crisprscore.R',
    'scripts/gtf_from_ribotish.py','scripts/gtf_from_price.py', 'scripts/replace_snps.py']
)
