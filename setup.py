
# setup.py
from setuptools import setup, find_packages

setup(
    name="SnakeHichipTF",
    version="0.1.4",
    description="A Snakemake pipeline for HiChIP data analysis",
    author="Yidan Sun",
    author_email="syidan@wustl.edu",
    url="https://github.com/YidanSunResearchLab/SnakeHichipTF",  # Updated username
    packages=find_packages(),  # Includes workflows/
    install_requires=[
        "snakemake>=8.28.0",  # Core dependency
    ],
    python_requires=">=3.8",  # Bumped to 3.8 for modern support
    package_data={
        "workflows": [
            "../scripts/*", 
            "../scripts/fithichip/*", 
            "../scripts/hicpro-3.1.0/*", 
            "../scripts/hicpro-3.1.0/bin/*", 
            "../scripts/hicpro-3.1.0/bin/utils/*", 
            "../scripts/hicpro-3.1.0/scripts/*", 
            "../rules/*.Snakefile",       # Snakemake workflows
            "../rules/envs/*.yaml",       # Conda envs
            "../organisms/*"              # Small organism files only
        ],
    },
    entry_points={
        "console_scripts": [
            "Genomesetup = workflows.Genomesetup:main",
            "Hichipsnake = workflows.Snakehichip:main",
            "TFtobias = workflows.TFtobias:main",
            "TFscprinter = workflows.TFscprinter:main",
            "Diffhichip = workflows.Diffhichip:main",
        ]
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",  # Matches python_requires
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
