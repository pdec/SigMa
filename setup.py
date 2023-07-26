import setuptools

with open("README.md") as inf:
    long_description = inf.read()

with open("VERSION") as inf:
    version = inf.read().strip()


def main():
    setuptools.setup(
        name="SigMa",
        version=version,
        scripts=["SigMa.py", "scripts/prepdb.py"],
        author="Przemyslaw Decewicz",
        author_email="p.decewicz@uw.edu.pl",
        description="SigMa - mapping phage signals on bacterial genomes for iterative prophage identification.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/pdec/SigMa",
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3.0",
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Unix',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'biopython>=1.78',
            'pandas>=1.4.3',
            'numpy>=1.22.3',
        ]
    )


if __name__ == "__main__":
    main()
