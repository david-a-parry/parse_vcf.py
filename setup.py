from setuptools import setup
setup(
    name = "parse_vcf",
    packages = [""],
    version = "0.1a1",
    description = "Variant Call Format parser and convenience methods",
    author = "David A. Parry",
    author_email = "david.parry@igmm.ed.ac.uk",
    url = "https://github.com/gantzgraf/parse_vcf",
    download_url = 'https://github.com/gantzgraf/parse_vcf/archive/0.1a1.tar.gz',
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=['pysam'],
    python_requires='>=3',
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
