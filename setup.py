from distutils.core import setup
setup(
    name = "parse_vcf",
    packages = [""],
    version = "0.0.1",
    description = "TO DO",
    author = "David A. Parry",
    author_email = "gantzgraf@github.com",
    url = "https://github.com/gantzgraf/parse_vcf",
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=[
          'pysam',
      ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        ],
)
