from distutils.core import setup
setup(
    name = "parse_vcf",
    packages = [""],
    version = "0.1",
    description = "Variant Call Format parser and convenience methods",
    author = "David A. Parry",
    author_email = "gantzgraf@github.com",
    url = "https://github.com/gantzgraf/parse_vcf",
    download_url = 'https://github.com/gantzgraf/parse_vcf/archive/0.1.tar.gz',
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=['pysam'],
    python_requires='>=3',
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        ],
)
