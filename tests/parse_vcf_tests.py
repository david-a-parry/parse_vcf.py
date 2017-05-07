from nose.tools import *
from parse_vcf import *
import subprocess

def test_exception_on_no_file():
    assert_raises(ParseError, VcfReader)

def test_header_error():
    assert_raises(HeaderError, VcfReader, 'tests/data/invalid_header.vcf')

def test_column_error():
    assert_raises(HeaderError, VcfReader, 'tests/data/invalid_col_header.vcf')
    assert_raises(HeaderError, VcfReader, 'tests/data/invalid_col_header2.vcf')

def test_open():
    vcf = VcfReader("tests/data/test1.vcf")
    assert(vcf)
    vgz = VcfReader("tests/data/test1.vcf.gz")
    assert(vgz)
    bcf = VcfReader("tests/data/test1.bcf")
    assert(bcf)

def test_read_header():
    v = VcfReader("tests/data/test1.vcf.gz")
    assert_equal(len(v.header.meta_header), 67)
    assert_equal(len(v.col_header), 247)
    assert_equal(len(v.metadata['INFO'].keys()), 32)
    assert_equal(len(v.metadata['FORMAT'].keys()), 9)
    assert_equal(len(v.metadata['ALT'].keys()), 1)
    assert_equal(len(v.metadata['FILTER'].keys()), 6)

def test_samples():
    v = VcfReader("tests/data/test1.vcf")
    assert_equal(len(v.header.samples), 238)
    for i in range(1, 239):
        sample = "Sample_" + str(i)
        assert_equal(v.header.sample_cols[sample], 8 + i)


