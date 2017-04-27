from nose.tools import *
from parse_vcf.parse_vcf import *
import subprocess

def test_exception_on_no_file():
    assert_raises(ParseError, VcfReader)

def test_header_error():
    assert_raises(HeaderError, VcfReader, 'tests/data/invalid_header.vcf')

def test_column_error():
    assert_raises(HeaderError, VcfReader, 'tests/data/invalid_col_header.vcf')
    assert_raises(HeaderError, VcfReader, 'tests/data/invalid_col_header2.vcf')

def test_open():
    v = VcfReader("tests/data/test1.vcf")
    assert(v)

def test_read_header():
    v = VcfReader("tests/data/test1.vcf")
    assert_equal(len(v.meta_header), 67)

def test_samples():
    v = VcfReader("tests/data/test1.vcf")
    assert_equal(len(v.samples), 238)
    for i in range(1, 239):
        sample = "Sample_" + str(i)
        assert_equal(v.sample_cols[sample], 8 + i)


