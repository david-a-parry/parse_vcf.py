from nose.tools import *
from parse_vcf.parse_vcf import *
import subprocess

def test_exception_on_no_file():
    assert_raises(TypeError, VcfReader)

def test_header_error():
    assert_raises(HeaderError, VcfReader, 'test_data/invalid_header.vcf')

def test_column_error():
    assert_raises(HeaderError, VcfReader, 'test_data/invalid_col_header.vcf')
    assert_raises(HeaderError, VcfReader, 'test_data/invalid_col_header2.vcf')

def test_open():
    vcf = VcfReader("test_data/test1.vcf")
    assert(vcf)
    vgz = VcfReader("test_data/test1.vcf.gz")
    assert(vgz)
    bcf = VcfReader("test_data/test1.bcf")
    assert(bcf)

def test_read_header():
    v = VcfReader("test_data/test1.vcf.gz")
    assert_equal(len(v.header.meta_header), 67)
    assert_equal(len(v.col_header), 247)
    assert_equal(len(v.metadata['INFO'].keys()), 32)
    assert_equal(len(v.metadata['FORMAT'].keys()), 9)
    assert_equal(len(v.metadata['ALT'].keys()), 1)
    assert_equal(len(v.metadata['FILTER'].keys()), 6)

def test_samples():
    v = VcfReader("test_data/test1.vcf")
    assert_equal(len(v.header.samples), 238)
    for i in range(1, 239):
        sample = "Sample_" + str(i)
        assert_equal(v.header.sample_cols[sample], 8 + i)

def test_read_variant():
    vcf = VcfReader("test_data/test1.vcf")
    record = next(vcf.parser)
    assert_equal(record.CHROM, '1')
    assert_equal(record.POS , 1025535)
    assert_equal(record.ID , 'rs113100937')
    assert_equal(record.REF, 'C')
    assert_equal(record.ALT, 'G')
    assert_equal(record.QUAL, 992.83)
    assert_equal(record.FILTER, 'PASS')
    assert_equal(record.INFO_FIELDS['AC'], '4')
    assert_equal(record.INFO_FIELDS['AF'], '8.439e-03')
    f = record.parsed_info_fields(['AF', 'AC', 'AN'])
    assert_equal(f['AC'][0], 4)
    assert_equal(f['AF'][0], 8.439e-03)
    assert_equal(f['AN'], 474)
    assert_equal(record.SPAN, 1025535)
#1       1025535 rs113100937     C       G       992.83  PASS    AC=4;AF=8.439e-03;AN=474;
    gts = record.parsed_gts(samples=['Sample_1', 'Sample_138'])
    assert_equal([gts[x]['Sample_1'] for x in gts], 
                 [(0, 0), (10, 0), 10, 30, (0, 30, 389)])
    assert_equal([gts[x]['Sample_138'] for x in gts], 
                 [(0, 0), (23, 0), 23, 66, (0, 66, 842)])
    gts2 = record.parsed_gts()
    assert_equal([gts[x]['Sample_1'] for x in gts], 
                 [gts2[x]['Sample_1'] for x in gts2])
    assert_equal([gts[x]['Sample_138'] for x in gts], 
                 [gts2[x]['Sample_138'] for x in gts2])
    assert('Sample_4' not in gts['GT'])
    assert_equal([gts2[x]['Sample_4'] for x in gts2], 
                 [(0, 0), (7, 0), 7, 21, (0, 21, 276)])
                 
