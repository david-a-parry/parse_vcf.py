import os
import gzip
import re
import warnings
from stat import S_ISREG

#common INFO fields and their types in case absent from header
COMMON_INFO = {
    '1000G' : {'Type' : 'Flag', 'Number' : 1},
    'AA' : {'Type' : 'String', 'Number' : 1},
    'AC' : {'Type' : 'Integer', 'Number' : 'A'},
    'AF' : {'Type' : 'Float', 'Number' : 'A'},
    'AN' : {'Type' : 'Integer', 'Number' : 1},
    'BQ' : {'Type' : 'Float', 'Number' : 1},
    'CIGAR' : {'Type' : 'String', 'Number' : 1},
    'DB' : {'Type' : 'Flag', 'Number' : 1},
    'DP' : {'Type' : 'Integer', 'Number' : 1},
    'END' : {'Type' : 'Integer', 'Number' : 1},
    'H2' : {'Type' : 'Flag', 'Number' : 1},
    'H3' : {'Type' : 'Flag', 'Number' : 1},
    'MQ' : {'Type' : 'Float', 'Number' : 1},
    'MQ0' : {'Type' : 'Integer', 'Number' : 1},
    'NS' : {'Type' : 'Integer', 'Number' : 1},
    'SB' : {'Type' : 'String', 'Number' : 1},
    'SOMATIC' : {'Type' : 'Flag', 'Number' : 1},
    'VALIDATED' : {'Type' : 'Flag', 'Number' : 1},
    # for SVs
    'BKPTID' : {'Type' : 'String', 'Number' : '.'},
    'CICN' : {'Type' : 'Integer', 'Number' : 2},
    'CICNADJ' : {'Type' : 'Integer', 'Number' : 1},
    'CIEND' : {'Type' : 'Integer', 'Number' : 2},
    'CILEN' : {'Type' : 'Integer', 'Number' : 2},
    'CIPOS' : {'Type' : 'Integer', 'Number' : 2},
    'CN' : {'Type' : 'Integer', 'Number' : 1},
    'CNADJ' : {'Type' : 'Integer', 'Number' : '.'},
    'DBRIPID' : {'Type' : 'String', 'Number' : 1},
    'DBVARID' : {'Type' : 'String', 'Number' : 1},
    'DGVID' : {'Type' : 'String', 'Number' : 1},
    'DPADJ' : {'Type' : 'Integer', 'Number' : 1},
    'EVENT' : {'Type' : 'String', 'Number' : 1},
    'HOMLEN' : {'Type' : 'Integer', 'Number' : '.'},
    'HOMSEQ' : {'Type' : 'String', 'Number' : '.'},
    'IMPRECISE' : {'Type' : 'Flag', 'Number' : 1},
    'MATEID' : {'Type' : 'String', 'Number' : 1},
    'MEINFO' : {'Type' : 'String', 'Number' : 4},
    'METRANS' : {'Type' : 'String', 'Number' : 4},
    'NOVEL' : {'Type' : 'Flag', 'Number' : 1},
    'PARID' : {'Type' : 'String', 'Number' : 1},
    'SVLEN' : {'Type' : 'Integer', 'Number' : 1},
    'SVTYPE' : {'Type' : 'String', 'Number' : 1},
}

#common genotype FORMAT Fields and their types in case absent from header
COMMON_FORMAT = {
    'DP'  : {'Type' : 'Integer', 'Number' : 1},
    'EC'  : {'Type' : 'Integer', 'Number' : 'A'},
    'FT'  : {'Type' : 'String', 'Number' : 1},
    'GL'  : {'Type' : 'Float', 'Number' : 'G'},
    'GLE' : {'Type' : 'String', 'Number' : '.'},
    'GP'  : {'Type' : 'Float', 'Number' : 'G'},
    'GQ'  : {'Type' : 'Integer', 'Number' : 1},
    'GT'  : {'Type' : 'String', 'Number' : 1},
    'HQ'  : {'Type' : 'Integer', 'Number' : 2},
    'MQ'  : {'Type' : 'Integer', 'Number' : 1},
    'PL'  : {'Type' : 'Integer', 'Number' : 1},
    'PQ'  : {'Type' : 'Integer', 'Number' : 1},
    'PS'  : {'Type' : 'Integer', 'Number' : 1},
    # for SVs
    'CN' : {'Type' : 'Integer', 'Number' : 1},
    'CNQ' : {'Type' : 'Float', 'Number' : 1},
    'CNL' : {'Type' : 'Float', 'Number' : '.'},
    'NQ' : {'Type' : 'Integer', 'Number' : 1},
    'HAP' : {'Type' : 'Integer', 'Number' : 1},
    'AHAP' : {'Type' : 'Integer', 'Number' : 1}
}

class VcfReader(object):
    """ 
        doc TODO
    """
    
    def __init__(self, filename=None, compressed=None, encoding='utf-8'):
        """ 
            Create a new VcfReader object
 
            Opens the given VCF file and stores the metaheader and 
            sample information
        """
        
        if not filename :
            raise ParseError('You must provide a filename')
        self.filename = filename
        if compressed is None:
            compressed = filename.endswith((".gz", ".bgz"))
        if compressed:
            self.file = gzip.open(filename, encoding=encoding, mode='rt')
        else:
            self.file = open(filename, encoding=encoding, mode='r')
        self.reader = (line.rstrip() for line in self.file if line.rstrip())
        self._mode = os.stat(self.filename).st_mode
        #read header information
        self.header      = self._read_header()
        #set some header values for convenience
        self.metadata    = self.header.metadata
        self.col_header  = self.header.col_header
        self.samples     = self.header.samples 
        self.sample_cols = self.header.sample_cols
        self.parser      = (VcfRecord(line, self.sample_cols, self.metadata) 
                                      for line in self.reader)
        self.var         = None
        

    def _read_header(self):
        """ 
            Called after opening VCF. This reads the meta header lines 
            into a list, gets columns names and sample names
        """

        meta_header = []
        coln_header  = []
        for line in self.reader:
            if line.startswith('##'):
                meta_header += [line]
            elif line.startswith('#CHROM'):
                coln_header = line.split("\t")
                break
            else:
                raise HeaderError('No column header found for VCF {}' 
                                  .format(self.filename))
        return VcfHeader(meta_header, coln_header)
     
    def search(self, chrom, start=None, end=None):
        """ 
            Retrieve lines by genomic location
            Sets self.reader to an iterator over the resulting lines
        """
        if not S_ISREG(self._mode):
            raise ParseError("Cannot run search on a non-regular file")
        pass
        #TODO
    
   
class VcfHeader(object):
    ''' Header class storing metadata and sample information for a vcf '''
    
    _meta_re = re.compile(r"""\#\#(\S+?)=(.*)""")#should capture any metadata
    _dict_re = re.compile(r"""                 #for capturing dict-key metadata
                          \#\#(\S+)            #captures field name (e.g. INFO)
                          =<ID=(\S+?)          #captures metadata ID
                          (,(.*))*             #capture remaining keys/values
                          >""",                #dict line should end with a >
                          re.X)
    _subd_re = re.compile(r"""                 #for extracting keys/values from
                                               #group(3) of _dict_re.match()
                          ,(\S+?)=             #get key
                          (".+?"|[^\s,]+)      #value can either be quoted or 
                                               #else should be all non-comma 
                                               #and non-whitespace chars
                         """, re.X)
    _required_keys = { 'INFO'   : ["Number", "Type", "Description"],
                       'FORMAT' : ["Number", "Type", "Description"],
                       'FILTER' : ["Description"],
                       'ALT'    : ["Description"],
                     }

    def __init__(self, meta_header, col_header):
        ''' Requires a list of metaheader lines and a list of column
            names
        '''

        self.meta_header = meta_header
        self.col_header  = col_header
        self.samples     = col_header[9:] or None
        self.sample_cols = dict([(x,i) for (i,x)
                            in enumerate(col_header[9:], start=9)]) or None
        for (i,c) in enumerate(
            "#CHROM POS ID REF ALT QUAL FILTER INFO".split()):
            #check essential column names and order
            if self.col_header[i] != c:
                raise HeaderError('Invalid column name. Expected {}, got {}'
                                  .format(c, self.col_header[i]))
        if (len(self.col_header) > 8):
            #9th column must be FORMAT if present, but is optional
            if self.col_header[8] != 'FORMAT':
                raise HeaderError('Invalid column name. Expected {}, got {}'
                                  .format('FORMAT', self.col_header[8]))
        self.metadata = {}
        self._parse_metadata()


    def _parse_metadata(self):
        ''' 
            Extract INFO, FORMAT, FILTER and contig information from VCF 
            meta header and store in dicts
        '''

        #check first line is essential fileformat line
        ff_match = self._meta_re.match(self.meta_header[0])
        if not ff_match or ff_match.group(1) != 'fileformat':
            raise ParseError("Error: First line of VCF should be fileformat" +
                             "metaheader (e.g. ##fileformat=VCFv4.2)")
        else:
            self.fileformat = ff_match.group(2)
        for h in self.meta_header:
            self._parse_header_line(h)

    def _parse_header_line(self, h):
        ''' 
            Parse a metaheader line and assign to self.metadata dict where
            keys are the type of metadata line and values are dicts of IDs to
            lists of either dicts of key-value pairs or string values.
        '''

        match_d = self._dict_re.match(h)
        match_m = self._meta_re.match(h)
        field = None
        fid   = None
        if match_d:
            #line is an e.g. INFO/FORMAT/FILTER/ALT/contig with multiple keys
            field = match_d.group(1)
            fid   = match_d.group(2)
            rest  = match_d.group(3)
            d = dict([(x, y) for (x, y) in self._subd_re.findall(rest)])
            if not field in self.metadata:
                self.metadata[field] = {fid : [d]}
            else:
                if fid in self.metadata[field]:
                    #multiple values - extend list
                    self.metadata[field][fid].append(d)
                else:
                    self.metadata[field][fid] = [d]
        elif match_m:
            field = match_m.group(1)
            fid   = match_m.group(2)
            if field in self.metadata:
                self.metadata[field].append(fid)
            else:
                self.metadata[field] = [fid]
        else:
            raise HeaderError("Invalid metaheader line {}".format(h))
        if field in self._required_keys:
            #ensure ALT/FORMAT/FILTER/INFO contain required keys
            last = self.metadata[field][fid][-1]#check entry we've just added
            for k in self._required_keys[field]:
                if not k in last:
                    raise HeaderError(
                            "Missing required key '{}' in metaheader line: {}" 
                            .format(k, h))
                

class VcfRecord(object):
    ''' 
        A single record from a Vcf created by parsing a non-header line 
        from a VCF file. May contain multiple alternate alleles.
    '''

    __slots__ = ['cols', '_metadata', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 
                 'QUAL', 'FILTER', 'INFO', 'FORMAT', '_sample_idx', '__CALLS', 
                 '__ALLELES', 'GT_FORMAT', '_SAMPLE_GTS', '_got_gts']

    def __init__(self, line, sample_idx, metadata):
        self.cols = line.split("\t", 9) #only collect first 9 columns initially
                                        #splitting whole line on a VCF with
                                        #lots of columns/samples is costly

        ( self.CHROM, self.POS, self.ID, self.REF, self.ALT,  
          self.QUAL, self.FILTER, self.INFO ) = self.cols[:8] 
        self.FORMAT     = None
        self.GT_FORMAT  = None
        self.CALLS      = None
        self.ALLELES    = None
        self._got_gts    = False         #flag indicating whether we've already 
                                        #retrieved GT dicts for every sample
        self._SAMPLE_GTS = {}
        self._sample_idx = sample_idx
        self._metadata = metadata
        if len(self.cols) > 8:
            self.FORMAT = self.cols[8]
            self.GT_FORMAT = self.FORMAT.split(':')

    @property
    def ALLELES(self):
        ''' list of REF and ALT alleles in order '''

        if self.__ALLELES is None:
            self.__ALLELES = [self.REF] + self.ALT.split(',')
        return self.__ALLELES
        
    @ALLELES.setter
    def ALLELES(self, alleles):
        self.__ALLELES = alleles

    @property
    def CALLS(self):
        '''
            split sample call fields and assign to self.CALLS dict of sample 
            id to call string. 

            self.CALLS does not get created in __init__ to save on overhead
            in VCF with many samples where we might not be interested in
            parsing sample calls 
            
            as of python 3.6 the CALLS dict will maintain sample order, but to 
            safely get a list of calls in the same order as the input VCF the 
            following syntax should be used:

            >>> v = VcfReader(my_vcf)
            >>> for record in v.parser:
            ...     calls = [ record.CALLS[x] for x in v.samples ]

        '''

        if self.__CALLS is None:
            if self._sample_idx:
                calls = self.cols.pop()
                self.cols.extend(calls.split("\t"))
                self.__CALLS = dict([(s, self.cols[self._sample_idx[s]]) 
                                      for s in self._sample_idx]) 
            else:
                self.__CALLS = {}
        return self.__CALLS
    
    @CALLS.setter
    def CALLS(self, calls):
        self.__CALLS = calls
     
    def sample_calls(self):
        ''' 
            Retrieve a dict of sample names to a dict of genotype field
            names to values. All returned values are strings.
            
            >>> d = record.sample_calls()
            >>> s1 = d['Sample_1']
            >>> gq = s1['GQ']

            ...or more concisely:

            >>> gq = record.sample_calls()['Sample_1']['GQ']
        '''

        if self._got_gts:
            return self._SAMPLE_GTS
        else:
            self._got_gts = True
            return dict([(s, self.get_sample_call(s)) 
                         for s in self._sample_idx ])

    def get_sample_call(self, sample):
        ''' 
            Retrieve a dict of genotype field names to values for a 
            single sample.
        
            This method creates dicts as needed for given samples, so 
            may be more efficient than using the 'sample_calls()' 
            method when parsing a VCF with many samples and you are 
            only interested in a information from a small number of 
            these samples.
            
            Args:
                sample: name of the sample to retrieve (as it appears 
                        in the VCF header

            Example: 
                >>> s1 = record.get_sample_call['Sample_1']
                >>> gq = s1['GQ']
        '''

        try:
            return self._SAMPLE_GTS[sample]
        except KeyError:
            if sample in self.CALLS:
                d = dict( [(f, v) for (f, v) in zip(self.GT_FORMAT, 
                                            (self.CALLS[sample].split(':')))] )
                self._SAMPLE_GTS[sample] = d
                return d
            else:
                raise ParseError("Sample {} is not in VCF" .format(sample))

    def get_parsed_gt_fields(self, field, values=[]):
        '''
            Retrieves values of genotype field parsed so that the 
            returned values are of the expected type (str, int or float)
            and are split into list format if appropriate. Fields are 
            handled according to the information present in the VCF 
            header metadata 

            Args:
                field:  genotype field to retrieve as it appears in the 
                        FORMAT field of the VCF record and in the VCF 
                        header.

                values: list of values from a genotype field
        '''
        
        f = self._get_field_translation('FORMAT', field)
        #f[0] is the class type of field, f[1] = True if values should be split
        pv = []
        for val in values:
            try: 
                if f[1]:
                    l = list(f[0](s) for s in val.split(','))
                    pv.append(l)
                else:
                    pv.append(f[0](val))
            except ValueError:
                if val == '.':
                    if f[1]:
                        pv.append([None])
                    else:
                        pv.append(None)
                else:
                    raise ParseError("Unexpected value type for " +
                                     "{} ".format(field) + "FORMAT field at " +
                                     " {}:{}" .format(self.CHROM, self.POS))
        return pv

    def _get_field_translation(self, field_type, field):
        '''
            returns a tuple of variable class type (int, float or str)
            and whether the value requires splitting
        '''

        try:
            f = self._metadata[field_type][field][0]
        except KeyError:
            try:
                f = COMMON_FORMAT[field]
            except KeyError:
                raise ParseError("Unrecognised FORMAT field '{}'"
                                 .format(field) + "at {}:{}. " 
                                 .format(self.CHROM, self.POS) + 
                                 "Non-standard  FORMAT fields should be " + 
                                 "represented in VCF header.")
        ctype = None
        if f['Type'] == 'String' or f['Type'] == 'Character':
            ctype = str
        elif f['Type'] == 'Float':
            ctype = float
        elif f['Type'] == 'Integer':
            ctype = int
        else:
            raise ParseError("Unrecognised FORMAT Type '{}' at {}:{}" 
                             .format(f['Type'], self.CHROM, self.POS))
        split = False
        try:
            if int(f['Number']) > 1:
                split = True
        except ValueError:
            split = True
        return (ctype, split)


class HeaderError(Exception):
    pass


class ParseError(Exception):
    pass
