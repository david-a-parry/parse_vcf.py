import os
import gzip
import re
import warnings
from stat import S_ISREG
try:
    import pysam
except ImportError:
    pysam = None


#common INFO fields and their types in case absent from header
COMMON_INFO = {
    '1000G':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'AA':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'AC':{'Type':'Integer', 'Class':int, 'Number':'A', 'Split':True},
    'AF':{'Type':'Float', 'Class':float, 'Number':'A', 'Split':True},
    'AN':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'BQ':{'Type':'Float', 'Class':float, 'Number':1, 'Split':False},
    'CIGAR':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'DB':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'DP':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'END':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'H2':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'H3':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'MQ':{'Type':'Float', 'Class':float, 'Number':1, 'Split':False},
    'MQ0':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'NS':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'SB':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'SOMATIC':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'VALIDATED':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    # for SVs
    'BKPTID':{'Type':'String', 'Class':str, 'Number':'.', 'Split':False},
    'CICN':{'Type':'Integer', 'Class':int, 'Number':2, 'Split':True},
    'CICNADJ':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'CIEND':{'Type':'Integer', 'Class':int, 'Number':2, 'Split':True},
    'CILEN':{'Type':'Integer', 'Class':int, 'Number':2, 'Split':True},
    'CIPOS':{'Type':'Integer', 'Class':int, 'Number':2, 'Split':True},
    'CN':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'CNADJ':{'Type':'Integer', 'Class':int, 'Number':'.', 'Split':False},
    'DBRIPID':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'DBVARID':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'DGVID':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'DPADJ':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'EVENT':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'HOMLEN':{'Type':'Integer', 'Class':int, 'Number':'.', 'Split':False},
    'HOMSEQ':{'Type':'String', 'Class':str, 'Number':'.', 'Split':False},
    'IMPRECISE':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'MATEID':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'MEINFO':{'Type':'String', 'Class':str, 'Number':4, 'Split':True},
    'METRANS':{'Type':'String', 'Class':str, 'Number':4, 'Split':True},
    'NOVEL':{'Type':'Flag', 'Class':None, 'Number':1, 'Split':False},
    'PARID':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'SVLEN':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'SVTYPE':{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
}

#common genotype FORMAT Fields and their types in case absent from header
COMMON_FORMAT = {
    'DP' :{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'EC' :{'Type':'Integer', 'Class':int, 'Number':'A', 'Split':True},
    'FT' :{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'GL' :{'Type':'Float', 'Class':float, 'Number':'G', 'Split':True},
    'GLE':{'Type':'String', 'Class':str, 'Number':'.', 'Split':False},
    'GP' :{'Type':'Float', 'Class':float, 'Number':'G', 'Split':True},
    'GQ' :{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'GT' :{'Type':'String', 'Class':str, 'Number':1, 'Split':False},
    'HQ' :{'Type':'Integer', 'Class':int, 'Number':2, 'Split':True},
    'MQ' :{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'PL' :{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'PQ' :{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'PS' :{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    # for SVs
    'CN':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'CNQ':{'Type':'Float', 'Class':float, 'Number':1, 'Split':False},
    'CNL':{'Type':'Float', 'Class':float, 'Number':'.', 'Split':False},
    'NQ':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'HAP':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False},
    'AHAP':{'Type':'Integer', 'Class':int, 'Number':1, 'Split':False}
}

class VcfReader(object):
    """ 
        doc TODO - individual methods documentation available
    """
    
    def __init__(self, filename=None, compressed=None, bcf=None, 
                 encoding='utf-8'):
        """ 
            Create a new VcfReader object
 
            Opens the given VCF file and stores the metaheader and 
            sample information

            Args:
                filename:   VCF file to open. Required.
                
                compressed: Boolean indicating whether input should be treated 
                            as bgzip compressed data. If not provided, this 
                            will be inferred from the file extension.

                bcf:        Boolean indicating whether the input is from a BCF 
                            file. If not provided, this will be inferred from 
                            the file extension.

                encoding:   Encoding of input data. Default = 'utf-8'.

        """
        
        if not filename :
            raise ParseError('You must provide a filename')
        self.filename = filename
        self.compressed = compressed
        self.bcf = bcf
        self.encoding = encoding
        self._tabix = None
        if self.bcf is None:
            self.bcf = filename.endswith(".bcf")
        if self.compressed is None:
            self.compressed = filename.endswith((".gz", ".bgz"))
        if self.bcf:
            if pysam is None:
                raise ParseError("pysam not available. Please install (e.g. " + 
                                 "via 'pip install pysam' to parse bcf files")
            self.file = pysam.VariantFile(filename)
            self.reader = (rec.__str__().rstrip() for rec 
                           in self.file.fetch() if rec.__str__().rstrip())
            head = self.file.header.__str__().split("\n")
            head.pop()
            cols = head.pop().split("\t")
            self.header = VcfHeader(head, cols)
        else:
            if self.compressed:
                self.file = gzip.open(filename, encoding=encoding, mode='rt')
            else:
                self.file = open(filename, encoding=encoding, mode='r')
            self.reader = (line.rstrip() for line in self.file if line.rstrip())
            self.header      = self._read_header()
        self._mode = os.stat(self.filename).st_mode
        #read header information
        #set some header values for convenience
        self.metadata    = self.header.metadata
        self.col_header  = self.header.col_header
        self.parser      = (VcfRecord(line, self,) for line in self.reader)
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
     
    def set_region(self, chrom, start=None, end=None):
        """ 
            Retrieve records by genomic location rather than reading 
            records line-by-line.

            Sets self.reader and self.parser to iterators over the 
            records retrieved.

            Args:
                chrom: chromosome or contig name. Required. You may 
                       instead provide a region in format 
                       'chr:start-end' if preferred instead of using 
                       start/end arguments below. 

                start: start position on chromosome/contig. 0-based
                       Default = None

                end:   end position on chromosome/contig. 
                       Default = None

            >>> v = VcfReader(my_vcf)
            >>> v.set_region('chr1') #get all variants on chr1
            >>> for record in v.parser:
            ... #do something with each record

            Because start coordinates are 0-based, to retrieve a variant
            at (or overlapping) chr1:1000000 the two examples below are
            equivalent:

            >>> v.set_region(chrom='chr1', start=999999 end=1000000)
            >>> v.set_region('chr1:1000000-1000000')

        """
        if not S_ISREG(self._mode):
            raise ParseError("Cannot run set_region() on a non-regular file")
        if (self.compressed):
            if not pysam:
                raise ParseError("pysam not available. Please install (e.g. " + 
                                 "via 'pip install pysam' to search by " + 
                                 "location on bgzip compressed VCFs.")
            idx = self.filename + '.tbi'
            if not os.path.isfile(idx):       #create index if it doesn't exist
                pysam.tabix_index(self.filename, preset="vcf")
            if not self._tabix:
                self._tabix = pysam.Tabixfile(self.filename, 
                                              encoding=self.encoding)
            try:
                self.reader = self._tabix.fetch(str(chrom), start, end)
                self.parser = (VcfRecord(line, self,) for line in self.reader)
            except ValueError as oops:
                contig = str(chrom).split(':')[0] #chrom could be a region or 
                                                  #just an int
                if contig not in self._tabix.contigs:
                    warnings.warn("Contig '{}' not found in VCF '{}'" 
                                  .format(contig, self.filename))
                else:
                    raise oops
        else:
            #easy solution - compress and index with pysam if not compressed,
            #but will be slow...
            #less easy solution, implement a custom index to seek to and from
            raise ParseError("Searching by location is not yet implemented " +
                             "for non-bgzip compressed VCFs.")


class VcfHeader(object):
    ''' Header class storing metadata and sample information for a vcf '''
    
    _meta_re = re.compile(r"""\#\#(\S+?)=(.*)""")  #should capture any metadata

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

    _csq_format_re = re.compile(r'''.*Format:\s*((\S+\|)*\S+)"''')
    #for capturing CSQ format in Description field of metaheader

    _required_keys = { 'INFO'   : ["Number", "Type", "Description"],
                       'FORMAT' : ["Number", "Type", "Description"],
                       'FILTER' : ["Description"],
                       'ALT'    : ["Description"],
                     }
    __slots__ = ['meta_header', 'col_header', 'samples', 'sample_cols',
                 'metadata', 'fileformat', '__csq_label', '__csq_fields', 
                 '_info_field_translater', '_format_field_translater'] 


    def __init__(self, meta_header, col_header):
        ''' 
            Requires a list of metaheader lines and a list of column
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
        self.csq_label = None
        self.csq_fields = None
        self._info_field_translater = {}
        self._format_field_translater = {}
        self._parse_metadata()

    @property
    def csq_label(self):
        ''' 
            String labelling the INFO field label of VEP consequence 
            annotations. Will raise a HeaderError if access is attempted 
            but no VEP CSQ or ANN field is present in the header.
        '''
        if self.__csq_label is None:
            self.csq_fields
        return self.__csq_label

    @csq_label.setter
    def csq_label(self, c):
        self.__csq_label = c

    @property
    def csq_fields(self):
        ''' 
            A list of CSQ field names in the order they are represented
            in CSQ INFO field entries. Set to None on initialization. 
            Will raise a HeaderError if access is attempted but no VEP 
            CSQ or ANN field is present in the header.
        '''

        if self.__csq_fields is None:
            try:
                csq_header = self.metadata['INFO']['CSQ'][-1]
                csq = 'CSQ'
            except KeyError:
                try:
                    csq_header = self.metadata['INFO']['ANN'][-1]
                    csq = 'ANN'
                except KeyError:
                    raise HeaderError("No CSQ or ANN field in INFO header - "+
                                      "unable to retrieve consequence fields.")
            self.csq_label = csq
            match = self._csq_format_re.match(csq_header['Description'])
            if match:
                self.__csq_fields = match.group(1).split('|')
            else:
                raise HeaderError("Could not parse {} Format in ".format(csq)
                                + "header. Unable to retrieve consequence "
                                + "annotations.")
        return self.__csq_fields
    
    @csq_fields.setter
    def csq_fields(self, csq):
        self.__csq_fields = csq

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

        for field_type in ['FORMAT', 'INFO']:
            try:
                for field in self.metadata[field_type]:
                    self._set_field_translation(field_type, field)
            except KeyError:
                warnings.warn("No '{}' field in header!".format(field_type))

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
            rest  = match_d.group(3) or ''
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

    def _set_field_translation(self, field_type, field):
        '''
            returns a tuple of variable class type (int, float or str)
            and whether the value requires splitting
        '''

        f = self.metadata[field_type][field][0]
        ctype = None
        if f['Type'] == 'String' or f['Type'] == 'Character':
            ctype = str
        elif f['Type'] == 'Float':
            ctype = float
        elif f['Type'] == 'Integer':
            ctype = int
        elif f['Type'] != 'Flag':
            raise ParseError("Unrecognised FORMAT Type '{}' in header" 
                             .format(f['Type']))
        split = False
        if f['Number'].isdigit():
            if int(f['Number']) > 1:
                split = True
        else:              #if not digit should be 'A', 'G', 'R' or '.' - split
            split = True
        if field_type == 'INFO':
            setter = self._info_field_translater    
        elif field_type == 'FORMAT':
            setter = self._format_field_translater
        else:
            raise ParseError("'{}' not recognised as a ".format(field_type) +
                             "field type for translation")
        setter[field] = (ctype, split)


                

class VcfRecord(object):
    ''' 
        A single record from a Vcf created by parsing a non-header line 
        from a VCF file. May contain multiple alternate alleles.
    '''

    _svalt_re = re.compile(r'''<(\w+)(:\w+)*>''') #group 1 gives SV type

    __slots__ = ['header', 'cols', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 
                 'FILTER', 'INFO', 'FORMAT', '__CSQ', 'samples', '_sample_idx',  
                 '__CALLS', '__ALLELES', '__DECOMPOSED_ALLELES', 
                 '__INFO_FIELDS',  'GT_FORMAT', '_SAMPLE_GTS', '_got_gts', 
                 '_vep_allele', ]

    def __init__(self, line, caller):
        ''' 
            VcfRecord objects require a line and a related VcfReader 
            object for initialization.

            Args:
                line:   a non-header line from a VCF file without 
                        newline character
                
                caller: a VcfReader object (normally the same VcfReader 
                        object that read the input line). Metadata and 
                        sample information will be read from this object 
                        in order to initialize the VcfRecord object.

        '''

        self.cols = line.split("\t", 9) #only collect first 9 columns initially
                                        #splitting whole line on a VCF with
                                        #lots of columns/samples is costly

        ( self.CHROM, self.POS, self.ID, self.REF, self.ALT,  
          self.QUAL, self.FILTER, self.INFO ) = self.cols[:8] 
        self.INFO_FIELDS        = None
        self.FORMAT             = None
        self.GT_FORMAT          = None
        self.CALLS              = None
        self.DECOMPOSED_ALLELES = None
        self.ALLELES            = None
        self.header             = caller.header
        self.CSQ                = None
        self._SAMPLE_GTS        = {}
        self._vep_allele        = {}
        self._got_gts           = False       #flag indicating whether we've already 
                                         #retrieved GT dicts for every sample

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
    def DECOMPOSED_ALLELES(self):
        ''' 
            list of AltAllele objects, one for each ALT allele in 
            order, after reducing them to their minimal representations
            (i.e. by trimming redundant nucleotides).
        '''

        if self.__DECOMPOSED_ALLELES is None:
            self._minimizeAlleles()
        return self.__DECOMPOSED_ALLELES
        
    @DECOMPOSED_ALLELES.setter
    def DECOMPOSED_ALLELES(self, alleles):
        self.__DECOMPOSED_ALLELES = alleles

    def _minimizeAlleles(self):
        self.DECOMPOSED_ALLELES = []
        for alt in self.ALLELES[1:]:
            ref = self.ALLELES[0]
            pos = self.POS
            while len(ref) > 1 and len(alt) > 1:
                if ref[-1] == alt[-1]:               #remove identical suffixes
                    ref = ref[0:-1]
                    alt = alt[0:-1]
                else:
                    break
            while len(ref) > 1 and len(alt) > 1:
                if ref[0] == alt[0]:                 #remove identical prefixes
                    ref = ref[1:]
                    alt = alt[1:]
                    pos += 1
                else:
                    break
            self.DECOMPOSED_ALLELES.append(AltAllele(chrom=self.CHROM, 
                                            pos=pos, ref=ref, alt=alt))
     
                

    @property 
    def INFO_FIELDS(self):
        ''' 
        A dict of INFO field names to values. All returned values are Strings 
        except for Flags which are assigned None.

        To obtain values parsed into the appropriate Type as defined by the 
        VCF header, use the 'parsed_info_fields()' method.
        '''
        if self.__INFO_FIELDS is None:
            self.__INFO_FIELDS = {}
            for i in self.INFO.split(';'):
                try:
                    (f, v) = i.split('=', 1)
                except ValueError:
                    (f, v) = (i, None)
                self.__INFO_FIELDS[f] = v
        return self.__INFO_FIELDS

    @INFO_FIELDS.setter
    def INFO_FIELDS(self, i):
        self.__INFO_FIELDS = i

    def parsed_info_fields(self, fields=None):
        if fields is not None:
            f_list = fields
        else:
            f_list = list(self.INFO_FIELDS)
        d = dict( (f, self._get_parsed_info_value(f, self.INFO_FIELDS[f])) 
                                                               for f in f_list)
        return d
            
    def _get_parsed_info_value(self, field, value):
        try:
            f = self.header._info_field_translater[field]
        except KeyError:
            try:
                f = (COMMON_INFO[field]['Class'], 
                     COMMON_INFO[field]['Split'])
                self.header._info_field_translater[field] = f
            except KeyError:
                raise ParseError("Unrecognised INFO field '{}'".format(field) 
                                 + "at {}:{}. ".format(self.CHROM, self.POS) + 
                                 "Non-standard  INFO fields should be " + 
                                 "represented in VCF header.")
        #f[0] is the class type of field, f[1] = True if values should be split
        try: 
            if (f[1]):
                pv = list(map(f[0], value.split(',')))
            else:
                pv = f[0](value)
        except (ValueError, TypeError, AttributeError):
            if value == '.' or value is None:
                pv = None
            else:
                raise ParseError("Unexpected value type for " +
                                 "{} ".format(field) + "FORMAT field at " +
                                 " {}:{}" .format(self.CHROM, self.POS))
        return pv
            
    @property
    def CALLS(self):
        '''
            split sample call fields and assign to self.CALLS dict of
            sample id to call string. 

            self.CALLS does not get created in __init__ to save on 
            overhead in VCF with many samples where we might not be 
            interested in parsing sample calls 
            
            As of python 3.6 the CALLS dict will maintain sample order,
            but to safely get a list of calls in the same order as the 
            input VCF the following syntax should be used:

            >>> v = VcfReader(my_vcf)
            >>> for record in v.parser:
            ...     calls = [record.CALLS[x] for x in record.samples]

        '''

        if self.__CALLS is None:
            if self.header.sample_cols:
                calls = self.cols.pop()
                self.cols.extend(calls.split("\t"))
                self.__CALLS = dict([(s, self.cols[self.header.sample_cols[s]]) 
                                      for s in self.header.samples]) 
            else:
                self.__CALLS = {}
        return self.__CALLS
    
    @CALLS.setter
    def CALLS(self, calls):
        self.__CALLS = calls
     
    def sample_calls(self):
        ''' 
            Retrieve a dict of sample names to a dict of genotype field
            names to values. All returned values are strings. For 
            values cast to the appropriate types (int/float/string/list)
            use the 'parsed_gts' function.
            
            >>> record.sample_calls()
            {'Sample_1': {'GT': '0/0', 'AD': '10,0', 'DP': '10', 'GQ': '30'},
            'Sample_2': {'GT': '0/1', 'AD': '6,6', 'DP': '12', 'GQ': '33'}}

            >>> d = record.sample_calls()
            >>> s1 = d['Sample_1']
            >>> s1['GQ']
            '30'
            
            ...or more concisely:

            >>> record.sample_calls()['Sample_1']['GQ']
            '30'
        '''

        if self._got_gts:
            return self._SAMPLE_GTS
        else:
            self._got_gts = True
            #get_sample_call() sets self._SAMPLE_GTS[s] for future use
            return dict([(s, self.get_sample_call(s)) 
                          for s in self.header.samples ]) 



    def get_sample_call(self, sample):
        ''' 
            Retrieve a dict of genotype field names to values for a 
            single sample.  
       
            This method creates dicts as needed for given samples, so 
            may be more efficient than using the 'sample_calls()' 
            method when parsing a VCF with many samples and you are 
            only interested in a information from a small number of 
            these samples.
            
            All values returned are strings. For values cast to an 
            appropriate type (int/float/string/list) use the 
            'parsed_gts(sample=[sample])' function.
           
            Args:
                sample: name of the sample to retrieve (as it appears 
                        in the VCF header

            Example: 
                >>> s1 = record.get_sample_call('Sample_1')
                >>> s1
                {'GT': '0/0', 'AD': '10,0', 'DP': '10', 'GQ': '30'}
                >>> s1['GQ']
                '30'
        '''

        try:
            return self._SAMPLE_GTS[sample]
        except KeyError:
            if sample in self.CALLS:
                d = dict( [(f, v) for (f, v) in zip(self.GT_FORMAT, 
                                            (self.CALLS[sample].split(':')))] )
                self._SAMPLE_GTS[sample] = d
                if len(self._SAMPLE_GTS) == len(self.header.samples):
                    #if we've now set self._SAMPLE_GTS for all samples set
                    #self._got_gts to True to prevent unnecessary looping in
                    #sample_calls() method
                    self._got_gts = True
                return d
            else:
                raise ParseError("Sample {} is not in VCF" .format(sample))

    def parsed_gts(self, samples=None, fields=None):
        ''' Returns a dict of GT field names to dicts of sample names 
            to values. By default, values for all samples and fields 
            will be retrieved, but a list of sample IDs and a list of 
            FORMAT fields to retrieve can be given.

            Missing values will be assigned to None.
    
            Because this is a relatively costly function, you are 
            advised to avoid calling this repeatedly for a single 
            record - you may speed things up by only calling for a 
            subset of samples and fields but in any case you probably 
            want to call this function once only per record, storing 
            the results in a variable.

            Args:
                samples: Optional list of sample names to retrieve 
                         values for. Default = None (values retrieved 
                         for all samples)
                                
                fields:  Optional list of field names (as they appear 
                         in the FORMAT field of the record) to retrieve 
                         values for. Default = None (values retrieved 
                         for all fields)
                                
            
            >>> record.parsed_gts()
            {'GT': {'Sample_1': '0/0', 'Sample_2': '0/1'},
            'AD': {'Sample_1': [10, 0], 'Sample_2': [6, 6]},
            'DP': {'Sample_1': 10, 'Sample_2': 12},
            'GQ': {'Sample_1': 30, 'Sample_2': 33}}
            
            >>> record.parsed_gts(samples=['Sample_2'])
            {'GT': {'Sample_2': '0/1'}, 'AD': {'Sample_2': [6, 6]},
            'DP': {'Sample_2': 12}, 'GQ': {'Sample_2': 33}}

            >>>  record.parsed_gts(fields=['GT', 'GQ'])
            {'GT': {'Sample_1': '0/0', 'Sample_2': '0/1'},
            'GQ': {'Sample_1': 30, 'Sample_2': 33}}

        '''
    
        d = {}
        if fields is not None:
            f_list = fields
        else:
            f_list = self.GT_FORMAT
        if samples is not None:
            s_list = samples
        else:
            s_list = self.header.samples
        for f in f_list:
            d[f] = dict(zip(s_list, self._get_parsed_gt_fields(f,
                        (self.sample_calls()[s][f] if f in 
                        self.sample_calls()[s] else None 
                        for s in s_list) ) ) )
        return d
        
    def _get_parsed_gt_fields(self, field, values=[]):
        #TODO - make this more efficient with a cython extension
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
        
        try:
            f = self.header._format_field_translater[field]
        except KeyError:
            try:
                f = (COMMON_FORMAT[field]['Class'], 
                     COMMON_FORMAT[field]['Split'])
                self.header._format_field_translater[field] = f
            except KeyError:
                raise ParseError("Unrecognised FORMAT field '{}'".format(field) 
                                 + "at {}:{}. ".format(self.CHROM, self.POS) + 
                                 "Non-standard  FORMAT fields should be " + 
                                 "represented in VCF header.")
        #f[0] is the class type of field, f[1] = True if values should be split
        pv = []
        for val in values:
            try:
                if f[1]:
                    pv.append(list(map(f[0], val.split(','))))
                else:
                    pv.append(f[0](val))
            except (ValueError, TypeError, AttributeError) as err:
                if val is None or val == '.':
                    if f[1]:
                        pv.append([None])
                    else:
                        pv.append(None)
                else:
                    raise err
                    raise ParseError("Unexpected value ('{}') ".format(val) 
                                    +"for {} " .format(field) + "FORMAT field "  
                                    +"at {}:{}" .format(self.CHROM, self.POS))
        return pv


    @property
    def CSQ(self):
        ''' 
            A list of dicts of CSQ/ANN annotations from VEP to values.
            Empty values are represented by empty Strings. Will raise 
            a HeaderError if the associated VCF header does not contain
            CSQ/ANN information and a ParseError if the record being 
            parsed does not contain a CSQ/ANN annotation in the INFO 
            field.
        '''
        if self.__CSQ is None:
            lbl = self.header.csq_label
            try:
                csqs = self.INFO_FIELDS[lbl].split(',')
            except KeyError:
                raise ParseError("Could not find '{}' label in ".format(lbl) +
                                 "INFO field of record at {}:{}"
                                 .format(self.CHROM, self.POS))
            self.__CSQ = []
            alleleToNum = {}
            for c in csqs:
                d = dict([(k,v) for (k, v) in zip(self.header.csq_fields, 
                                                              c.split('|'))]) 
                if 'ALLELE_NUM' in d:
                    d['alt_index'] = d['ALLELE_NUM']
                else:
                    d['alt_index'] = self._vepToAlt(d['Allele'])
                self.__CSQ.append(d)
        return self.__CSQ

    @CSQ.setter
    def CSQ(self, c):
        self.__CSQ = c
        
    def _vepToAlt(self, allele):
        #figure out how alleles will be handled by looking at the REF vs ALTs
        if self._vep_allele:
            return self._vep_allele[allele]
        is_sv = False
        is_snv = False
        is_indel = False
        is_mnv = False
        ref = self.ALLELES[0]
        for i in range(1, len(self.ALLELES)):
            alt = self.ALLELES[i]
            if alt == '*':
                self._vep_allele[alt] = i
            else:
                matches_sv = self._svalt_re.match(alt)
                if matches_sv:
                    is_sv = True
                    sv_type = matches_sv.group(1)
                    if sv_type == 'DUP':
                        self._vep_allele['duplication'] = i
                    elif sv_type == 'INS':
                        self._vep_allele['insertion'] = i
                    elif sv_type == 'DEL':
                        self._vep_allele['deletion'] = i
                    else:
                        self._vep_allele[sv_type] = i #should catch CNVs
                else:
                    if len(alt) == 1 and len(ref) == 1:
                        is_snv = True
                    elif len(alt) == len(ref):
                        is_mnv = True
                    else:
                        is_indel = True
                    self._vep_allele[alt] = i
                
        if is_sv:
            #no more editing required as long as 
            #not at the same site as short variant
            if is_snv or is_mnv or is_indel:
                raise ParseError("Unable to parse structural variants at the "
                               + "same site as a non-structural variant")
        else:
            #if is_snv or (is_mnv and not is_indel):
                #do nothing - VEP does not trim alleles in this situation
            if not is_snv and is_indel:
                #VEP trims first base unless REF and ALT differ at first base
                first_base_differs = False
                ref_start = ref[:1]
                for alt in self.ALLELES[1:]:
                    if alt != '*':
                        alt_start = alt[:1]
                        if alt_start != ref_start:
                            first_base_differs = True 
                            break
                if not first_base_differs:
                    #no trimming if first base differs for any ALT, 
                    #otherwise first base is trimmed
                    trimmed = {}
                    pop = []
                    for alt in self._vep_allele:
                        if alt != '*':
                            i = self._vep_allele[alt]
                            pop.append(alt)
                            if len(alt) > 1:
                                alt = alt[1:]
                            else:
                                alt = '-'
                            trimmed[alt] = i
                    for p in pop:
                        self._vep_allele.pop(p, None)
                    self._vep_allele.update(trimmed)
        return self._vep_allele[allele]


class AltAllele(object):
    '''
        Represents basic genomic features of a single alternative 
        allele call. Features are 'CHROM', 'POS', 'REF' and 'ALT'.
    '''
    
    __slots__ = ['CHROM', 'POS', 'REF', 'ALT']

    def __init__(self, record=None, allele_index=1, chrom=None, pos=None,
                 ref=None, alt=None):
        '''
            Either created from a given VcfRecord and the index of the 
            allele to be represented or from chrom, pos, ref and alt 
            arguments.

            Args:
                record:       VcfRecord object containing the allele of     
                              interest. Uses 'allele_index' argument to
                              determine the allele to represent.

                allele_index: index of the allele to represent (e.g. 1 
                              for the first ALT allele, 2 for the 
                              second or 0 for the REF allele). 
                              Default = 1.

                chrom:        chromosome/contig (required if not using 
                              a VcfRecord for construction, ignored 
                              otherwise).

                pos:          position of ref allele (required if not 
                              using a VcfRecord for construction, 
                              ignored otherwise).

                ref:          reference allele (required if not using a
                              VcfRecord for construction, ignored 
                              otherwise).

                alt:          alternative allele (required if not using 
                              a VcfRecord for construction, ignored 
                              otherwise).

        '''

        if record is not None:
            self.CHROM = record.CHROM
            self.POS   = record.POS    
            self.REF   = record.REF    
            self.ALT   = record.ALT
        else:
            self.CHROM = chrom
            self.POS   = pos    
            self.REF   = ref    
            self.ALT   = alt


class HeaderError(Exception):
    pass


class ParseError(Exception):
    pass
