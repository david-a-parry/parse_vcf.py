import os
import gzip
import re
import warnings
from stat import S_ISREG


class VcfReader(object):
    """ doc TODO
    """
    
    def __init__(self, filename=None, compressed=None, encoding='utf-8'):
        """ Create a new VcfReader object
            Opens the given VCF file and stores the metaheader and sample 
            information
        """
        
        if not filename :
            raise ParseError('You must provide a filename')
        self.filename = filename
        if compressed is None:
            compressed = filename.endswith((".gz", ".bgz"))
        if compressed:
            self.file = gzip.open(filename, encoding=encoding, mode='r')
        else:
            self.file = open(filename, encoding=encoding, mode='r')
        self.reader = (line.rstrip() for line in self.file if line.rstrip())
        self._mode = os.stat(self.filename).st_mode
        #read header information
        self.header      = self._read_header()
        #set some header values for convenience
        self.metadata    = self.header.meta_header
        self.col_header  = self.header.col_header
        self.samples     = self.header.samples 
        self.sample_cols = self.header.sample_cols
        self.var         = None
        

    def _read_header(self):
        """ called after opening VCF. This reads the meta header lines into a
            list, gets columns names and sample names
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
        """ Retrieve lines by genomic location
            Sets self.reader to an iterator over the resulting lines
        """
        if not S_ISREG(self._mode):
            raise ParseError("Cannot run search on a non-regular file")
        pass
        #TODO
    
   
class VcfHeader(object):
    ''' Header class storing metadata and sample information for a vcf
    '''
    
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
    #_alt_re = re.compile(r"""\#\#ALT=<ID=(\w+),Description=(.*)
    #                     (,(\w+)=(.*))*""", re.X)
    #_info_re = re.compile(r"""\#\#INFO=<ID=(\w+),Number=([\w\d\.]),Type=(\w+),
    #                      Description=(.*)(,(\w+)=(.*))*""", re.X)
    #_format_re = re.compile(r"""\#\#FORMAT=<ID=(\w+),Number=([\w\d\.]),
    #                        Type=(\w+),Description=(.*)(,(\w+)=(.*))*""", re.X)
    #_filter_re = re.compile(r"""\#\#FILTER=<ID=(\w+),Description=(.*)
    #                        (,(\w+)=(.*))*""", re.X)
    #_contig_re = re.compile(r"""\#\#contig=<ID=(\w+)(,(\w+)=(.*))*""", re.X)

    def __init__(self, meta_header, col_header):
        ''' Requires a list of metaheader lines and a list of column names
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
        ''' Extract INFO, FORMAT, FILTER and contig information from VCF 
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
        ''' Parse a metaheader line and assign to self.metadata dict where
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
                
class HeaderError(Exception):
    pass

class ParseError(Exception):
    pass
