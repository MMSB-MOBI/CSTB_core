import hashlib, gzip, os, re, shutil
from Bio import SeqIO
# Expected format

#fasta   taxid   name    gcf     accession  ftp


def tsvReader(tsvFilePath, _min=None, _max=None):
    with open(tsvFilePath) as f : 
        f.readline()
        i = 0
        for l in f:            
            l_split = list( filter(lambda x:  x != '', l.strip("\n").split("\t") ) )
            if len(l_split) < 5:
                raise ValueError(f"Current tsv record length is less than 5 ({len(l_split)})\n=>{l_split}") 
            fasta = l_split[0]
            if not fasta:
                raise FormatError("genomes list, first column (fasta) is empty")  
            taxid = l_split[1]
            if  taxid == '-':
                taxid = None
            else:
                taxid = int(taxid)
            name = l_split[2]
            if not name:
                raise FormatError("genomes list, third column (name) is empty")
            gcf = l_split[3]
            if  gcf == '-':
                gcf = None
            acc = l_split[4]
            if acc == '-':
                acc = None

            t = (fasta, name, taxid, gcf, acc)           
            if not _min is None:
                if not _max is None:
                    if i >= _min and i <= _max:                       
                        yield (t) 
            elif not _max is None:
                if i <= _max:
                    yield (t) 
            else : 
                yield (t)
            i += 1

''' ZIP and straight  OPEN delivers SIMILAR MD5
>>> import gzip
>>> import hashlib
>>> fp = gzip.open('test.gz', 'rb')
>>> hasher = hashlib.md5()
>>> buf = fp.read()
>>> hasher.update(buf)
>>> hasher.hexdigest()
'a02b540693255caec7cf9412da86e62f'

>>> import hashlib
>>> fp = open('test.txt', 'rb')
>>> hasher = hashlib.md5()
>>> buf = fp.read()
>>> hasher.update(buf)
>>> hasher.hexdigest()
'a02b540693255caec7cf9412da86e62f'
'''
MAX_CHAR = 100000
## Compute hash on fasta striping: header-line, spaces, and new-line
def fileHash(filePath, noHeader=True, stripinSpace=True):
    global MAX_CHAR
    hasher = hashlib.md5()
    with Zfile(filePath, 'rb') as f:
        if noHeader:
            f.readline()# discard header
        while True:
            buf = f.read(MAX_CHAR)
            if not buf:
                break
            if stripinSpace:
                buf = stripBytes(buf)
            hasher.update(buf)

    return hasher.hexdigest()

def stripBytes(byte):
    _str = str(byte, 'utf-8')# convert byte to string for striping               
    _strStriped = re.sub( r'\s+', '',  _str)
    return _strStriped.encode('utf-8')# encode back to byte

def hashStripedBytes(byte):
    hasher = hashlib.md5()
    buf = stripBytes(byte)
    hasher.update(buf)

    return hasher.hexdigest()

def hashStripedString(string):
    hasher = hashlib.md5()
    _strStriped = re.sub( r'\s+', '',  string)
    buf = stripBytes(_strStriped.encode('utf-8'))
    hasher.update(buf)

    return hasher.hexdigest()

class Zfile(object):
    def __init__(self, filePath, mode='r'):
        self.file_obj = zOpen(filePath, mode)
    def __enter__(self):
        return self.file_obj
    def __exit__(self, type, value, traceback):
        self.file_obj.close()

# return a stream using base name of gzip extension:
def zOpen(filePath, mode='r'):
    m = mode
    mz = 'rt' if mode == 'r' else 'rb'
    if filePath.endswith(".gz"):
        fp = gzip.open(filePath, mz)
        return fp

    try:         
        fp = open(filePath, m)
        return fp
    except (OSError, IOError) as e:
        fp = gzip.open(filePath + '.gz', mz)
        return fp


def zExists(filePath):
    return os.path.isfile(filePath) or os.path.isfile(filePath + '.gz')


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def fileToGunzip(flatFile):
    targetFile = f"{flatFile}.gz"
    with open(flatFile, 'rb') as f_in, gzip.open(targetFile, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return targetFile

def gunzipToFile(gzipedFile):
    targetFile = gzipedFile.replace('.gz', '')
    with gzip.open(gzipedFile, 'rb') as f:
        file_content = f.read()
        with open(targetFile, 'wb') as fp:
            fp.write(file_content)
    return targetFile

# yield (full_header, sequence, header_id)
def zFastaReader(filePath):
     with Zfile(filePath) as handle:
        for genome_seqrecord in SeqIO.parse(handle, "fasta"):
            genome_seq = genome_seqrecord.seq
            ref = genome_seqrecord.id
            header = genome_seqrecord.description
            yield(str(header), str(genome_seq), str(ref))



def sgRNAIndexWriter(data, fname, wLen, mode='w'):
    bOcc = None
    with open(fname, mode) as fp:
        fp.write(f"# {len(data)} {wLen}\n")
        for datum in data:
            if bOcc is None:
                bOcc = guessOccurenceFormat(datum)
            if bOcc:
                fp.write( ' '.join([str(d) for d in datum]) + "\n")
            else:
                fp.write(str(datum) + "\n")
    return len(data)

def sgRNAIndexReader(indexFilePath):
    skipFirst = True
    with open(indexFilePath, 'r') as fp:    
        for l in fp:
            _ = l.split()
            if skipFirst:
                if not l.startswith("#"):
                    raise IOError(f"Irregular header line in index file \"{l}\"")
                skipFirst = False
                yield ( int(_[2]), None )
                continue
            if len(_) != 1 and  len(_) != 2:
                raise IOError(f"Irregular line in index file \"{l}\"")
            yield ( int(_[0]), None if len(_) == 1 else int(_[1]) )

def sgRNAplainWriter(data, fname, mode= 'w'):
    bOcc = None
    with open(fname, mode) as fp:
        fp.write(str(len(data)) + "\n")

        for datum in data:
            if bOcc is None:
                bOcc = guessOccurenceFormat(datum)
            if bOcc:
                fp.write( ' '.join([str(d) for d in datum]) + "\n")
            else:
                fp.write(str(datum) + "\n")
    return len(data) 


def guessOccurenceFormat(datum):
    if datum[1] is None:
        print("Writing index in no occurence format")
        return False
    
    print("Writing index in occurences format")
    return True
