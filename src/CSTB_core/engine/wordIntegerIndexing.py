"""
    Translate a dictionary of constant-length CRISPR motifs into an ordered set of integer, using base4 encoding.
    Usage:
        wordIntegerIndexing.py code <pickledDictionary> [--dbase] [--out=<outFile> --occ]
        wordIntegerIndexing.py decode <indexedFile> [--dbase] [--out=<outFile>]
        wordIntegerIndexing.py translate <indexedFile> [--dbase] [--out=<outFile>]

    Options:
        -h --help                               Show this screen
        --occ                                   notify the number of occurence of each word           
        -o <outFile>, --out <outFile>          Name of the output file [default : ./pickledDictionary.index]
        --dbase                                 Switch to power of two encoder, default is twobits
"""

import os
import pickle
import math
import twobits
from CSTB_core.utils.io import sgRNAplainWriter, sgRNAIndexReader, sgRNAIndexWriter
from docopt import docopt

ENCODER=twobits.encode
DECODER=twobits.decode
CODEC='twobits'


def encode(word, codec='twobits'):
    if codec == 'twobits':
        fn = twobits.encode
    elif codec == 'pow2':
        fn = pow2encoderWrapper
    else:
        raise TypeError(f"Unknown codec \"{codec}\"")

    if len(word) > 32:
        raise ValueError(f"Word can't be encoded, too long ({len(word)}> 32).")

    return fn(word)

def decode(value, wLen, codec='twobits'):
    if codec == 'twobits':
        fn = twobits.decode
    elif codec == 'pow2':
        fn = pow2decoderWrapper
    else:
        raise TypeError(f"Unknown codec \"{codec}\"")

    if wLen > 32:
        raise ValueError(f"Word can't be decoded, too long ({wLen}> 32).")
    return fn(value, wLen)

def occWeight(k,datum):
    n = 0
    for o in datum:
        for _o in datum[o]:
            n += len(datum[o][_o])
    #print(k,n)
    return n

# we may add a second field to each wordCode line, the occurence number
def indexAndMayOccurence(data, occ=True):
    global ENCODER
    word_list = list(data.keys())
    wLen = len(word_list[0])
    print(f"Encoding what seems like {wLen} length words")
    
    indexData = []
    for w in word_list:
        if len(w) != wLen:
            raise IOError(f"Not even size word to encode {len(w)} != {wLen}")
        _ = ( ENCODER(w), occWeight(w, data[w]) ) if occ else ENCODER(w)
        indexData.append( _ )
        
    indexData =  sorted(indexData , key=lambda x: x[0] ) if occ else sorted(indexData)
    return indexData, wLen

# Read a pickle sgRNA dictionary encodes its keys   
def indexFromPickle(file_path, occ=True):
    p_data = pickle.load(open(file_path, "rb"))
    data, wLen = indexAndMayOccurence(p_data,  occ=occ)

    return data, wLen

def translate(data):
    global ENCODER   
    return [ (  ENCODER(word), mayOcc ) for word, mayOcc in data ]


def pow2encoderWrapper(word):
    return weightWord( word, "ATCG", length=len(word) )

def weightWord(word, alphabet, length=None):
    """
    Code a word by base len(alphabet) and return this int
    """
    rank = 0
    if length:
        if length != len(word):
            raise ValueError("Irregular word length " + str(len(word)) +
                             " (expected " + str(length) + ")")
    for i, letter in enumerate(word[::-1]):
        wei = alphabet.index(letter)
        base = len(alphabet)
        rank += wei * pow(base, i)
    return rank

def project(value, lenFrom, lenTo, alphabet="ATCG"):
    base = len(alphabet)
    _value = value
    offset = 0
    for i in range(lenFrom, lenTo-1, -1):
        w = math.trunc(value / pow(base, i))
        offset += w * pow(base, i)
        value = value % pow(base, i)
    
    return _value - offset

def pow2decoderWrapper(code, wLen):
    return pow2Decoder(code, "ATCG", length=wLen)

def pow2Decoder(rank, alphabet, length=20):
    """
    Decode the rank (int) to a word according to an alphabet given
    """
    word = ""
    base = len(alphabet)
    for i in range(length - 1, -1, -1):
        index = math.trunc(rank / pow(base, i))
        assert index < len(alphabet)            
        word += alphabet[index]
        rank = rank % pow(base, i)
    return word

def reverse(indexFilePath):
    data = []
    motifLength=None
    #motifLength = sgRNAIndexReader(indexFilePath)
   
    for interCode, mayOccurence in sgRNAIndexReader(indexFilePath):
            if motifLength is None:
                motifLength = interCode
                print (f"Detected word Length is {motifLength}")
                continue
            try :
                s = DECODER( interCode, motifLength )
            except AssertionError:
                print(f"Can't decode {interCode}. Specified motif length {motifLength} is probably too short")
                return
            data.append( (s, mayOccurence) )
    return data, motifLength

def toggleEncoding():
    global ENCODER
    global DECODER
    global CODEC

    if CODEC == 'twobits':
        ENCODER = pow2encoderWrapper
        DECODER = pow2decoderWrapper
        CODEC = 'pow2'
    else :      
        ENCODER = twobits.encode
        DECODER = twobits.decode
        CODEC = 'twobits'

    print(f"Toggling to {CODEC} encoding")
    return CODEC

if __name__ == "__main__":
    
    ARGUMENTS = docopt(__doc__, version='wordIntegerIndexing 1.0')

    if ARGUMENTS['--dbase']:
        toggleEncoding()

    targetFile =  ARGUMENTS['--out'] if ARGUMENTS['--out'] else None
            
    # Reading sgRNA words, writing integers
    if ARGUMENTS['code']:
        dataEncoded, wLen = indexFromPickle(ARGUMENTS['<pickledDictionary>'], ARGUMENTS['--occ'])

        if not targetFile:
            targetFile = ('.'.join(os.path.basename(ARGUMENTS['<pickledDictionary>']).split('.')[0:-1])
                    + '.index')
        sgRNAIndexWriter(dataEncoded, targetFile, wLen)    

        print(f"Successfully indexed {len(dataEncoded)} words of size {wLen}\nfrom:",
              f"{ARGUMENTS['<pickledDictionary>']}\ninto:{targetFile}"
              )
        exit(0)

    # Reading integers words, writing sgRNA
    if ARGUMENTS['decode']:
        dataWords, wLen = reverse(ARGUMENTS['<indexedFile>'])
        if not targetFile:
            targetFile = ('.'.join(os.path.basename(ARGUMENTS['<indexedFile>']).split('.')[0:-1])
                    + '_decoded.motif')
        sgRNAplainWriter(dataWords, targetFile)

        print("Successfully decoded", len(dataWords), "words\nfrom:",
            ARGUMENTS['<indexedFile>'], "\ninto:", targetFile)
        exit(0)

     # Reading integers words, writing integers code
    if ARGUMENTS['translate']:
        codecFrom=CODEC
        sgRNAdecoded, wLen = reverse(ARGUMENTS['<indexedFile>'])
        toggleEncoding()
        translatedEncoding = translate(sgRNAdecoded)
        if not targetFile:
            targetFile = ('.'.join(os.path.basename(ARGUMENTS['<indexedFile>']).split('.')[0:-1])
                    + '_translated.index')
        sgRNAIndexWriter(translatedEncoding, targetFile, wLen)

        print( (f"Successfully translated{len(translatedEncoding)} words\n"
                f"from:\"{ARGUMENTS['<indexedFile>']}\"({codecFrom})\n"
                f"into:\"{targetFile}\"({CODEC})")
        )
        exit(0)