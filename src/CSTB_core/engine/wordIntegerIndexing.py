"""
    Translate a dictionary of constant-length CRISPR motifs into an ordered set of integer, using base4 encoding.
    Usage:
        wordIntegerIndexing.py <pickledDictionary> [--dbase] [--out=<outFile> --occ]
        wordIntegerIndexing.py reverse <motifLength> <indexedFile> [--dbase]

    Options:
        -h --help                               Show this screen
        -o <outFile>, --out <outFile>           Name of the output file [default : ./pickledDictionary.index]
        --dbase,                                Switch to power of two encoder, default is twobits
"""

import os
import pickle
import math
import twobits
from docopt import docopt

ENCODER=twobits.encode
DECODER=twobits.decode

def occWeight(data):
    k,datum = data
    n = 0
    for o in datum:
        for _o in datum[o]:
            n += len(datum[o][_o])
    #print(k,n)
    return n

# same as index pickle, coding and order wise, 
# we just add a second field to each wordCode line, the occurence number
def indexAndOccurencePickle(file_path, target_file):
    global ENCODER
    p_data = pickle.load(open(file_path, "rb"))
    word_list = list(p_data.keys())
    print(f"Encoding what seems like {len(word_list[0])} length words")
    data = sorted( [ ( ENCODER(w), occWeight((w, p_data[w])) ) for w in word_list], key=lambda x: x[0])
    with open(target_file, "w") as filout:
        filout.write(str(len(data)) + "\n")
        for datum in data:
            filout.write( ' '.join([str(d) for d in datum]) + "\n")

    return len(data)

# same as index pickle, coding and order wise, 
# we just add a second field to each wordCode line, the occurence number
def indexAndOccurence(data):
    word_list = list(data.keys())
    data = sorted( [ ( ENCODER(w), occWeight((w, data[w])) ) for w in word_list], key=lambda x: x[0])
    return data
    
def indexPickle(file_path, target_file):
    """
    Take a pickle file, code it and write it in a file
    """
    global ENCODER
    p_data = pickle.load(open(file_path, "rb"))
    word_list = list(p_data.keys())
    data = sorted([ ENCODER(w) for w in word_list])
    with open(target_file, "w") as filout:
        filout.write(str(len(data)) + "\n")
        for coding_int in data:
            filout.write(str(coding_int) + "\n")

    return len(data)

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
    return decode(code, "ATCG", length=wLen)

def decode(rank, alphabet, length=20):
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

def writeIndexes(indexData, output):
    with open(output, "w") as o:
        o.write(f"{len(indexData)}\n")
        for datum in indexData:
            o.write(f"{' '.join([str(d) for d in datum])}\n")

def reverse(indexFilePath, motifLength, skipFirst=True):
    with open(indexFilePath, 'r') as fp:    
        for l in fp:
            if skipFirst:
                skipFirst = False
                continue
            _ = l.split()           
            try :
                s = DECODER( int(_[0]), motifLength )
            except AssertionError:
                print(f"Can't decode {_[0]}. Specified motif length {motifLength} is probably too short")
                return
            print(s)

def toggleEncoding():
    global ENCODER
    global DECODER
    ENCODER = pow2encoderWrapper
    DECODER = pow2decoderWrapper
    
if __name__ == "__main__":
    
    ARGUMENTS = docopt(__doc__, version='wordIntegerIndexing 1.0')

    if ARGUMENTS['--dbase']:
        print("Toggling to pow2 encoding")
        toggleEncoding()

    if ARGUMENTS['reverse']:
        reverse(ARGUMENTS['<indexedFile>'], int(ARGUMENTS["<motifLength>"]))
        exit(1)

    TARGET_FILE = ('.'.join(os.path.basename(ARGUMENTS['<pickledDictionary>']).split('.')[0:-1])
                   + '.index')
    if ARGUMENTS['--out']:
        TARGET_FILE = ARGUMENTS['--out']
    
    indexFn = indexPickle
    if ARGUMENTS['--occ']:
        indexFn = indexAndOccurencePickle
    TOTAL = indexFn(ARGUMENTS['<pickledDictionary>'], TARGET_FILE)
    print("Successfully indexed", TOTAL, "words\nfrom:",
          ARGUMENTS['<pickledDictionary>'], "\ninto:", TARGET_FILE)
