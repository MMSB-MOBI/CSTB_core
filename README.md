# ENCODING DECODING sgRNA motifs

This librairy performs the encoding of sgRNA words  into integer and the decoding of integer into sgRNA words.
The two strategies used to convert word into integer are the **power of two multipliers** or the **two bits per base** encoding.
These general schemes of encoding can be applied to any *"ATCG"* k-mers up to 32 characters. This limit is due to the implementation of the two-bits per base encoder which stores word bits representation as unsigned 64bits integers.

## Using CLI

### Encoding from pickled sgRNAs

By **default**, the **two bits per base** encoder is used to convert word into integer,

```sh
python wordIntegerIndexing.py code 1a26fb6786e323a436d12439f42b0afa.p \
                            --out test_twobits.index --occ
```

The `--occ` flags is to account for the number of occurences of each word, displayed in the second column of the output file.
Here is a sample of the ouput file,

```text
# 704200 23 twobits
233379311 1
3579170655 1
4245510159 1
4652947967 1
````

The header line shows the total number of words encoded, the nucleotide-length of the words and code used.

Toggle to **power of two encoding**, with the `--dbase` flag,

```sh
python wordIntegerIndexing.py code 1a26fb6786e323a436d12439f42b0afa.p \
                            --out test_pow2.index --occ \
                            --dbase
```

Which will output encoded words in a similar format.
In all outputs, words/integer are sorted numerically ascending by their integer representation.

### Decoding from an encoded sgRNAs index file

Decoding is triggered by the use of the **decode** first positional argument.
Codec will be automatically ead from the provided sgRNA index file.

```sh
python wordIntegerIndexing.py decode test_twobits.index /
                            --out decoded_from_twobits.motifs
```

Here is a sample of the output file,

```text
704200
AAAAAAAAAGCGTTCACCCGTGG 1
AAAAAAAGCCCCCCCGAGGCCGG 1
AAAAAAAGGGCAAGCCCTAAAGG 1
AAAAAACACCCCCCTCCTCGGGG 1
AAAAAACCGCTGGAAAGCGTTGG 1
```

The header line shows the total number of decoded words. Following lines shows the sgRNA words and thei occurences.

### Translating sgRNA index from one codec to another

The conversion between the integer coding schemes is triggered by the use of the **translate** first positional argument.
Input codec will be automatically read from the provided sgRNA index file and the outpout codec guessed accordingly.

```sh
python wordIntegerIndexing.py translate test_twobits.index /
                            --out pow2_from_twobits.motifs
```

The output is similar to a `python wordIntegerIndexing.py code` call, with a toggled codec. Here is a sample,

```text
# 704200 23 pow2
179924927 1
165118735 1
717157295 1
660474943 1
```

## Using module functions
Most features of the CLI are available through the import of the ``CSTB_core.engine.wordIntegerIndexing`` module.

The default codec is **two-bits per base**.
The two codec can produce identical integer representations.
```python
import CSTB_core.engine.wordIntegerIndexing as kmer
kmer.encode('GGGGGGGGGGGGGGGGGGGGGGG')
#70368744177663
kmer.encode('GGGGGGGGGGGGGGGGGGGGGGG', codec='pow2')
#70368744177663
kmer.decode(70368744177663, 23)
#'GGGGGGGGGGGGGGGGGGGGGGG'
kmer.decode(70368744177663, 23, codec='pow2')
#'GGGGGGGGGGGGGGGGGGGGGGG'
```

But not always,

```python
import CSTB_core.engine.wordIntegerIndexing as kmer
kmer.encode('AAAAAAAAAGCGTTCACCCGTGG')
#233379311
kmer.encode('AAAAAAAAAGCGTTCACCCGTGG', codec='pow2')
#248916703
kmer.decode(233379311, 23)
#'AAAAAAAAAGCGTTCACCCGTGG'
kmer.decode(248916703, 23, codec='pow2')
#'AAAAAAAAAGCGTTCACCCGTGG'
```