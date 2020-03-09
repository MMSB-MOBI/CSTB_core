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
# 704200 23
233379311 1
3579170655 1
4245510159 1
4652947967 1
````

The header line shows the total number of words encoded and the nucleotide-length of the words.

Toggle to **power of two encoding**, with the `--dbase` flag,

```sh
python wordIntegerIndexing.py code 1a26fb6786e323a436d12439f42b0afa.p \
                            --out test_pow2.index --occ \
                            --dbase
```

Which will output encoded words in a similar format.
In all outputs, words/integer are sorted numerically ascending according to their integer representation.

### Decoding from an encoded sgRNAs index file

Decoding is triggered by the use of the **decode** first positional argument.
Similarly to encoding, the programms assumes that the **two bits per base encoding** was used to encode the data.

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

The header line shows the total number of decoded words.

Use the `--dbase` flag to toggle to the **pow2 decoder**.

```sh
python wordIntegerIndexing.py decode test_pow2.index /
                            --out decoded_from_pow2.motifs
                            --dbase
```

## Using module function



## TO rework BELOW


### By default in twobits

```sh
python wordIntegerIndexing.py <myPickle> --occ -o test_2bits.index
``` 

### Using the pow of 2 encoder

```sh
python wordIntegerIndexing.py <myPickle> --occ -o test_pow2.index --dbase
```

## Decoding from unsigned 64 integers

### By default in twobits

```sh
python wordIntegerIndexing.py reverse 23 test_2bits.index >! motif_2bits.bak
```

### Using the pow of 2 decoder

```sh
python wordIntegerIndexing.py reverse 23 test_pow2.index --dbase >! motif_pow2.bak
```

## Assess operation reciprocity & egality

Following command should be mute

```sh
diff <(sort motif_2bits.bak) <(sort motif_pow2.bak)
```
