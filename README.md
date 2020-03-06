# ENCODING DECODING kmers

## Encoding from pickled sgRNAs

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
