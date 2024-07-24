# Availability

This repository contains the source code developed during the master thesis "Implementing a sensitive and long-read optimized DNA alignment tool for tree of life scale sequence research". This thesis is based on the Diamond software by Benjamin Buchfinks, therefore the repository includes all of its source code.

## Compilation for Optimal Performance

To compile the code for optimal performance, use the following commands:

```bash
mkdir build 
cd build
cmake -DWITH_DNA=ON -DEXTRA=ON -DCMAKE_BUILD_MARCH=native ..
```

## Database Creation

To create a database, use the `makedb` command as shown below:

```bash
makedb --in db.fasta -d db.dmnd --dbtype nucl
```

## Basic Run Command

The basic command to run the tool is shown below. For controlling other parameters, refer to the manual of DIAMOND. Parameters are optimized for long reads. To align short reads, leave out the `--align-long-reads` parameter.

```bash
./diamond blastn -d db.dmnd -q query.fasta --align-long-reads  -o output.txt 
```

## Sensitivity Modes

The tool supports the following sensitivity modes:

- `--faster`
- `--fast`
- `--default`
- `--sensitive`
- `--very-sensitive`
- `--ultra-sensitive`

## Specific BLASTN Options

- `--zdrop <value>`: zdrop for gapped dna alignment. Default is 40.
- `--repetition-cutoff <value>`: Filter out top FLOAT fraction of repetitive minimizers. Default is 0.0002.
- `--extension <value>`: Extension algorithm (wfa, ksw=default). Default is "ksw".
- `--chaining-out`: Use chaining without extension.
- `--align-long-reads`: Use chaining with extension.
- `--best-hsp-only`: Only show the best HSP for any single query-subject pair.
- `--chain-pen-gap-scale <value>`: Scaling factor for the chaining gap penalty. Default is 0.8.
- `--chain-pen-skip-scale <value>`: Scaling factor for the chaining skip penalty. Default is 0.0.
- `--penalty <value>`: blastn mismatch penalty. Default is -3.
- `--reward <value>`: blastn match reward. Default is 2.
- `--chain-align-cutoff`: Percentage of chains not mapped/aligned in comparison to the best chaining score.
- `--min-chain-score`: Minimum chaining score to be considered for an alignment.
- `--max-overlap-extension`: Minimum chaining score to be considered for an alignment.
- `--zdrop-extension <value>`: zdrop for extension at chain ends. Default is 80.
- `--zdrop-global <value>`: zdrop for global extension between anchors. Default is 150.
- `--band-extension <value>`: min band for extension at chain ends. Default is 40.
- `--band-global <value>`: min band for global extension between anchors. Default is 40.
