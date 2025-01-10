from __future__ import annotations
import argparse
from collections import abc
from concurrent.futures import ProcessPoolExecutor, as_completed
import io
import os
import pathlib
import subprocess
import sys
import tempfile
from typing import Literal
import pandas as pd
from Bio.Seq import Seq
import json

def makeblastdb(input_file: pathlib.Path, type: Literal["nucl"] | Literal["prot"]):
    """Run `makeblastdb` to create a blast database for the given FASTA file. Specify `type="nucl"` for nucleotide FASTA files and `type="prot"` for amino acids."""
    cmd = ["makeblastdb", "-in", str(input_file), "-dbtype", type]
    process = subprocess.run(cmd, capture_output=True)
    if len(process.stderr) > 0:
        print(process.stderr, file=sys.stderr)
        raise RuntimeError("makeblastdb failed.")

def _blastnp(executable: Literal["blastn"] | Literal["blastp"], query_filename: pathlib.Path, db_name: pathlib.Path):
    query_filename = query_filename.absolute()

    # We use -outfmt 13, which creates a json file per query sequence entry.
    # These files are being written to this temporary directory by blastn/blastp, read by us and then the directory is deleted again.
    with tempfile.TemporaryDirectory() as directory:
        cmd = [
            executable,
            "-task", executable,                # Either "blastn" or "blastp"
            "-query", str(query_filename),
            "-db", db_name.name,
            "-outfmt", "13",
            "-out", directory + "/out",
            "-num_threads", "4",  # Number of available CPUs
        ]
        process = subprocess.run(cmd, capture_output=True, cwd=db_name.parent)
        if len(process.stderr) > 0:
            raise ValueError(process.stderr[:100])
        
        # Read all blast output files
        data = {}
        for filename in pathlib.Path(directory).glob("out*.json"):
            with open(filename, "r") as f:
                # All relevant results are under BlastOutput2 - report - results - search
                entry = json.load(f)['BlastOutput2']["report"]["results"]["search"]
            # Gather all hits, and all hsps for each hit, saving the target title.
            hits = [(hit['description'][0]['title'].split('|')[0].strip(), hsp) for hit in entry["hits"] for hsp in hit["hsps"]]
            if hits:
                # If we have any hits, save the one with the lowest evalue
                name, hit = min(hits, key=lambda e: e[1]["evalue"])
                hit["b"] = name
                data[entry["query_title"]] = hit

    # Return a DataFrame with all the data for the lowest-evalue hit for each query sequence.
    return pd.DataFrame.from_dict(data, orient="index").rename_axis(index="a")

def blastp(query_filename: pathlib.Path, db_name: pathlib.Path):
    """Run blastp with the given query FASTA file on the given DB."""
    return _blastnp("blastp", query_filename, db_name)

def blastn(query_filename: pathlib.Path, db_name: pathlib.Path):
    """Run blastn with the given query FASTA file on the given DB."""
    return _blastnp("blastn", query_filename, db_name)

class FASTAFile(abc.Mapping):
    """Class to read and subset sequences in a FASTA file."""
    def __init__(self, filename: pathlib.Path, type="nucl"):
        """Instantiate a FASTAFile with the given file name and type (either `nucl` or `prot`)."""
        self.filename = pathlib.Path(filename)
        self.type = type

        self._indexed = False
        self._sequences = None

    def _index(self):
        # Make a blast database if not already done.
        if self._indexed:
            return
        makeblastdb(self.filename, self.type)
        self._indexed = True

    def _read(self):
        # Lazy read the sequences in the file.
        if self._sequences is not None:
            return
        sequences = {}
        with open(self.filename, "r") as f:
            name = None
            s = ""
            for line in f:
                if line.startswith('>'):
                    if name is not None:
                        sequences[name] = s
                    s = ""
                    name = line[1:].split("|")[0].strip()
                else:
                    s += line.strip()
            sequences[name] = s
        self._sequences = sequences

    def __len__(self):
        """Number of sequences in the file."""
        self._read()
        return len(self._sequences)
    
    def __getitem__(self, key: abc.Sequence[str]|str):
        """Slice the FASTA file to get a subset of sequences as dict."""
        self._read()
        if (not isinstance(key, str)) and (isinstance(key, abc.Sequence) or isinstance(key, pd.Series)):
            # Return a dict with the sequences if key is a list
            return {name: self._sequences[name] for name in key}
        # Return a single sequence if key is a str
        return self._sequences[key]
    
    def __iter__(self):
        """Iterate over all sequence names."""
        self._read()
        return iter(self._sequences)

    def query(self, other: FASTAFile):
        """Run blastn or blastp with this FASTA file as the reference and the given FASTA file as the query."""
        self._index()
        if self.type == "nucl":
            output = blastn(other.filename, self.filename)
        elif self.type == "prot":
            output = blastp(other.filename, self.filename)
        else:
            raise ValueError(f"Unknown alphabet type {self.type}")
        
        return output
    
    def translate(self):
        """Translate the nucleotide sequences in this file to amino acid sequences and return them as a new temporary FASTA file object."""
        if self.type != "nucl":
            raise ValueError("Cannot translate FASTA files with a type other than 'nucl'.")
        protein_sequences = {name: Seq(_pad_Ns(sequence)).translate() for name, sequence in self.items()}
        return TemporaryFASTAFile(protein_sequences, type="prot")

    @classmethod
    def from_dict(cls, sequences: dict[str, str], filename: None | pathlib.Path=None, type="nucl"):
        new = cls(filename, type=type)
        new._sequences = dict(sequences)
        if filename is not None:
            new.save(filename)
        return new

    def save(self, file: io.TextIOWrapper | pathlib.Path):
        if isinstance(file, pathlib.Path) or isinstance(file, str):
            file = open(file, "w")
        with file:
            for name, sequence in self.items():
                file.write(f">{name}\n{sequence}\n")

def _pad_Ns(s: str):
    # Pad the given sequence with "N"s to make it a multiple of three.
    return s + (((3 - len(s) % 3) % 3)*"N")

class TemporaryFASTAFile(FASTAFile):
    """Class describing a temporary FASTA file."""
    def __init__(self, sequences: dict[str, str], type="nucl"):
        """Instantiate a temporary FASTA file, saving the given sequences to a temporary known location. This is deleted during object destruction."""
        tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
        filename = pathlib.Path(tmp.name)
        FASTAFile.__init__(self, filename, type=type)
        self._sequences = dict(sequences)
        self.save(tmp)
        
    def __del__(self):
        # Remove the temporary file.
        os.remove(str(self.filename))

def reciprocal_best_blast(a: FASTAFile, b: FASTAFile):
    """Identify the reciprocal best blast hits between the two given FASTA files."""
    # Get the best blast hits from a in b (there...).
    df_there = b.query(a)
    # Get the sequences of the hits in b as a temporary FASTA file.
    reciprocal_query = TemporaryFASTAFile(b[df_there.loc[:, "b"]])
    # Get the reverse best blast hit from the identified sequences in a (...and back again).
    df_back = a.query(reciprocal_query)
    # Join the two DataFrames on the "b" column (since the query names are from "b" in the back step)
    df: pd.DataFrame = df_there.join(df_back, on="b", how="inner", lsuffix="there", rsuffix="back")
    # Only retain entries in which the back step identified the same entry as we started with in the there step.
    return df[df.index == df.bback].reset_index()

def _insertion(alignedpseq: str, trimmednuclseq: str) -> str:
    # Take the given aligned amino acid sequence and insert --- for each inserted amino acid.
    finalnuclseq: str = ""
    # Pad in case the nucleotide sequence length is not divisible by 3.
    trimmednuclseq = _pad_Ns(trimmednuclseq)
    nuclseqpos = 0
    for ppos, aa in enumerate(alignedpseq):
        if aa == "-":
            finalnuclseq += "---"
        else:
            finalnuclseq += trimmednuclseq[nuclseqpos:nuclseqpos+3]
            nuclseqpos += 3
    return finalnuclseq

def _nr_different_kmers(a: str, b: str, k: int=3):
    # Count the number of different kmers between strings a and b
    if len(a) != len(b):
        raise ValueError("Strings need to be of same length")
    differences = 0
    for i in range(0, len(a), k):
        if a[i:i+k] != b[i:i+k]:
            differences += 1
    return differences

def _calculate_Ka(r: pd.Series) -> int:
    # Count the number of "+" and " " in the blast midline from the back step, this gives us the number of non-synonymous mutations.
    midline: str = r.midlineback
    return midline.replace("+", " ").count(" ")

def _calculate_KsKa(r: pd.Series):
    # Correct the nucleotide sequences with any insertions from the alignment of the amino acid sequences
    # (i.e. insert "---" into the nucl sequences every time we encounter a "-" in the prot sequences).
    a = _insertion(r.hseqback, r.anuclseq[(r.hit_fromback-1)*3:(r.hit_toback-1)*3+3])
    b = _insertion(r.qseqback, r.bnuclseq[(r.query_fromback-1)*3:(r.query_toback-1)*3+3])
    # Count the number of different kmers (this is Ka+Ks)
    return _nr_different_kmers(a, b)

def _row_KaKs(r: pd.Series):
    try:
        # Get the number of non-synonymous mutations Ka
        Ka = _calculate_Ka(r)
        # Get the number of synonymous mutations by identifying all kmer changes and then subtracting Ka
        Ks = _calculate_KsKa(r) - Ka
        return pd.Series({
            "Ka": Ka,
            "Ks": Ks,
            "KaKs": Ka/Ks if Ks > 0 else float("NaN"),
        })
    except Exception as e:
        # Fail gracefully if there is a problem.
        print(
            f"Could not calculate Ka/Ks for {r.name} and {r.b}:\n" + str(e),
            file=sys.stderr,
        )
        return pd.Series({
            "Ka": float("NaN"),
            "Ks": float("NaN"),
            "KaKs": float("NaN"),
        })

def kaks_single(an: FASTAFile, bn: FASTAFile) -> pd.Series:
    """Calculate Ka/Ks for the given pair of FASTA files."""
    # Translate nucleotide sequences to amino acid sequences.
    ap = an.translate()
    bp = bn.translate()
    # Identify the reciprocal best blast hits
    rbb = reciprocal_best_blast(ap, bp)
    # Insert the nucleotide sequences into the DataFrame
    rbb.loc[:, "anuclseq"] = [an[name] for name in rbb.a]
    rbb.loc[:, "bnuclseq"] = [bn[name] for name in rbb.b]
    # Set the index as the query sequence names 
    rbb = rbb.set_index("a")
    # Calculate Ka/Ks for each entry
    return rbb.apply(_row_KaKs, axis=1).KaKs.rename(bn.filename.stem)

def kaks_multiple(query: FASTAFile, references: list[FASTAFile]) -> pd.DataFrame:
    """Calculate Ka/Ks of the query sequences against each of the given reference sequence files."""
    # Display a progress bar if tqdm is available
    try:
        from tqdm.auto import tqdm
    except ImportError:
        tqdm = lambda x, *args, **kwargs: x

    # Perform the single comparisons and concatenate the DataFrames.
    results = []
    # We allow blast to use 4 threads, so use num_cpu / 4 worker processes.
    with ProcessPoolExecutor(len(os.sched_getaffinity(0)) // 4) as executor:
        futures = [
            executor.submit(kaks_single, query, reference)
            for reference in references
        ]
        for future in tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())
    
    return pd.concat(results, axis=1)

def cli_main():
    parser = argparse.ArgumentParser(
        description="Perform a Ka/Ks analysis for the given nucleotide query and reference FASTA files. The result is given in CSV format to stdout."
    )
    parser.add_argument("query", type=FASTAFile, help="The filename of the query FASTA.")
    parser.add_argument("reference", nargs="+", type=FASTAFile, help="The filename(s) of the reference FASTA.")

    args = parser.parse_args()
    kaks_multiple(args.query, args.reference).to_csv("/dev/stdout")
