#!/usr/bin/env python3
"""
This script compares a bacterial whole genome assembly to a reference sequence.
It assumes that all replicons are circular, the reference sequence is complete 
(i.e. one contig per replicon) and that the reference sequences are cleanly circularised 
(i.e. no missing or extra bases at the start/end of the sequence).

It produces tab-delimited output with the following columns:
1. Reference filename
2. Reference sequence count
3. Reference sequence total length
4. Assembly filename
5. Assembly sequence count
6. Assembly sequence total length
7. Number of errors in assembly-to-reference alignments
8. The size of the largest indel in assembly-to-reference alignments
9. Number of additional bases in the assembly
10. Number of missing bases in the assembly
11. Fragmentation count: how many times the assembly fragmented a sequence into multiple contigs
12. Missing reference sequences: number of reference sequences with no alignments
13. Junk contig count: number of assembly contigs with no alignments
14. Redundant contig count: number of contigs that had alignments in raw alignments but disappeared after filtering
"""

import argparse
import gzip
import os
import random
import re
import subprocess
import sys
import tempfile


def get_arguments():
    parser = argparse.ArgumentParser(description='Assess a bacterial genome assembly against a '
                                                 'reference genome sequence')

    parser.add_argument('-r', '--reference', type=str,
                        help='reference FASTA file (can be gzipped)')
    parser.add_argument('-a', '--assembly', type=str,
                        help='FASTA file to be assessed (can be gzipped)')

    parser.add_argument('--header', action='store_true',
                        help='Print the header line')
    parser.add_argument('--alignment_preset', type=str, default='asm5',
                        help='Minimap2 alignment preset')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of minimap2 alignment threads')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    if args.header:
        print_header()
    if args.reference is None or args.assembly is None:
        exit()

    ref_seqs = load_fasta(args.reference)
    assembly_seqs = load_fasta(args.assembly)
    
    total_ref_size = sum(len(seq) for _, seq in ref_seqs)
    total_assembly_size = sum(len(seq) for _, seq in assembly_seqs)

    raw_alignments = get_alignments(ref_seqs, args.assembly, args.threads, args.alignment_preset)
    raw_alignments = sorted(raw_alignments, key=lambda a: (-a.ref_align_length, a.ref_start))
    alignments = filter_alignments(raw_alignments, ref_seqs)

    errors = sum(a.errors for a in alignments)
    substitutions = sum(a.total_substitutions() for a in alignments)
    insertions = sum(a.total_insertions() for a in alignments)
    deletions = sum(a.total_deletions() for a in alignments)
    assert errors == substitutions + insertions + deletions
    extra_bases, missing_bases = count_extra_missing_bases(ref_seqs, assembly_seqs, alignments)
    max_indel = max(a.max_indel() for a in alignments)

    fragmentation_count = count_fragmentations(alignments)
    missing_refs = count_missing_references(raw_alignments, ref_seqs)
    junk_contigs = count_junk_contigs(raw_alignments, assembly_seqs)
    redundant_contigs = count_redundant_contigs(raw_alignments, alignments)

    print(f'{args.reference}\t{len(ref_seqs)}\t{total_ref_size}\t'
          f'{args.assembly}\t{len(assembly_seqs)}\t{total_assembly_size}\t'
          f'{errors}\t{substitutions}\t{insertions}\t{deletions}\t'
          f'{max_indel}\t{extra_bases}\t{missing_bases}\t'
          f'{fragmentation_count}\t{missing_refs}\t{junk_contigs}\t{redundant_contigs}')


def print_header():
    print('reference_filename\treference_sequence_count\treference_size\t'
          'assembly_filename\tassembly_sequence_count\tassembly_size\t'
          'error_count\tsubstitutions\tinsertions\tdeletions\t'
          'largest_indel\textra_assembly_bases\tmissing_assembly_bases\t'
          'assembly_fragmentation_count\tmissing_reference_sequences\t'
          'junk_contig_count\tredundant_contig_count')


def get_alignments(ref_seqs, assembly_filename, threads, alignment_preset):
    with tempfile.TemporaryDirectory() as temp_dir:
        multipled_ref_filename = multiply_reference(ref_seqs, temp_dir)
        p = subprocess.run(['minimap2', '-c', '-t', str(threads), '-x', alignment_preset, '--eqx', '-H',
                             multipled_ref_filename, assembly_filename],
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        output = p.stdout.decode()
    return [Alignment(line) for line in output.splitlines()]


def multiply_reference(ref_seqs, temp_dir):
    """
    Write a reference FASTA file with each sequence is represented multiple times, allowing a
    single alignment to span the circular boundary. For larger sequences, tripling the reference
    is enough, but for smaller sequences (where assemblers can make many copies in one contig),
    more copies are necessary.
    """
    multipled_ref_filename = os.path.join(temp_dir, 'multipled_ref.fasta')
    with open(multipled_ref_filename, 'wt') as f:
        for name, seq in ref_seqs:
            if len(seq) > 10000:
                multipled_seq = seq + seq + seq
            else:
                multipled_seq = seq + seq + seq + seq + seq
            f.write(f'>{name}\n{multipled_seq}\n')
    return multipled_ref_filename


def filter_alignments(alignments, ref_seqs):
    """
    Filter alignments to avoid redundant coverage of the circular reference.
    """
    ref_coverage = {name: [False] * len(seq) for name, seq in ref_seqs}
    clean_alignments = []
    for a in alignments:
        ref_len = len(ref_coverage[a.ref_name])
        ranges = get_normalised_ranges(ref_len, a.ref_start, a.ref_end)
        coverage_sum = 0
        for s, e in ranges:
            coverage_sum += sum(ref_coverage[a.ref_name][s:e])
        if coverage_sum > a.ref_align_length / 2:
            continue
        clean_alignments.append(a)
        for s, e in ranges:
            ref_coverage[a.ref_name][s:e] = [True] * (e - s)
    return clean_alignments


def get_normalised_ranges(ref_len, start, end):
    """
    Convert an alignment defined by a start and end (which may span multiple copies)
    into a list of intervals on the single-copy reference.
    """
    total = end - start
    offset = start % ref_len
    intervals = []
    first_length = min(ref_len - offset, total)
    intervals.append((offset, offset + first_length))
    total -= first_length
    while total >= ref_len:
        intervals.append((0, ref_len))
        total -= ref_len
    if total > 0:
        intervals.append((0, total))
    return intervals


def count_extra_missing_bases(ref_seqs, assembly_seqs, alignments):
    """
    Computes:
      - extra bases: portions of the reference with multiple alignments and assembly parts not covered
      - missing bases: portions of the reference not covered at all.
    For the reference, alignment positions are normalised to [0, L).
    """
    ref_coverage = {name: [0] * len(seq) for name, seq in ref_seqs}
    assembly_coverage = {name: [0] * len(seq) for name, seq in assembly_seqs}

    for a in sorted(alignments, key=lambda a: a.ref_align_length, reverse=True):
        ref_len = len(ref_coverage[a.ref_name])
        ranges = get_normalised_ranges(ref_len, a.ref_start, a.ref_end)
        for s, e in ranges:
            for i in range(s, e):
                ref_coverage[a.ref_name][i] += 1
        for i in range(a.query_start, a.query_end):
            assembly_coverage[a.query_name][i] += 1

    ref_covered_multi, ref_not_covered = tally_coverage(ref_seqs, ref_coverage)
    _, assembly_not_covered = tally_coverage(assembly_seqs, assembly_coverage)
    return ref_covered_multi + assembly_not_covered, ref_not_covered


def tally_coverage(seqs, coverage):
    covered_multi, not_covered = 0, 0
    for name, _ in seqs:
        cov = coverage[name]
        covered_multi += sum(count-1 for count in cov if count > 1)
        not_covered += sum(1 for count in cov if count == 0)
    return covered_multi, not_covered


def count_fragmentations(alignments):
    ref_to_contigs = {}
    for a in alignments:
        ref_to_contigs.setdefault(a.ref_name, set()).add(a.query_name)
    return sum(len(contigs) - 1 for contigs in ref_to_contigs.values())


def count_missing_references(alignments, ref_seqs):
    aligned_refs = set(a.ref_name for a in alignments)
    return sum(1 for name, _ in ref_seqs if name not in aligned_refs)


def count_junk_contigs(alignments, assembly_seqs):
    aligned_contigs = set(a.query_name for a in alignments)
    return sum(1 for name, _ in assembly_seqs if name not in aligned_contigs)


def count_redundant_contigs(raw_alignments, filtered_alignments):
    raw_contigs = set(a.query_name for a in raw_alignments)
    filtered_contigs = set(a.query_name for a in filtered_alignments)
    return sum(1 for name in raw_contigs if name not in filtered_contigs)


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict.values())
    with open(filename, 'rb') as unknown_file:
        file_start = unknown_file.read(max_len)
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:
        return open


def load_fasta(filename):
    fasta_seqs = []
    current_name, current_seq_lines = None, []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name is not None:
                    fasta_seqs.append((current_name.split()[0], ''.join(current_seq_lines)))
                current_name, current_seq_lines = line[1:], []
            else:
                current_seq_lines.append(line)
        if current_name is not None:
            fasta_seqs.append((current_name.split()[0], ''.join(current_seq_lines)))
    return fasta_seqs


class Alignment(object):
    def __init__(self, paf_line):
        parts = paf_line.strip().split('\t')
        if len(parts) < 11:
            sys.exit('Error: alignment file does not seem to be in PAF format')

        self.query_name = parts[0]
        self.query_start = int(parts[2])
        self.query_end = int(parts[3])
        self.ref_name = parts[5]
        self.ref_start = int(parts[7])
        self.ref_end = int(parts[8])
        self.errors = int(parts[10]) - int(parts[9])
        self.ref_align_length = self.ref_end - self.ref_start

        cigar_entry = next((p for p in parts if p.startswith('cg:Z:')), None)
        if cigar_entry is None:
            sys.exit('Error: no CIGAR string found')
        self.cigar = cigar_entry[5:]

    def __repr__(self):
        return f'{self.query_name}: {self.query_start}-{self.query_end}; {self.ref_name}: {self.ref_start}-{self.ref_end}; {self.cigar}'

    def max_indel(self):
        return max((int(n) for n, op in re.findall(r'(\d+)([ID])', self.cigar)), default=0)

    def total_substitutions(self):
        return sum(int(n) for n in re.findall(r'(\d+)X', self.cigar))

    def total_insertions(self):
        return sum(int(n) for n in re.findall(r'(\d+)I', self.cigar))

    def total_deletions(self):
        return sum(int(n) for n in re.findall(r'(\d+)D', self.cigar))


if __name__ == '__main__':
    main()


# Unit tests for Pytest
# =====================

def test_get_normalised_ranges():
    assert get_normalised_ranges(10, 2, 8) == [(2, 8)]
    assert get_normalised_ranges(10, 12, 18) == [(2, 8)]
    assert get_normalised_ranges(10, 22, 28) == [(2, 8)]
    assert get_normalised_ranges(10, 0, 10) == [(0, 10)]
    assert get_normalised_ranges(10, 6, 14) == [(6, 10), (0, 4)]
    assert get_normalised_ranges(10, 5, 15) == [(5, 10), (0, 5)]
    assert get_normalised_ranges(10, 2, 18) == [(2, 10), (0, 8)]
    assert get_normalised_ranges(10, 0, 20) == [(0, 10), (0, 10)]
    assert get_normalised_ranges(10, 2, 28) == [(2, 10), (0, 10), (0, 8)]
    assert get_normalised_ranges(10, 2, 38) == [(2, 10), (0, 10), (0, 10), (0, 8)]
    assert get_normalised_ranges(10, 12, 48) == [(2, 10), (0, 10), (0, 10), (0, 8)]
