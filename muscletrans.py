#! /usr/bin/env python
# muscletrans.py

"""
Given a fasta file of dna sequences, performs a translation alignment using muscle with default settings.
Saves the interleaved phylip format alignment to a file with the suffix _aln.phy.
Uses the bacterial translation table by default (see configuration section).
Requires muscle in your $PATH
Usage: muscletrans.py <file with multiple dna sequences in fasta format>
Aaron M Duffy aduffy70{at}gmail.com
July 2010
"""

from sys import argv
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq, SeqRecord
import subprocess
import sys
import fileinput
import random

def main():
    # Configuration
    #Select the desired NCBI translation table
    translationTable = 11

    # Open the DNA sequence file and read the fasta sequences into a dictionary
    if (len(argv) > 1):
        dnaFileName = argv[1]
    else:
        dnaFileName = None
    dnaSeqFile = fileinput.input(dnaFileName)
    dnaSeqDict = SeqIO.to_dict(SeqIO.parse(dnaSeqFile, "fasta"))

    # Translate the sequences
    aaSeqRecords = []
    for key in dnaSeqDict:
        aaSeq = SeqRecord(dnaSeqDict[key].seq.translate(table=translationTable), id=key)
        aaSeqRecords.append(aaSeq)
    dnaSeqFile.close()

    # Replace stop codons with X (unknown aa) so muscle doesn't drop them
    for aaSeq in aaSeqRecords:
        noStopCodonSeq = str(aaSeq.seq).replace('*', 'X')
        aaSeq.seq = Seq(noStopCodonSeq)

    # Align the aa sequences
    commandLine = str(MuscleCommandline(seqtype='protein'))
    childProcess = subprocess.Popen(commandLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=(sys.platform!="win32")) #don't pipe stderr or muscle hangs
    SeqIO.write(aaSeqRecords, childProcess.stdin, "fasta")
    childProcess.stdin.close()
    aaAlignment = AlignIO.read(childProcess.stdout, "fasta")

    # Convert the aa alignment into a dna alignment
    dnaAlignment = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
    for taxon in aaAlignment:
        aaCount = 0
        dnaSeq = ''
        for aaResidue in taxon.seq:
            if (aaResidue == '-'):
                dnaSeq = dnaSeq + '---'
            else:
                dnaSeq = dnaSeq + dnaSeqDict[taxon.id].seq[aaCount*3:aaCount*3+3]
                aaCount+=1
        # As we add the sequences to the alignment remove gene name from the sequence id so they taxon match the PAML constraint tree
        dnaAlignment.add_sequence(taxon.id.split('_')[0], str(dnaSeq))
    if (dnaFileName):
        outFileName = dnaFileName.split('.')[0] + '_aln.phy'
    else:
        outFileName = 'out_aln.phy'
    outFile = open(outFileName, 'w+')
    AlignIO.write([dnaAlignment], outFile, "phylip")

    # Biopython doesn't tag Interleaved phylip files and PAML requires it so...
    outFile.seek(0,0)
    modifiedAlignmentText = outFile.readlines()
    modifiedAlignmentText[0] = modifiedAlignmentText[0].rstrip() + ' I\n'
    outFile.seek(0,0)
    outFile.writelines(modifiedAlignmentText)
    outFile.close()

    #We could remove the alignment columns that contain stop codons but maybe it is best to leave them until we've made manual adjustments to the alignment


if __name__=='__main__':
    main()
