#! /usr/bin/env python
# assae.py

"""
Aaron's Super Simple Alignment Editor
An interactive console-based alignment editor.
Requires BioPython 1.54
Aaron M Duffy aduffy70{at}gmail.com
July 2010
"""

from sys import argv
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import fileinput
import re

class DisplayedAlignment(object):
    """
    Provides tools for displaying and manipulating an alignment and storing all previous versions
    """

    def __init__(self, alignment):
        self.displayedColumn = 0
        self.alignment = alignment
        self.alignmentHistory = [alignment[:,:]]
        self.changed = False
        self.translated = False
        self.translationTable = 1

    def ParseIndex(self, text):
        """
        Parses a text string specifying a range of rows (taxa) and columns.  Expects the text to be in the format used to specify a range from a Bio.Align.MultipleSeqAlignment.  Returns indices for the start and stop taxon and the start and stop columns
        """
        taxonStart = 0
        taxonStop = len(self.alignment) - 1
        columnStart = 0
        columnStop = self.alignment.get_alignment_length() - 1
        if (',' not in text):
            self.AlertMessage('Invalid index format.  (taxa or columns missing)', 'high')
            return (-1,-1,-1,-1)
        else:
            text = text.strip()
            indices = text.split(',')
            if (len(indices) > 2):
                self.AlertMessage('Invalid index format.  (too many fields)', 'high')
                return (-1,-1,-1,-1)
            else:
                if (':' in indices[0]): #there is a range specified in the taxon index
                    taxonIndices = indices[0].split(':')
                    if (taxonIndices[0]): #a start taxon is specified
                        try:
                            taxonStart = int(taxonIndices[0].strip())
                        except:
                            self.AlertMessage('Invalid index format. (taxon start index not an integer)', 'high')
                            return (-1, -1, -1, -1)
                    if (taxonIndices[1]): #a stop taxon is specified
                        try:
                            taxonStop = int(taxonIndices[1].strip())
                        except:
                            self.AlertMessage('Invalid index format. (taxon stop index not an integer)', 'high')
                            return (-1, -1, -1, -1)
                elif (indices[0]): #a single taxon is specified
                    try:
                        taxonStart = int(indices[0].strip())
                        taxonStop = int(indices[0].strip())
                    except:
                        self.AlertMessage('Invalid index format. (taxon start or stop index not an integer)', 'high')
                        return (-1, -1, -1, -1)
                if (':' in indices[1]): #there is a range specified in the taxon index
                    columnIndices = indices[1].split(':')
                    if (columnIndices[0]): #a start taxon is specified
                        try:
                            columnStart = int(columnIndices[0].strip())
                        except:
                            self.AlertMessage('Invalid index format. (column start index not an integer)', 'high')
                            return (-1, -1, -1, -1)
                    if (columnIndices[1]): #a stop taxon is specified
                        try:
                            columnStop = int(columnIndices[1].strip())
                        except:
                            self.AlertMessage('Invalid index format. (column stop index not an integer)', 'high')
                            return (-1, -1, -1, -1)
                elif (indices[1]): #a single taxon is specified
                    try:
                        columnStart = int(indices[1].strip())
                        columnStop = int(indices[1].strip())
                    except:
                        self.AlertMessage('Invalid index format. (column start or stop index not an integer)', 'high')
                        return (-1, -1, -1, -1)
                if ((0 <= taxonStart <= taxonStop) & (0 <= columnStart <= columnStop)):
                    return (taxonStart, taxonStop, columnStart, columnStop)
                else:
                    self.AlertMessage('Invalid index range. (start > stop or index < 0)', 'high')
                    return (-1,-1,-1,-1)

    def ColorizeDNA(self, text):
        """
        Colorizes output based on nucleotide
        """
        if (text == 'A'):
            escape = '\033[92m' # Green
        elif (text == 'G'):
            escape = '\033[93m' # Yellow
        elif (text == 'T'):
            escape = '\033[91m' # Red
        elif (text == 'C'):
            escape = '\033[96m' # Blue
        else:
            return text
        return escape + text + '\033[0m'

    def ColorizeAA(self, text):
        """
        Colorize output based on amino acid polarity or nonpolarity
        """
        if (text in ['A', 'F', 'H', 'I', 'K', 'L', 'M', 'P', 'R', 'V', 'W']):
            escape = '\033[91m' # Red
        elif (text in ['C', 'G', 'N', 'Q', 'S', 'T', 'Y', 'B', 'Z']):
            escape = '\033[96m' # Blue
        elif (text in ['D', 'E']):
            escape = '\033[92m' # Green
        elif (text in ['X', '*']):
            escape = '\033[93m' # Yellow
        else:
            return text
        return escape + text + '\033[0m'

    def AlertMessage(self, text, severity='low'):
        """
        Display an alert message with a tag and color corresponding to the severity of the alert ('low', 'medium', 'high')
        """
        if (severity == 'high'):
            escape = '\033[91m' # Red
            tag = '!!!'
        elif (severity == 'medium'):
            escape = '\033[93m' # Yellow
            tag = '***'
        else:
            escape = '\033[92m' # Green
            tag = '   '
        print escape + tag, text, tag + '\033[0m'

    def Show(self, column=0):
        """
        Displays 100 columns of the alignment beginning at 'column'
        """
        if column < 0:
            column = 0
        row = 0
        marker = '|    :    ' * 10
        spacer = ' ' * 15
        markerRow = spacer + marker
        if (self.translated == False):
            indexRow = spacer
            for index in range(column, column + 100, 10):
                indexRow = indexRow + str(index).ljust(10)
            print indexRow
            print markerRow
            for sequence in self.alignment[:,column:column + 100]:
                print '%2d) %10s' % (row, sequence.id),
                dnaSequence = ''
                for nucleotide in str(sequence.seq):
                    dnaSequence += self.ColorizeDNA(nucleotide)
                print dnaSequence,
                if (column + 100 < self.alignment.get_alignment_length()):
                    print '...'
                else:
                    print
                row += 1
            print markerRow
            print indexRow
        else:
            indexRow = spacer
            for index in range(column / 3, (column / 3) + 100, 10):
                indexRow = indexRow + str(index).ljust(10)
            print indexRow
            print markerRow
            for sequence in self.alignment[:, column:column + 300]:
                proteinSequence = ''
                for codonPosition in range(0, len(sequence), 3):
                    codon = sequence.seq[codonPosition:codonPosition + 3]
                    if (str(codon) == '---'):
                        proteinSequence += '-'
                    elif ('-' in codon):
                        proteinSequence += '?'
                    else:
                        proteinSequence += self.ColorizeAA(str(codon.translate(table = self.translationTable)))
                print '%2d) %10s %s' % (row, sequence.id, proteinSequence),
                if (column + 300 < self.alignment.get_alignment_length()):
                    print '...'
                else:
                    print
                row += 1
            print markerRow
            print indexRow
        self.displayedColumn = column

    def BackupAlignment(self):
        """
        Stores the current alignment state to the alignment change history
        """
        self.alignmentHistory.append(self.alignment[:,:])

    def UndoChanges(self):
        """
        Reverts to the previous state in the alignment change history.  Does not effect which column index is displayed or whether the sequence is displayed as translated or not since those are not changes to the alignment.
        """
        if (len(self.alignmentHistory) > 1):
            self.alignmentHistory.pop()
            self.alignment = self.alignmentHistory[-1][:,:]
            self.Show(self.displayedColumn)
        else:
            self.AlertMessage('Nothing to undo.', 'low')

    def DeleteRange(self, rangeText, silent=False):
        """
        Removes a row and column range from the alignment
        """
        startTaxon, stopTaxon, startColumn, stopColumn = self.ParseIndex(rangeText)
        if (self.translated == True):
            startColumn = startColumn * 3
            stopColumn = (stopColumn * 3) + 2
        if (startTaxon >= 0): #Make sure we had a valid range
            changeLength = 0
            deleteTaxon = False
            if ((startColumn == 0) & (stopColumn == len(self.alignment[0]) - 1)):
                deleteTaxon = True
            if ((startTaxon > 0) | (stopTaxon < len(self.alignment) - 1)):
                changeLength = (stopColumn - startColumn) + 1
            taxon = 0
            newSequences = []
            for Sequence in self.alignment:
                if (taxon in range(startTaxon, stopTaxon + 1)):
                    if (not deleteTaxon):
                        if (startColumn > 0):
                            Sequence.seq = Sequence.seq[:startColumn] + Sequence.seq[stopColumn + 1:]
                        else:
                            Sequence.seq = Sequence.seq[stopColumn + 1:]
                        if (changeLength):
                            Sequence.seq = Sequence.seq + Seq('-' * changeLength)
                        newSequences.append(Sequence)
                else:
                    newSequences.append(Sequence)
                taxon += 1
            self.alignment = MultipleSeqAlignment(newSequences)
            if (not silent):
                self.Show(self.displayedColumn)
                self.BackupAlignment()

    def ModifyRange(self, rangeText, nucleotide='-'):
        """
        Changes the nucleotides in a row and column range to a specified nucleotide.  Has no effect when the alignment is translated since the corresponding change to the underlying nucleotide alignment would be ambiguous at best.
        """
        nucleotide = nucleotide.upper()
        if (self.translated == True):
            self.AlertMessage("Can't modify protein sequences.", 'medium')
        elif (nucleotide not in ['A', 'G', 'C', 'T', 'R', 'K', 'S', 'W', 'M', 'Y', 'D', 'V', 'B', 'H', 'N', '-']):
            self.AlertMessage('Invalid nucleotide.  (only AGTC- and IUB nucleotide codes are permitted)', 'high')
        else:
            startTaxon, stopTaxon, startColumn, stopColumn = self.ParseIndex(rangeText)
            if (startTaxon >= 0): #Make sure we have a valid range
                taxon = 0
                newSequences = []
                modificationLength = (stopColumn - startColumn) + 1
                for Sequence in self.alignment:
                    if (taxon in range(startTaxon, stopTaxon + 1)):
                        if (startColumn > 0):
                            Sequence.seq = Sequence.seq[:startColumn] + Seq(nucleotide * modificationLength) + Sequence.seq[stopColumn + 1:]
                        else:
                            Sequence.seq = Seq(nucleotide * modificationLength) + Sequence.seq[stopColumn + 1:]
                    newSequences.append(Sequence)
                    taxon += 1
                self.alignment = MultipleSeqAlignment(newSequences)
                self.Show(self.displayedColumn)
                self.BackupAlignment()

    def InsertRange(self, rangeText):
        """
        Inserts a row and column range into the alignment and fills it with gaps ('-' for nucleotides or '---' for translated alignments)
        """
        startTaxon, stopTaxon, startColumn, stopColumn = self.ParseIndex(rangeText)
        if (self.translated == True):
            startColumn = startColumn * 3
            stopColumn = (stopColumn * 3) + 2
        if (startTaxon >= 0): #Make sure we had a valid range
            changeLength = (stopColumn - startColumn) + 1
            taxon = 0
            newSequences = []
            for Sequence in self.alignment:
                if (taxon in range(startTaxon, stopTaxon + 1)):
                    if (startColumn > 0):
                        Sequence.seq = Sequence.seq[:startColumn] + Seq('-' * changeLength) + Sequence.seq[startColumn:]
                    else:
                        Sequence.seq = Seq('-' * changeLength) + Sequence.seq[:]
                else:
                    Sequence.seq = Sequence.seq + Seq('-' * changeLength)
                newSequences.append(Sequence)
                taxon +=1
            self.alignment = MultipleSeqAlignment(newSequences)
            self.Show(self.displayedColumn)
            self.BackupAlignment()

    def Jump(self, column):
        """
        Moves the displayed column to a specified column index
        """
        if (self.translated == True):
            column = column * 3
        self.Show(column)

    def ScrollRight(self, offset=100):
        """
        Scroll the display 'offset' columns to the right
        """
        if (self.translated == True):
            offset = offset * 3
        self.Show(self.displayedColumn + offset)

    def ScrollLeft(self, offset=100):
        """
        Scroll the display 'offset' columns to the left
        """
        if (self.translated == True):
            offset = offset * 3
        self.Show(self.displayedColumn - offset)

    def Reverse(self):
        """
        Reverses the order of the columns in the alignment.  Has no effect on translated sequences.
        """
        if (self.translated == False):
            self.alignment = self.alignment[:,::-1]
            self.Show(self.displayedColumn)
            self.BackupAlignment()
        else:
            self.AlertMessage("Can't reverse protein sequences.", 'medium')

    def Complement(self):
        """
        Give the complement of the alignment.  Has no effect on translated sequences.
        """
        if (self.translated == False):
            for i in range(len(self.alignment)):
                self.alignment[i].seq = self.alignment[i].seq.complement()
            self.Show(self.displayedColumn)
            self.BackupAlignment
        else:
            self.AlertMessage("Can't complement protein sequences.", 'medium')

    def ReverseComplement(self):
        """
        Reverse and complement the alignment.  Has no effect on translated sequences.
        """
        if (self.translated == False):
            for i in range(len(self.alignment)):
                self.alignment[i].seq = self.alignment[i].seq.reverse_complement()
            self.Show(self.displayedColumn)
            self.BackupAlignment()
        else:
            self.AlertMessage("Can't reverse-complement protein sequences.", 'medium')

    def Translate(self, translationTable=11):
        """
        Switch to displaying and manipulating the sequence as a protein sequence.
        Still works on translated sequences if a different translation table is specified, otherwise it backtranslated translated sequences.
        """
        if ((self.translated == False) | ((self.translated == True) & (self.translationTable != translationTable))):
            self.translated = True
            self.translationTable = translationTable
            self.displayedColumn = self.displayedColumn - (self.displayedColumn % 3)
            self.Show(self.displayedColumn)
        else:
            self.BackTranslate()

    def BackTranslate(self):
        """
        Revert to displaying and manipulating the sequence as a dna sequence.  Has no effect if the sequence is already dna.
        """
        if (self.translated == True):
            self.translated = False
            self.Show(self.displayedColumn)
        else:
            self.AlertMessage("Can't back-translate.  Alignment contains DNA sequences", 'medium')

    def Save(self, fileName='alignment.phy', alignmentFormat='phylip'):
        """
        Write alignment to disk
        """
        AlignIO.write(self.alignment, fileName, alignmentFormat)
        self.AlertMessage('Saved alignment to ' + fileName + ' in ' + alignmentFormat + ' format.', 'low')

    def CleanUp(self):
        """
        Condense the alignment by removing any columns that contain spaces in all taxa.
        """
        blankColumnPattern = re.compile('^-*$')
        blankColumns = []
        for columnIndex in range(self.alignment.get_alignment_length() - 1):
            columnValues = self.alignment[:,columnIndex]
            match = blankColumnPattern.search(columnValues)
            if (match):
                blankColumns.append(str(columnIndex))
        for column in blankColumns[::-1]:
            self.DeleteRange(',' + str(column), True)
        self.Show(self.displayedColumn)
        self.BackupAlignment()


if __name__ == '__main__':
    if (len(argv) > 1):
        alignmentFileName = argv[1]
    else:
        alignmentFileName = None
    alignmentFile = fileinput.input(alignmentFileName)
    alignmentFromFile = AlignIO.read(alignmentFile, "phylip")
    alignment = DisplayedAlignment(alignmentFromFile)
    alignment.Show(0)
    choice = ''
    while (choice != 'X'):
        print '(#)Jump', '(<>)Scroll', '(t)Translate', '(r)Reverse', '(c)Complement', '(rc)Rev-comp', '(d)Delete', '(i)Insert', '(m)Modify', '(cl)Cleanup'
        print '(z)Undo', '(s)Save', '(x)Exit'
        choice = raw_input('Command: ').upper()
        if (choice == 'Z'):
            alignment.UndoChanges()
        elif (choice == 'R'):
            alignment.Reverse()
        elif (choice == 'C'):
            alignment.Complement()
        elif (choice == 'RC'):
            alignment.ReverseComplement()
        elif (choice == '>'):
            alignment.ScrollRight(100)
        elif (choice == '<'):
            alignment.ScrollLeft(100)
        elif (choice == 'T'):
            alignment.Translate()
        elif (choice == ''):
            alignment.ScrollRight(1)
        elif (choice == 'S'):
            alignment.Save()
        elif (choice == 'CL'):
            alignment.CleanUp()
        else:
            displayPattern = re.compile('^(\d+)')
            match = displayPattern.search(choice)
            if (match):
                alignment.Jump(int(match.group(1)))
            else:
                translatePattern = re.compile('^T(\d+)')
                match = translatePattern.search(choice)
                if (match):
                    alignment.Translate(int(match.group(1)))
                else:
                    scrollLeftPattern = re.compile('^<(\d+)')
                    match = scrollLeftPattern.search(choice)
                    if (match):
                        alignment.ScrollLeft(int(match.group(1)))
                    else:
                        scrollRightPattern = re.compile('^>(\d+)')
                        match = scrollRightPattern.search(choice)
                        if (match):
                            alignment.ScrollRight(int(match.group(1)))
                        else:
                            savePattern = re.compile('^S(\S+),(\S+)')
                            match = savePattern.search(choice)
                            if (match):
                                alignment.Save(match.group(1), match.group(2))
                            else:
                                deleteRangePattern = re.compile('^D(.*)')
                                match = deleteRangePattern.search(choice)
                                if (match):
                                    alignment.DeleteRange(match.group(1))
                                else:
                                    modifyRangePattern = re.compile('^M(.*)(\S)')
                                    match = modifyRangePattern.search(choice)
                                    if (match):
                                        alignment.ModifyRange(match.group(1), match.group(2))
                                    else:
                                        insertRangePattern = re.compile('^I(.*)')
                                        match = insertRangePattern.search(choice)
                                        if (match):
                                            alignment.InsertRange(match.group(1))

