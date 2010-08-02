#! /usr/bin/env python
# alignedit.py

"""
An interactive alignment editor.
Requires BioPython 1.54
Aaron M Duffy aduffy70{at}gmail.com
July 2010
"""

from sys import argv
from Bio import AlignIO
import fileinput
import re
import curses

class DisplayedAlignment(object):
    """
    Provides tools for displaying and manipulating an alignment and storing all previous versions
    """

    def __init__(self, alignment):
        self.displayedColumn = 0
        self.alignment = alignment
        self.alignmentHistory = [alignment]
        self.changed = False
        self.translated = False
        self.translationTable = 1

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
        Colorize output based on Amini Acid polarity
        """
        if (text in ['A', 'F', 'H', 'I', 'K', 'L', 'M', 'P', 'R', 'V', 'W']):
            escape = '\033[91m'
        elif (text in ['C', 'G', 'N', 'Q', 'S', 'T', 'Y', 'B', 'Z']):
            escape = '\033[96m'
        elif (text in ['D', 'E']):
            escape = '\033[92m'
        elif (text in ['X', '*']):
            escape = '\033[93m'
        else:
            return text
        return escape + text + '\033[0m'

    def Show(self, column=0):
        """
        Displays 100 columns of the alignment beginning at 'column'
        """
        if column < 0:
            column = 0
        if (self.translated == False):
            marker = '|         ' * 10
            print '           ' + str(column) + '\n           ' + marker
            for sequence in self.alignment[:,column:column + 100]:
                #print '%10s %s' % (sequence.id, sequence.seq),
                print '%10s' % sequence.id,
                dnaSequence = ''
                for nucleotide in str(sequence.seq):
                    dnaSequence += self.ColorizeDNA(nucleotide)
                print dnaSequence,
                if (column + 100 < self.alignment.get_alignment_length()):
                    print '...'
                else:
                    print
        else:
            marker = '|         ' * 10
            print '           ' + str(column) + '\n           ' + marker
            column = column * 3
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
                print '%10s %s' % (sequence.id, proteinSequence),
                if (column + 300 < self.alignment.get_alignment_length()):
                    print '...'
                else:
                    print
        self.displayedColumn = column

    def MarkVersion(self):
        """
        Stores the current alignment to the alignment history
        """
        if (self.changed):
            self.alignmentHistory.append(self.alignment)
            self.changed = False

    def UndoChanges(self):
        """
        Reverts the current alignment to the most recently saved alignment in the alignment history.
        """
        if (self.changed):
            self.alignment = self.alignmentHistory[-1]
            self.changed = False
        if (self.translated == True):
            self.Show(self.displayedColumn / 3)
        else:
            self.Show(self.displayedColumn)

    def PreviousVersion(self):
        """
        Restores a previously saved alignment from the alignment history
        """
        if (self.changed):
            self.alignment = self.alignmentHistory[-1]
        else:
            if (len(self.alignmentHistory) > 1):
                self.alignmentHistory.pop()
                self.alignment = self.alignmentHistory[-1]
        self.changed = False
        if (self.translated == True):
            self.Show(self.displayedColumn / 3)
        else:
            self.Show(self.displayedColumn)

    def DeleteColumns(self, start, extend=0):
        """
        Removes columns start to start + extend (inclusive) from the alignment
        """
        if (self.translated == True):
            start = start * 3
            extend = (extend * 3) + 2
        if (start > 0):
            self.alignment = self.alignment[:,:start] + self.alignment[:, start + extend + 1:]
            self.Show(self.displayedColumn / 3)
        else:
            self.alignment = self.alignment[:, start + extend + 1:]
            self.Show(self.displayedColumn)
        self.changed = True

    def ScrollRight(self, offset=100):
        """
        Scroll the display 'offset' columns to the right
        """
        if (self.translated == True):
            self.Show((self.displayedColumn / 3) + offset)
        else:
            self.Show(self.displayedColumn + offset)

    def ScrollLeft(self, offset=100):
        """
        Scroll the display 'offset' columns to the left
        """
        if (self.translated == True):
            self.Show((self.displayedColumn / 3) - offset)
        else:
            self.Show(self.displayedColumn - offset)

    def Reverse(self):
        """
        Reverses the order of the columns in the alignment.
        Has no effect on translated sequences.
        """
        if (self.translated == False):
            self.alignment = self.alignment[:,::-1]
            self.changed = True
            self.Show(self.displayedColumn)
        else:
            print "Can't reverse protein sequences."

    def Complement(self):
        """
        Give the complement of the alignment.
        Has no effect on translated sequences.
        """
        if (self.translated == False):
            for i in range(len(self.alignment)):
                self.alignment[i].seq = self.alignment[i].seq.complement()
            self.changed = True
            self.Show(self.displayedColumn)
        else:
            print "Can't complement protein sequences."

    def ReverseComplement(self):
        """
        Reverse and complement the alignment.
        Has no effect on translated sequences.
        """
        if (self.translated == False):
            for i in range(len(self.alignment)):
                self.alignment[i].seq = self.alignment[i].seq.reverse_complement()
            self.changed = True
            self.Show(self.displayedColumn)
        else:
            print "Can't reverse-complement protein sequences."

    def Translate(self, translationTable=11):
        """
        Switch to displaying and manipulating the sequence as a protein sequence.
        Still works on translated sequences (so a different translation table can be called).
        """
        self.translated = True
        self.translationTable = translationTable
        self.displayedColumn = self.displayedColumn - (self.displayedColumn % 3)
        self.Show(self.displayedColumn / 3)

    def BackTranslate(self):
        """
        Revert to displaying and manipulating the sequence as a dna sequence.
        Has no effect if the sequence is already dna.
        """
        if (self.translated == True):
            self.translated = False
            self.Show(self.displayedColumn)
        else:
            print "Can't back-translate.  Alignment contains DNA sequences"

    def Save(self, fileName='alignment.phy', alignmentFormat='phylip'):
        """
        Write alignment to disk
        """
        AlignIO.write(self.alignment, fileName, alignmentFormat)
        print 'Saved alignment to ' + fileName + ' in ' + alignmentFormat + ' format.'

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
    while (choice == ''):
        print '(##)Jump'.ljust(11), '(<>)Scroll'.ljust(11), '(t)Translate'.ljust(11), '(b)Backtrans'.ljust(11), '(r)Reverse'.ljust(11), '(c)Complement'.ljust(11), '(rc)Rev-comp'.ljust(11), '(dc)Delete col)'.ljust(11)
        print '(mv)Mark version'.ljust(11),'(z)Undo'.ljust(11), '(pv)Prev version'.ljust(11), '(s)Save'.ljust(11), '(x)Exit'.ljust(11)
        choice = raw_input(': ').upper()
        if (choice == 'MV'):
            alignment.MarkVersion()
        elif (choice == 'Z'):
            alignment.UndoChanges()
        elif (choice == 'PV'):
            alignment.PreviousVersion()
        elif (choice == 'R'):
            alignment.Reverse()
        elif (choice == 'C'):
            alignment.Complement()
        elif (choice == 'RC'):
            alignment.ReverseComplement()
        elif (choice == 'B'):
            alignment.BackTranslate()
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
        else:
            displayPattern = re.compile('^(\d+)')
            match = displayPattern.search(choice)
            if (match):
                alignment.Show(int(match.group(1)))
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
                            deleteColumnPattern = re.compile('^DC(\d+)-(\d+)')
                            match = deleteColumnPattern.search(choice)
                            if (match):
                                alignment.DeleteColumns(int(match.group(1)), int(match.group(2)) - int(match.group(1)))
                            else:
                                deleteColumnPattern = re.compile('^DC(\d+)')
                                match = deleteColumnPattern.search(choice)
                                if (match):
                                    alignment.DeleteColumns(int(match.group(1)))
                                else:
                                    savePattern = re.compile('^S(\S+),(\S+)')
                                    match = savePattern.search(choice)
                                    if (match):
                                        alignment.Save(match.group(1), match.group(2))
        if (choice != 'X'):
            choice = ''





