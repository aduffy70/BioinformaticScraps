#! /usr/bin/env python
# pamlize.py

"""
Adds an 'I' for 'Interleaved' to the end of the first line of a phylip formatted alignment.  The package paml requires that interleaved phylip alignnments be tagged with an 'I' but biopython doesn't place one and will not read an alignment that has one.  I don't know if Biopython is right or if paml is right, but either way this script solves the problem - though perhaps this script should be called unbiopythonize
Usage: pamlize.py <phylip format alignment file>
Aaron M Duffy aduffy70{at}gmail.com
August 2010
"""

from sys import argv

def main():
    # Biopython doesn't tag Interleaved phylip files and PAML requires it so...
    alignmentFile = open(argv[1], 'r+')
    alignmentFile.seek(0,0)
    modifiedAlignmentText = alignmentFile.readlines()
    modifiedAlignmentText[0] = modifiedAlignmentText[0].rstrip() + ' I\n'
    alignmentFile.seek(0,0)
    alignmentFile.writelines(modifiedAlignmentText)
    alignmentFile.close()

    #We could remove the alignment columns that contain stop codons but maybe it is best to leave them until we've made manual adjustments to the alignment


if __name__=='__main__':
    main()
