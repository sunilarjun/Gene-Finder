import sys

# one letter amino acid table converter
aminoAcids = {
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "AAU":"N","AAC": "N",
    "GAU":"D","GAC":"D",
    "UGU":"C","UGC":"C",
    "GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
    "CAA":"Q","CAG":"Q", 
    "CAU":"H","CAC":"H",
    "AUU":"I","AUC":"I","AUA":"I",
    "AAA":"K","AAG":"K",
    "UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "AUG":"Met",
    "UUU":"F","UUC":"F",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "UAA":"Stop","UAG":"Stop","UGA":"Stop",
    "AGU":"S","AGC":"S","UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UGG":"W",
    "UAU":"Y","UAC":"Y",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "GUU":"V","GUG":"V","GUA":"V","GUG":"V"
}

# DNA to Complement
DNAComplement = {
    "T":"A",
    "G":"C",
    "C":"G",
    "A":"T"
}

# DNA to mRNA
mRNA = {
    "T":"A",
    "G":"C",
    "C":"G",
    "A":"U"
}

# mRNA complement mapping
mRNAComplement = {
    "A":"U",
    "U":"A",
    "C":"G",
    "G":"C"
}

def longestORF(frame, direction, data):
    
    localDNA = data.upper()
    startFrame = frame-1
    
    mRNASeq = transmRNA(localDNA)
    tRNASeq = transtRNA(mRNASeq)
    
    #print ''.join(localDNA) testing if input is correct

    if direction is "natural":
        print "mRNA:"
        print "5'->3'"
        print ''.join(mRNASeq)
        print "\n"
        print "Longest ORF:"
        print "5'3' Frame "+str(frame)
    elif direction is "reverse":
        print "Complementary mRNA:"
        print "3'->5'"
        print ''.join(tRNASeq)
        print "\n"
        print "Longest ORF:"
        print "3'5' Frame "+str(frame)
        tRNASeq = reverseStrand(tRNASeq)

    aminoAcidMap(tRNASeq[startFrame:])
    print "\n"

def transmRNA(DNASequence):
    
    print "Input DNA from 5'->3':"
    print ''.join(DNASequence)
    print "\n"

    print "DNA Complement:"
    print "3'->5'"
    print ''.join([DNAComplement[nucleotide] for nucleotide in DNASequence])
    print "\n"

    print "DNA to mRNA:"
    print "5'->3'"
    print ''.join([mRNA[nucleotide] for nucleotide in DNASequence])
    print "\n"

    return ''.join([mRNA[nucleotide] for nucleotide in DNASequence])


def reverseStrand(sequence):
    
    sequence = sequence[::-1]
    return ''.join([mRNAComplement[nucleotide] for nucleotide in sequence])
    

def aminoAcidMap(sequence):
    
    sequenceLen = len(sequence)
    current=0
    protein = []
    while (current < sequenceLen):
        
        codon = sequence[current:current+3]

        if codon in aminoAcids:
            aminoAcid = aminoAcids[codon]
            protein.append(aminoAcid)
        current = current + 3
    
    print ' '.join(protein)
    
def transtRNA(mRNASeq):

    #print ''.join([mRNAComplement[nucleotide] for nucleotide in mRNASeq])

    return ''.join([mRNAComplement[nucleotide] for nucleotide in mRNASeq])
    

def transMac(data):

    #print ''.join([mRNAComplement[nucleotide] for nucleotide in mRNASeq])
    #print ''.join([mRNA[nucleotide] for nucleotide in DNASequence])

    longestORF(3, "reverse", data)

def main():
    
    if(len(sys.argv) <= 1):
        print "Usage: $ python GeneFinder2.py gene.txt"
        sys.exit()
    
    with open(sys.argv[1], 'r') as geneFile:
        data = geneFile.readlines()
        data = [line.replace(' ','').strip() for line in data]
        data = ''.join(data)
    
    transMac(data)

main()