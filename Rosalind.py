import itertools
import os
string="CTAAGCCCCATTGTGGTCGACAGGATTGGACCGTGTAGCCCAAATATATCACTCATTTTTTTCCAATCTGCCGATGAAGGTACACATCTTGGCCTATGCACCTACCCGAGGCGCACGGCGACCTTTAGCACCATTCTACGTCACCGGAGACTCTCGAGAATGATGGAGTGTCAGAGGTAGACATCTGTAGAGCGATTTCGTGCGCCGACTCACTCAAGGGCATCTATTGACCCTTTGACAGAGGACATCCCGCGAGAGGATTTTAGATCTTGCGTTTCTTATCGGACACCGTTCCTTTAACGGGAAGTCACTAAGAGCCCCAACCGAAGATGATCTTGGAGGTTACTTGGGGACCAGAATATGGGTTTTCCGCTTGGCTCTTTCTCGATTGTTTACGGGGGCTGGGCAGCAGTACTCCCTCTATGTAGAGGTGTAGCAAAAATATGTAGAATCTGCCATCGGGCGTACGTAGCGGCATGAACCCGCAAGTCAGAAAGGGACTATAATGGCATTAAGTGCCTCGCATTAAACCTGTAGTAAAAACCCTGCCGGGGGAATGCGCTCTGTAGCAGACAACTTGCTAGCAGCAAGCGGGTAATGAAGACCTCTGCATGCCTTATCCTACTCTTGTGAAACTTGCCGGGCTATGTTACAGCTCTTCGGCAACCTGGGGCCTAGGAAGGGACATCCCCCCTGGACTCCGCCATCTCTGATTGCTGATACGCCGTCGGCACAGCCGGTCACTAAAAAAGTAAAGAGAATGTTATCCCTGTTGCCAGGTTAGAGCCTACGCATCGGTCCACAGGATACTGCTCACAGCGACGCCAAGGCAATAGTAGCATCCCACACTCTCGGACAGTT"

print(string.count("A"),string.count("C"),string.count("G"),string.count("T"))

reversed=''.join(reversed(string))
array = []
for i in range(0,len(reversed)):
    if reversed[i] == "G":
        array.append("C")
    elif reversed[i] == "C":
        array.append("G")
    elif reversed[i] == "A":
        array.append("T")
    else:
        array.append("A")


print("".join(array))

"""
Next Problem
"""

def Calculate(data):
    counter = 0
    for i in data:
        if i == 'G' or i =='C':
            counter += 1
        else:
            pass

    percentage = (counter/len(data)) * 100

    return percentage

id_percentage = {}
id_percentage['Rosalind_5135'] = [Calculate('CGGATCAGAAGCGCTATGCAATCTATCTAATCTGAATCTGGAACTTTAGGGAGAGTACCA\
GGTCTCCAACGGAGCCGGAGACATCACAGACTAGTGTTTATTTCTCCACTCCCCATGACT\
CTAGATACGTTGAGGACTCTCGGGAGACCGATTTCACACAACCGCCGAGCCTGGATATGA\
ACCATTAGAATCCTGTGGTCCTTCTCGCTCATGTCAGTGAGGCGGAACCTAGAGGATTCC\
CACTGCAATCCCATACTTATAATCATATAGCTCCCCGTGATAATCGAAAGCTCCCTAGTC\
TGCCGTCTCTGCTGGGAATGAAGAGGGGAGGTGAATATATCATTTTGAAATGTAGCGAAC\
AGCGCGACTATATGATCTTGAGTCTCTGCTGACCACCTATTCCGTCCAGTTACCAGCCAA\
AACGTGCTATGAAGGGTATCTCAGACTAGCCCTAAACTAAGTGAATCGTTACATCTTGTG\
AAGATATCTCGAAACTAACCGGTCCTATAGGCCTAGCGTCCAGTTTTTTGTTAAGGTCGC\
ATGTGGACGGTTGAGAGAGTGAGTAACCATGCATACCCTACAGGTCCTCGTCGACCCGGC\
ACAAATGAGCCATGCTGGCGGTTCCTTGCGCATTCACTGCTGTAATAGCTATTCCGCCCT\
CTAATGGTGGGCGAGATTACGTCGTAGCAGTGGAACCGTGGCCCTCGCAACCTACAGGCC\
ACGACAAGCTCTGAGGAATCATCACCTGGGGGTCACCAGTTCATTTACACCTACTGATCC\
TTAGATTAAGCCAGTGTTCA')]
id_percentage['Rosalind_6463'] = [Calculate('ACTTACCCTGAGTACACTCATTCGCGGATTAGCAGGCCAATATTCAACATAGCCAAGTAA\
AGAGTGCGCCCTATGGCCCACTCCGCATGGAGCTCAAATTTGACAGTCGCCGCGAAATCG\
GAAGGTCTTTAATGGTTGTGAAGCGCAGTCATAGCACCGAGCCACCGATAGATGGGACCG\
CGGAATAGGACGCGCCGGCCCCCTGAATAAACTCCAGCATGCTTCGGATCTGTAGGCACG\
TCTTATCACCGAAGCAATTCCATAGGGAGTATACGATACTGGACTAATCTATCACATCTC\
CTATGCGACGCGTATGAGCCTGCTCACAAACCCGAGGAGGTCACCTTTTCTAAGCAGGCT\
ACCAGCCGTCAGTCAATGGTCGAGACCTCACCCGAACATACTATCGTTTAGCCAATTCAG\
CAAGTGACTGGCGTAATAAAAGTAGCACTCTCTGGGGGATTACCCGGTTCTTACCTGCCA\
ACAGGCGTGGCCGATACATAATCAGAATCTATATACAAGCAGTTAGCATGTGATCCCGAT\
AATTATTCCCACAATGCTGATGGTACGAAGCGTTTTCCCTTCTCTTGCAATTATTGGTCT\
CCCGACGGGTGAAATTGCTGTTTAATTTTATTTTGGCAGTCTTCTGATCTGTCTTACGTA\
GGCCCGAGTAGGAGGTTGCAAAGTTCTACCGAGACGGGGTACCAGTATCCGTCGCGCGAT\
CCATTATCATACCACGCCACGCGGACACTGCGAATCGAGCTCAACGCGCCTCTGCCACAA\
GTTGCAATGATATTGCCGTCACACGAACGCGGAATAGCAGACCCTGATGACGCCCCTCCG\
ATCTCCGATGCCTTTACGGGCACTATCAAGGCCAGTCATAGTAGTATTTAATATTTGCGG\
AGTGGTCACGCAAGGGCTATATTACTCTACGATTGTATGGGGTCGCATGAACACACACAT\
ATGCGATATTTGGTCTAACCCATGCGATGGTGCACGACG')]
id_percentage['Rosalind_6363'] = [Calculate('TAGTTATACCCGCAGGTTAGAGTGGAGAGTAAGATATGCTGCCACGAGACTTAGGCTAGA\
GAGTTGCCAGCCGACCGACGACTATCATTATATGGCTTGCTGTTAGCACAAGGGCCGGCC\
AGATAGAATAGTGACTGCGCTTGTTTAATAACGACCCCTAATCAAGCAAACTTCTTGTTT\
GTTGAGTTTCGGTCGGCAAGCCTGTATTTCTCCGAGACCAAATCCTCGAAACGTCCGTCT\
TTACAAGCCAGTTATGGAGGGGGTCCGGAGGATGACGATGCCTAAGTTGGCCTTGTGCCA\
CGTCTCTACGGTGGACCAAACCCGAACCAATTTACTACTTCTTAAACGCGGGTATATGAA\
AAACGGCCGGGACGGAGTGCAGGTGCCAGGGCGGGATGAACAGAACGGATCAATCCCTTC\
CTTATTTTGCAATGGCGTGACGCGTTCGAACCATTATATTACGATACGATACAGTATGCC\
TCTGGTAGCTCCGATTTGCATATAGACGGGCGGTCGGGAGTTTACGTTTAGGCACCTTTT\
AACTCTGTGCTTGGTAGCGCCACCGTGTAAGCATTCCGAAGAAATTAACAGCTGTAATCC\
AGACAGCGGGTACCATACCGGAAGAGTGCAAAACCTACAGGAAGCTAGAGGAATGGATGT\
CTCGTTGTACTAGAATTGGTAGATACTAACAGGAAGACTGCTAATATCGATACCGAAGTC\
CAAGCGGACCCGTGAAGTGACTACCCCGGTGTTAGGGCGAAAGGGAATGGGACTAACTAT\
TTTCGTTCGGTGCGGAGTAATCTTTATCGTCCCGATGATACCAGGCTTCACATTTCCGTT\
GCCAACGTCAGATGTGAGAACCGTGTCCAGGTGCCAGGGGATTATTGCATCCTCCATGGT\
CCTAAAATAAGCACCGTACTAATCTTTATCACCGCTCAACGGGGCGTGGTCAGCGAATTT\
CGTTTGCCTTCCCGACTCATCATTCTC')]
id_percentage['Rosalind_7494'] = [Calculate('ATATCTTTTGCTAACCCCGCCCCTACCTTCAAGGGGCGTTTAGTCTGCGTATCGTCCCTT\
TCTCGCTAAACGTAATTCAAGCCCAAATTATAAAAAACTACCAGCGAGAGATTTCGAGCT\
GTTGTTGTAGCTTATACCTTACCCAAGCCATTCCCAAAGTCAGTAGAGAGAAGGTAGCCA\
GGCAAGTATTCCACCGTGAGGCCTTACGGAACCTTAATTCCCTTATTTCTTTCTGTGAAA\
TGAGGGATTACCCCCGTAGAGCGGTCCACGGTGGTTACCGCACTGCGAGTAGCATATCTG\
GCATGTGCTGTCGCATAAAGAGCATCCTGTACAACGGTATATCTTATGCCTGGTATCGGT\
GCGATAGAACGTGGCAGCTTTGTTACGGGATGTCATATCCTGAGCCTCCGTTCATATTCT\
ATTAAAGTATACCCTCTACACTATGTGGGAGGTTCGGCCTCTTAAGAGATATTACATGCG\
CAGGGAACGCCAATCGGGCAGAATTACTCAGCGTGGAGACCCTCCATCATAGCGTCAAGC\
ATCTCCGCCCTACTCCAGTTGACGCTCTCGCAGGGTAATCTACGATGTGCGGCCGATACT\
GCGGTTATGCGACTAGTACGGTTCCGACCTTTATTCACGTCTTTTAGCCACCGCAAGTCT\
CCATGCGAGGGTTTTTCAGACCGTACACAAGGATCGGTCCGTCGCCTACCCAAGGGTGAC\
AGCCCTATAGGGCCGCAGGTGATCTTTCATAACCTCAGGTACACTTGGCACGAACTCCAT\
GCCTTCTGTTGGAGAGGAGGCGTTAAGCTCCACCGGTCCACAACAGGCGTTAAGGCTTGG\
TCGCAGGACTATC')]
id_percentage['Rosalind_7658'] = [Calculate('AAAAACGCGCCCTATACCAAACAACTCAAGTGTATGTGCAAAGCCCATTGTCGAATTAGA\
ACGTCGGCGCTTTACGGCATTGGTGGGGGTGCACTGTTAATTACACGGACTAGCTTCCGC\
CAGAGCCAGCCGTCTTATCGAGACTACTACTTTCAGAGTTGAAACTACAGCGATAATAGT\
AACCAAGCACTATAGCCCCCCTAAAAATCCTTCACTATCCAGAACGCAACAGTCTGTCCA\
GCGAAGGGACCGCACTCCAGTTAGTCTTAGATCCTTCTTCGTACTGCATTCAGAGGAAAC\
TTACAATGTTACTAACCTCACTCTATAGAAGCGAGATTAAAGCGATGAACATATGCTCGT\
CCGCCATCCTACCCCGGGCTAGCATTACAACTCGATAGACTGTTGCCCCGAGAAGAATAT\
AGCACTAACAGGACGCTTCGTTTTCCTCACGGTCTAAAATCCTCCCATCGCCTTCGAGAC\
AGCATTAACAGCGAGAAGTGTGAAGAGCACATATCCCGGAGTGTTGACCGATGTTTATAT\
AAAGTTAAAGCGTACATTCTCTGTGGTAGGGCGTTATTTGACTAGTGCAACCAAGCCCTA\
TTCTGTGAATTCCGATGTGTCACTACTTGCACCTCTCATGGCGGGGATCCCCTGGGACGC\
GCTGCTGTGGATGTCGGCGTACCGCCGGTCCACGGTGGAATTGAGGGGTACTACTCTCGG\
AACCAAGCGACTCAAAAAAGATCGAGAGTTACTTAGACTGTCTCCAAGCAAAACCGAGCA\
AGTGGGGTCTGCCGATGCATCGTACTCCGCCAATTGCAATCGTACAAGGTAAGTATCAAC\
GCCGGACGGGGCCAGATTGTCGACCGCGTTATTTAATGCCGCGCGCGAAATTGACAATCC\
AGTTTCCTAGTTGATCTCATTAGTTGACGGTCTGAGCGCCTCATAA')]
id_percentage['Rosalind_7840'] = [Calculate('GGTATCTTCCAGAGCTTCAAGAACAACACAATAGTTAGCCTCTGCGATACTCTCGCGCGT\
TTGGCTGTCGTTAGATGAACGGCCGAGACTCTTCAACATCTATGCTGGCTTGCTCTTATT\
CTTGACGTTATTATCAGCGTACAGATATCTTTCCTCGCAGTTTGATCGTACTTTGTGCGA\
CTCCGCGTATCGAGGGCCTTAAAGGTCTCAGCGCATCTAGGCTCATTCTGTGAAGGACTA\
ATGGCCTCTATGAACAATGATCCCGTGATTGCGTGACTCACAATAAGCCCATAGCACCAG\
TGCTAGTGTTATTATTCGACCCAACAATACTATTACGGGCACAGGGGACTTGAAGGACAG\
AGCTCGAACCGGTTCCGAGCCGAACAATGTACTCGCCATTTGTTCGAATCAGCTGGTTAG\
AATCTCGGCTATCACACCGTCGCGGGTTCTAAATAGAAGGTGCCCTGAGTGCATAGAGGC\
GTGGGGCCGCATGCTGGACCAATATTACCGAGGCATTTACTCCACAGGAAACCATAATGT\
TATAGGGCGAACGGGTTTCGTCCGGTGCCGTGTCAACAAGGGGTTTCTCAGGTCCAGCAT\
TTTTGAACAACTCTCCGTGAAAGGAAATCGCTGAGTAACCAGCAGATTGAATTTGGTTAA\
AAAGCAGTTGCTGAATATTAAGTGCTAGTCAGTACTATACCGGATATGGCTGATCGCAGT\
GCGTTGAGATGACCGCAAAGTACGAGGGTCATATGTACAGTCACCTCGTATAGCTCGTGG\
TTCCCTATACCCTCACGCAAGGAGTGACATACGATAACCACAATGTGGACGGAGATGATC\
GTAAGTTTTCACCTATTTAGACTCGTTGCCACAGCAATTCAAAGCCAAACATTGTACTCA\
GTCATCGCGTGAGGGGTAGAAGGCATTCTTAACCTTTGACCTAATACTAAGT')]
id_percentage['Rosalind_6170'] = [Calculate('AGAAGGCGGATAGGGGGCTTCATCGCGCGAATTTACTACCGTGGGATCTCGACCGCTGGG\
ATGGCCTTACCGACGGGTATAGCATCCTAGCATAGGGTCATCCGCTGGTTGAATTGGTGA\
TTAGGGGCACCAAGACTTTCACCCAACTTGTTCCTACATTGTGAAACAAAGCTTTCCCTT\
CCGAGAGGGTTACCGGCCTTAGTGAGTAACACGCGAGTTGTACGCTAAGAAGACAGTCTT\
GACTAACCTATGCATGGATGGAGCTACCGCTCATGAACCTCCTTCTATGCAAGTTCCCGC\
TATACGTGATTGAGGAAAACCCGTTCGGAAGCCGAGGAATGTGATCACTCATTACCCCGA\
TACGCTTTAAAAGAAAGTACGGATCTAGCGGGACGGATGGTGGTAGAATCGGCCCTAGAG\
TTCGTCGCCGGAGAAGTGGCGGCAGGAATCGTGTACCAACCGTAATTCCCGCTTACACAT\
TTCAAATGCGGAGGGACTGGTCCTGTAACACGCAACCGGCCGGGTGAAGCTGCTCGTTAG\
CTCGCAACCCTGTCTGTAGAGTTAACCCAGGACCGTTCAGTGAAATGGTAAAATAAGATT\
CAGTATGTCGGACCTACTCGACAGAATTAACGATTGCTCCGACACATAAGGTCGGCTTCT\
TAGAAGAATTATCAAGATGTGAGAACTTCTTGTCAGCGGAACATTTGCATTCGGTAAGGC\
TGATTGGGGGGGAGTTGCTCAAGTCGGGCTCAGCTAAGCCGATACGACTCGCTCCGGCGG\
GCGTACCCGGCCAGCTGGTCGCCCGTAGCAACCCATATCAGCACAACAAATCACATAAGC\
GGGGGCGGTCACCATTTCATCAGATTTTAAGAATCTCTTCCTAGCCGGAAGGTTGATGGC\
ACACTGTGTGCTATTTCTACGCGACGTCGCTTG')]
id_percentage['Rosalind_2981'] = [Calculate('GTACGTTTAGCCAAGCACCACTGTGGCAGCCAGCAGTAGTCCGATCTGAACTCAGGGCCT\
ATTCCATATTTGCCGGCGGGAACCTGAGGCGTTCTATGCGATCATCTATCTGGGAATCGG\
ACGCAGTTCGAGACCGCAAACAGTTGTTCTCGTCGGTTACCGGGAAAAAACGACCCCGGT\
AGGCGGCCCTGTTGGCGCCTCAAAGACGTTCGTGCGACCAATTATAATTGATCGCAATAT\
CGAATGCAGCTGCCACTATTGCAATGATTCATAACAGGGATATCTAATAGACAGTTATAC\
GCGTGTCCCATGCGGTACTGTCCTAGTCGTGCGATAGTCCCTCCAACTCCTCGAGGTCTA\
TTACGTCCGTGGGTTCTACTCACGTATCAGCAAACGTCTACCGCTCTAGAACGCTGACTT\
CCAGCTTGCGAAAATTGTCGTTCGCCACTCATCGTCACTAGTGGGAGCGGTTAAATAACA\
CTCACTTACGGCTCCGTGTGTTTGAACATGGGAATGAAACATAATATCCCGGCTTCGGTC\
GTGGCTGCCCGAGCTAACTTTGCAATTTCGCCCTTTGTTGTGTCGAAAGTAAGGAGTTGG\
GGCGGCAAGCAGACGCGAGTTATTCGAGTACCTTACTATAGGTCCTAATGTGGACATTAC\
TGGCACTGCCGGCACAAGGCAAACGTACGCATCCGTATCTGGAAGATAAAAGGAACGCAC\
CACTACCCCATTGAGCCCACGTTTTGGGGCGATTGAAACGTGGTAAGCGGGTAGCACCAC\
CATTTAATACGGATAACAAATGTCGGGGCGTTACCGACCAGTAGGTGGATTAGAT')]
id_percentage['Rosalind_0109'] = [Calculate('TCTTATTCTTGTAGGGTGGCCGGAGGTAGTTCCTCAGTGGGCTGAAAGGAGGGCAGTCCG\
TGAGGCGTGCTTGGCTAGAAATACGCGCTCCTATGTCAAAGATCACGTCACCTTGTCCCT\
TGCGTATGATAGGTGACCCTGATCTACGAGCCGGATTTCCATTAAGTATACCGTCTTCAC\
GAGAAGTGTGTGTTACCGACTGAGTGGCCTCTAGCCGTCATCGTTCTAGAAGAGCCGGCG\
TGCTCGGATAATTCCCTTTAGTAATTTGGATAAGGCTACAAGTGGGGACCATCTTGCTAT\
AGGCCCCAAGCGTGGGGTAATAGGTATGAGGGCCAGACGGCTCAAACAATGAGCGATACT\
AGGGCACGTACAACGGAGTTTGAGTACAACAACAGAGTGCACTGGTTAGTACGCTGAAAG\
CCAACGGTACGTCGGAGCTAGAGCAAACGCTGGAGAGTATCGCCAACACCCTTAGAGTGC\
ATTCCTCGCAGTCCGCCAATGTGTCCTCCCATGGACGGGCGGTGGGTGGATTTGGCAATA\
ATATCCGTAAGCAAATTTACACAAACGTCACGCCAGGCAAGTAGACGACATCATGACGTA\
GCTATTTTGTAACCTAATCACGTTGCGGGGGGTATTCCGCGGTATGCAACAGGCGCCCTC\
AATAAGCATGCGCAAAATGAATCTCTCATCATTTGGATTCCCAGATATGTTATACCGACT\
TGTCCTGTTTTTGAATGACGCTACCGTAACAGTATGATTCGTTATATACTATGTCAACTG\
AGATGCTACATTGCCCTTGTCGTGGATAGTGGACACGCACAGTAATAAGGGATGAGTGTT\
TGTGACCAGTATGAAGGTATACGGGTAACACAGGTATACGGCCTGGAATCTGCTAGCTTT\
CGTCGCAACGGGAATTAAAGGTGTCCAAGGTGACCCACTGCGAACAATATAAAGGGACCA\
TCTAGTTAGCGAAGTTCAGTTTTG')]
id_percentage['Rosalind_6803'] = [Calculate('TATTAGAATCACCCGCTTAGACCATTAAACGAGATATAGCTTAGCATGCTCCCATCGGTT\
CCATTAGAGTCGGCCCCACAGCCTCCTCAGGTTAGATCCAAGTAGCGTAGCACATTGGGA\
TCTCATGCGCACTAAGCGGATCGGGGTGCAGCGTGAAGCAGGAATTAGATAGCATTCCCC\
ACTGACTCAGCAATGTCGAGGTGGCTGTGGGATCATTTAACCATAAGTATAGCAATAATG\
TTGCCAGCTAAACGTTAGTCAGTTTCGTGCAGTCGGCGCGCTGGTTGGACCCCCAGGTTT\
CAATAGTATTTGGTGCGGGAGCACGTCATTTTAATTGATTTTTGCGGTCCTATACGTACA\
CTACTAGTGGGCACGGTGGGGCGCTTGACCTCTCTGCAGTCCGTCGCCCGGCAATCTATC\
GTACGTCCGCCAAAGAGCTCAGGCGGAGACACTGTCTATATTGGCGGCCCAGATCGAGTT\
GTGCGTCTATAACCGTCTACCCGCCAGTGAGTCAATAGCACGGAACCATGCCGAAGCTAT\
ACAAGACAAACATATCATATAATATATCTGGCCCCCGTCTTCTTACCACGGATGGTACTG\
CGCTCACGGGAAACACAACTGCGCTCAGGAATCTTAACTAGAAATAGAGCGGCTACGATC\
TATCAAGAGCCCCCTATGTTGTCCTAGCCGTGAGTATGCAGTGCCGTCAGACACTTCATG\
TTGATGGTCTAGTCCCGCACGCCTAAGACCCTAACCAGGACGATCCCAGAATGTCTACTG\
AAGATGGCGCGTCATTGAGATTGCTAACTAGCTGGCATTCGACCGGCCTTTGGCGAAGCC\
GATAATCGGCCTCTTTCGGAGCCCCTGTATCAAGCAGGGTATCCACTTCCTTTCAGCTGC\
ACCTGAATAGAAAGTCAGGATAAGAACACTGCT\
')]
print(id_percentage)
"""
Next Problem
"""

def CheckMatching(s,t):
    s_arr = []
    t_arr = []
    match = 0
    for i in s:
        s_arr.append(i)
    for x in t:
        t_arr.append(x)

    for letter in range(len(s_arr)):
        if s_arr[letter] == t_arr[letter]:
            match += 1
        else:
            pass

    return len(s_arr)-match


s = "CACCCTACGTCAAAAGGCAGTGATCTCACTTACTTTGCCAGCAATAACGTATGGTGAGGACTAATGGCGGTTGTAGGCCCAAGCCCAAATCACCTATGCAGTGTTCTACTCTGTCGCACGATCTTTAGGAGCCGTAGACAAAATGCCAAACCAGCATCAAGAATGGGTGACTGTTCAATCTACGGTCCCAATGCAAGCAGCTCGTTAAGTTACCGTTGCCACCCGGGTGTACAGAACAGAAGTCCACGTCCTTACTGCGATGTTAGCTGCGGAACTCCAAAGTTAAAAGACCTAGGGCCTGCTCTTAAGCTGCGATAGAGATACCTCTGCACGTACGTGGCATTGCCAGATCCGCTTATGCTTGCATGGTTGCGATCAGTAAAGTTGATGTGAATATTTTCCCTGGGTTTATATGGACGGCCCCACCTTTCTTATGCTATGTTCTGGCTCACTACGTCGCCAAATTCGTAGGTGCCTCATAGATTCTCCAACGAATCCTTGGACTGAACCAGATTACCAGGTAGTAAGTGAATCCGTCGGAGTGACGCTTAGGCGTTTGAGCAACGGAGGTGGTGCCCAGAGGGCGTGGAACGCTCTTAATTGATATAGACAGTCTATTGTTTGATTCGAATAAGGCTGACTTTATCCTTCCCGGTATTATTTCTTCACTGTTTGTGTGTGATTCAGGAACCAGGCTTGTGTATGGAATCCTCTGCAACTAGTTCAAGGGTCTGAAACCTCTGACATATCTCTCTCCACGTTCCTGGCGCGGCTGCGTGTACAACGCTTAGAAGTCTAACAGTCCTCACGACATGTATCGTGACCTTATACCCTACCTTCTGGTTTTCGGGAAGCTGGAGTCGGCCAGGTGTTACTCGCTTTGTGACTTGCCAATTGCCCACCTGGTTCCCACCATTCAGAAAATGTGCGAGAATACCCTCTGTCCATAACAAATCT"
t = "CTCTCGGAGTTTATGAGCACTGCGTTAATTTACTATCCCATAAATAACCCTCAGCTTGGCGTATCGAGCGAAGGAGACGGCAGCATATATAAGTCAACAACGTCTCTACGCTACCCTAAGTTCTACATATCCCCACGGGAAAATCTCACCACTTTGATAAGATAGTCAGCAGATTCAACAGACGAACCCTATTAATGCGGATCCTTACTAAACCTATACAACAAGCAGTTGCGCAGCTAAAGTCCGTTTTAGTAACGTTATGCTCGCTGGAAAGTCTAAAAGTGTAAAGGTCGTGGGCATGCACTTAAACTACCTTAGACGGACAGCTACCAGGTGGGCGCAATTCTAGGTCCGCTTATACTAGCACGCGTGAGGTCAGAGAGGTGAATGGGCTTTGTTCCCTTTGGTTTATGTCTGGTCCCCAACCCCGGTAGTGCAACGGTCAACCTCATTCCTTCACCCAATTCGTCTGGGTCGTATTGATTGTTCGCCGTCTCCTGGGACTGAAGTCATTGTTTATATAACAATCGAGTATTGCCGTTTGAAATTTCCGCACATTAACAGCCGAGCTGACACCCAAAGGACTTGGAGCGTGCTTCAGTGGTATTACCAGTGCGATGTGTGAACAGTCTAAGAGCATCATTCACCGTAAATCTCATAATTACTCACTTCTTAATGGTGTTTTGTTAGCAGGGATGCTGTAATACATCTAGGCGAACTTGTCTATGGTAGTGAAACCCCTTAAATACCATCATCCAAGATTGGCCGCATGACCCATCCACAACCCGAAGAATTATAACGGTATTAACACAATCCGTAGTAAAATGTCGGATTTACTTCGGGAGGAGGCGACCTTGTAGGCTGCTACGTATGTATCGCTTTTTTACGTACACAGAACAGAGGTGGGATGAACCTTTCCGCTGACTGTCTTCAGTGCTCTATGAACGTAAGCAAACA"
print(CheckMatching(s,t))

"""
Next Problem - Converting RNA to Protein
"""

def RNAtoProtein(data):

    data_arr = []
    protein_arr = []
    for x in data:
        data_arr.append(x)

    codonTable = {"UUU":"F","UUC":"F","UUA":"L",
                  "UUG":"L","UCU":"S","UCC":"S",
                  "UCA":"S","UCG":"S","UAU":"Y",
                  "UAC":"Y","UAA":"Stop","UAG":"Stop",
                  "UGU":"C","UGC":"C","UGA":"Stop",
                  "UGG":"W","CUU":"L","CUC":"L",
                  "CUA":"L","CUG":"L","CCU":"P",
                  "CCC":"P","CCA":"P","CCG":"P",
                  "CAU":"H","CAC":"H","CAA":"Q",
                  "CAG":"Q","CGU":"R","CGC":"R",
                  "CGA":"R","CGG":"R","AUU":"I",
                  "AUC":"I","AUA":"I","AUG":"M",
                  "ACU":"T","ACC":"T","ACA":"T",
                  "ACG":"T","AAU":"N","AAC":"N",
                  "AAA":"K","AAG":"K","AGU":"S",
                  "AGC":"S","AGA":"R","AGG":"R",
                  "GUU":"V","GUC":"V","GUA":"V",
                  "GUG":"V","GCU":"A","GCC":"A",
                  "GCA":"A","GCG":"A","GAU":"D",
                  "GAC":"D","GAA":"E","GAG":"E",
                  "GGU":"G","GGC":"G","GGA":"G",
                  "GGG":"G"}

    for i in range(0,len(data_arr),3):
        print(i)
        print(data_arr[i:i+3])
        if "".join(data_arr[i:i+3]) in codonTable.keys():
            if codonTable["".join(data_arr[i:i+3])] != "Stop":
                protein_arr.append(codonTable["".join(data_arr[i:i+3])])
        else:
            print("not there!")


    return protein_arr

protein_string = "".join(RNAtoProtein("AUGAUCAGAAUGGUCUUCGCCAGCGGAUGUCCGUUCGGUCUAAGAAAGAGCUCACUACGGGCCCAUGUUGGGCACUUUCAUGUCGCAAGCAGCAUGUUUCGGCUCCUAUUAUUAUGCCCGGAAGUGGGGUCACAUCCUCUACACGCCAACAGUAAUGUGUACAACUUGUCUAAUCUAGUACGCAACUACAACAGAGUGGUAGCAGUGGUGGCCCGCGUCGAAAGGCGAUGUGCGAGCUCAUGCCUAGACAUGGCUCUCGCAUAUGCACAUGGUAAUAGUCGCAGGUUAUGUUCCGACGCGCUAAGCGCGUCGAUUUGCCGAGAAUUCGUAUCCAAUUAUCAAUCUUAUUCUCCUACGUAUCAAGAGGCCGGCGGAAGCUCAGACUAUAAAGCCACGAAGAUUGAAGGAGCGACGCAUGAGGUGCAAGGAGUAGAAUUCGACGGAUCCAACUAUUAUCCGGGUUUUAAGCUGCAAAUAUGCAUGAUGCCUGGGAGCAGGAGGUCAUCCUGCAUGGCCCGCCGCAAGCGUCAGUGUGAGUGUGUCUUUAACUGUGGCAUUCGGCACAUAACAAAUAAUCUGCCGAACGUUACAAGGACAUCCAUACCAAAUAGACCCCAGGCGAUUCAUAUACUGUGCCUUCCAAUGGGCCAGGCACACCUCAUGGCCGGUCUCACGUUCUCGAUCACAGGUUUUCCACUAAUAUCAGUCUCCGGAUUACCUGCUAUGGCUACAGCUUCAGUCUAUGCAAAAAGAUUGUCUUCGCGAUGUUCAUGCUUAGAGAUUCCAUACGGUAUGAUUUUGGGUAGAGAUGUCUGCUCCUCUGUUUGUAAUAUCAAUGGUUUGUGCUUUGUGCUGCGGCGAUGUUCGAGUCACGUCCGCGCGUCAUGUGGGGCGCCCGAAGCUUGUCAGCCCCGUUGGAACGAGUUUGUGGUGAUCCGCACGCGGCACCCUUUUUGUUGUGCACUAACUAUUUUAUCCGUGUAUGGCUCCCCACCGGAGAGGAAGGCUGUUUUCUUCGCCGGUAACCGAUUCCCUAUUGGUAUAUCGCGACUAAGCUAUUGUGUGAACCAAACGUCCGGGCCAUCGGAUAUGCUAUGUCGAAGGACUCACGGAUUUAACGGCUAUACAUGCAAGCGAUCUAUCCUCGCUAGGGAGCCUGAACCCAGCCCGAGGGCUAACGAGAGAGCACUUGCAAGCGCCCUUAGUCGCUACCGGUGGAGAGCUAAUCGUUCGUGGAUCCUCGUAUCUGAGGAGACAGAAGGUCCUAUCUCAGACGACUCCCAUGGCCUGGAGUUAAAAUUAAUCCGAGUUCUGACGAACGACUAUGUACCCGCUUUUGAGCCAGAAGUUCUACGCGUAUCAGCGAGUGGAACGACACUUCAAGUCCUGGUCUUAUUCCUAGAUGCCACUCUCCGAGGGUAUCCGUUUUUGUGCCGUCCGAAAGUUGGAGGUACAGAAUGUAUGUCCUUGGAGCGUGUAUGUACUAGGGCAAGAUGUAUAGAACAGAUUGCGGUAGCGGCAUGUCAUCCGGACGUCGCCAUGUCGUGCGUUUCCGCCUCAUUAAAGGGGCGCCCAGCUAUUACACGCAUCGCAUGCUGGGGGAGUCAGCCCGCAGGUCCGCAUAUGUACGAGCUUUCGGCGGACAAAGGUGUACUAUUGUCACCGAUACAUCGCAACCAAAGGGUGACCCUACUUGAUGGACGACCUUCCUUCGGCAGUUGUGAGGAGCUCAUAGCAGGUUUUAUACUUCUCAUGAGACGUCCUGAAAUUUACGACGUCGCGGUCGAUGUCGUAAUGGGUAGCUGGUGGGGGCGAUACACUGAAACUAACGCCGGGAGGGUUAUUAUGCGAUUAGGACGACCAAAUACUGGCCCCCAGAUACUACGCGUGGAGCGGGGACACGGAGGCCCUAGCAAGCCAAUUCCGCACGUUGCUAGCACCUCUAUAAGAGCUCCGCGCGAGAUCUACGUGAGAGGCCAUGCACGUAACCGUUCGUGGCAGUUAAAAAGAUCCGUCCCUAUGACACGUAAUGGGUGCAGACCAACGAUUUCAAUAGUUAUCGGUCGACGCCCGGAAUAUGCGGCUUUGACCAAUGGUACGAGCUGUGGGGGAGUCGAUGUGUCUCGUACCCUUUUGGAAGGGGAGCAAACACAGAUUGCACUUGCCAACUGCCAGAGCGCGAAUGUUUUACUAACUAAAAGUUGUCUUUCGAGCAGUUUUCUUAGGCCAGAUUCAGAAGCAACUACUUCAAGACCAAUUAGCGGCACUUAUAUAAUAAUAAGUUAUACAUUUGAGAGGGACUACAUAAGUCUACACCUCAGCUGUUCGUCCUGGUAUAUGGUAACUCAGCGUACCUUGCGUGCUUAUACCCCUUUGCAUAUACCCCGUUGUCGCCUCGUUCCAUUGUACUUAGGAGGGAUAUAUCGUCGGAGACCGACAAUACGCCCUCCGUCGCGGAUCACGGAAGUCCAUAGUGUUGUAUCCGGUCGUAGCUCGGGUGUGGAAGACCGAAUUCCUGUUAUUAAGAGUGGCACUAAUUUUACUGGGGACGGGGCGGUAUCUCUAAAGUCCGCGCUCGGGACUGCAGUGGCGCUAGGUAUGCACGCAAAGACGCUAUAUCUGCCUCUGGUGGCCUAUGAGAUUCCGUAUACAGUACGUCCUUCUCACUCUUGUCCUCAAUCAGAGUUGUAUAAAGAGUCCAUCAUCAAGUGUCGGCUGAGUGCUUCACCGUUGUGGUUAAUUGGCGCAGUUCUUGGAUUGCACGAAAUUAAAACGCUGAAGUCUACAUGUCUGAACAUUUCGGUUCUCAGCUUAGGAGCGCGCAUUUGCUGCCCCGAACCCAUUCUAGUCCGCAUGGUGAAACACAGUGAUCACAAUGUUCCCCUACAUUGCUGCGUGAUGCCUGCGCGGCCGAUGAAUACCCCCGCUUGUCAUCUGGAGUCUAAGCAUAGGAACACCGCACUCCUCGCACACCUGAGUCGGUUGCCUGAAAGAGUUACGAGAGGUCUUUCGAAGAAUGGGGCGUGUACCUAUGCUUUCCGCAGGUCGGCAAUCCAGGACGUCCCGCCAUUUCGGAUGCCCCUUGCACCGCUAUUAACCUUGACCUUGAAAGAGCUACGAUCUAGACUGCCACGUGUAAUUAACCCCCGGGUACGGGAUGGCGGCGGCCGGUUCGGCUACAGCAGGCCUGACGCCGUACCGAGUAUCCUUUUUAAAACGUCAGCUGCCUCGGAUUCUAAGGAUAUUAUUUGGCUGCCAUCCGGGGAAAAUUCAUCAGCGUCAUCAAACGGAUCAGUGACGGUCAACGGGUCAGGCUUGUGCUUGUUAUCACAACUCUUUGUAGACGUUAUGGAACGCCCGGUGAAUGCCGGAAAACACCCGCGCGACCAACGCUAUGGUAGGGACAACUGCUCUUCUCCCGUGCGUCAAGAGCUUAGCCGCCUUUACUCCAUCGAAAUAUAUCGAGGCUACGUUGAGAACCGUACACAACGUUGUGAGAAGGGGCAGAUGCGGCUGGCGAAAGUCCGUGUGUUCCACCAUCGCAUUCGUGGCUGGUCCGUGGCGUUAGUUCAAGAUGUCGCUUAUAACAACGUACGAGGGGCUCUUAUAUACGCGGUGGUCUCCGUGGGCGUACUAUGCAAUAACCCGAUUCGCUUUAGACGGCAUUAUUUUGAUGCAGACCACCAGUUGCCAGAAGUACGGUAUCGUCGGCCCACGGGAUUUACGUUGAAUCUGCAAUACAACGUAACAAUUAAAGCUUCCAUCCCCUGUAAAUUUAUACUUAUACAAGUACAGAUUCUCAUGUCGAAUCCUAGAAUACGAACCCCCUCUGACAGACCUCAGGAUGACUCUGUGGUGUCCAUUGAAGGCUACGCGAACUACCGUGUAAUCAAUCUCUUCGUGCCUGUAUUGAUAGGCAUACUCCAGAUAUCACUGGUCUCAAUCAAUUACCGUCAGUUCCGUUGGACGCGGCUACUAGGACCACGGCACCGCUCUCUUUACAUGAACUGGUACCCGCCCUGUUUCGCAAUUACAGGGAGUUUUAUCCACAAAGCACAUUUUGGCUGUGGGGCUCCUCACGGCGACCUCAUCGCUCGGACUAUUGCUUCCCUGGAUCUACGUUCAAUCCCAGUGUGCUAUCCACCUUGUAUUGUCGGUUUGGGGGCUGCAUGCCUAGCGAUGAACCUGCCUAGGAUCAUGCGUCAGUGCCCACCUCUUCGUAGUGUAUCUUUAGUGCCAGACCCAAAUCCACGGCCGAAUUAUGACACUGAGUUUGAUGUGGAUCUGAAACUAUUCCCGACUUUGACUGAGUACCAACCCGGGCAGGUGCCAAUUGUUCAAGGGUACAUCCACCUCACCCGCAGGGGUCGUAACCCUACUCGCGGCGAGCAGGCAUACGAAGAACCUUGUCGCAGUAUACGCGAUUUUUGUUCACUAUUACGAAACAUCCCUUACGCGUUCAUUUUCCGUGAUGAGUUUGCUUGCGGCCGGCUAGGAAAAACUGUAUCGUCCAUGUUCAUUGGUCACCACAUUGUGUUUACCGGGCCCCCUUACCUCAACUUAGUGCUAGACGUCUUGGCAGCGAGGGUUGUAUUGCAGUCGCUGGUGACGUGGGUUUCUUGCGCAGAGCGCGGUCCGUCCCCGAAUACAACAUACCUAGAGCGUGGUCACAUCCGGUUAUACGCCGGAGCGAAUGUUUCACUUCGAUCCAAAUGUCUAUCACUCGUAUACUGGUGGGGUCACAAGGGGACGGCAUGGGGUGCCUCCUCGCAGUGUUCGGUUCGCGCCCUCGUCCACGUUCUACGUGAGGCUCUUUCCAGCUAUAAGUACGUCACUAAACUACACCAACCCGUCCAGCAACAACGCGCGGUGUUAGUCCUCACGAACUACAGAAUGUUAAAUAUUUGGCGCGAACGGCGCGUACACAGCGCUUCGAUUACUCGGCCUUCCACUCAAAAGACACUUUUCGUGUGGCCAUUUGCUAGUCCCGCAUUGGUGGCUUCGCCGGAACCUGGUGCGCACUUUUAUUUGACCACUCUCCAUAAUAAGUGUCUACUUGUGCUAGUAUUUAAAAAGAUCCGACAUAUGCUCGGUCGGCGUUUAGGCUACGGUCUCCGGGGUGCACGAUUACAAAGUAAGUCAUUUGCAACUACUCCGCUUCAAUCGGUUCGUGCAUCUACUCUGGUAAACGUUCUACGGUAUCUGUGUGUGAGUUACGCCGAAUCGAUUGACGCUACACCACGAACCCCCCGUGGUAUUGCACAAUCGUCUGGAGUACUCAAGGGAGUCAACACUUCCUUGUUACAGAUUCUAGGGGACGGGUUCCACGAUACUGCUCCAAAUCUUGUUUGGCUGAAUACGCGAUCGCCGCGUCUGACGAUACUUCCUAACUCUUGCGUAGUAUCCGCUGGCAAAGAAUAUGGCGGGCACUGUGGCCCUGGAAAGGCGGUCUCCUGUGAGCAGGUGUGCACUUGGCCAUACCCGGAACUAGACGACAUCUUACGGAUCGAGCUCGCCCUACACAGAACGUUAGAUACACACCGUGAGUGUAACUAUCGGAGCGGAACGGGCAACGGGCGCAUACUCAACAUCUCACUAAUUGCUCAAAAAGAGAAGAUAGGUAGACCAGAAUCGUUGUAUUACCGUCUCUCGUGUAAUCGCGUGUUCAUCAUCUGGGUACCGAGAGAUAUCCCUGGGAUCCGUUCAGACGUUUUUGUAGGGAGCGAGACAUUCUCGUGUAAUUCCAUAAGGAACUUUGGUUCUACUUCCCCAUGCGACUGCAACAUAAAAUUGUUAAAACCUUCCUCAUGGGGAAUCUCUUGGAGGCAACCGGACCCCGACAUUGAAAGCGUUCUUACUUGGUUAACUCAUGAAGCAGCCGUAGCAGGGAGCCACGCGCCGGUUAACGGCUGCUCCAGCUUGUGGCCCAUGUGCGAUAUAGUCUGUGUCUGUCUUUCCUUGGGACGUAAAUACAUCUUCACACUACGUAAUUUUGAUCUGUACCCAGACCGGAAAGUGACUGACCGUCUGAUGAGGGAGUUGAAAAGAAACAAAAUAUGCUGGCUAGUAUUAAUUAGUUUCAGCCACUACCGUCUUUAUUCGAACUUCAUUUCAAACCAUCUACGCGGUCGAUUUUAUAGGUUAGAGAGUCCCAAUAACGGUCUGUGCGCUAUGAUUGAAGGAACCGUGAACACGAUGUCGGGGCUAAAUUGGGUCAUAAUCGAUGGAAGGCACCCGGCGCACGGUUAUCAGACACUGCAGAAUUACGACACCAAGCAAGUUAUAAGAAGUAGUAAACGGCGCCGACACUCGCCGGAGGCUCCAAUUCCGCAGUGCAAAUCGACUUUGAUCCACUUCGCCAGGAAUCAUCCACCCGCGUAUGACGAUCUACGGAGCGAUAUACCUAGCUUUUACUCAAUACGGGUAAAAAGGCCAGCGACGUCACGGCACGAACGAACUACCGUGAUGUAUAUCCCAAGGGCGUUGUACGAAACCUUAAGGUUGUAUAUGGAUACAUCCGAGGCUUUAGGCUAUCAUAUAUCGCGUAGUGGAUUGAUACUUAUCACGAAAUUGGGCAAUAAUAAUCGCGACCGCAGCCUCCGUGCUAUCAGAGUUUCAUCGUACGGGUCUAUCCCUGCUGGGAUUGUUUCCUAUCCGGUCGUGUCCGAAGAUCGACCUACCUUAACUGGUAGAUGUAAAGAUUUGGAGGUCCUGCGCAAAGCGGAGCGUAAGCUGGCUUUCAUUGCGUACACUCAACGGAGUUACCGUACACACGGUAAUCUGACGCCACGUUCAUGCCGGCCUUUUUCUUCAUUGCUCCAUUAUGGUAGCAAGUAUAUCCGCCUUUACUGGUGGGGUGUCGACACGAGAGGCCUCAUGCGACUUCGGUCAGUCUCUGACCCUCGAAGACGAAUCCGGUGUCCGAGGUGUCGUCGCAUGGCGGAAGGUAUUGUUACCGUUGGGAGAGGAAAUCCGGCCGACUCCGAUAGAGUGCAAACGCUAAAUUGUACCGGAAAGAAUCAGGAGCCCCGGCAUCCCCGUCUUAUCCGAUUUGUUACAACUGUCCGGACCAGCUCGACGCCUCGAGUAGGCAUGCUCGCUUGCGAACGUGCUGGGUCUUGGGGAUGCCGGACACUGAAGGAUUACUCCGCGGGGUACGUGCCUGGCGAACAUAUCUCACCGUUGCGUCACAGUCUGGCAGAUCUCGGCAGGUCUGCGUCACGCUACUUUGGAACACGCUUAAACCGUCAUGCUCCUCCGCUUAGGAAGCAAGAAAUAUCCGUUGCCGGCACACAGAACCGAAUAGCGACGGAUCUAAGCGAGCGGUAUGCGUUUGCGCCUCAUGCGUUGUCAGUUUGUCGGGGCGAAAACGGGAAGGAAAUAGUGUGUGAAGGUGAUCCCUGCGUUGCAGGCCGUUGCAGGGUAGAACCUCAUAAGAUUGAAUUGAAAUCACGAGACUGUCUUCUGUGCCACCGUAGGUUUUCAUCGAACUUGCGUUUAUGCCCGGAAGCGAGAAUGCGCAUUGACACCACGACCCAAUCAGCCUGGCGGUUGUACUGUCUGUUAACGCUGAGUAGGGUACUUGAGAGGUGGGUGACAUCCCGCGAGCGUGUUUAUGUAGGCUACAUGUCGGUCAAGCAGACCACAGAGGGAACUGUAUCGUCUUGCCCGUGUACCGAUCUCCGGGGUAGCAGCCUUUGGAGCUCGGCUGACACCGAGGUCGCGGCUUGUCGGCAAUCGGUACAGACCAAUAGCUUUUUAUACGGACCAUACCGACCGCCUUUCUCGCACGAAAAACGCGCGUGCCAUGGCAAGUCCGUGUCCAAAACGAACCAUUGGCAUGUCGCAAUCUCAAGAUUGACAUUAGCGAUUACACGUAUUAUAUACUCAGUAUGUACUAAUAGACGGAGGGUCGUCAUCCUUGCUCCUGGAUUGGACGAAUUACAUCAAAGUACUAACAAAGGACACAAUCCAAGCAAUGGACACUGCGUCCUUCAGUGUUCGAUGACGUUUCACACGCGUAGACUCACAUAUGGGGCAUUGUAA"))
print(protein_string)


"""
Next Problem - Calculating Protein Mass

"""

def CalculateProteinMass(data):
    data_protein_arr = []
    masses = 0
    for x in data:
        data_protein_arr.append(x)

    ProteinMassTable = {"A":71.03711,"C":103.00919,
                        "D":115.02694,"E":129.04259,
                        "F":147.06841,"G":57.02146,
                        "H":137.05891, "I":113.08406,
                        "K":128.09496, "L":113.08406,
                        "M":131.04049, "N":114.04293,
                        "P":97.05276, "Q":128.05858,
                        "R":156.10111, "S":87.03203,
                        "T":101.04768, "V":99.06841,
                        "W":186.07931, "Y":163.06333,
                        }
    for i in range(0,len(data_protein_arr)):
        if "".join(data_protein_arr[i]) in ProteinMassTable.keys():
            masses += ProteinMassTable["".join(data_protein_arr[i])]
        else:
            pass

    return masses

totalmass = CalculateProteinMass("EMGDMLPASTMKQLDKFQNCLPRFMSKCRVCIATMVCTHHWTMHEKKWWKHWCQDAAKGYKHMRIQDSSFLRGSYGGLFSIEFTHYFWAFDCYEMWWPVAINCCSPNAPTLWMGLTCHVLDYHWTEIKAMDVIELVMRYCGQNSDEVLAIQILGWARSQSQECWEAPRIPVWWGGNHPLWTWNCLSQQDNRLKMFKHRPHSSPFNYEHGMEKCYHGMGFLNTAMFHACHKRTYVCVFFKMNKDCKTDNNCRTIRYWVHVEYVKIWQAHMNEDDIYCETQMAKGKRYNYLTQHIDPTPIGYPCLDMVYTNEMGRWLFRVPYVLCHTWQHDFNQRNIRHLEWNAMHHFQYFQAEQTREHMSSNVHFFYNHKDTSPSVTHHEYAAVCHFCNEYKLSHFDRWTETRFNKMINDPMYWYSAPLTWFNGPLWDNIRTTYLHPRSLQTYMIYDMCLYMRDARCRSHVLWFYTGGGAVNGWECMVLYNAMYEKQSQDSTQWDPCQMKYEELGLQVTSPVMIWGACKHKWRSDPQMWVKDYVRNLKMLFFSQHKVVSTGQPVFQRCKVRKDCQWDEDAPGCVRSFVTCRSIGRLYIRGLDETMCMLQASLDHSFEDTFRGACRIGANPLSPMRGEAMIAPWLGLGMANIACTGNRNRIMKCSYEPAPYWVQSCISTFPKEKAWMTATFVYFCCFFWVQWWRITMLRFMNKIYTVQPAWLSLQSSFWLTLMWSYSDHVDRTITRIDNVSYWQGTPNEQTRMWVCLNWNSVSAVTSFMFHRKHDYISAWKQTTNVSAYDCSYEYNAYLTNCQNQPHQLHHMKDRHYIRMIEITSSANIDSVRMWNTNWKHSLHSQFTHLRIHAWINCKSNFEEMFTNGHKRERCVMCELCLWAVCHNLMCSWPCSIWDFNGWYQPKMMRHGTDEKNYFDKQAALVWRDRAVVCIWALIFKVMCQYWHLMYPDKQELGQTLDPGKNWCEQVAKEF")
print(totalmass)


"""
Next Problem - Finding Motif in DNA
"""

def FindMotif(string,substring):
    arr_string = []
    arr_substring = []
    results_arr = []
    for x in string:
        arr_string.append(x)
    for x_sub in substring:
        arr_substring.append(x_sub)
    for i in range(0, len(arr_string)):
        if "".join(arr_string[i:i+len(substring)]) == substring:
            results_arr.append(i+1)
            print("Yes")
        else:
            print("No")
            pass

    return results_arr


#motifindex = FindMotif("AGGGGGTACAGGGGGTAAACAGGGGGTATACGGGGGTAGGGGGTACTTGGGGGTAGTCCGGGGGTAGGGGGTATAGGGGGTACCGGGGGTATAGGGGGGTAGCGGGGGTAGGGGGTAGGGGGGTAAGGGGGTAGGGGGTAGGGGGTAAGGGGGTACGCCCACGCTCTCGGGGGTAATAATGGGGGTAGGGGGTAACCGGGGGTAGGGGGTATTTTGGGGGTAGGGGGTAGGGGGTATTGTATGGGGGTAGGGGGTAGGGGGTACCGAGGGGGTAGTTTTTCTCCCAGGGGGTAGGGGGTAGGGGGTACGACAAGGCCGGGGGGTATTCGGGGGGTAGTGGGGGTAGGGGGTAGGGGGTAGAGGGGGTAGGGGGTATCGGGGGTAGTATAACAGGGGGGTAATCGGTACATGGGGGTAGGGGGTACTAGGGGGTAGGGGGTAGGGGGTAGGGGGTAAAAGAAGGGGGTAAGGGGGTAGGGGGTAATTGGGGGTACGTTTCGAGTGGGGGTATTGGGGGTAGGGGGTATGGGGGTAAATATGGGGGTAGGGGGTACGGGGGTAAGCGGGGGTATAAAGAGGGGGTACTAAGTGGGGGTAAGGGGGTAGGGGGTACGGGGGGTATGGGGGTAGTTCGGGGGTAAGCAGGGGGTACGGGGGTAGGGGGTAGGGGGTAGGGGGGTAGAGGGGGTAGGGGGGTAACAAGGGGGGTAGGGGGTAAGGGGGGTAGGGGGTAGGGGGTAGGGGGTAGGGGGTACACTGTGAGGGGGGTATCCGGGGGTAGGGGGTAGGGGGGTAGGTTCAGGGGGTATGGGGGTAGGGGGGTATCAGGGGGTAGGGGGGTACGGGGGTAAAAACGGGGGGTAGGTATGGGGGTAAAAGGGGGTAGTGGAGGGGGTAGGGGGTAGGGGGTATCCCTCGACGGGGGTAGGGGGTAACGGGGGTAGGGGGTATGTGGGGGGGTAGGGGGTA","GGGGGTAGG")
#print(motifindex)

"""
Next Problem - Locating Restriction Sites
"""

def ReverseCompliment(data):
    reverse_data_compliment = data.replace("A", "t").replace("G", "c").replace("T", "a").replace("C","g").upper()
    return reverse_data_compliment

def LocateRestrictionSites(data):
    org_arr = []
    rev_arr = []
    locations = {}
    list1 = list()
    list2 = list()
    for x in data:
        org_arr.append(x)
    reverse = ReverseCompliment(data)
    for y in reverse:
        rev_arr.append(y)
    print(org_arr)
    print(rev_arr)
    print(len(data)-1)

    for index in range(0, len(data)):
        y = (len(data) - 1)
        for x_range in range(3,13):
            if index+x_range <= len(data):
                print(index,org_arr[index:index+x_range])
                print(rev_arr[index:index+x_range][::-1])
                if org_arr[index:index+x_range] == rev_arr[index:index+x_range][::-1]:
                    locations[index+1] = x_range
                    list1.append(index+1)
                    list2.append(x_range)
                else:
                    pass

            elif org_arr[y-4:y] == rev_arr[0:4][::-1]:
                list1.append(index + 1)
                list2.append(x_range)


    return locations, list1,list2


x,l1,l2= LocateRestrictionSites("CGCGAGAGGCATGCGTATTATCGCAAAACGTTAAGCCTCAAAAGGTGGCCACATCGTCGA\
ATGACCACACACGAGTGTCCCCGACGCTCCCCCAGCGCAAGGTCGCCGGGCCAGAGGGGA\
CGTGAGGGATCTTTTAAGACAATAACGCATGTGGAACCAAGTTTCGATTGAGGCTAGTAC\
TCAACCAATCGTCGGCGCTAAGGGGACGCATACTGGCGGAAGGAGAGTCCCAGAAGATCG\
AGCGCCAGCCCGTGTCCGATGCTCCGGTTCCTCCAAAAGGAAAAGCGCTCCACTTATAAG\
CACGCACCTATACGGTACTAAGAATAACATTGGGCTGCCCTCGCTTTTCATGTACGGCCA\
AACATATCGGACTCCATCCGTTATGTAAACGTAGTGTCGATAGAGTACGCCTCGGGTTTT\
AAGAAGTCGGCCTATTTATCGCCCTACACAGGCCCTCGATTAATCGACGACGGGGTTTTA\
CAGTCAACGTATCTACCTCAAAGTAGTCCGCGTCTAGTTCCCCCTTGGAAACCCTAGGCC\
GGTTCCTCTATCGGTCATCGATCAATCATAATGATACCAGCCTACGGTCAGTAAGGTTCG\
TAAGTGGGAAGGTGGTTTTTACGTGTGGCGGGCTTCACAACACTCTAAATTCAAAGTGCA\
TACAAACCCTACGGGTGGTTTTCTTTGGCCGTGAGAAGCACATTACCTGGGCAACTGGCC\
GAGCCAGATTAGTAGCGACTATGTAAGCTTACAAACCCCCGGTTGAGACTATTCTGAGAT\
GTCTGCCGACAGGCACAGGCCCTTTTAAGCACGTTATCTGCCAATGGTGTCCCTGGGACA\
GAAGACTAACGGTACCGCATCTCTTCTCGGCAAGATTGATGATCAAACGATCCTATGGAC\
GTGCCTCACAGACATTGGATAT")
print(x)
for m in range(0,len(l1)):
    print(str(l1[m])+"\t"+str(l2[m]))