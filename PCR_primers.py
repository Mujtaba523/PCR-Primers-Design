#Write first any import statements that you need in your program
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

#Define here any functions you will be using
def readin_fasta(filename):
    """Given a *string* with a file name/path of a FASTA file, the function `readin_fasta` returns one string with the DNA sequence (no spaces or gaps) contained in the file."""
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        return sequence

def CG_edge(sequence):
    """Given a string with a DNA sequence as an argument, the function `CG_edge` returns a boolean indicating if the argument string starts and ends with a 'G' or a 'C'."""
    sequence = sequence.upper()
    start_G = sequence.startswith("G")
    start_C = sequence.startswith("C")
    end_G = sequence.endswith("G")
    end_C = sequence.endswith("C")
    if (start_G or start_C) and (end_G or end_C):
        return True
    return False

def CG_count(sequence):
    """Given a string with a DNA sequence as an argument, the function `CG_count` returns an integer with the number of 'C' and 'G' present in the string."""
    sequence = sequence.upper()
    return int(sequence.count("C")+sequence.count("G"))

def ReverseComplement(sequence):
    """Given a string with a DNA sequence as an argument, the function `ReverseComplement` returns a string consisting of a DNA sequence that is complementary of, and 
    in reverse order with respect to the input string. For example, for an input string 'CATGG', the output of the function should be 'CCATG'."""
    seq = Seq(sequence)
    return str(seq.reverse_complement())

def nucleotide_frequency(sequence):
    """Count nucleotides in a given sequence. Return the integer"""
    seq = sequence.upper()
    number = dict(Counter(seq))
    return sum(number.values())

def CG_content_subsec(sub_sequence, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
    res = []
    for i in range(0, len(sub_sequence) - k + 1, k):
        subseq = sub_sequence[i:i + k]
        res.append(
            round(CG_count(subseq) / len(subseq) * 100))
    return res[0]

def complement5(sub_sequence1, sub_sequence2):
    """Given 2 strings as arguments, each containing a DNA sequence, the function `complement5` 
    returns a boolean indicating if there exists a contiguous sub-fragment of length 5 or above in sub_sequence1
    that can pair up with a fragment of sub_sequence2."""
    for i in range(len(sub_sequence1)-4):
        newString = sub_sequence1[i:i+5]
        sr = newString[::-1]
        for i in range(len(sr)):
            if (sr[i] == "A"):
                sr = sr[:i] + "T" + sr[i+1:]
            elif (sr[i] == "T"):
                sr = sr[:i] + "A" + sr[i + 1:]
            elif (sr[i] == "G"):
                sr = sr[:i] + "C" + sr[i + 1:]
            elif (sr[i] == "C"):
                sr = sr[:i] + "G" + sr[i + 1:]
        newString = sr
        res = sub_sequence2.find(newString)
        if(res>0):
            return True
        return False

if __name__=="__main__":
    #Remove the pass statement below and write the rest of the
    #program inside this if block
    #Do not forget indentation
    filename = input("Enter filename: ")                                                                
    position = input("Enter starting and ending postion number with space between them: ")
    position = position.split(" ")
    start = int(position[0])                                                                      #Starting position specified by user
    end = int(position[1])                                                                        #Ending position specified by user
    sequence = readin_fasta(filename)
    if start>end:
       raise Exception("Enter the position correctly")                           # if starting position is greater than ending position, exception is raised hence terminating the execution of programs.
    seq = sequence[start:end]                                                    # Selecting the desired part from the whole sequence
    if start<0 or end<0:                                                         # If start, end are negative then slicing cannot be done properly.
        print("Invalid starting position or ending postions.")                   
    elif len(seq)>=20:                                                           # Condition is each primer must be 20 nucleotides long, therefore the sequence should be greater then 20.
        for i in reversed(range(start)):
        # Checking each condition for start primer sequence.
            start_primer = ReverseComplement(sequence[i:i+20])
            start_primer_substring_1 = start_primer[:len(start_primer)//2]
            start_primer_substring_2 = start_primer[len(start_primer)//2:]
            if CG_edge(start_primer)==True:  
                if nucleotide_frequency(start_primer)==20:                                                        
                    if CG_content_subsec(start_primer)>=40 and CG_content_subsec(start_primer)<=60:
                        if complement5(start_primer_substring_1,start_primer_substring_2) == False:
                            start_primer_sequence = True
                            print("Start primer sequence: " + str(start_primer))
                            break
                        else:
                            start_primer_sequence = ""
                    else:
                            start_primer_sequence = ""
                else:
                            start_primer_sequence = ""
            else:
                            start_primer_sequence = ""
        if start_primer_sequence == "":
            print("Start primer sequence: ")
        if end>len(sequence):
            reverse_primer_sequence = ""
        else:
            counter = 0
            end_point = int(end-20)
            for j in range(end):
            # Checking each condition for reverse primer sequence.
                reverse_primer = sequence[end_point+counter:end+counter]
                reverse_primer_substring_1 = reverse_primer[:len(reverse_primer)//2]
                reverse_primer_substring_2 = reverse_primer[len(reverse_primer)//2:]
                counter+=1
                if CG_edge(reverse_primer)==True:
                    if nucleotide_frequency(reverse_primer)==20:
                        if CG_content_subsec(reverse_primer)>=40 and CG_content_subsec(reverse_primer)<=60:
                            if complement5(reverse_primer_substring_1,reverse_primer_substring_2) == False:
                                reverse_primer_sequence = True
                                print("Reverse primer sequence: " + str(reverse_primer))
                                break  
                            else:
                                reverse_primer_sequence = ""
                        else:
                                reverse_primer_sequence = ""
                    else:
                                reverse_primer_sequence = ""
                else:
                                reverse_primer_sequence = ""
                
                if (start_primer_sequence == True) and (reverse_primer_sequence==True):
                    print("") 
                    break
        if reverse_primer_sequence == "":
            print("Reverse primer sequence: ")            
    else:
        print("The sequence given have short length, it should be greater than 20 nucleotides.")

    
