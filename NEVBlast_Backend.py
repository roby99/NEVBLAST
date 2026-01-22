from Bio.Blast import NCBIXML
import csv
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

def blastparser(slist, file):
    import ast
    seqlists = []
    sbjctqueryLists = []
    
    # Use ast.literal_eval to properly parse the signature string
    try:
        parsed_sigs = ast.literal_eval(slist)
    except:
        # Fallback to old parsing if format is different
        trim = slist[1: len(slist) - 1]
        sigs = trim.split("], [")
        parsed_sigs = []
        for a in sigs:
            aalist = []
            templist = a.split(", ")
            for b in templist:
                blist = []
                blist.append(b[1:])
                blist.append(b[0])
                aalist.append(blist)
            parsed_sigs.append(aalist)
    
    seqlists = parsed_sigs



    # Opens XML storage file and parses it resulting in an output
    result_handle = open(file + ".xml")
    rawTextFile = open(r"" + file + ".txt", "w")
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:

            # Begins parsing through signature sequences
            for y in seqlists:
                Tempscorelist = []
                for x in y:
                    # Clean the position string of any non-digit characters
                    position_str = ''.join(c for c in x[0] if c.isdigit())
                    r = int(position_str) if position_str else 0
                    downcount = r - hsp.query_start + 1
                    e = 0
                    while e < downcount:
                        if hsp.query[e:e + 1] == "-":
                            downcount = downcount + 1
                        e = e + 1

                    # Case 1: The signature is before the start of the sequence
                    if downcount <= 0:
                        Tempscorelist.append("0,-")
                        Tempscorelist.append("0," + x[1])

                    # Case 2: The signature is within the sequence
                    else:
                        qh = str(downcount - hsp.query[:downcount].count("-"))
                        sh = str(downcount - hsp.sbjct[:downcount].count("-"))
                        Tempscorelist.append(str(sh) + "," + str(hsp.sbjct[downcount - 1:downcount]))
                        Tempscorelist.append(str(qh) + "," + x[1]) #needed for scoring - cannot comment out
                sbjctqueryLists.append(Tempscorelist)

            # Writing to the text file
            rawTextFile.write(alignment.title + ' ~ ')
            rawTextFile.write(str(hsp.sbjct) + ' ~ ')
            rawTextFile.write(str(hsp.expect) + ' ~ ')
            for a in sbjctqueryLists:
                rawTextFile.write(str(a) + ' ~ ')
            rawTextFile.write("\n")
            sbjctqueryLists = []
    rawTextFile.close()


# Hash takes in the text file output from blastParser. The text file has
# outputs of [Name~Sequence~E-Value~X expected signature~X actual signature].
# The X is a variable and the number of expected and actual signatures vary
# by what was given by the user. Each line is parsed and scored against the
# Blosum62 matrix for a score that is output to a CSV file.

def hash(file):
    blosumHash = [
        [['C'], ['-', '-5'], ['C', '9'], ['S', '-2'], ['T', '-1'], ['P', '-4'], ['A', '-1'], ['G', '-4'], ['N', '-3'],
         ['D', '-4'], ['E', '-5'], ['Q', '-4'], ['H', '-4'], ['R', '-4'], ['K', '-4'], ['M', '-2'], ['I', '-2'],
         ['L', '-2'], ['V', '-1'], ['F', '-3'], ['W', '-3'], ['Y', '-3']],
        [['S'], ['-', '-5'], ['C', '-2'], ['S', '5'], ['T', '1'], ['P', '-1'], ['A', '1'], ['G', '-1'], ['N', '0'],
         ['D', '-1'], ['E', '0'], ['Q', '0'], ['H', '-1'], ['R', '-1'], ['K', '-1'], ['M', '-2'], ['I', '-3'],
         ['L', '-3'],
         ['V', '-2'], ['F', '-3'], ['W', '-4'], ['Y', '-2']],
        [['T'], ['-', '-5'], ['C', '-1'], ['S', '1'], ['T', '5'], ['P', '-2'], ['A', '0'], ['G', '-2'], ['N', '0'],
         ['D', '-1'], ['E', '-1'], ['Q', '-1'], ['H', '-2'], ['R', '-1'], ['K', '-1'], ['M', '-1'], ['I', '-1'],
         ['L', '-2'], ['V', '0'], ['F', '-2'], ['W', '-4'], ['Y', '-2']],
        [['P'], ['-', '-5'], ['C', '-4'], ['S', '-1'], ['T', '-2'], ['P', '8'], ['A', '-1'], ['G', '-3'], ['N', '-3'],
         ['D', '-2'], ['E', '-2'], ['Q', '-2'], ['H', '-3'], ['R', '-2'], ['K', '-1'], ['M', '-3'], ['I', '-4'],
         ['L', '-3'], ['V', '-3'], ['F', '-4'], ['W', '-5'], ['Y', '-4']],
        [['A'], ['-', '-5'], ['C', '-1'], ['S', '1'], ['T', '0'], ['P', '-1'], ['A', '5'], ['G', '0'], ['N', '-2'],
         ['D', '-2'], ['E', '-1'], ['Q', '-1'], ['H', '-2'], ['R', '-2'], ['K', '-1'], ['M', '-1'], ['I', '-2'],
         ['L', '-2'], ['V', '0'], ['F', '-3'], ['W', '-3'], ['Y', '-2']],
        [['G'], ['-', '-5'], ['C', '-4'], ['S', '-1'], ['T', '-2'], ['P', '-3'], ['A', '0'], ['G', '6'], ['N', '-1'],
         ['D', '-2'], ['E', '-3'], ['Q', '-2'], ['H', '-3'], ['R', '-3'], ['K', '-2'], ['M', '-4'], ['I', '-5'],
         ['L', '-4'], ['V', '-4'], ['F', '-4'], ['W', '-4'], ['Y', '-4']],
        [['N'], ['-', '-5'], ['C', '-3'], ['S', '0'], ['T', '0'], ['P', '-3'], ['A', '-2'], ['G', '-1'], ['N', '6'],
         ['D', '1'], ['E', '-1'], ['Q', '0'], ['H', '0'], ['R', '-1'], ['K', '0'], ['M', '-3'], ['I', '-4'],
         ['L', '-4'],
         ['V', '-4'], ['F', '-4'], ['W', '-4'], ['Y', '-3']],
        [['D'], ['-', '-5'], ['C', '-4'], ['S', '-1'], ['T', '-1'], ['P', '-2'], ['A', '-2'], ['G', '-2'], ['N', '1'],
         ['D', '6'], ['E', '1'], ['Q', '-1'], ['H', '-2'], ['R', '-2'], ['K', '-1'], ['M', '-4'], ['I', '-4'],
         ['L', '-5'],
         ['V', '-4'], ['F', '-4'], ['W', '-6'], ['Y', '-4']],
        [['E'], ['-', '-5'], ['C', '-5'], ['S', '0'], ['T', '-1'], ['P', '-2'], ['A', '-1'], ['G', '-3'], ['N', '-1'],
         ['D', '1'], ['E', '6'], ['Q', '2'], ['H', '0'], ['R', '-1'], ['K', '1'], ['M', '-2'], ['I', '-4'], ['L', '-4'],
         ['V', '-3'], ['F', '-4'], ['W', '-4'], ['Y', '-3']],
        [['Q'], ['-', '-5'], ['C', '-4'], ['S', '0'], ['T', '-1'], ['P', '-2'], ['A', '-1'], ['G', '-2'], ['N', '0'],
         ['D', '-1'], ['E', '2'], ['Q', '6'], ['H', '1'], ['R', '1'], ['K', '1'], ['M', '0'], ['I', '-3'], ['L', '-3'],
         ['V', '-3'], ['F', '-4'], ['W', '-3'], ['Y', '-2']],
        [['H'], ['-', '-5'], ['C', '-4'], ['S', '-1'], ['T', '-2'], ['P', '-3'], ['A', '-2'], ['G', '-3'], ['N', '0'],
         ['D', '-2'], ['E', '0'], ['Q', '1'], ['H', '8'], ['R', '0'], ['K', '-1'], ['M', '-2'], ['I', '-4'],
         ['L', '-3'],
         ['V', '-4'], ['F', '-2'], ['W', '-3'], ['Y', '2']],
        [['R'], ['-', '-5'], ['C', '-4'], ['S', '-1'], ['T', '-1'], ['P', '-2'], ['A', '-2'], ['G', '-3'], ['N', '-1'],
         ['D', '-2'], ['E', '-1'], ['Q', '1'], ['H', '0'], ['R', '6'], ['K', '2'], ['M', '-2'], ['I', '-3'],
         ['L', '-3'],
         ['V', '-3'], ['F', '-4'], ['W', '-4'], ['Y', '-3']],
        [['K'], ['-', '-5'], ['C', '-4'], ['S', '-1'], ['T', '-1'], ['P', '-1'], ['A', '-1'], ['G', '-2'], ['N', '0'],
         ['D', '-1'], ['E', '1'], ['Q', '1'], ['H', '-1'], ['R', '2'], ['K', '5'], ['M', '-2'], ['I', '-3'],
         ['L', '-3'],
         ['V', '-3'], ['F', '-4'], ['W', '-4'], ['Y', '-3']],
        [['M'], ['-', '-5'], ['C', '-2'], ['S', '-2'], ['T', '-1'], ['P', '-3'], ['A', '-1'], ['G', '-4'], ['N', '-3'],
         ['D', '-4'], ['E', '-2'], ['Q', '0'], ['H', '-2'], ['R', '-2'], ['K', '-2'], ['M', '6'], ['I', '1'],
         ['L', '2'],
         ['V', '1'], ['F', '0'], ['W', '-2'], ['Y', '-2']],
        [['I'], ['-', '-5'], ['C', '-2'], ['S', '-3'], ['T', '-1'], ['P', '-4'], ['A', '-2'], ['G', '-5'], ['N', '-4'],
         ['D', '-4'], ['E', '-4'], ['Q', '-3'], ['H', '-4'], ['R', '-3'], ['K', '-3'], ['M', '1'], ['I', '5'],
         ['L', '1'],
         ['V', '3'], ['F', '-1'], ['W', '-3'], ['Y', '-2']],
        [['L'], ['-', '-5'], ['C', '-2'], ['S', '-3'], ['T', '-2'], ['P', '-3'], ['A', '-2'], ['G', '-4'], ['N', '-4'],
         ['D', '-5'], ['E', '-4'], ['Q', '-3'], ['H', '-3'], ['R', '-3'], ['K', '-3'], ['M', '2'], ['I', '1'],
         ['L', '4'],
         ['V', '1'], ['F', '0'], ['W', '-2'], ['Y', '-2']],
        [['V'], ['-', '-5'], ['C', '-1'], ['S', '-2'], ['T', '0'], ['P', '-3'], ['A', '0'], ['G', '-4'], ['N', '-4'],
         ['D', '-4'], ['E', '-3'], ['Q', '-3'], ['H', '-4'], ['R', '-3'], ['K', '-3'], ['M', '1'], ['I', '3'],
         ['L', '1'],
         ['V', '4'], ['F', '-1'], ['W', '-3'], ['Y', '-2']],
        [['F'], ['-', '-5'], ['C', '-3'], ['S', '-3'], ['T', '-2'], ['P', '-4'], ['A', '-3'], ['G', '-4'], ['N', '-4'],
         ['D', '-4'], ['E', '-4'], ['Q', '-4'], ['H', '-2'], ['R', '-4'], ['K', '-4'], ['M', '0'], ['I', '-1'],
         ['L', '0'],
         ['V', '-1'], ['F', '6'], ['W', '0'], ['Y', '3']],
        [['W'], ['-', '-5'], ['C', '-3'], ['S', '-4'], ['T', '-4'], ['P', '-5'], ['A', '-3'], ['G', '-4'], ['N', '-4'],
         ['D', '-6'], ['E', '-4'], ['Q', '3'], ['H', '3'], ['R', '-4'], ['K', '-4'], ['M', '-2'], ['I', '-3'],
         ['L', '-2'],
         ['V', '-3'], ['F', '0'], ['W', '11'], ['Y', '2']],
        [['Y'], ['-', '-5'], ['C', '-3'], ['S', '-2'], ['T', '-2'], ['P', '-4'], ['A', '-2'], ['G', '-4'], ['N', '-3'],
         ['D', '-4'], ['E', '-3'], ['Q', '-2'], ['H', '2'], ['R', '-3'], ['K', '-3'], ['M', '-2'], ['I', '-2'],
         ['L', '-2'],
         ['V', '-2'], ['F', '3'], ['W', '2'], ['Y', '7']]]
    
    blosum_dict = {}
    for row in blosumHash:
        aa1 = row[0][0]  # First amino acid
        for pair in row[1:]:  # All the scoring pairs
            aa2 = pair[0]  # Second amino acid
            score = int(pair[1])  # Score value
            blosum_dict[(aa1, aa2)] = score

    file_open = open(r"" + file + ".txt", "r")
    hitArray = []
    maxVal = 0

    # This loop parses the file input placing the values in a 2-D list.
    for line in file_open:
        lineList = line.split("~")
        lineList[-1] = lineList[-1].strip("\n")
        hitArray.append(lineList)

    # This loop parses across each hit in the 2-D list for scoring.
    for a in hitArray:
        num = 0
        loops = int((len(a) - 4))
        for x in range(3, 3 + loops):
            sigsList = a[x].split("', '")
            sigsList[0] = sigsList[0][3:]
            sigsList[-1] = sigsList[-1][:-3]

            for l in range(0, len(sigsList) - 1, 2):
                letter1 = sigsList[l + 1]
                letter2 = sigsList[l]
                # Get maximum score using dictionary
                maxVal += blosum_dict.get((letter1[-1], letter1[-1]), 0)
                # Get actual score using dictionary
                num += blosum_dict.get((letter1[-1], letter2[-1]), 0)

            # Normalizes the score before adding it to the list and resetting the variables.
            if maxVal != 0:
                a.append(num / maxVal)
            else:
                a.append(0.0)
            
            num = 0
            maxVal = 0


    # Places the lists into the CSV output file
    with open(r"" + file + ".csv", "w") as out_file:
        nevBlast_writer = csv.writer(out_file)
        nevBlast_writer.writerows(hitArray)
