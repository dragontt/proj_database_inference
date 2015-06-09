import sys
import os
import argparse
import numpy as nmp

"""
I am a docstring.  I am your friend.
"""

def parse_args(argv):
    parser = argparse.ArgumentParser(description="this prints before the usage string")
    # examples with default values and types, and short and long args
    parser.add_argument('input_file', metavar='input_file',
                        type=str, help='if specified as \'-\', takes input from standard input')
    parser.add_argument('-n', '--name', type=str, dest="name", default=None,
                        help="the TF name of the pfm, used if passed in stdin")
    parser.add_argument('--num_processes', dest="num_processes",
                        type=int, default=6)
    parsed = parser.parse_args(argv[1:])
    return parsed
#end function

def parse_scer_pfm(filename, handle_passed = False):
    '''
    Given the filename of a PFM file from ScerTF, parse it into a nmp array.
    The returned array represents the PFM where each row is a base (A, C, G, T),
    and each column value is the base's relative frequency in that position in
    the motif.
    '''
    retval = None
    if handle_passed:
        reader = filename
    else:
        reader = open(filename)

    for n, line in enumerate(reader):
        line = line.strip()
        base = line.split('|')[0].strip()
        # print line.split('|')[1].strip().split()
        frequencies = [float(x.strip()) for x in line.split('|')[1].strip().split()]
        # each row holds the frequencies for one base
        # order is ACGT
        # 
        if retval == None:
            retval = nmp.zeros((4, len(frequencies)))
        retval[n, :] = nmp.array(frequencies)
    return retval
#end function

def fire_motif_to_meme(fire_motif_str, name):
    ''' A light convenience function to convert a fire motif str to meme format. '''
    return pfm_as_meme_str(fire_motif_to_pfm(fire_motif_str), name)

def fire_motif_to_pfm(fire_motif_str):
    ''' Given a string representing a FIRE motif, transform it into a PFM where
    each row holds the freqs for a base, and each column holds the freqs for 
    one position. '''

    # initialize an empty pfm
    pfm = nmp.zeros((4, 0))
    waiting_for_close = False
    cur_freqs = nmp.zeros((4, 1))
    pos_mapper = {'A':0, 'C':1, 'G':2, 'T':3}

    # read through the characters of the motif string
    for char in fire_motif_str:
        if char == '[':
            waiting_for_close = True
        elif char == '.':
            cur_freqs = nmp.ones((4,1))
        elif char == ']':
            if not waiting_for_close:
                raise ValueError('] present without preceding [')
            waiting_for_close = False
        else:
            cur_freqs[pos_mapper[char]] += 1
        # close out this loop, reset/account for variables
        if not waiting_for_close:
            cur_freqs /= nmp.sum(cur_freqs)
            pfm = nmp.hstack([pfm, cur_freqs])
            cur_freqs[:] = 0
    return pfm
    

def pfm_as_meme_str(pfm, pwm_name):
    outstr_lines = []
    outstr_lines += ["MEME version 4", ""]
    outstr_lines += ["ALPHABET= ACGT", ""]
    outstr_lines += ["strands: + -", ""]
    outstr_lines += ["Background letter frequencies", 
                     "A 0.30 C 0.20 G 0.20 T 0.30", ""]
    outstr_lines += ["MOTIF %s"%(pwm_name)]
    outstr_lines += ["letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0" % (int(nmp.shape(pfm)[1]))]

    # each row in meme format is one base position, with the cols being ACGT freqs
    for ix in range(nmp.shape(pfm)[1]):
        stb = map(str, pfm[:, ix])
        st = ' '.join(stb)
        outstr_lines += [st]
    return "\n".join(outstr_lines)
#end function

def main(argv):
    """ The main module should take in inputs from the command line, 
    carry out the 'CLI' functionality of the script, then write the
    results to stdout. """
    parsed = parse_args(argv)
    instream = sys.stdin
    name = parsed.name
    if parsed.input_file != "-":
        instream = open(parsed.input_file, 'r')
        name = parsed.input_file.split('.')[1]
    print pfm_as_meme_str(parse_scer_pfm(instream, handle_passed=True), name)
#end function

if __name__ == "__main__":
    main(sys.argv)
