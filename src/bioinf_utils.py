import csv
import sys
import os


"""
A non-runnable module with utils for doing bioinformatics-y stuff.

Provided utilities:
Yeast:
Functions:
- translate between sgd names and systematic names
"""

### Locations of data sources to draw on
database_cross_reference = "data/dbxref.tab"
sgd_to_systematic = "data/sgd-names-systematic.tab"

def sgd_translate_names(names):
    ''' Given a list of gene names, returns a list of their systematic names. '''
    reader = csv.reader(open(sgd_to_systematic, 'r'), delimiter='\t')
    translator = dict([(line[0], line[1]) for line in reader])
    retval = []
    for n in names:
        try:
            retval.append(translator[n])
        except KeyError:
            sys.stderr.write(n+": translation not found, using original name\n")
            retval.append(n)
    return retval
# end function


if __name__ == "__main__":
    raise NotImplemented("This module can't be run from the command line.")
