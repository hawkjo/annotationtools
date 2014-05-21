"""
GO_tools.py

A set of simple tools for working with the Gene Ontology.
"""
import sys
from collections import defaultdict

__author__      = 'John Hawkins'
__copyright__   = 'Copyright 2014 John Hawkins'
__license__     = 'Apache 2.0'

def processGOTerm(goTerm):
    """
    In an object representing a GO term, replace id, name, def lists as strings.
    Returns the modified object as a dictionary.
    """
    #-----------------------------------------------------------------------------------------------
    # This function is taken, with modification, from Uli Koehler, Copyright 2013.
    # http://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/
    #-----------------------------------------------------------------------------------------------
    ret = dict(goTerm) #Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.iteritems():
        if key in ['id','name','def','namespace']:
            ret[key] = value[0]
    return ret

def parseGOOBO(filename):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each 
    Keyword arguments:
        filename: The filename to read
    """
    #-----------------------------------------------------------------------------------------------
    # This function is taken, with modification, from Uli Koehler, Copyright 2013.
    # http://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/
    #-----------------------------------------------------------------------------------------------
    with open(filename, "r") as infile:
        currentGOTerm = None
        for line in infile:
            line = line.strip()
            if not line: continue #Skip empty
            if line == "[Term]":
                if currentGOTerm: yield processGOTerm(currentGOTerm)
                currentGOTerm = defaultdict(list)
            elif line == "[Typedef]":
                #Skip [Typedef sections]
                if currentGOTerm: yield processGOTerm(currentGOTerm)
                currentGOTerm = None
            else: #Not [Term]
                #Only process if we're inside a [Term] environment
                if currentGOTerm is None: continue
                key, sep, val = line.partition(":")
                currentGOTerm[key].append(val.partition('!')[0].strip())
        #Add last term
        if currentGOTerm is not None:
            yield processGOTerm(currentGOTerm)

def GO_OBO_to_dict_with_id_as_key(filename, add_alt_id_keys=True):
    """
    Returns dictionary where keys are GO ids and values are dictionaries containing all the
    properties of the given GO term.
    Keyword arguments:
        filename:           Filename of GO OBO file.
        add_alt_id_keys:    (Default True) Bool to create dict entries using alt_ids as keys.
    """
    GO_dict = {}
    for goTerm in parseGOOBO(filename):
        GO_dict[goTerm['id']] = goTerm
        if add_alt_id_keys and 'alt_id' in goTerm: 
            for alt_id in goTerm['alt_id']:
                 GO_dict[alt_id] = goTerm
    return GO_dict

def slim_id_set_given_go_id_dict(GO_dict, slim_terms_set, include_part_of=True, include_all_relationships=False):
    """
    Returns a dict which takes as keys all possible GO ids, including alt_ids, and returns as values
    the corresponding list of equidistant parents to slim to.
    Keyword arguments:
        GO_dict:                    Dict created with GO_OBO_to_dict_with_id_as_key.
        slim_terms_set:             Set of GO ids of all key terms to slim to.
        include_part_of:            (Default True) Bool to use 'part_of' relationships for slimming.
        include_all_relationships:  (Default False) Bool to use all relationships for slimming.
    """
    slim_info_given_go_id = {}
    for go_id in GO_dict.keys():
        slim_info_given_go_id[go_id] = {'slims_to_set':set(), 'dist':None}

    def set_slim_info(go_id):
        # Recursive function to dynamically find closest slim key terms
        if slim_info_given_go_id[go_id]['slims_to_set'] and slim_info_given_go_id[go_id]['dist'] is not None: 
            # Check if term has already been taken care of
            return

        elif go_id in slim_terms_set or GO_dict[go_id]['name'] in \
                ['biological_process', 'molecular_function', 'cellular_component']:
            # Check if go_id is one of the key words
            slim_info_given_go_id[go_id]['slims_to_set'] = { GO_dict[go_id]['id'] }
            slim_info_given_go_id[go_id]['dist'] = 0

        else:
            # Create parent list
            # Start with all 'is_a' relationships
            if 'is_a' not in GO_dict[go_id]:
                if 'is_obsolete' not in GO_dict[go_id]:
                    sys.exit('Non-obsolete non-root term with no parents: ' + GO_dict[go_id]['name'])
                return
            else:
                parent_list = GO_dict[go_id]['is_a']
            # Add 'part_of' and/or all relationships
            if include_part_of or include_all_relationships:
                if 'relationship' in GO_dict[go_id]:
                    for relationship in GO_dict[go_id]['relationship']:
                        rel_type, rel_id = relationship.split()
                        if rel_type == 'part_of' or include_all_relationships:
                            parent_list.append( rel_id )

            for parent in parent_list:
                set_slim_info(parent)
            slim_info_given_go_id[go_id]['dist'] = \
                    min(slim_info_given_go_id[parent]['dist'] + 1 for parent in parent_list)
            for parent in parent_list:
                if slim_info_given_go_id[go_id]['dist'] == slim_info_given_go_id[parent]['dist'] + 1:
                    slim_info_given_go_id[go_id]['slims_to_set'].update(
                            slim_info_given_go_id[parent]['slims_to_set'] )

    for go_id in GO_dict.keys():
        set_slim_info(go_id)

    slim_id_set_given_go_id = {}
    for go_id in GO_dict.keys():
        slim_id_set_given_go_id[go_id] = slim_info_given_go_id[go_id]['slims_to_set']

    return slim_id_set_given_go_id
    
def map2slim(GO_filename, slim_filename, gaf_filename, out_filename, include_part_of=True,
        include_all_relationships=False):
    """
    Takes a gaf file and replaces GO ids with slimmed GO ids from slim file.
    Where GO id slims to more than one slimmed GO id, a line for each GO id is output.
    Duplicate lines are removed.
    Keyword arguments:
        GO_filename:                Filename of GO OBO file.
        slim_filename:              Filename of slim OBO file.
        gaf_filename:               Filename for input gaf file.
        out_filename:               Filename for output gaf file.
        include_part_of:            (Default True) Bool to use 'part_of' relationships for slimming.
        include_all_relationships:  (Default False) Bool to use all relationships for slimming.
    """
    GO_dict = GO_OBO_to_dict_with_id_as_key(GO_filename)
    slim_terms_set = set( goTerm['id'] for goTerm in parseGOOBO(slim_filename) )
    slim_id_set_given_go_id = slim_id_set_given_go_id_dict(GO_dict, slim_terms_set, include_part_of,
            include_all_relationships)

    # First collect all desired output lines in a set to de-duplicate slimmed results.
    lines_to_write = set()
    with open(gaf_filename) as f:
        # Check first line for header lines.
        headerlines = []
        line = f.readline()
        while line[0] == '!': 
            headerlines.append( line )
            line = f.readline() 
        while line:
            var = line.strip().split('\t')
            for slim_id in slim_id_set_given_go_id[ var[4] ]:
                lines_to_write.add( '\t'.join(var[:4] + [ slim_id ] + var[5:]) + '\n' )
            line = f.readline()

    with open(out_filename,'w') as outfile:
        for line in headerlines:
            outfile.write(line)
        for line in sorted(lines_to_write):
            # First two columns are db and db id, so regular string sorting results in desired output
            outfile.write(line)

if __name__ == "__main__":
    GO_filename = '/home/hawkjo/scratch/dbs/GO/gene_ontology_ext.obo'
    slim_filename = '/home/hawkjo/scratch/dbs/GO/goslim_generic.obo'

    GO_dict = GO_OBO_to_dict_with_id_as_key(GO_filename)
    slim_terms_set = set( goTerm['id'] for goTerm in parseGOOBO(slim_filename) )

    print 'GO OBO file:',GO_filename
    print 'slim OBO file:',slim_filename
    print 'Number of slim terms:',len(slim_terms_set)

    slim_id_set_given_go_id = slim_id_set_given_go_id_dict(GO_dict, slim_terms_set,
            include_part_of=True, include_all_relationships=False)

    print 'Len of slim-go set:',len(slim_id_set_given_go_id)

    from collections import Counter
    slim_id_list_hist = Counter(len(lst) for lst in slim_id_set_given_go_id.values() )
    print 'slim list length hist:'
    for k,v in sorted(slim_id_list_hist.items()):
        print '\t',k,v

    print 'Example key-value pairs:'
    maps2self = 0
    for i,(k,v) in enumerate(slim_id_set_given_go_id.items()):
        if i <= 5:
            print '\t',k,v
        if len(v) == 5:
            print '\t',k,v
        if k in v: maps2self += 1

    print 'Ids which map to self:',maps2self

    print 'Terms without parents:', sum(1 for go_id in GO_dict.keys() if 'is_a' not in GO_dict[go_id])

    print 'Root ids:'
    for go_id in GO_dict.keys():
        if 'is_a' not in GO_dict[go_id] and len(slim_id_set_given_go_id[go_id]) != 0:
            print '\t',go_id,':',GO_dict[go_id]['name']
