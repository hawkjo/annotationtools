"""
GO_tools.py
Version 0.0
Author: John Hawkins

A set of simple tools for working with the Gene Ontology.
"""
import sys
from simple_go_obo_parser import parseGOOBO

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
        if add_alt_id_keys:
            if 'alt_id' not in goTerm: continue
            elif type(goTerm['alt_id']) is str:
                GO_dict[goTerm['alt_id']] = goTerm
            elif type(goTerm['alt_id']) is list:
                for alt_id in goTerm['alt_id']:
                    GO_dict[alt_id] = goTerm
            else: sys.exit('Type error of alt_id in ' + goTerm['id'])
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
        slim_info_given_go_id[go_id] = {'slims_to_list':set(), 'dist':None}

    def set_dist(go_id):
        if slim_info_given_go_id[go_id]['slims_to_list'] and slim_info_given_go_id[go_id]['dist'] is not None: 
            # Check if term has already been taken care of
            return

        elif go_id in slim_terms_set or GO_dict[go_id]['name'] in \
                ['biological_process', 'molecular_function', 'cellular_component']:
            # Check if go_id is one of the key words
            slim_info_given_go_id[go_id]['slims_to_list'] = { GO_dict[go_id]['id'] }
            slim_info_given_go_id[go_id]['dist'] = 0

        else:
            # Create parent list
            # Start with all 'is_a' relationships
            if 'is_a' not in GO_dict[go_id]:
                if 'is_obsolete' not in GO_dict[go_id]:
                    sys.exit('Non-obsolete non-root term with no parents: ' + GO_dict[go_id]['name'])
                return
            elif type(GO_dict[go_id]['is_a']) is str:
                parent_list = [ GO_dict[go_id]['is_a'] ]
            elif type(GO_dict[go_id]['is_a']) is list:
                parent_list = GO_dict[go_id]['is_a']
            else:
                sys.exit('Type error of is_a: ' + go_id)
            # Add 'part_of' and/or all relationships
            if include_part_of or include_all_relationships:
                if 'relationship' in GO_dict[go_id]:
                    if type(GO_dict[go_id]['relationship']) == str:
                        relationship_list = [ GO_dict[go_id]['relationship'] ]
                    if type(GO_dict[go_id]['relationship']) == list:
                        relationship_list = GO_dict[go_id]['relationship']
                    for relationship in relationship_list:
                        rel_type, rel_id = relationship.split()
                        if rel_type == 'part_of' or include_all_relationships:
                            parent_list.append( rel_id )

            for parent in parent_list:
                set_dist(parent)
            slim_info_given_go_id[go_id]['dist'] = \
                    min(slim_info_given_go_id[parent]['dist'] + 1 for parent in parent_list)
            for parent in parent_list:
                if slim_info_given_go_id[go_id]['dist'] == slim_info_given_go_id[parent]['dist'] + 1:
                    slim_info_given_go_id[go_id]['slims_to_list'].update(
                            slim_info_given_go_id[parent]['slims_to_list'] )

    for go_id in GO_dict.keys():
        set_dist(go_id)

    slim_id_set_given_go_id = {}
    for go_id in GO_dict.keys():
        slim_id_set_given_go_id[go_id] = slim_info_given_go_id[go_id]['slims_to_list']

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
    slim_terms_set = set( GO_OBO_to_dict_with_id_as_key(slim_filename, add_alt_id_keys=False).keys() )
    slim_id_set_given_go_id = slim_id_set_given_go_id_dict(GO_dict, slim_terms_set)

    with open(out_filename,'w') as outfile:
        with open(gaf_filename) as f:
            outfile.write( f.readline() )
            for line in f:
                var = line.strip().split('\t')
                try:
                    for slim_id in slim_id_set_given_go_id[ var[4] ]:
                        outfile.write( '\t'.join(var[:4] + [ slim_id ] + var[5:]) + '\n' )
                except IndexError:
                    print var
                    sys.exit()

if __name__ == "__main__":
    GO_filename = '/home/hawkjo/scratch/dbs/GO/gene_ontology_ext.obo'
    slim_filename = '/home/hawkjo/scratch/dbs/GO/goslim_generic.obo'

    GO_dict = GO_OBO_to_dict_with_id_as_key(GO_filename)
    slim_terms_set = set( GO_OBO_to_dict_with_id_as_key(slim_filename, add_alt_id_keys=False).keys() )
    print 'Number of slim terms:',len(slim_terms_set)

    slim_id_set_given_go_id = slim_id_set_given_go_id_dict(GO_dict, slim_terms_set)

    print 'Len of slim-go dict:',len(slim_id_set_given_go_id)

    from collections import Counter
    slim_id_list_hist = Counter()
    for lst in slim_id_set_given_go_id.values():
        slim_id_list_hist[ len(lst) ] += 1
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

    no_defs = 0
    def_str = 0
    def_list = 0
    for go_id in GO_dict.keys():
        if 'def' not in GO_dict[go_id]:
            no_defs += 1
        elif type(GO_dict[go_id]['def']) is str:
            def_str += 1
        elif type(GO_dict[go_id]['def']) is list:
            def_list += 1
    print 'Def summary (None, One, More): %d, %d, %d' % (no_defs, def_str, def_list)
