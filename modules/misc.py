#!/usr/bin/python
import re
import string
import libsbml
import numpy
import scipy
import scipy.linalg
import scipy.optimize
import random
import copy
import math
try:
    from . import SBtab
except:
    import SBtab


def table_type(sbtab):
    '''
    determines table_type of SBtab file
    '''
    for row in sbtab.split('\n'):
        if row.startswith('!!'):
            try:
                row = row.replace('"', "'")
                tabletype = re.search("TableType='([^']*)'", row).group(1)
                return tabletype
            except: pass
    return False


def count_tabs(sbtab_string):
    '''
    count how many SBtabs are in in this string
    '''
    counter = 0
    for row in sbtab_string:
        if row.startswith('!!SBtab'):
            counter += 1
    return counter


    

def validate_file_extension(file_name, file_type):
    '''
    returns Boolean to evaluate if the file has the correct extension:
    sbml => xml
    sbtab => tsv
    '''
    # check extension for sbml file
    if file_type == 'sbml' and file_name[-3:] == 'xml': return True
    elif file_type == 'sbml': return False
    else: pass

    # check extension for sbtab file
    if file_type == 'sbtab' and file_name[-3:] == 'tsv': return True
    elif file_type == 'sbtab': return False
    else: pass

    # if something is completely off, return False
    return False


def check_delimiter(sbtab_file):
    '''
    determine the delimiter of the tabular file
    '''
    sep = False

    for row in sbtab_file.split('\n'):
        if row.startswith('!!'): continue
        if row.startswith('!'):
            s = re.search('(.)(!)', row[1:])
            # if there is only 1 column, we have to define a default separator
            # let's use a tab.
            try: sep = s.group(1)
            except: sep = '\t'

    return sep


def valid_prior(sbtab_prior):
    '''
    if the given SBtab file is a prior for parameter balancing, it needs to be
    checked thorougly for the validity of several features
    '''
    validity = []

    # check table type
    if sbtab_prior.table_type != 'Quantity':
        validity.append('Error: The TableType of the prior file is not '\
                        'correct: %s. '\
                        'It should be Quantity' % sbtab_prior.table_type)

    # check for required columns
    required_columns = ['!QuantityType', '!Unit', '!MathematicalType',
                        '!PriorMedian', '!PriorStd', '!PriorGeometricStd',
                        '!DataStd', '!Dependence',
                        '!UseAsPriorInformation', '!MatrixInfo']
    for column in required_columns:
        if column not in sbtab_prior.columns_dict:
            validity.append('Error: The crucial column %s is missing in'\
                            ' the prior file.' % column)

    # check for required row entries
    required_rows = ['standard chemical potential',
                     'catalytic rate constant geometric mean',
                     'concentration', 'concentration of enzyme',
                     'Michaelis constant', 'inhibitory constant',
                     'activation constant', 'chemical potential',
                     'product catalytic rate constant',
                     'substrate catalytic rate constant',
                     'equilibrium constant', 'forward maximal velocity',
                     'reverse maximal velocity', 'reaction affinity']
    for row in sbtab_prior.value_rows:
        try: required_rows.remove(row[sbtab_prior.columns_dict['!QuantityType']])
        except: pass

    for row in required_rows:
        validity.append('Error: The prior file is missing an entry for th'\
                        'crucial value %s.' % row)

    return validity


def extract_pseudos_priors(sbtab_prior):
    '''
    extracts the priors and pseudos of a given SBtab prior table
    '''
    pseudo_list = ['chemical potential', 'product catalytic rate constant',
                   'substrate catalytic rate constant',
                   'equilibrium constant', 'forward maximal velocity',
                   'reverse maximal velocity', 'reaction affinity']
    pmin = {}
    pmax = {}
    pseudos = {}
    priors = {}
    
    for row in sbtab_prior.value_rows:
        pmin[row[sbtab_prior.columns_dict['!QuantityType']]] = float(row[sbtab_prior.columns_dict['!LowerBound']])
        pmax[row[sbtab_prior.columns_dict['!QuantityType']]] = float(row[sbtab_prior.columns_dict['!UpperBound']])
        if row[sbtab_prior.columns_dict['!MathematicalType']] == 'Additive':
            std = row[sbtab_prior.columns_dict['!PriorStd']]
        else:
            std = row[sbtab_prior.columns_dict['!PriorGeometricStd']]
        median = row[sbtab_prior.columns_dict['!PriorMedian']]

        if row[sbtab_prior.columns_dict['!QuantityType']] in pseudo_list:
            pseudos[row[sbtab_prior.columns_dict['!QuantityType']]] = [float(median),
                                                                       float(std)]
        else:
            priors[row[sbtab_prior.columns_dict['!QuantityType']]] = [float(median),
                                                                      float(std)]
       
    return pseudos, priors, pmin, pmax


def id_checker(sbtab, sbml):
    '''
    this function checks, whether all the entries of the SBML ID columns of
    the SBtab file can also be found in the SBML file. If not, these are
    omitted during the balancing. But there should be a warning to raise user
    awareness.
    '''
    sbtabid2sbmlid = []
    reaction_ids_sbml = []
    species_ids_sbml = []
    s_id = None
    r_id = None

    for reaction in sbml.getListOfReactions():
        reaction_ids_sbml.append(reaction.getId())
    for species in sbml.getListOfSpecies():
        species_ids_sbml.append(species.getId())

    for row in sbtab.split('\n'):
        splitrow = row.split('\t')
        if len(splitrow) < 3: continue
        if row.startswith('!!'): continue
        elif row.startswith('!'):
            for i, element in enumerate(splitrow):
                if element == '!Compound:SBML:species:id': s_id = i
                elif element == '!Reaction:SBML:reaction:id': r_id = i
            if s_id is None:
                sbtabid2sbmlid.append('''Error: The SBtab file lacks the obliga
                tory column "'"!Compound:SBML:species:id"'" to link the paramet
                er entries to the SBML model species.''')
            if r_id is None:
                sbtabid2sbmlid.append('''Error: The SBtab file lacks the obliga
                tory column "'"!Reaction:SBML:reaction:id"'" to link the parame
                ter entries to the SBML model reactions.''')
        else:
            try:
                if splitrow[s_id] != '' \
                   and splitrow[s_id] not in species_ids_sbml \
                   and splitrow[s_id] != 'nan' and splitrow[s_id] != 'None':
                    sbtabid2sbmlid.append('''Warning: The SBtab file holds a sp
                    ecies ID which does not comply to any species ID in the SBM
                    L file: %s''' % (splitrow[s_id]))
            except: pass
            try:
                if splitrow[r_id] != '' \
                   and splitrow[r_id] not in reaction_ids_sbml \
                   and splitrow[r_id] != 'nan' and splitrow[r_id] != 'None':
                    sbtabid2sbmlid.append('''Warning: The SBtab file holds a re
                    action ID which does not comply to any reaction ID in the S
                    BML file: %s''' % (splitrow[r_id]))
            except: pass

    return sbtabid2sbmlid





def xml2html(sbml_file):
    '''
    generates html view out of xml file
    '''
    old_sbml = str(sbml_file).split('\\n')
    new_sbml = '<xmp>'

    for row in old_sbml:
        new_sbml += row + '\n'

    new_sbml += '</xmp>'
        
    return new_sbml


def id_checker(sbtab, sbml):
    '''
    this function checks, whether all the entries of the SBML ID columns of the SBtab file can also be
    found in the SBML file. If not, these are omitted during the balancing. But there should be a warning
    to raise user awareness.
    '''
    sbtabid2sbmlid = []

    reaction_ids_sbml = []
    species_ids_sbml  = []

    s_id = None
    r_id = None

    for reaction in sbml.getListOfReactions():
        reaction_ids_sbml.append(reaction.getId())
    for species in sbml.getListOfSpecies():
        species_ids_sbml.append(species.getId())

    for row in sbtab.value_rows:
        if len(row) < 3: continue
        try:
            s_id = sbtab.columns_dict['!Compound:SBML:species:id']
            r_id = sbtab.columns_dict['!Reaction:SBML:reaction:id']
        except:
            sbtabid2sbmlid.append('Error: The SBtab file lacks either of the obligatory columns "'"!Compound:SBML:species:id"'" or "'"!Reaction:SBML:reaction:id"'" to link the parameter entries to the SBML model species.')
        try:
            if row[s_id] != '' and row[s_id] not in species_ids_sbml and row[s_id] != 'nan' and row[s_id] != 'None':
                sbtabid2sbmlid.append('Warning: The SBtab file holds a species ID which does not comply to any species ID in the SBML file: %s'%(row[s_id]))
        except: pass
        try:
            if row[r_id] != '' and row[r_id] not in reaction_ids_sbml and row[r_id] != 'nan' and row[r_id] != 'None':
                sbtabid2sbmlid.append('Warning: The SBtab file holds a reaction ID which does not comply to any reaction ID in the SBML file: %s'%(row[r_id]))
        except: pass
            
    return sbtabid2sbmlid


def cut_tabs(sbtab_strings):
    '''
    cuts an SBtab document in single SBtab files
    '''
    sbtabs = []
    sbtab_string = ''
    counter = 1

    for row in sbtab_strings.split('\n'):
        if row.startswith('!!'):
            if sbtab_string == '':
                sbtab_string = row + '\n'
                continue
            else:
                try:
                    sbtabs.append(sbtab_string)
                    sbtab_string = row + '\n'
                    counter += 1
                except:
                    print('Warning: Could not write SBtab %s' % counter)
                    counter += 1
        else:
            sbtab_string += row + '\n'
    sbtabs.append(sbtab_string) 
                    
    return sbtabs


    
    '''
    # THIS HERE DOESN'T HELP US WITHOUT THE SBTAB DOCUMENT CLASS
    sbtabs = []
    sbtab_string = ''
    counter = 1

    for row in sbtab_strings.split('\n'):
        print(row)
        if row.startswith('!!'):
            if sbtab_string == '':
                sbtab_string = row + '\n'
                continue
            else:
                try:
                    sbtab = SBtab.SBtabTable(sbtab_string, 'sbtab_%s.tsv' % counter)
                    sbtabs.append(sbtab)
                    sbtab_string = row + '\n'
                    print('Tis written')
                    counter += 1
                except:
                    print('Warning: Could not write SBtab %s' % counter)
                    counter += 1
        else:
            sbtab_string += row + '\n'

    sbtab = SBtab.SBtabTable(sbtab_string, 'sbtab_%s.tsv' % counter)
    sbtabs.append(sbtab) 
                    
    return sbtabs
    '''

def tsv_to_html(sbtab, filename=None):
    '''
    generates html view out of tsv file
    '''
    if type(sbtab) == str and filename:
        ugly_sbtab = sbtab.split('\n')
        nice_sbtab = '<p><h2><b>%s</b></h2></p>' % filename
        delimiter = check_delimiter(sbtab)
    else:
        ugly_sbtab = sbtab.return_table_string().split('\n')
        nice_sbtab = '<p><h2><b>'+sbtab.filename+'</b></h2></p>'
        delimiter = sbtab.delimiter

    first = True
    for row in ugly_sbtab:
        if row.startswith('!!') and first:
            nice_sbtab += '<a style="background-color:#00BFFF">'+row+'</a><br>'
            nice_sbtab += '<table>'
            first = False
        elif row.startswith('!!'):
            nice_sbtab += '</table>'
            nice_sbtab += '<a style="background-color:#00BFFF">'+row+'</a><br>'
            nice_sbtab += '<table>'
        elif row.startswith('!'):
            nice_sbtab += '<tr bgcolor="#87CEFA">'
            splitrow = row.split(delimiter)
        elif row.startswith('%'):
            nice_sbtab += '<tr bgcolor="#C0C0C0">'
        elif row.startswith('Parameter balancing log file'):
            nice_sbtab += '<table><tr>'
        else: nice_sbtab += '<tr>'
        
        for i,thing in enumerate(row.split(delimiter)):
            if thing.startswith('!!'): continue
            new_row = '<td>'+str(thing)+'</td>'
            nice_sbtab += new_row
        nice_sbtab += '</tr>'
    nice_sbtab += '</table>'     

    return nice_sbtab  
