"""
SBtab2SBML
==========

Python script that converts SBtab file/s to SBML.
"""
#!/usr/bin/env python
import re
import libsbml
try: from . import SBtab
except: import SBtab
try: from . import tablibIO
except: import tablibIO
import string
import random
import sys

# all allowed secondary SBtab table types
sbtab_types = ['Quantity', 'Event', 'Rule']
urns = ['obo.chebi', 'kegg.compound', 'kegg.reaction', 'obo.go', 'obo.sgd',
        'biomodels.sbo', 'ec-code', 'kegg.orthology', 'uniprot']

class ConversionError(Exception):
    '''
    Base class for errors in the SBtab conversion class.
    '''    
    def __init__(self,message):
        self.message = message
        
    def __str__(self):
        return self.message

class SBtabDocument:
    '''
    SBtab document to be converted to SBML model
    '''
    def __init__(self, sbtab_doc):
        '''
        Initalizes SBtab document, checks it for SBtab count.
        If there are more than 1 SBtab file to be converted, provide a "tabs" parameter higher than 1.

        Parameters
        ----------
        sbtab : SBtab Document object
           SBtab Document Class.
        '''
        self.sbtab_doc = sbtab_doc
        self.filename = sbtab_doc.name
        self.warnings = []

    def convert_to_sbml(self, sbml_version):
        '''
        Generates the SBML file using the provided SBtab file/s.
        '''
        # initialize new model
        self.new_document = libsbml.SBMLDocument()
        self.new_model = self.new_document.createModel()
        self.new_model.setId('default_id')
        self.new_model.setName('default_name')
        if sbml_version == '24':
            self.new_document.setLevelAndVersion(2,4)
        elif sbml_version == '31':
            self.new_document.setLevelAndVersion(3,1)
        else:
            self.warnings.append('The given SBML version %s could not be gene'\
                                 'rated.' % sbml_version)
            return (False, self.warnings)

        # initialize some required variables for conversion
        self.reaction_list = []
        self.species_list = []
        self.compartment_list = []
        self.modifier_list = []
        self.id2sbmlid = {}

        # 1. build compartment
        try:
            compartment = self.check_compartments()
            if not compartment:
                self.warnings.append('No compartment could be initialised for'\
                                     'the model. Please check provided'\
                                     'compartment information.')
                return (False, self.warnings)
        except:
            self.warnings.append('Error: The compartment initialisation crash'\
                                 'ed. Please check for valid compartment info'\
                                 'rmation.')
            return (False, self.warnings)

        # 2. build compounds
        #self.compound_sbtab()

        try:
            if 'Compound' in self.sbtab_doc.types:
                self.compound_sbtab()
        except:
            self.warnings.append('Warning: The provided compounds could not b'\
                                 'e initialised properly. Please check for va'\
                                 'lid compound information.')
        # 3. build reactions
        try:
            if 'Reaction' in self.sbtab_doc.types:
                self.reaction_sbtab()
        except:
            self.warnings.append('Error: The provided reaction information co'\
                                 'uld not be converted. Please check for vali'\
                                 'd reaction information.')

        # 4. check for secondary SBtab table types
        for table_type in sbtab_types:
            if table_type in self.sbtab_doc.types:
                try:
                    name = 'self.' + table_type.lower() + '_sbtab()'
                    eval(name)
                except:
                    self.warnings.append('Warning: Could not process informat'\
                                         'ion from SBtab %s.' % table_type)

        # write generated information to SBML model
        new_sbml_model = libsbml.writeSBMLToString(self.new_document)
        return (new_sbml_model, self.warnings)

    def return_warnings(self):
        '''
        return warnings from the SBML conversion.
        '''
        return self.warnings
    
    def check_compartments(self):
        '''
        compartment build up is neccessary and tricky:
        either we have a compartment SBtab, or if not, we have to check a
        possible reaction and/or compound SBtab for compartments; if all of
        these do not work, we need one default compartment
        '''
        def_comp_set = False
        #1. check for compartment SBtab
        if 'Compartment' in self.sbtab_doc.types:
            try:
                self.compartment_sbtab()
                return True
            except:
                self.warnings.append('There was a compartment SBtab but it '\
                                     'could not be used for SBML compartment '\
                                     'initialisation.')

        #2. if there was no compartment SBtab given, check whether it is given
        # in the other SBtabs; if we find a compartment, we return True since
        # the compartment will be built upon SBtab Compound procession
        if 'Compound' in self.sbtab_doc.types:
            sbtab_compound = self.sbtab_doc.type_to_sbtab['Compound']
            if '!Location' in sbtab_compound.columns:
                for row in sbtab_compound.value_rows:
                    if row[sbtab_compound.columns_dict['!Location']] != '':
                        return True

        #3. if there was no compartment SBtab given and no Location found in
        # compound SBtab, check whether it is given
        # in the other SBtabs; if we find a compartment, we return True since
        # the compartment will be built upon SBtab Reaction procession
        if 'Reaction' in self.sbtab_doc.types:
            sbtab_reaction = self.sbtab_doc.type_to_sbtab['Reaction']
            if '!Location' in sbtab_reaction.columns:
                for row in sbtab_reaction.value_rows:
                    if row[sbtab_reaction.columns_dict['!Location']] != '':
                        return True

        #4. Nothing yet? Then create a default compartment
        self.def_comp_set = True
        default_compartment = self.new_model.createCompartment()
        default_compartment.setId('Default_Compartment')
        default_compartment.setName('Default_Compartment')
        default_compartment.setSize(1)
        self.compartment_list.append('Default_Compartment')
        return True

    def set_annotation(self, element, annotation, urn, elementtype):
        '''
        Set an annotation for a given SBML element.

        Parameters
        ----------
        element : libsbml object
           Element that needs to be annotated.
        annotation : str
           The identifier part of the annotation string.
        urn : str
           URN that links to the external web resource.
        elementtype : str
           What kind of element needs to be annotated? Model or Biological?
        '''
        element.setMetaId(element.getId() + "_meta")
        cv_term = libsbml.CVTerm()
        if elementtype == 'Model':
            cv_term.setQualifierType(0)
            cv_term.setModelQualifierType(libsbml.BQB_IS)
        else:
            cv_term.setQualifierType(1)
            cv_term.setBiologicalQualifierType(libsbml.BQB_IS)

        resource_term = "http://identifiers.org/" + urn + '/' + annotation
        cv_term.addResource(resource_term)

        return cv_term

    def compartment_sbtab(self):
        '''
        extract information from the Compartment SBtab
        '''
        sbtab_compartment = self.sbtab_doc.type_to_sbtab['Compartment']

        # build compartments
        for row in sbtab_compartment.value_rows:
            # name and id of compartment (optional SBML id)
            if row[sbtab_compartment.columns_dict['!ID']] not in self.compartment_list:
                compartment = self.new_model.createCompartment()
                if '!SBML:compartment:id' in sbtab_compartment.columns and \
                   row[sbtab_compartment.columns_dict['!SBML:compartment:id']] != '':
                    compartment.setId(str(row[sbtab_compartment.columns_dict['!SBML:compartment:id']]))
                else:
                    compartment.setId(str(row[sbtab_compartment.columns_dict['!ID']]))
                if '!Name' in sbtab_compartment.columns and \
                   row[sbtab_compartment.columns_dict['!Name']] != '':
                    compartment.setName(str(row[sbtab_compartment.columns_dict['!Name']]))
                else:
                    compartment.setName(str(row[sbtab_compartment.columns_dict['!Compartment']]))
            self.compartment_list.append(row[sbtab_compartment.columns_dict['!ID']])

            # set the compartment size and SBOterm if given
            if '!Size' in sbtab_compartment.columns and \
               row[sbtab_compartment.columns_dict['!Size']] != '':
                try: compartment.setSize(float(row[sbtab_compartment.columns_dict['!Size']]))
                except: pass

            if '!SBOTerm' in sbtab_compartment.columns and \
               row[sbtab_compartment.columns_dict['!SBOTerm']] != '':
                try: compartment.setSBOTerm(int(row[sbtab_compartment.columns_dict['!SBOTerm']][4:]))
                except: pass

            # search for identifiers and annotations
            for column in sbtab_compartment.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab_compartment.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    try:
                        cv_term = self.set_annotation(compartment, annot,
                                                      urn, 'Model')
                        compartment.addCVTerm(cv_term)
                    except:
                        print('There was an annotation that could not be assi'\
                              'gned properly: ',compartment.getId(), annot)
           
    def compound_sbtab(self):
        '''
        extract information from the Compound SBtab and writes it to the model
        '''
        sbtab_compound = self.sbtab_doc.type_to_sbtab['Compound']
        # build compounds
        for row in sbtab_compound.value_rows:
            if row[sbtab_compound.columns_dict['!ID']] not in self.species_list:
                species = self.new_model.createSpecies()

                # name and id of compartment (optional SBML id)
                if '!SBML:species:id' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!SBML:species:id']] != '':
                    species.setId(str(row[sbtab_compound.columns_dict['!Compound:SBML:species:id']]))
                    self.id2sbmlid[row[sbtab_compound.columns_dict['!ID']]] = row[sbtab_compound.columns_dict['!Compound:SBML:species:id']]
                else:
                    species.setId(str(row[sbtab_compound.columns_dict['!ID']]))
                    self.id2sbmlid[row[sbtab_compound.columns_dict['!ID']]] = None
                if '!Name' in sbtab_compound.columns and \
                   not row[sbtab_compound.columns_dict['!Name']] == '':
                    if '|' in row[sbtab_compound.columns_dict['!Name']]:
                        species.setName(str(row[sbtab_compound.columns_dict['!Name']].split('|')[0]))
                    else: species.setName(str(row[sbtab_compound.columns_dict['!Name']]))
                self.species_list.append(species.getId())

                # speciestype (if given)
                if '!SBML:speciestype:id' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!SBML:speciestype:id']] != '':
                    species_type = self.new_model.createSpeciesType()
                    species_type.setId(str(row[sbtab_compound.columns_dict['!SBML:speciestype:id']]))
                    species.setSpeciesType(row[sbtab_compound.columns_dict['!SBML:speciestype:id']])

                # if compartments are given, add them
                if '!Location' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!Location']] != '':
                    if not row[sbtab_compound.columns_dict['!Location']] in self.compartment_list:
                        new_comp = self.new_model.createCompartment()
                        new_comp.setId(str(row[sbtab_compound.columns_dict['!Location']]))
                        self.compartment_list.append(row[sbtab_compound.columns_dict['!Location']])
                    species.setCompartment(row[sbtab_compound.columns_dict['!Location']])
                elif self.def_comp_set:
                    species.setCompartment('Default_Compartment')

                # some more options
                if '!InitialConcentration' in sbtab_compound.columns \
                   and row[sbtab_compound.columns_dict['!InitialConcentration']] != '':
                    species.setInitialConcentration(float(row[sbtab_compound.columns_dict['!InitialConcentration']]))
                elif '!InitialValue' in sbtab_compound.columns and \
                     row[sbtab_compound.columns_dict['!InitialValue']] != '':
                    species.setInitialConcentration(float(row[sbtab_compound.columns_dict['!InitialValue']]))

                if '!IsConstant' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!IsConstant']] != '':
                    if row[sbtab_compound.columns_dict['!IsConstant']].lower() == 'false':
                        try:
                            species.setConstant(0)
                            species.setBoundaryCondition(0)
                        except: pass
                    else:
                        try:
                            species.setConstant(1)
                            species.setBoundaryCondition(1)
                        except: pass

                if '!Comment' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!Comment']] != '':
                    try:
                        note = '<body xmlns="http://www.w3.org/1999/xhtml"><p>%s</p></body>' % row[sbtab_compound.columns_dict['!Comment']]
                        species.setNotes(note)
                    except: pass

                if '!ReferenceName' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!ReferenceName']] != '':
                    try:
                        note = '<body xmlns="http://www.w3.org/1999/xhtml"><p>%s</p></body>' % row[sbtab_compound.columns_dict['!ReferenceName']]
                        species.setNotes(note)
                    except: pass
                    
                if '!SBOTerm' in sbtab_compound.columns and \
                   row[sbtab_compound.columns_dict['!SBOTerm']] != '':
                    try: species.setSBOTerm(int(row[sbtab_compound.columns_dict['!SBOTerm']][4:]))
                    except: pass

                if '!Unit' in sbtab_compound.columns and self.unit_mM == False:
                    if row[sbtab_compound.columns_dict['!Unit']] == 'mM':
                        self.unit_def_mm()
                        self.unit_mM = True
                    elif row[sbtab_compound.columns_dict['!Unit']].lower().startswith('molecules'):
                        self.unit_def_mpdw()
                        self.unit_mpdw = True

                # search for identifiers and annotations
                for column in sbtab_compound.columns_dict.keys():
                    if 'Identifiers' in column:
                        annot = row[sbtab_compound.columns_dict[column]]
                        if annot == '': continue
                        for pattern in urns:
                            if pattern in column:
                                urn = pattern
                        try:
                            cv_term = self.set_annotation(species, annot,
                                                          urn, 'Biological')
                            species.addCVTerm(cv_term)
                        except:
                            print('There was an annotation that I could not a'\
                                  'ssign properly: ',species.getId(), annot)

        #species without compartments yield errors --> set them to the first available compartment
        for species in self.new_model.getListOfSpecies():
            if not species.isSetCompartment():
                species.setCompartment(self.compartment_list[0])

    def is_number(self, s):
        '''
        test if a given string is a number masked as string; this is mainly
        important while setting SBML IDs (which must NOT be numbers)
        '''
        try:
            float(s)
            return True
        except:
            return False
                
    def reaction_sbtab(self):
        '''
        extract information from the Reaction SBtab and write it to the model
        '''
        sbtab_reaction = self.sbtab_doc.type_to_sbtab['Reaction']

        # preprocessing: if there are species in the reaction formulas, which
        #                  we have not yet created for the model, create them
        if '!ReactionFormula' in sbtab_reaction.columns:
            self.get_reactants(sbtab_reaction)
            for reaction in self.reaction2reactants:
                try: compartment = self.reaction2compartment[reaction]
                except: compartment = False
                educts = self.reaction2reactants[reaction][0]
                for educt in educts:
                    if educt == '': continue
                    if educt not in self.id2sbmlid.keys() and \
                       not educt in self.species_list:
                        sp = self.new_model.createSpecies()
                        sp.setId(str(educt))
                        sp.setName(str(educt))
                        sp.setInitialConcentration(1)
                        if compartment: sp.setCompartment(compartment)
                        elif self.def_comp_set:
                            sp.setCompartment('Default_Compartment')
                        self.species_list.append(educt)
                products = self.reaction2reactants[reaction][1]
                for product in products:
                    if product == '': continue
                    if not product in self.id2sbmlid.keys() and \
                       not product in self.species_list:
                        sp = self.new_model.createSpecies()
                        sp.setId(str(product))
                        sp.setName(str(product))
                        sp.setInitialConcentration(1)
                        if compartment: sp.setCompartment(compartment)
                        elif self.def_comp_set:
                            sp.setCompartment('Default_Compartment')
                        self.species_list.append(product)

        #if compartments are given for the reactions and these compartments are not built yet:
        if '!Location' in sbtab_reaction.columns:
            for row in sbtab_reaction.value_rows:
                if row[sbtab_reaction.columns_dict['!Location']] == '':
                    continue
                if row[sbtab_reaction.columns_dict['!Location']] not in self.compartment_list:
                    compartment = self.new_model.createCompartment()
                    compartment.setId(row[sbtab_reaction.columns_dict['!Location']])
                    compartment.setName(row[sbtab_reaction.columns_dict['!Location']])
                    compartment.setSize(1)
                    self.compartment_list.append(row[sbtab_reaction.columns_dict['!Location']])

        try:
            sbtab_reaction.columns_dict['!KineticLaw']
            self.warnings.append('Warning: Please be aware that the SBtab -> SBML conversion does not include a validation of the provided kinetic rate laws. Thus, invalid SBML code may be produced which cannot be simulated. Please check the correctness of your kinetic rate laws manually.')
        except: pass

        #creating the reactions
        for row in sbtab_reaction.value_rows:
            # if the reaction must not be included in the model: continue
            if '!BuildReaction' in sbtab_reaction.columns and \
               row[sbtab_reaction.columns_dict['!BuildReaction']] == 'False':
                continue
            react = self.new_model.createReaction()

            # set id and name
            if '!SBML:reaction:id' in sbtab_reaction.columns and \
               row[sbtab_reaction.columns_dict['!SBML:reaction:id']] != '' and \
                                                                        not self.is_number(row[sbtab_reaction.columns_dict['!SBML:reaction:id']]):
                react.setId(str(row[sbtab_reaction.columns_dict['!SBML:reaction:id']]))
            else: react.setId(str(row[sbtab_reaction.columns_dict['!ID']]))

            if '!Name' in sbtab_reaction.columns:
               if row[sbtab_reaction.columns_dict['!Name']] != '':
                if '|' in row[sbtab_reaction.columns_dict['!Name']]:
                    react.setName(str(row[sbtab_reaction.columns_dict['!Name']].split('|')[0]))
                else: react.setName(str(row[sbtab_reaction.columns_dict['!Name']]))
            else: react.setName(str(row[sbtab_reaction.columns_dict['!ID']]))

            # some more options
            if '!SBOTerm' in sbtab_reaction.columns:
                if row[sbtab_reaction.columns_dict['!SBOTerm']] != '':
                    try: react.setSBOTerm(int(row[sbtab_reaction.columns_dict['!SBOTerm']][4:]))
                    except: pass

            if '!IsReversible' in sbtab_reaction.columns and \
               row[sbtab_reaction.columns_dict['!IsReversible']] != '':
                if string.capwords(row[sbtab_reaction.columns_dict['!IsReversible']]) == 'False':
                    try: react.setReversible(0)
                    except: pass
                elif string.capwords(row[sbtab_reaction.columns_dict['!IsReversible']]) == 'True':
                    try: react.setReversible(1)
                    except: pass                   

            #if sumformula is at hand: generate reactants and products
            if '!ReactionFormula' in sbtab_reaction.columns and \
               row[sbtab_reaction.columns_dict['!ReactionFormula']] != '':
                for educt in self.reaction2reactants[row[sbtab_reaction.columns_dict['!ID']]][0]:
                    if educt == '': continue
                    reactant = react.createReactant()
                    if educt in self.id2sbmlid:
                        if self.id2sbmlid[educt] != None:
                            reactant.setSpecies(self.id2sbmlid[educt])
                        else: reactant.setSpecies(educt)
                    else: reactant.setSpecies(educt)
                    reactant.setStoichiometry(self.rrps2stoichiometry[row[sbtab_reaction.columns_dict['!ID']],educt])
                for product in self.reaction2reactants[row[sbtab_reaction.columns_dict['!ID']]][1]:
                    if product == '': continue
                    reactant = react.createProduct()
                    if product in self.id2sbmlid.keys():
                        if self.id2sbmlid[product] != None: reactant.setSpecies(self.id2sbmlid[product])
                        else: reactant.setSpecies(product)
                    else: reactant.setSpecies(product)
                    reactant.setStoichiometry(self.rrps2stoichiometry[row[sbtab_reaction.columns_dict['!ID']],product])
            '''
            #Uncomment, if we start using SBML Level 3, Version 1, or higher
            #if location is at hand: link reaction to the compartment
            if sbtab_reaction.columns_dict['!Location'] and row[sbtab_reaction.columns_dict['!Location']] != '':
                if not row[sbtab_reaction.columns_dict['!Location']] in self.compartment_list:
                    new_comp = self.new_model.createCompartment()
                    new_comp.setId(row[sbtab_reaction.columns_dict['!Location']])
                    react.setCompartment(new_comp)
                    self.compartment_list.append(row[sbtab_reaction.columns_dict['!Location']])
                else:
                    react.setCompartment(row[loc_column])
            '''
            #if an enzyme is given, mark it as modifier to the reaction
            try:
                sbtab_reaction.columns_dict['!Regulator']
                if row[sbtab_reaction.columns_dict['!Regulator']] != '':
                    if "|" in row[sbtab_reaction.columns_dict['!Regulator']]:
                        splits = row[sbtab_reaction.columns_dict['!Regulator']].split('|')
                        for element in splits:
                            #element = element.strip()
                            if element.startswith('+') and element[1:] in self.species_list:
                                try:
                                    mod = react.createModifier()
                                    mod.setSpecies(element[1:])
                                    mod.setSBOTerm(459)
                                    self.modifier_list.append(element)
                                except: pass
                            elif element.startswith('-') and element[1:] in self.species_list:
                                try:
                                    mod = react.createModifier()
                                    mod.setSpecies(element[1:])
                                    mod.setSBOTerm(20)
                                    self.modifier_list.append(element)
                                except: pass
                            elif element[1:] in self.species_list:
                                try:
                                    mod = react.createModifier()
                                    mod.setSpecies(element[1:])
                                    self.modifier_list.append(element)
                                    letring = str('Warning: The reaction modifier '+element+' could not be identified as either stimulator or inhibitor. Please add an SBO Term.')
                                    self.warnings.append(letring)
                                except: pass
                    else:
                        if (row[sbtab_reaction.columns_dict['!Regulator']]).startswith('+') and row[sbtab_reaction.columns_dict['!Regulator']] in self.species_list:
                            try:
                                mod = react.createModifier()
                                mod.setSpecies(row[sbtab_reaction.columns_dict['!Regulator']][1:])
                                mod.setSBOTerm(459)
                                self.modifier_list.append(row[sbtab_reaction.columns_dict['!Regulator']])
                            except: pass
                        elif (row[sbtab_reaction.columns_dict['!Regulator']]).startswith('-') and row[sbtab_reaction.columns_dict['!Regulator']] in self.species_list:
                            try:
                                mod = react.createModifier()
                                mod.setSpecies(row[sbtab_reaction.columns_dict['!Regulator']][1:])
                                mod.setSBOTerm(20)
                                self.modifier_list.append(row[sbtab_reaction.columns_dict['!Regulator']])
                            except: pass
                        elif row[sbtab_reaction.columns_dict['!Regulator']] in self.species_list:
                            try:
                                mod = react.createModifier()
                                mod.setSpecies(row[sbtab_reaction.columns_dict['!Regulator']])
                                self.modifier_list.append(row[sbtab_reaction.columns_dict['!Regulator']])
                                letring = str('Warning: The reaction modifier '+row[sbtab_reaction.columns_dict['!Regulator']]+' could not be identified as either stimulator or inhibitor. Please add an SBO Term.')
                                self.warnings.append(letring)
                            except: pass
                        self.modifier_list.append(row[sbtab_reaction.columns_dict['!Regulator']])

            except: pass
            '''
            #if metabolic regulators are given: extract them and create them
            try:
                sbtab_reaction.columns_dict['!MetabolicRegulators']
                if row[sbtab_reaction.columns_dict['!MetabolicRegulators']] != '':
                    acts,inhs = self.extractRegulators(row[sbtab_reaction.columns_dict['!MetabolicRegulators']])
                    for activator in acts:
                        if activator in self.species_list and activator not in self.modifier_list:
                            acti = react.createModifier()
                            acti.setSpecies(activator)
                            acti.setSBOTerm(459)
                            #react.addModifier(acti)
                            self.modifier_list.append(acti)
                    for inhibitor in inhs:
                        if inhibitor in self.species_list and inhibitor not in self.modifier_list:
                            inhi = react.createModifier()
                            inhi.setSpecies(inhibitor)
                            inhi.setSBOTerm(20)
                            #react.addModifier(inhi)
                            self.modifier_list.append(inhi)
            except: pass
            '''
            #set annotations if given:
            ####
            ####commented out, because by some reason this crashes the reaction (!)
            ####
            '''
            for column in sbtab_reaction.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab_reaction.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    print urn
                    try:
                        cv_term = self.set_annotation(react,annot,urn,'Biological')
                        react.addCVTerm(cv_term)
                    except:
                        print 'There was an annotation that I could not assign properly: ',react.getId(),annot #,urn

            '''
            #since local parameters need to be entered *after* reaction creation, but *before* setting
            try:
                sbtab_reaction.columns_dict['!KineticLaw']
                if row[sbtab_reaction.columns_dict['!KineticLaw']] != '':
                    kl = react.createKineticLaw()
                    formula = row[sbtab_reaction.columns_dict['!KineticLaw']]
                    kl.setFormula(formula)
                    react.setKineticLaw(kl)
                    #for erraneous laws: remove them
                    if react.getKineticLaw().getFormula() == '':
                        react.unsetKineticLaw()
            except: pass

            
    def unit_def_mm(self):
        '''
        build unit definition
        '''
        ud = self.new_model.createUnitDefinition()
        ud.setId('mM')
        ud.setName('mM')

        mole = ud.createUnit()
        mole.setScale(-3)
        mole.setKind(libsbml.UNIT_KIND_MOLE)

        litre = ud.createUnit()
        litre.setExponent(-1)
        litre.setKind(libsbml.UNIT_KIND_LITRE)

    def unit_def_mpdw(self):
        '''
        build unit definition
        '''
        ud = self.new_model.createUnitDefinition()
        ud.setId('mmol_per_gDW_per_hr')
        ud.setName('mmol_per_gDW_per_hr')

        mole = ud.createUnit()
        mole.setScale(-3)
        mole.setExponent(1)
        mole.setKind(libsbml.UNIT_KIND_MOLE)

        litre = ud.createUnit()
        litre.setScale(0)
        litre.setExponent(-1)
        litre.setKind(libsbml.UNIT_KIND_GRAM)
       
        second = ud.createUnit()
        second.setScale(0)
        second.setExponent(-1)
        second.setMultiplier(0.000277777777777778)
        second.setKind(libsbml.UNIT_KIND_SECOND)

    def extractRegulators(self,mods):
        '''
        Extracts the regulators from the column "Regulator".

        Parameters
        ----------
        mods : str
           The modifiers of a reaction.
        '''
        activators = []
        inhibitors = []

        splits = mods.split(' ')

        for i,element in enumerate(splits):
            if element == '+': activators.append(splits[i+1])
            elif element == '-': inhibitors.append(splits[i+1])

        return activators,inhibitors

        
    def get_reactants(self, sbtab):
        '''
        Extracts the reactants from the a reaction formula.

        Parameters
        ----------
        sbtab : SBtab object
           SBtab file as SBtab object.
        '''
        self.reaction2reactants = {}
        self.rrps2stoichiometry = {}
        self.reaction2compartment = {}
        self.specrect2compartment = {}
        educts = []
        products = []
        
        for reaction in sbtab.value_rows:
            r_id = reaction[sbtab.columns_dict['!ID']]
            if '!Location' in sbtab.columns:
                self.reaction2compartment[r_id] = reaction[sbtab.columns_dict['!Location']]
            sum_formula  = reaction[sbtab.columns_dict['!ReactionFormula']]

            #is a compartment given for the reaction? (nice, but we cannot set it (only in SBML version 3))
            if sum_formula.startswith('['):
                self.reaction2compartment[r_id] = re.search('[([^"]*)]',sum_formula).group(1)

            #check the educts
            try:
                educt_list = re.search('([^"]*)<=>',sum_formula).group(1)
                educts = []
                for educt in educt_list.split('+'):
                    try:
                        float(educt.lstrip().rstrip().split(' ')[0])
                        self.rrps2stoichiometry[r_id,educt.lstrip().rstrip().split(' ')[1]] = float(educt.lstrip().rstrip().split(' ')[0])
                        educts.append(educt.lstrip().rstrip().split(' ')[1])
                    except:
                        self.rrps2stoichiometry[r_id,educt.lstrip().rstrip()] = 1
                        educts.append(educt.lstrip().rstrip())
            except: pass

            #check the products
            try:
                product_list = re.search('<=>([^"]*)',sum_formula).group(1)
                products = []
                for product in product_list.split('+'):
                    try:
                        float(product.lstrip().rstrip().split(' ')[0])
                        self.rrps2stoichiometry[r_id,product.lstrip().rstrip().split(' ')[1]] = float(product.lstrip().rstrip().split(' ')[0])
                        products.append(product.lstrip().rstrip().split(' ')[1])
                    except:
                        self.rrps2stoichiometry[r_id,product.lstrip().rstrip()] = 1
                        products.append(product.lstrip().rstrip())
            except: pass

            self.reaction2reactants[r_id] = [educts,products]

    def quantitySBtab(self):
        '''
        Extracts the information from the Quantity SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Quantity']

        for row in sbtab.value_rows:
            try:
                row[sbtab.columns_dict['!Description']]
                if row[sbtab.columns_dict['!Description']] == 'local parameter':
                    for reaction in self.new_model.getListOfReactions():
                        kl      = reaction.getKineticLaw()
                        formula = kl.getFormula()
                        if row[sbtab.columns_dict['!SBML:reaction:parameter:id']] in formula:
                            lp = kl.createParameter()
                            lp.setId(row[sbtab.columns_dict['!SBML:reaction:parameter:id']])
                            try: lp.setValue(float(row[sbtab.columns_dict['!Value']]))
                            except: lp.setValue(1.0)
                            try: lp.setUnits(row[sbtab.columns_dict['!Unit']])
                            except: pass
                            if '!Unit' in sbtab.columns and self.unit_mM == False:
                                if row[sbtab.columns_dict['!Unit']] == 'mM':
                                    self.unit_def_mm()
                                    self.unit_mM = True
                                elif row[sbtab.columns_dict['!Unit']].lower().startswith('molecules'):
                                    self.unit_def_mpdw()
                                    self.unit_mpdw = True
                            if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                                try: lp.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                                except: pass
                else:
                    parameter = self.new_model.createParameter()
                    parameter.setId(row[sbtab.columns_dict['!SBML:reaction:parameter:id']])
                    parameter.setUnits(row[sbtab.columns_dict['!Unit']])
                    parameter.setValue(float(row[sbtab.columns_dict['!Value']]))
                    if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                        try: parameter.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                        except: pass
            except:
                parameter = self.new_model.createParameter()
                parameter.setId(row[sbtab.columns_dict['!SBML:reaction:parameter:id']])
                try: parameter.setValue(float(row[sbtab.columns_dict['!Value']]))
                except: parameter.setValue(1.0)
                try: parameter.setUnits(row[sbtab.columns_dict['!Unit']])
                except: pass
                if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                    try: parameter.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                    except: pass

    def eventSBtab(self):
        '''
        Extracts the information from the Event SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Event']

        for row in sbtab.value_rows:
            event = self.new_model.createEvent()
            event.setId(row[sbtab.columns_dict['!Event']])
            try: event.setName(row[sbtab.columns_dict['!Name']])
            except: pass
            try:
                if row[sbtab.columns_dict['!Assignments']] != '':
                    asses = row[sbtab.columns_dict['!Assignments']].split('|')
                    if len(asses) > 1:
                        for ass in asses:
                            ea  = event.createEventAssignment()
                            var = ass.split('=')[0].strip()
                            val = ass.split('=')[1].strip()
                            ea.setMath(libsbml.parseL3Formula(val))
                            ea.setVariable(var)
                    else:
                        ea  = event.createEventAssignment()
                        var = asses[0].split('=')[0].strip()
                        val = asses[0].split('=')[1].strip()
                        ea.setMath(libsbml.parseL3Formula(val))
                        ea.setVariable(var)
            except: pass            
            try:
                if row[sbtab.columns_dict['!Trigger']] != '' and row[sbtab.columns_dict['!Trigger']] != 'None':
                    trig = event.createTrigger()
                    trig.setMetaId(row[sbtab.columns_dict['!Event']]+'_meta')
                    trig.setMath(libsbml.parseL3Formula(row[sbtab.columns_dict['!Trigger']]))
            except: pass
            if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                try: event.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                except: pass
            try:
                if row[sbtab.columns_dict['!Delay']] != '' and row[sbtab.columns_dict['!Delay']] != 'None':
                    dl = event.createDelay()
                    dl.setMath(libsbml.parseL3Formula(row[sbtab.columns_dict['!Delay']]))
            except: pass
            try:
                if row[sbtab.columns_dict['!UseValuesFromTriggerTime']] == 'False':
                    event.setUseValuesFromTriggerTime(0)
                else:
                    event.setUseValuesFromTriggerTime(1)
            except: pass

            for column in sbtab.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    try:
                        cv_term = self.set_annotation(event,annot,urn,'Biological')
                        event.addCVTerm(cv_term)
                    except:
                        print('There was an annotation that I could not assign properly: ',event.getId(),annot)

    def ruleSBtab(self):
        '''
        Extracts the information from the Rule SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Rule']

        for row in sbtab.value_rows:
            if row[sbtab.columns_dict['!Name']] == 'assignmentRule':
                rule = self.new_model.createAssignmentRule()
            elif row[sbtab.columns_dict['!Name']] == 'algebraicRule':
                rule = self.new_model.createAlgebraicRule()
            elif row[sbtab.columns_dict['!Name']] == 'rateRule':
                rule = self.new_model.createRateRule()
            else: continue
            rule.setMetaId(row[sbtab.columns_dict['!Rule']]+'_meta')
            try: rule.setName(row[sbtab.columns_dict['!Name']])
            except: pass
            try: rule.setUnits(row[sbtab.columns_dict['!Unit']])
            except: pass
            try:
                if row[sbtab.columns_dict['!Formula']] != '':
                    asses = row[sbtab.columns_dict['!Formula']]
                    var = asses.split('=')[0].strip()
                    val = asses.split('=')[1].strip()
                    rule.setMath(libsbml.parseL3Formula(val))
                    rule.setVariable(var)
            except:
                pass
            for column in sbtab.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    try:
                        cv_term = self.set_annotation(event,annot,urn,'Biological')
                        rule.addCVTerm(cv_term)
                    except:
                        print('There was an annotation that I could not assign properly: ',rule.getId(),annot)

if __name__ == '__main__':

    try: sys.argv[1]
    except:
        print('You have not provided input arguments. Please start the script by also providing an SBtab file and an optional SBML output filename: >python sbtab2sbml.py SBtabfile.csv Output')
        sys.exit()
        
    file_name    = sys.argv[1]
    sbtab_file_o = open(file_name,'r')
    sbtab_file   = sbtab_file_o.read()

    try: output_name = sys.argv[2]+'.xml'
    except: output_name = file_name[:-4]+'.xml'

    Converter_class = SBtabDocument(sbtab_file,file_name)
    SBML_output     = Converter_class.makeSBML()
    new_SBML_file   = open(output_name,'w')
    new_SBML_file.write(SBML_output[0])
    new_SBML_file.close()

    print('The SBML file has been successfully written to your working directory or chosen output path.')
