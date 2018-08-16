"""
SBML2SBtab Converter
====================

Python script that converts an SBML file to SBtab file/s.

See specification for further information.
"""
#!/usr/bin/env python
import re, libsbml, numpy
try: import SBtab
except: from . import SBtab
import sys

#Rule and Event not updated yet (later. lack of time atm.)
#supported_table_types = ['Compartment', 'Compound', 'Reaction', 'Rule',
#                         'Quantity', 'Event']
supported_table_types = ['Compartment', 'Compound', 'Reaction', 'Quantity']


class ConversionError(Exception):
    '''
    Base class for errors in the SBtab conversion class.
    '''
    def __init__(self,message):
        self.message = message
        
    def __str__(self):
        return self.message

class SBMLDocument:
    '''
    SBML model to be converted to SBtab file/s
    '''
    def __init__(self, sbml_model, filename):
        '''
        Initalizes SBtab document, checks it for SBtabs

        Parameters
        ----------
        sbml_model : libsbml model object
            SBML model as libsbml object.
        filename : str
            Filename with extension.
        '''
        self.model = sbml_model
        self.fbc = False
        if filename.endswith('.xml') or filename.endswith('.sbml'):
            cut = re.search('(.*)\.', filename)
            self.filename = cut.group(1)
        else: raise ConversionError('Wrong extension of file %s.' % filename)

    def convert_to_sbtab(self):
        '''
        Generates the SBtab files.
        '''
        self.warnings = []
        sbtab_doc = SBtab.SBtabDocument(self.filename)

        try:
            fbc = self.model.getPlugin('fbc')
            self.fbc = True
        except: pass

        for table_type in supported_table_types:
            try:
                function_name = 'self.'+table_type.lower()+'_sbtab()'
                sbtab = eval(function_name)
                if sbtab != False:
                    sbtab_doc.add_sbtab(sbtab)
            except:
                self.warnings.append('Could not generate SBtab %s.' % table_type)

        if self.fbc:
            try:
                sbtab = self.fbc_objective()
                if sbtab != False:
                    sbtab_doc.add_sbtab(sbtab)
            except:
                self.warnings.append('Could not generate SBtab FBC Objective Function.')
                
        return (sbtab_doc, self.warnings)

    def compartment_sbtab(self):
        '''
        build a compartment SBtab
        '''
        # header row
        sbtab_compartment  = '!!SBtab SBtabVersion="1.0" Document="%s" TableT'\
                             'ype="Compartment" TableName="Compartment"'\
                             '\n' % self.filename
        # columns
        columns = ['!ID', '!Name', '!Size', '!Unit', '!SBOTerm']
        sbtab_compartment += '\t'.join(columns) + '\n'

        # value rows
        for compartment in self.model.getListOfCompartments():
            value_row = [''] * len(columns)
            value_row[0] = compartment.getId()
            if compartment.getName() != '':
                value_row[1] = compartment.getName()
            else: value_row[1] = compartment.getId()
            try: value_row[2] = str(compartment.getSize())
            except: pass
            try: value_row[3] = str(species.getUnits())
            except: pass
            if str(compartment.getSBOTerm()) != '-1':
                value_row[4] ='SBO:%.7d'%compartment.getSBOTerm()
            sbtab_compartment += '\t'.join(value_row) + '\n'

        # build SBtab Table object from basic information
        sbtab_compartment = SBtab.SBtabTable(sbtab_compartment,
                                             self.filename + '_compartment.tsv')

        # test and extend for annotations
        for i, compartment in enumerate(self.model.getListOfCompartments()):
            try:
                annotations = self.get_annotations(compartment)
                for j, annotation in enumerate(annotations):
                    col_name = '!Identifiers:' + annotation[1]
                    if col_name not in columns:
                        new_column = [''] * (self.model.getNumCompartments() + 1)
                        new_column[0] = col_name
                        new_column[i + 1] = annotation[0]
                        columns.append(col_name)
                        sbtab_compartment.add_column(new_column)
                    else:
                        sbtab_compartment.change_value_by_name(compartment.getId(),
                                                               col_name,
                                                               annotation[0])
            except:
                self.warnings.append('Could not add the annotation %s for'\
                                     'compartment %s' % annotations[1],
                                     compartment.getId())           

        return sbtab_compartment
        
    def compound_sbtab(self):
        '''
        build a compound SBtab
        '''
        # header row
        sbtab_compound = '!!SBtab SBtabVersion="1.0" Document="%s" TableType='\
                         '"Compound" TableName="Compound"\n' % self.filename

        # columns
        columns = ['!ID', '!Name', '!Location', '!Charge', '!IsConstant',
                   '!SBOTerm', '!InitialConcentration', '!hasOnlySubstanceUnits']
        if self.fbc:
            columns = columns + ['!SBML:fbc:chemicalFormula', '!SBML:fbc:charge']
            
        sbtab_compound += '\t'.join(columns) + '\n'

        # value rows
        for species in self.model.getListOfSpecies():
            value_row = [''] * len(columns)
            value_row[0] = species.getId()
            try: value_row[1] = species.getName()
            except: pass
            try: value_row[2] = species.getCompartment()
            except: pass
            try: value_row[3] = str(species.getCharge())
            except: pass
            try: value_row[4] = str(species.getConstant())
            except: pass
            if str(species.getSBOTerm()) != '-1':
                value_row[5] ='SBO:%.7d'%species.getSBOTerm()
            try: value_row[6] = str(species.getInitialConcentration())
            except: pass
            try: value_row[7] = str(species.getHasOnlySubstanceUnits())
            except: pass

            if self.fbc:
                try:
                    fbc_plugin = species.getPlugin('fbc')
                    value_row[8] = str(fbc_plugin.getChemicalFormula())
                    value_row[9] = str(fbc_plugin.getCharge())
                except:
                    self.warnings.append('FBC Species information could not be read.')

            sbtab_compound += '\t'.join(value_row) + '\n'


        sbtab_compound = SBtab.SBtabTable(sbtab_compound,
                                          self.filename + '_compound.tsv')
            
        # test and extend for annotations
        for i, species in enumerate(self.model.getListOfSpecies()):
            try:
                annotations = self.get_annotations(species)
                for j, annotation in enumerate(annotations):
                    col_name = '!Identifiers:' + annotation[1]
                    if col_name not in columns:
                        new_column = [''] * (self.model.getNumSpecies() + 1)
                        new_column[0] = col_name
                        new_column[i + 1] = annotation[0]
                        columns.append(col_name)
                        sbtab_compound.add_column(new_column)
                    else:
                        sbtab_compound.change_value_by_name(species.getId(),
                                                            col_name,
                                                            annotation[0])
            except:
                self.warnings.append('Could not add the annotation %s for'\
                                     'compound %s' % annotations[1],
                                     species.getId())

        return sbtab_compound

    def event_sbtab(self):
        '''
        Builds an Event SBtab.
        '''
        if len(self.model.getListOfEvents()) == 0:
            return False
            
        event    = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Event" TableName="Event"'],['']]
        header   = ['!Event','!Name','!Assignments','!Trigger','!SBOterm','!Delay','!UseValuesFromTriggerTime']
        identifiers  = []
        column2ident = {}

        for eve in self.model.getListOfEvents():
            value_row = ['']*len(header)
            value_row[0] = eve.getId()
            value_row[1] = eve.getName()
            if eve.getNumEventAssignments() > 1:
                try:
                    eas = eve.getListOfEventAssignments()
                    ea_entry = ''
                    for ea in eas:
                        var       = ea.getVariable()
                        ea_e      = ea.getMath()
                        ea_entry += var+' = '+libsbml.formulaToL3String(ea_e)+' | '
                    value_row[2] = ea_entry[:-3]
                except: pass
            else:
                try:
                    eas = eve.getListOfEventAssignments()
                    for ea in eas:
                        var          = ea.getVariable()
                        vr           = ea.getMath()
                        value_row[2] = var+' = '+libsbml.formulaToL3String(vr)
                except: pass
            try:
                trigger      = eve.getTrigger().getMath()
                value_row[3] = libsbml.formulaToL3String(trigger)
            except: pass
            if str(eve.getSBOTerm()) != '-1': value_row[4] ='SBO:%.7d'%eve.getSBOTerm()
            try: value_row[5] = str(eve.getDelay())
            except: pass
            try: value_row[6] = str(eve.getUseValuesFromTriggerTime())
            except: pass
            try:
                annot_tuples = self.get_annotations(eve)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header)+1)
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]-1] = annotation[0]
            except: pass
            event.append(value_row)

        event[1] = header
        event_SB = event[0]
        for row in event[1:]:
            event_SB.append('\t'.join(row))
        event_SBtab = '\n'.join(event_SB)
            
        return [event_SBtab,'event']

    def rule_sbtab(self):
        '''
        Builds a Rule SBtab.
        '''
        if len(self.model.getListOfRules()) == 0:
            return False
            
        rule     = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Rule" TableName="Rule"'],['']]
        header   = ['!Rule','!Name','!Formula','!Unit']
        identifiers  = []
        column2ident = {}

        for ar in self.model.getListOfRules():
            value_row = ['']*len(header)
            value_row[0] = ar.getId()
            value_row[1] = ar.getElementName()
            try:
                var          = ar.getVariable()
                vr           = ar.getMath()
                value_row[2] = var+' = '+libsbml.formulaToL3String(vr)
            except: pass
            try: value_row[3] = ar.getUnits()
            except: pass            
            try:
                annot_tuples = self.get_annotations(ar)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header)+1)
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]-1] = annotation[0]
            except: pass
            rule.append(value_row)

        rule[1] = header
        rule_SB = rule[0]
        for row in rule[1:]:
            rule_SB.append('\t'.join(row))
        rule_SBtab = '\n'.join(rule_SB)
            
        return [rule_SBtab,'rule']

    def fbc_objective(self):
        '''
        builds a (preliminary?) SBtab of the (not established) TableType FbcObjective
        '''
        if self.model.getPlugin('fbc') == None:
            return False

        fbc_plugin = self.model.getPlugin('fbc')
        active_obj = fbc_plugin.getActiveObjectiveId()

        # header row
        sbtab_fbc = '!!SBtab SBtabVersion="1.0" Document="%s" TableType='\
                    '"FbcObjective" TableName="FBC Objective"\n' % self.filename

        # columns
        columns = ['!ID', '!Name', '!Type', '!Active']
            
        sbtab_fbc += '\t'.join(columns) + '\n'

        # value rows
        for obj in self.model.getListOfObjectives():
            value_row = [''] * len(columns)
            value_row[0] = obj.getId()
            try: value_row[1] = obj.getName()
            except: pass
            try: value_row[2] = obj.getType()
            except: pass
            if obj.getId() == active_obj: value[3] = 'True'
            else: value[3] = 'False'




        
    def get_annotations(self, element):
        '''
        Tries to extract an annotation from an SBML element.
        '''
        annotation = False
        urn = False
        annot_tuples = []

        pattern2urn = {"CHEBI:\d+$": "obo.chebi",
                       "C\d+$": "kegg.compound",
                       "GO:\d{7}$": "obo.go",
                       "((S\d+$)|(Y[A-Z]{2}\d{3}[a-zA-Z](\-[A-Z])?))$": "sgd",
                       "SBO:\d{7}$": "biomodels.sbo",
                       "\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$": "ec-code",
                       "K\d+$": "kegg.orthology",
                       "([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$": "uniprot"}

        if element.getCVTerms() == None: return annot_tuples
        
        for cvterm in element.getCVTerms():
            for j in range(cvterm.getNumResources()):
                resource = cvterm.getResourceURI(j)
                # this is required since URL encoding of a colon is %3A and crashes the RE
                resource = resource.replace('%3A',':')
                for pattern in pattern2urn.keys():
                    try:
                        search_annot = re.search(pattern, resource)
                        annotation = search_annot.group(0)
                        urn = pattern2urn[pattern]
                        annot_tuples.append([annotation, urn])
                        break
                    except: pass
                    try:
                        search_annot = re.search('identifiers.org/(.*)/(.*)',
                                                 resource)
                        annotation = search_annot.group(2)
                        urn = search_annot.group(1)
                        annot_tuples.append([annotation, urn])
                        break
                    except: pass

        return annot_tuples

    def reaction_sbtab(self):
        '''
        build a reaction SBtab
        '''
        # header row
        sbtab_reaction = '!!SBtab SBtabVersion="1.0" Document="%s" TableType='\
                         '"Reaction" TableName="Reaction"\n' % self.filename

        # columns
        columns = ['!ID', '!Name', '!ReactionFormula', '!Location',
                   '!Regulator', '!KineticLaw', '!SBOTerm', '!IsReversible']
        if self.fbc:
            columns = columns + ['!SBML:fbc:GeneAssociation', '!SBML:fbc:LowerBound',
                                 '!SBML:fbc:UpperBound']
            
        sbtab_reaction += '\t'.join(columns) + '\n'

        for reaction in self.model.getListOfReactions():
            value_row = [''] * len(columns)
            value_row[0] = reaction.getId()
            try: value_row[1] = reaction.getName()
            except: pass
            try: value_row[2] = self.make_sum_formula(reaction)
            except: pass
            try: value_row[3] = str(reaction.getCompartment())
            except: pass
            try:
                modifiers = reaction.getListOfModifiers()
                if len(modifiers) > 1:
                    modifier_list = ''
                    for i, modifier in enumerate(modifiers):
                        if i != len(reaction.getListOfModifiers()) - 1:
                            modifier_list += modifier.getSpecies() + '|'
                        else: modifier_list += modifier.getSpecies()
                    value_row[4] = modifier_list
                elif len(modifiers) == 1:
                    for modifier in modifiers:
                        value_row[4] = modifier.getSpecies()
            except:
                pass
            try:
                fm = reaction.getKineticLaw().getFormula()
                value_row[5] = fm.replace('\n', '')
            except: pass
            if str(reaction.getSBOTerm()) != '-1':
                value_row[6] ='SBO:%.7d' % reaction.getSBOTerm()
            try: value_row[7] = str(reaction.getReversible())
            except: pass

            if self.fbc:
                try:
                    fbc_plugin = reaction.getPlugin('fbc')
                    try:
                        ga = fbc_plugin.getGeneProductAssociation()
                        fbc_object = libsbml.FbcExtension()
                        try:
                            type_code = ga.getAssociation().getTypeCode()
                            ass_name = fbc_object.getStringFromTypeCode(type_code)
                        except: ass_name = '|'
                        if ass_name == 'FbcOr' or ass_name == 'FbcAnd':
                            k = ga.getAssociation()
                            associations = k.getListOfAssociations()
                            ass_list = ''
                            for a in associations:
                                ass_list += a.getGeneProduct() + ' ' + ass_name + ' '
                            value_row[8] = ' '.join(ass_list.split(' ')[:-2])
                        elif ass_name == 'GeneProductRef':
                            gpr = ga.getAssociation()
                            value_row[8] = gpr.getGeneProduct()
                    except: pass
                    value_row[9] = str(fbc_plugin.getLowerFluxBound())
                    value_row[10] = str(fbc_plugin.getUpperFluxBound())
                except:
                    self.warnings.append('FBC Reaction information could not be read.')
            
            sbtab_reaction += '\t'.join(value_row) + '\n'

        sbtab_reaction = SBtab.SBtabTable(sbtab_reaction,
                                          self.filename + '_reaction.tsv')

        # test and extend for annotations
        for i, reaction in enumerate(self.model.getListOfReactions()):
            try:
                annotations = self.get_annotations(reaction)
                for j, annotation in enumerate(annotations):
                    col_name = '!Identifiers:' + annotation[1]
                    if col_name not in columns:
                        new_column = [''] * (self.model.getNumReactions() + 1)
                        new_column[0] = col_name
                        new_column[i + 1] = annotation[0]
                        columns.append(col_name)
                        sbtab_reaction.add_column(new_column)
                    else:
                        sbtab_reaction.change_value_by_name(reaction.getId(),
                                                            col_name,
                                                            annotation[0])
            except:
                self.warnings.append('Could not add the annotation %s for'\
                                     'reaction %s' % annotations[1],
                                     reaction.getId())

        return sbtab_reaction

    def quantity_sbtab(self):
        '''
        Builds a Quantity SBtab.
        '''
        # if there are no parameters/quantities in the model, return False
        parameters = True
        if len(self.model.getListOfParameters()) == 0:
            parameters = False
            for reaction in self.model.getListOfReactions():
                kinetic_law = reaction.getKineticLaw()
                if len(kinetic_law.getListOfParameters()) != 0:
                    parameters = True
        if not parameters: return False

        # header row        
        sbtab_quantity = '!!SBtab SBtabVersion="1.0" Document="%s" TableType='\
                         '"Quantity" TableName="Quantity"\n' % self.filename
        # columns
        columns = ['!ID', '!Parameter:SBML:parameter:id', '!Value',
                   '!Unit', '!Type']
        sbtab_quantity += '\t'.join(columns) + '\n'

        # required for later iteration of parameter annotations
        local_parameters = []
        
        # value rows for local parameters
        for reaction in self.model.getListOfReactions():
            kinetic_law = reaction.getKineticLaw()
            if kinetic_law:
                value_row = [''] * len(columns)
                for parameter in kinetic_law.getListOfParameters():
                    value_row[0] = parameter.getId() + '_' + reaction.getId()
                    value_row[1] = parameter.getId()
                    value_row[2] = str(parameter.getValue())
                    try: value_row[3] = parameter.getUnits()
                    except: pass
                    value_row[4] = 'local parameter'
                    local_parameters.append(parameter)
                sbtab_quantity += '\t'.join(value_row) + '\n'
                

        sbtab_quantity = SBtab.SBtabTable(sbtab_quantity,
                                          self.filename + '_quantity.tsv')
                    

        # test and extend for annotations of local parameters
        for i, quantity in enumerate(local_parameters):
            try:
                annotations = self.get_annotations(quantity)
                for j, annotation in enumerate(annotations):
                    col_name = '!Identifiers:' + annotation[1]
                    if col_name not in columns:
                        new_column = [''] * (len(local_parameters) + 1)
                        new_column[0] = col_name
                        new_column[i + 1] = annotation[0]
                        columns.append(col_name)
                        sbtab_quantity.add_column(new_column)
                    else:
                        sbtab_quantity.change_value_by_name(quantity.getId(),
                                                            col_name,
                                                            annotation[0])
            except:
                self.warnings.append('Could not add the annotation %s for'\
                                     'quantity %s' % annotations[1],
                                     quantity.getId())

        # value_rows for global parameters 
        for parameter in self.model.getListOfParameters():
            value_row = [''] * len(sbtab_quantity.columns)
            value_row[0] = parameter.getId()
            value_row[1] = parameter.getId()
            value_row[2] = str(parameter.getValue())
            try: value_row[3] += parameter.getUnits()
            except: pass
            value_row[4] = 'global parameter'
            sbtab_quantity.add_row(value_row)

        # test and extend for annotations of global parameters
        for i, parameter in enumerate(self.model.getListOfParameters()):
            try:
                annotations = self.get_annotations(parameter)
                for j, annotation in enumerate(annotations):
                    col_name = '!Identifiers:' + annotation[1]
                    if col_name not in columns:
                        new_column = [''] * (self.model.getNumParameters() + 1)
                        new_column[0] = col_name
                        new_column[i + 1] = annotation[0]
                        columns.append(col_name)
                        sbtab_quantity.add_column(new_column)
                    else:
                        sbtab_quantity.change_value_by_name(parameter.getId(),
                                                            col_name,
                                                            annotation[0])
            except:
                self.warnings.append('Could not add the annotation %s for'\
                                     'quantity %s' % annotations[1],
                                     parameter.getId())
                
        return sbtab_quantity

    def make_sum_formula(self, reaction):
        '''
        Generates the reaction formula of a reaction from the list of products and list of reactants.

        Parameters
        ----------
        reaction : libsbml object reaction
           Single reaction object from the SBML file.
        '''
        sumformula = ''
        id2name = {}
        
        for species in self.model.getListOfSpecies():
            id2name[species.getId()] = species.getName()

        for i, reactant in enumerate(reaction.getListOfReactants()):
            if i != len(reaction.getListOfReactants()) - 1:
                if reactant.getStoichiometry() != 1.0:
                    sumformula += str(float(reactant.getStoichiometry())) + ' ' + reactant.getSpecies() + ' + '
                else:
                    sumformula += reactant.getSpecies() + ' + '
            else:
                if numpy.isnan(reactant.getStoichiometry()):
                    sumformula += '1 ' + reactant.getSpecies() + ' <=> '
                elif reactant.getStoichiometry() != 1.0:
                    sumformula += str(float(reactant.getStoichiometry())) + ' ' + reactant.getSpecies() + ' <=> '
                else:
                    sumformula += reactant.getSpecies() + ' <=> '
        
        if sumformula == '': sumformula += '<=> '
        
        for i, product in enumerate(reaction.getListOfProducts()):
            if i != len(reaction.getListOfProducts()) - 1:
                if product.getStoichiometry() != 1.0:
                    sumformula += str(float(product.getStoichiometry())) + ' ' + product.getSpecies() + ' + '
                else:
                    sumformula += product.getSpecies() + ' + '
            else:
                if numpy.isnan(product.getStoichiometry()):
                    sumformula += '1 ' + product.getSpecies() + ' <=> '
                elif product.getStoichiometry() != 1.0:
                    sumformula += str(float(product.getStoichiometry())) + ' ' + product.getSpecies()
                else:
                    sumformula += product.getSpecies()

        return sumformula
        

if __name__ == '__main__':

    try: sys.argv[1]
    except:
        print('You have not provided input arguments. Please start the script by also providing an SBML file and an optional SBtab output filename: >python sbml2sbtab.py SBMLfile.xml Output')
        sys.exit()

    file_name  = sys.argv[1]
    try: output_name = sys.argv[2]+'.csv'
    except: output_name = file_name[:-4]+'.csv'

    reader     = libsbml.SBMLReader()
    sbml       = reader.readSBML(file_name)
    model      = sbml.getModel()
    Sbml_class = SBMLDocument(model,file_name)

    (sbtabs,warnings) = Sbml_class.makeSBtabs()

    #print warnings

    for sbtab in sbtabs:
        sbtab_name = output_name[:-4]+'_'+sbtab[1]+output_name[-4:]
        sbtab_file = open(sbtab_name,'w')
        sbtab_file.write(sbtab[0])
        sbtab_file.close()

    print('The SBtab file/s have been successfully written to your working directory or chosen output path.')
