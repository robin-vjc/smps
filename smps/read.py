from __future__ import print_function

import os
import re
import copy
from collections import defaultdict, OrderedDict

import numpy as np
import networkx as nx
from scipy.sparse import dok_matrix
try:
    import gurobipy as gb
except ImportError as e:
    print('Gurobipy required to parse the smps files.')
    raise

# Supported stochastic model types
MODEL_TYPE_SCENARIOS_DISCRETE = 1
ALLOWED_MODEL_TYPES = [MODEL_TYPE_SCENARIOS_DISCRETE]

# Sections in the .sto file
SECTION_SCENARIO_DATA = 10
# Section in the .cor file
SECTION_ROWS_IN_COR_FILE = 11
SECTION_COLUMNS_IN_COR_FILE = 12


class StochasticModel(object):
    """
    Main class containing the information of the stochastic model. It expects the three SMPS files .cor, .tim and .sto
    to be present at the location specified in filename (e.g., filename='./data/invest' requires the files
    './data/invest.cor', './data/invest.tim' and './data/invest.sto').
    """
    def __init__(self, model_filename):
        # The following dictionaries will be filled with the matrix blocks defining the problem; we handle according to
        # the following notation standard set in the SMPS documentation:
        # min/max   1/2 (x_1, x_2,...,x_T)' |  Q_21  Q_22 ...  Q_2T  |(x_1, x_2,...,x_T)
        #                                   |                        |
        #                                   |  Q_T1  Q_T2 ...  Q_TT  |
        #
        #           +   (c_1, c_2,...,c_T)' (x_1, x_2,...,x_T)
        #
        # s.t.
        #          A_11 x_1 + A_12 x_2 + ... + A_1T x_T <> b_1
        #          A_21 x_1 + A_22 x_2 + ... + A_2T x_T <> b_2
        #             ...        ...              ...      ...
        #          A_T1 x_1 + A_T2 x_2 + ... + A_TT x_T <> b_T
        #
        #          lb_i <= x_i <= ub_i,  i=1,...,T
        #
        #          x_i in R^{n^1_i} x Z^{n^2_i}.
        #
        # It will parse
        # - the matrices A_11, A_12, ... into entries of self.A
        # - the RHS vectors b_1, b_2, ... into self.b
        # - the cost vectors c_1, c_2, ... into self.c
        # - the constraint direction (<=, ==, >=) into self.b_sense
        # - upper and lower bounds for the variables in lb and ub
        # - the optimization variable type (continuous, binary, integer)
        self.A = {}
        self.b = {}
        self.b_sense = {}
        self.c = {}
        self.lb = {}
        self.ub = {}
        self.vtype = {}

        # Additional auxiliary attributes
        self.model_filename = model_filename
        self.scenarios = OrderedDict()  # OrderedDict of instances of Scenario
        self.periods = []  # list of strings
        self.var_names = []  # var_names = ['Var1','Var2',...]
        self.constr_names = [] # constr_names = ['Constr1','Constr2',...]
        self.objective_name = ''  # label of the objective (string)
        # Record of which variable/constr belongs to which stage. Example:
        # stage_vars = { 1: ['Var1'], 2: ['Var2','Var3'], 3: ['Var4']}
        # means that Var1 belongs to stage 1, 'Var2' and 'Var3' belong to stage 2, etc. Similar for constraints:
        # stage_constrs = {1: ['Consrt1','Constr2'], 2: ['Constr3'], 3: ['Consrt4','Constr5']}
        self.stage_vars = {}
        self.stage_constrs = {}
        self.mode_of_modification = 'REPLACE'  # default, possible vaulues supported: ['REPLACE', 'ADD']

        # Parse nominal model
        self._parse_nominal_model()

        # Parse time information: using .tim file, updates value of constr_names, var_names, stage_vars, stage_constrs,
        # as well as A, b, c, lb, ub, vtype.
        self._parse_time_information()

        # Parse stochastic model: fills the self.scenarios{} dictionary
        self._parse_stochastic_information()

    def _parse_nominal_model(self):
        print('Parsing nominal model information from ' + self.model_filename + '.cor and .tim ...')
        # rename the .cor file into .mps to have Gurobi read it, as it MUST have the correct ending for it to be parsed
        os.rename(self.model_filename+'.cor',self.model_filename+'.mps')
        self.nominal_model = gb.read(self.model_filename+'.mps')
        os.rename(self.model_filename+'.mps',self.model_filename+'.cor')

        # We have to manually parse the .cor file too, because gurobi does not register the label of the objective.. >:(
        with open(self.model_filename+'.cor', 'r') as f:
            current_section = None
            for line in f:
                # Extract and clean line words
                line_word = re.split(' |\t', line)
                line_word = [x.strip() for x in line_word]
                line_word = filter(None, line_word)

                if line_word:
                    if line_word[0] == 'ROWS':
                        current_section = SECTION_ROWS_IN_COR_FILE
                    elif line_word[0] == 'COLUMNS':
                        current_section = SECTION_COLUMNS_IN_COR_FILE
                    if current_section == SECTION_ROWS_IN_COR_FILE:
                        if line_word[0] == 'N':
                            self.objective_name = line_word[1]
                    else:
                        # not interested.
                        pass

        if self.nominal_model.ModelSense != 1:
            # TODO Instead of outputting a warning, transform problem into minimization instead, and make them
            # aware that all the results have a '-' in front when it comes to objectives.
            print('WARNING: The nominal model in the .cor file is defined as a MAXIMIZATION.'
                   'The decomposition methods assume a MINIMIZATION (the rest should work).')

    def _parse_time_information(self):
        # First open .tim file and get the data
        periods_start_position = []
        with open(self.model_filename+'.tim', 'r') as f:
            for line in f:
                # Extract and clean line words
                line_word = re.split(' |\t', line)
                line_word = [x.strip() for x in line_word]
                line_word = filter(None, line_word)

                if line_word[0] in ['TIME','PERIODS','ENDATA']:
                    pass  # just padding lines
                else:
                    self.periods.append(line_word[2])
                    periods_start_position.append((line_word[0],line_word[1]))

        # Then deduce which variables belong to which stage, and which constraints belong to which stage.
        # First, fill list of variables and constraints names
        for i in self.nominal_model.getConstrs():
            self.constr_names.append(i.constrName)

        for i in self.nominal_model.getVars():
            self.var_names.append(i.varName)

        # Uopdate self.stage_vars and self.stage_constrs
        index_of_first_var_of_this_stage = 0
        index_of_first_cstr_of_this_stage = 0
        for i, stage in enumerate(self.periods):
            # SMPS convention in .tim file is to have the first entry point at the first var/constr of the model; the
            # first stage is however defined by the var/constr up to the secon entry; the first entry can thus be
            # discarded as it is redundant
            if i > 0:
                index_of_first_var_of_next_stage = self.var_names.index(periods_start_position[i][0])
                index_of_first_cstr_of_next_stage = self.constr_names.index(periods_start_position[i][1])

                self.stage_vars[i] = self.var_names[index_of_first_var_of_this_stage : index_of_first_var_of_next_stage]
                self.stage_constrs[i] = self.constr_names[index_of_first_cstr_of_this_stage : index_of_first_cstr_of_next_stage]

                # move forward the indices, to get appropriate slicing
                index_of_first_var_of_this_stage = index_of_first_var_of_next_stage
                index_of_first_cstr_of_this_stage = index_of_first_cstr_of_next_stage

        # finally append last period
        self.stage_vars[len(self.periods)] = self.var_names[index_of_first_var_of_this_stage : ]
        self.stage_constrs[len(self.periods)] = self.constr_names[index_of_first_cstr_of_this_stage : ]

        # Recover matrix blocks A[i,j] and b.
        self.A = {}
        for t, period in enumerate(self.periods):
            # pad indices so that A[1,1] (and not A[0,0]) is the first block, conform to the convention in the SMPS
            # documentation.
            tp = t+1
            # while the A matrix will typically be staircase (i.e., A[2,1] = 0 matrix), here we assume that every block
            # could contain values, again to be fully compatible with SMPS; otherwise this would have been a
            # ``for tau in range(t+1):''
            for tau, period_ in enumerate(self.periods):
                taup = tau+1
                self.A[tp,taup] = dok_matrix((len(self.stage_constrs[tp]), len(self.stage_vars[taup])), dtype=np.float32)
                # fill b[tp] entry - a vector. This should be done outside of this for loop (efficiency), but we
                # keep it here for code readability (it would be messier to swap the order of the loops).
                # self.b[tp] = dok_matrix((len(self.stage_constrs[tp]),1), dtype=np.float32)
                self.b[tp] = np.zeros(len(self.stage_constrs[tp]))
                self.b_sense[tp] = np.empty((len(self.stage_constrs[tp]),1), dtype=str)
                self.b_sense[tp][:] = ''
                # now fill row and columns of the (sparse) matrix contained in A[tp, taup] and b[tp]
                for row_local_index, row in enumerate(self.stage_constrs[tp]):
                    # get reference to gurobi constraint object
                    gb_cstr_reference = self.nominal_model.getConstrByName(row)
                    # get object with sparse representation of the constraint
                    # (awkward Gurobi feature that constr and getRow are different)
                    gb_row = self.nominal_model.getRow(gb_cstr_reference)
                    self.b[tp][row_local_index] = gb_cstr_reference.RHS
                    self.b_sense[tp][row_local_index] = gb_cstr_reference.Sense
                    # OLD: LinExpr._vars has been removed in new versions of gurobipy.
                    # for col_index, col in enumerate(gb_row._vars):
                    # NEW:
                    for col_index in range(gb_row.size()):
                        col = gb_row.getVar(col_index)
                        if col.varName in self.stage_vars[taup]:
                            # we have to create a new local index for the submatrix
                            col_local_index = self.stage_vars[taup].index(col.varName)
                            # OLD: _coeffs was removed in enw version of gurobipy
                            # self.A[tp,taup][row_local_index,col_local_index] = gb_row._coeffs[col_index]
                            # NEW:
                            self.A[tp,taup][row_local_index,col_local_index] = gb_row.getCoeff(col_index)

        # fill c, ub, lb, vtype
        for t, period in enumerate(self.periods):
            tp = t+1
            self.c[tp] = np.zeros(len(self.stage_vars[tp]))
            self.lb[tp] = np.zeros(len(self.stage_vars[tp]))
            self.ub[tp] = np.zeros(len(self.stage_vars[tp]))
            self.vtype[tp] = np.empty(len(self.stage_vars[tp]), dtype=str)
            self.vtype[tp][:] = ''
            for col_local_index, variable in enumerate(self.stage_vars[tp]):
                self.c[tp][col_local_index] = self.nominal_model.getVarByName(variable).Obj
                self.lb[tp][col_local_index] = self.nominal_model.getVarByName(variable).LB
                self.ub[tp][col_local_index] = self.nominal_model.getVarByName(variable).UB
                self.vtype[tp][col_local_index] = self.nominal_model.getVarByName(variable).VType

    def _parse_stochastic_information(self):
        print('Parsing stochastic information from '+self.model_filename+'.sto ...')
        current_section = None

        with open(self.model_filename+'.sto', 'r') as f:
            for l, line in enumerate(f):
                # Extract and clean line words
                line_word = re.split(' |\t', line)  # split by whitespace or tab
                line_word = [x.strip() for x in line_word]  # remove "\n" etc.
                line_word = filter(None, line_word)  # remove empty elements in the line_word list

                # parse to see if it's a section header
                if line_word[0] == '*':
                    pass  # it's a comment
                elif line_word[0] == 'STOCH':
                    # parse name if given
                    if len(line_word) > 1:
                        self.model_name = line_word[1]
                elif line_word[0] == 'SCENARIOS' and line_word[1] == 'DISCRETE':
                    print('Stochastic model is of type SCENARIOS DISCRETE')
                    self.model_type = MODEL_TYPE_SCENARIOS_DISCRETE
                    if len(line_word) > 2:
                        if line_word[2] in ['ADD','MULTIPLY','REPLACE']:
                            self.mode_of_modification = line_word[2]
                        else:
                            print('Non-understood parameter(s) on line {}: {}'.format(l, line_word[3:]))
                elif line_word[0] == 'SC':
                    current_section = SECTION_SCENARIO_DATA
                    current_scenario = line_word[1]
                    if line_word[2] == r"'ROOT'":
                        parent = 'ROOT'
                    else:
                        parent = line_word[2]
                    self._add_scenario(line_word[1], parent, float(line_word[3]), line_word[4])
                elif line_word[0] == 'ENDATA':
                    pass

                # otherwise we are within some section
                else:
                    if current_section == SECTION_SCENARIO_DATA:
                        # add information on modification to the nominal problem data related to this scenario
                        self.scenarios[current_scenario].add_data_modification((line_word[1], line_word[0]),
                                                                               float(line_word[2]))

        if self.model_type not in ALLOWED_MODEL_TYPES:
            raise ValueError('Only models of type '+str(ALLOWED_MODEL_TYPES)+' are supported.')

    def _add_scenario(self, scenario_id, parent, probability, branch_period):
        """
        :param id: (string) id of the scenario
        :param parent: (string) id of the parent scenario; 'ROOT' if it is an independent scenario
        :param probability: (float) probability of the scenario
        :param branch_period: (string) period starting from which the scenario differs from its parent
        :return: None
        """
        if self.model_type != MODEL_TYPE_SCENARIOS_DISCRETE:
            raise ValueError('This method can only be called on Stochastic Models of type SCENARIOS DISCRETE.')

        self.scenarios[scenario_id] = Scenario(scenario_id, parent, probability, branch_period, self)

    def generate_deterministic_equivalent(self):

        print('Generating deterministic equivalent for the model {}...'.format(self.nominal_model.ModelName))
        deterministic_equivalent = gb.Model('Deterministic Equivalent')

        #####################
        # TWO STAGE VERSION #
        #####################

        x = {}  # x[stage,'NODE']
        # generate first-stage variables
        for i in range(len(self.c[1])):
            if (1,'ROOT') not in x.keys():
                x[1,'ROOT'] = {}
            x[1,'ROOT'][i] = deterministic_equivalent.addVar(obj=self.c[1][i],
                                                             lb=self.lb[1][i],
                                                             ub=self.ub[1][i],
                                                             vtype=self.vtype[1][i],
                                                             name='x_t{}_n{}_i{}'.format(1,'ROOT',i))
        deterministic_equivalent.update()

        # first stage constraints
        lhs = defaultdict(gb.LinExpr)
        for entry in self.A[1,1].items():
            row, column = entry[0]
            value = entry[1]
            lhs[row] = lhs[row] + value*x[1,'ROOT'][column]

        for row, value in enumerate(self.b[1]):
            deterministic_equivalent.addConstr(lhs[row], self.b_sense[1][row], self.b[1][row],
                                               name='stage1_row{}'.format(row))

        # second stage constraints
        for scn in self.scenarios:
            # Contribution to lhs from A[1,1]
            lhs_1 = defaultdict(gb.LinExpr)

            for entry in self.scenarios[scn].A[2,1].items():
                row, column = entry[0]
                value = entry[1]
                lhs_1[row] = lhs_1[row] + value*x[1,'ROOT'][column]

            # Contribution from A[2,2]
            lhs_2 = defaultdict(gb.LinExpr)

            # generate second stage variables
            if (2,scn) not in x.keys():
                x[2,scn] = {}
            for i in range(len(self.c[2])):
                x[2,scn][i] = deterministic_equivalent.addVar(obj=self.scenarios[scn].c[2][i]*self.scenarios[scn].probability,
                                                              lb=self.scenarios[scn].lb[2][i],
                                                              ub=self.scenarios[scn].ub[2][i],
                                                              vtype=self.scenarios[scn].vtype[2][i],
                                                              name='x_stage{}_node{}_i{}'.format(2,scn,i))
            deterministic_equivalent.update()

            for entry in self.scenarios[scn].A[2,2].items():
                row, column = entry[0]
                value = entry[1]
                lhs_2[row] = lhs_2[row] + value*x[2,scn][column]

            # add constraints to stochastic model
            for row, value in enumerate(self.b[2]):
                deterministic_equivalent.addConstr(lhs_1[row] + lhs_2[row],
                                                   self.scenarios[scn].b_sense[2][row],
                                                   self.scenarios[scn].b[2][row],
                                                   name='stage2_sc{}_row{}'.format(scn,row))
            deterministic_equivalent.update()

        # Set objective to minimize or maximize
        deterministic_equivalent.modelSense = self.nominal_model.modelSense

        return deterministic_equivalent

    def plot_scenario_tree(self):
        G = nx.DiGraph()
        G.add_node('ROOT')

        # build tree
        for scn in self.scenarios:
            genealogy = self.scenarios[scn].genealogy
            for i, node in enumerate(genealogy):
                if node != 'ROOT':
                    node_name = '{}_{}'.format(str(node), str(i))
                    if node_name not in G.nodes():
                        G.add_node(node_name)
                        if genealogy[i-1] == 'ROOT':
                            parent = 'ROOT'
                        else:
                            parent = '{}_{}'.format(genealogy[i-1],i-1)
                        G.add_edge(parent,node_name)

        #nx.draw_graphviz(G, with_labels=True, arrows=True, prog='dot', node_color='black')
        # node_position = nx.graphviz_layout(G, prog='dot')
        node_position = nx.pygraphviz_layout(G, prog='dot')
        label_position = copy.deepcopy(node_position)
        for p in label_position:
            # Now sure why networkx stores position as a tuple...
            label_position[p] = list(label_position[p])
            # label_position[p][1] += 15
            label_position[p] = tuple(label_position[p])
        # nx.draw_graphviz(G, position, with_labels=False, arrows=True, prog='dot', node_color='black')
        nx.draw(G, node_position, with_labels=False, arrows=True, node_color='red')
        nx.draw_networkx_labels(G, pos=label_position)
        # nx.write_dot(G,'tree_plot.dot')

    def scenario_tree_reduction(self):
        # TODO to be implemented
        pass

    # This 'workaround' is necessary to enable parallel computations.
    def __getstate__(self):
        d = dict(self.__dict__)
        for key in d.keys():
            if key == 'nominal_model':
                del d[key]
        return d

    # def __setstate__(self, d):
    #    self.__dict__.update(d)


class Scenario(object):
    """
    Stores information related to a single scenario. Manages (lazy) generation of single scenario data matrices.
    """
    def __init__(self, scenario_id, parent, probability, branch_period, stochastic_model):
        self.scenario_id = scenario_id
        self.parent = parent
        self.probability = probability
        self.branch_period = branch_period
        self.stochastic_model = stochastic_model  # a reference to the stochastic model to which the scenario belongs
        self.data_modification = {}

        # Problem's data modified according to the scenario
        self._A = {}
        self._b = {}
        self._c = {}
        self._lb = {}
        self._ub = {}
        self._b_sense = {}
        self._vtype = {}

        if parent != 'ROOT':
            # if its parent is not the root, we have to inherit the modifications to the model from its parent
            self.data_modification = copy.deepcopy(stochastic_model.scenarios[parent].data_modification)
            # the remainder of data_modification is carried over while parsing the sto. file in
            # StochasticModel._parse_stochastic_information, through usage of the self.add_data_modification method

        # figure out genealogy (i.e., sequence of nodes in the tree)
        if self.parent == 'ROOT':
            # if our parent is root, then the entire tree branch is our own!
            self.genealogy = [self.scenario_id for x in range(len(stochastic_model.periods))]
            self.genealogy[0] = 'ROOT'
        else:
            # otherwise we have to inherit genealogy from our parent, and append our own modification
            start = stochastic_model.periods.index(branch_period)
            self.genealogy = copy.deepcopy(stochastic_model.scenarios[parent].genealogy)
            self.genealogy[start:] = [self.scenario_id for x in range(len(self.genealogy) - start)]

    def add_data_modification(self, position, value):
        """
        :param position: ('ROW','COL') tuple of strings determine the position in the matrix where the modification
        takes place
        :param value: value of modification
        :return: None
        """
        self.data_modification[position] = value

    def __str__(self):
        return 'Scenario ID: {}, probability: {}, parent: {}'.format(self.scenario_id, self.probability, self.parent)

    def _update_model_matrices_according_to_scenario_modifications(self):
        # I cannot put this in init because when I run the init I still do not have all the modifications in place;
        # these come as I parse the .sto file.
        # TODO implement modifications to lb and ub (unclear how they would be specified in the .sto file, not
        # explained in the documentation).
        self._A = copy.deepcopy(self.stochastic_model.A)
        self._b = copy.deepcopy(self.stochastic_model.b)
        self._c = copy.deepcopy(self.stochastic_model.c)
        self._lb = copy.deepcopy(self.stochastic_model.lb)
        self._ub = copy.deepcopy(self.stochastic_model.ub)
        self._b_sense = copy.deepcopy(self.stochastic_model.b_sense)
        self._v_type = copy.deepcopy(self.stochastic_model.vtype)

        for position, value in self.data_modification.items():
            row_label = position[0]
            col_label = position[1]

            # Check if it's a RHS modification
            if col_label == 'RHS':
                # find stage i of the modification
                i = (key for key,value in self.stochastic_model.stage_constrs.items() if row_label in value).next()
                local_row = self.stochastic_model.stage_constrs[i].index(row_label)
                if self.stochastic_model.mode_of_modification == 'REPLACE':
                    self._b[i][local_row] = value
                elif self.stochastic_model.mode_of_modification == 'ADD':
                    self._b[i][local_row] = self.stochastic_model.b[i][local_row] + value
                else:
                    self._b[i][local_row] = self.stochastic_model.b[i][local_row] * value
            # Then check if it's a modification in the objective vector
            elif row_label == self.stochastic_model.objective_name:
                # find stage i of the modification
                i = (key for key,value in self.stochastic_model.stage_vars.items() if col_label in value).next()
                local_col = self.stochastic_model.stage_vars[i].index(col_label)
                if self.stochastic_model.mode_of_modification == 'REPLACE':
                    self._c[i][local_col] = value
                elif self.stochastic_model.mode_of_modification == 'ADD':
                    self._c[i][local_col] = self.stochastic_model.c[i][local_col] + value
                else:  # MULTIPLY case
                    self._c[i][local_col] = self.stochastic_model.c[i][local_col] * value
            # Otherwise it's a modification of A
            else:
                try:
                    # We have to figure out which (i,j) we have to query from A[i,j] to extract the appropriate
                    # sub matrix to modify
                    i = (key for key,value in self.stochastic_model.stage_constrs.items() if row_label in value).next()
                    j = (key for key,value in self.stochastic_model.stage_vars.items() if col_label in value).next()
                except StopIteration:
                    print('The .sto file is inconsistent with the .cor file: the given label for the row or column \n' \
                          'at which the scenario modification should take place is not in the list of labels \n' \
                          'determined in the .cor file.')
                local_row = self.stochastic_model.stage_constrs[i].index(row_label)
                local_col = self.stochastic_model.stage_vars[j].index(col_label)
                if self.stochastic_model.mode_of_modification == 'REPLACE':
                    self._A[i,j][local_row, local_col] = value
                elif self.stochastic_model.mode_of_modification == 'ADD':
                    self._A[i,j][local_row, local_col] = self.stochastic_model.A[i,j][local_row, local_col] + value
                else:
                    self._A[i,j][local_row, local_col] = self.stochastic_model.A[i,j][local_row, local_col] * value

    # Lazy generation of modified matrices, only upon request.
    @property
    def A(self):
        # if it hasn't been updated yet, we first update the local A matrix with the scenario modifications
        if not self._A:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._A

    @property
    def b(self):
        if not self._b:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._b

    @property
    def c(self):
        if not self._c:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._c

    @property
    def lb(self):
        if not self._lb:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._lb

    @property
    def ub(self):
        if not self._ub:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._ub

    @property
    def b_sense(self):
        if not self._b_sense:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._b_sense

    @property
    def vtype(self):
        if not self._v_type:
            self._update_model_matrices_according_to_scenario_modifications()
        return self._v_type
