import os
import pytest
from smps.read import StochasticModel

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
TEST_PATH = THIS_PATH+'./test_datasets/'


@pytest.fixture
def small_sm_instance():
    return StochasticModel(TEST_PATH + 'sizes10')


def test_stochastic_model(small_sm_instance):
    sm = small_sm_instance

    # basic checks of the stochastic model
    assert len(sm.scenarios) == 10
    assert len(sm.periods) == 2
    assert sm.periods[0] == 'STAGE-1'
    assert sm.periods[1] == 'STAGE-2'
    assert max(sm.b[1]) == 200
    assert min(sm.b[1]) == 0
    # checks of scenario
    assert sm.scenarios['SCEN01'].scenario_id == 'SCEN01'
    assert sm.scenarios['SCEN01'].probability == 0.1
    assert sm.scenarios['SCEN01'].branch_period == 'STAGE-2'
    assert sm.scenarios['SCEN01'].parent == 'ROOT'

    # print sm.var_names
    # print sm.constr_names
    # sm.plot_scenario_tree()  # try the plotting function


def test_generate_det_eq(small_sm_instance):
    # Generate and solve the deterministic equivalent
    sm = small_sm_instance
    det_eq = sm.generate_deterministic_equivalent()
    # print('Solving...')
    det_eq.params.MIPGap = 0.01
    det_eq.optimize()
    assert det_eq.status == 2, 'sizes10 model should be solved to optimality without problems'