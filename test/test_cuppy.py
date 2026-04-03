from coinor.cuppy.cuttingPlanes import solve, MILPInstance, gomoryMixedIntegerCut

def test_cuppy():
    m = MILPInstance(module_name = 'coinor.cuppy.examples.MIP6')
    solve(m, whichCuts = [(gomoryMixedIntegerCut, {})],
          display = False, debug_print = True, use_cglp = False)
    assert(-0.9 >= m.lp.objectiveValue >= -1.1)

test_cuppy()