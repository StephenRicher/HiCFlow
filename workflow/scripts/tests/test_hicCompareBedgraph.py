#!/usr/bin/env python

import os
import sys
import pytest
import filecmp
from testUtils import datadir, validateOutput
from hicCompareBedgraph import *

hicCompareBedgraph_params = (
    [('testMat1.homer', 'testMat1-all.exp', 'testMat1-up.exp',  'testMat1-down.exp')])
@pytest.mark.parametrize(
    'infile, expectAll, expectUp, expectDown', hicCompareBedgraph_params)
def test_hicCompareBedgraph(infile, expectAll, expectUp, expectDown, datadir):
    infile = datadir.join(infile)
    # Retrieve paths to expected outputs
    expectAll = datadir.join(expectAll)
    expectUp = datadir.join(expectUp)
    expectDown = datadir.join(expectDown)
    # Define paths to observed outputs
    obsAll = datadir.join('testMat1-all.obs')
    obsUp = datadir.join('testMat1-up.obs')
    obsDown = datadir.join('testMat1-down.obs')

    hicCompareBedgraph(infile, obsAll, 100, obsUp, obsDown)
    assert filecmp.cmp(expectAll, obsAll)
    assert filecmp.cmp(expectDown, obsDown)
    assert filecmp.cmp(expectUp, obsUp)
