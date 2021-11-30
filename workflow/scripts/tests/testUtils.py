#!/usr/bin/env python

import os
import sys
import pytest
from distutils import dir_util
sys.path.insert(1, os.path.join(sys.path[0], '..'))


@pytest.fixture
def datadir(tmpdir, request):
    """ Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely. datadir can be used just like tmpdir. """

    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


def validateOutput(captured, datadir, expectedStdout=None, expectedStderr=None):
    """ Compare observed capsys with expected outputs line by line """

    if expectedStdout:
        assertLine(datadir.join(expectedStdout), captured.out.splitlines(True))
    if expectedStderr:
        assertLine(datadir.join(expectedStderr), captured.err.splitlines(True))


def assertLine(expected, observed):
    """ Read expectedFile line by line and check same as observed output """

    with open(expected) as expected:
        for i, (observeLine, expectLine) in enumerate(zip(observed, expected)):
            assert observeLine == expectLine, f'Mismatch on line {i+1}.'
