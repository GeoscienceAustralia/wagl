#!/bin/env python

from __future__ import absolute_import
import unittest
from ULA3.fc import EndMemberFactory, EndMember

class TestEndMembers(unittest.TestCase):
    def test_getAllIds(self):
        ids = EndMemberFactory.getAllIds()
        assert ids is not None

    def test_getId(self):
        ids = sorted(EndMemberFactory.getAllIds())
        id = ids[0]
        assert id == '2009_08_10' 

    def test_getEndMember(self):
        ids = sorted(EndMemberFactory.getAllIds())
        endMember =  EndMemberFactory.getEndMemberById(ids[0])
        assert endMember is not None
        assert endMember.id == ids[0]
        assert endMember.sumWeight == 0.01 
        assert endMember.values.shape  == (56, 4)

    def test_get_bad_id(self):
        try:
            endMember =  EndMemberFactory.getEndMemberById('silly_key')
            fail("should have thrown exception")
        except KeyError:
            pass   # expected this



if __name__ == '__main__':
    unittest.main()    
