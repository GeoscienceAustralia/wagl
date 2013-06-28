import unittest, pprint, os
from ULA3.utils import execute

class test_Subproc(unittest.TestCase):

	def print_result(self, result):
		print
		pprint.pprint(result)

	def test_true(self):
		result = execute('/bin/true')
		#self.print_result(result)
		self.assertEqual(result['returncode'], 0)
		self.assertEqual(result['stdout'], '')
		self.assertEqual(result['stderr'], '')

	def test_false(self):
		result = execute('/bin/false')
		#self.print_result(result)
		self.assertEqual(result['returncode'], 1)
		self.assertEqual(result['stdout'], '')
		self.assertEqual(result['stderr'], '')

	def test_ls(self):
		result = execute('ls -l')
		#self.print_result(result)
		self.assertEqual(result['returncode'], 0)
		self.assertNotEqual(result['stdout'], '')
		self.assertEqual(result['stderr'], '')

	def test_ls_wdir(self):
		tdir = os.getcwd()
		xdir = '/tmp'
		assert os.path.isdir(xdir)
		result = execute('ls -l', cwd=xdir)
		self.assertEqual(result['returncode'], 0)
		self.assertNotEqual(result['stdout'], '')
		self.assertEqual(result['stderr'], '')
		self.assertEqual(result['caller_wd'], tdir)
		self.assertNotEqual(result['stdout'], execute('ls -l')['stdout'])
