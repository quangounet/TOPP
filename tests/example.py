import unittest


class TestExample(unittest.TestCase):
    def setUp(self):
        self.x = 12

    def test_x(self):
        self.assertEqual(self.x, 12)

    def test_ok(self):
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
