import os
from compress_bin2tar import bin2targz
import unittest

class TestBin2TarGz(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.test_dir = "test_directory"
        os.makedirs(self.test_dir, exist_ok=True)
        os.makedirs(os.path.join(self.test_dir, "subdir1", "outputs"), exist_ok=True)
        os.makedirs(os.path.join(self.test_dir, "subdir2", "outputs"), exist_ok=True)
        os.makedirs(os.path.join(self.test_dir, "subdir3", "outputs"), exist_ok=True)

        # Create dummy .bin files
        with open(os.path.join(self.test_dir, "subdir1", "outputs", "file1.bin"), "w") as f:
            f.write("Dummy file 1")
        with open(os.path.join(self.test_dir, "subdir2", "outputs", "file2.bin"), "w") as f:
            f.write("Dummy file 2")
        with open(os.path.join(self.test_dir, "subdir3", "outputs", "file3.txt"), "w") as f:
            f.write("Dummy file 3")
    
    def tearDown(self):
        # Remove the temporary directory after testing
        os.system("rm -rf {}".format(self.test_dir))

    def test_bin2targz(self):
        bin2targz(self.test_dir)

        # Check if the .bin files are compressed as .tar.gz files
        self.assertTrue(os.path.exists(os.path.join(self.test_dir, "subdir1", "outputs", "file1.tar.gz")))
        self.assertTrue(os.path.exists(os.path.join(self.test_dir, "subdir2", "outputs", "file2.tar.gz")))
        self.assertFalse(os.path.exists(os.path.join(self.test_dir, "subdir3", "outputs", "file3.tar.gz")))

        # Check if the original .bin files are removed
        self.assertFalse(os.path.exists(os.path.join(self.test_dir, "subdir1", "outputs", "file1.bin")))
        self.assertFalse(os.path.exists(os.path.join(self.test_dir, "subdir2", "outputs", "file2.bin")))
        self.assertTrue(os.path.exists(os.path.join(self.test_dir, "subdir3", "outputs", "file3.txt")))

        # Cleaning by removing the created directories
        self.tearDown()

if __name__ == "__main__":
    unittest.main()
