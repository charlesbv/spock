import os.path
import time

file_path = 'test_wait_file.txt'
while not os.path.exists(file_path):
    time.sleep(1)
