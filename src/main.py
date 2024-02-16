## python script to run main function 
## Add any argument parsing, function calling, here

import argparse
from functions import * 




## Parsing from command line, and running the script. 
if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  args = parser.parse_args()
