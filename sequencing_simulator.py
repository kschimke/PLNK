import sys
import os
import time
import glob

def commandline():
    import argparse
    parser = argparse.ArgumentParser(description='Reorients/demuxes/trims consensus reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--input','-i', type=str, help='directory containing fast5 files')
    parser.add_argument('--output','-o', type=str, help='path to directory to copy fast5s')

    args = parser.parse_args()
    return args

def main():
    args = commandline()

    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    file_list = glob.glob(inputFolder + '/*.fast5')
    file_list.sort(key=lambda x:os.path.getmtime(x))

    previous=os.path.getmtime(file_list[0])
    for f5_path in file_list:
        current=os.path.getmtime(f5_path)

        wait_time=current-previous
        print('waiting %f seconds to move file %s' % (wait_time,f5_path.split('/')[-1]))
        time.sleep(wait_time)
        os.system('scp %s %s' %(f5_path,outputFolder+fast5))

        previous=current

if __name__ == '__main__':
    main()
