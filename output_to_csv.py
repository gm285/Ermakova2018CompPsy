import glob
import sys

def main():
  for output_filename in sorted(glob.iglob('runs/110.*/output')):
    output_number = output_filename.split('/')[1]
    with open(output_filename, 'r') as f:
      last_line = lastLine(f)
      items = last_line.split()
      if len(items) == 9 and items[0] == '25' and items[-1] == ']':
        sys.stdout.write(
            ','.join([output_number] + items[1:6]))
        sys.stdout.write('\n')
      else:
        sys.stderr.write('Skipped output %s, %s\n' % (output_filename, items))

def lastLine(f):
  line = None
  for line in f:
    pass
  return line

if __name__ == '__main__':
  main()
