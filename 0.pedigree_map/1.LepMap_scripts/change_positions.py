import sys
import os.path

### CHECKS ###
if len(sys.argv) != 3: #check input
    exit('Argument 1 should be the input map file, argument 2 the file with the positions to change')

if os.path.isfile(sys.argv[1]) == False | os.path.isfile(sys.argv[2]) == False : #check if file exists
    exit('input file(s) not found')

# read positions 
pos = open(sys.argv[2],'r')
poslist = []

for line in pos:
    poslist.append( int ( line.rstrip() ) )
pos.close()
# open map
map = open(sys.argv[1],'r')

lines = map.readlines()

map.close()

sum1 = 0
for element in poslist:
    sum1 += int(lines[element])
    lines[element] = '0\n'

print(sum1)

out = sys.argv[1] + ".replaced"
outf = open(out, 'w')
outf.writelines(lines)
outf.close()


