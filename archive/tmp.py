import sys 


for line in sys.stdin:
    arg = line.strip()



    div = arg[1:3]

    file = arg.replace('_','-') + '.pdb'

    file = f'lib/rcsb/{div}/{file}'

    print(file)
