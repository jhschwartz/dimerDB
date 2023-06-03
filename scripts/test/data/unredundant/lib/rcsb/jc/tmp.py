import os

files = ['fake-a1-m1-cA2.pdb',
'fake-a1-m1-cA3.pdb',
'fake-a1-m1-cA6.pdb',
'fake-a1-m1-cAy.pdb']

for fn in files:
    with open('tmp', 'w') as tf, open(fn, 'r') as f:
        maxnum = -1
        for line in f:
            tf.write(line)
            try:
                n = int(line.split()[5])
            except IndexError:
                continue
            if n > maxnum:
                maxnum=n
        f.seek(0)
        for line in f:
            try:
                n = maxnum + int(line.split()[5])
            except IndexError:
                continue
            n = str(n)
            line = list(line)
            if len(n) == 3:
                line[23:26] = n
            elif len(n) == 2:
                line[24:26] = n
            elif len(n) == 1:
                line[25] = n
            tf.write(''.join(line))

    os.remove(fn)
    os.rename('tmp', fn)
