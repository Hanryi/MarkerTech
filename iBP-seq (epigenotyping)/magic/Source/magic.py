import sys

with open(sys.argv[1], 'r') as f1:
    with open(sys.argv[2], 'r') as f2:
        with open(sys.argv[3], 'w') as f3:
            f11 = f1.readlines()
            f22 = f2.readlines()
            f3.write('chr\tpos\tstrand\tcontext\n')
            for i in f11:
                i = i.strip('\n').split('\t')
                for l in f22:
                    l = l.strip('\n').split('\t')
                    if l[0] == i[0] and l[1] == i[1]:
                        f3.write(l[0] + '\t' + l[1] + '\t' + l[2] + '\t' + l[5] + '\n')
