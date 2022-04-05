import sys

def compare_vcf(my_vcf, bcf_vcf):
        '''
        This function takes in two vcfs and compares them to see if they are the same.
        If they are, it returns True (vv)
        '''

        with open(my_vcf, 'rb') as f1, open(bcf_vcf, 'rb') as f2:
                i = 0
                identical = True
                for line1 in f1:
                        i += 1
                        for line2 in f2:
                                if i == 2:
                                        continue
                                elif line1 != line2:
                                        identical = False
                                        print("Line ", i, ": DIFFERENT")
                                        print("\tbcf:", line1)
                                        print("\tpython:", line2)
                                break
                print(identical)
                f1.close()
                f2.close()

compare_vcf(sys.argv[1], sys.argv[2])