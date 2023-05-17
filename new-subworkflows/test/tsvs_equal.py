def tsvs_equal(tsv1, tsv2):
    with open(tsv1, 'r') as f1, open(tsv2, 'r') as f2:
        for line1, line2 in zip(f1, f2):
            for item1, item2 in zip(line1.split(), line2.split()):
                if item1 != item2:
                    return False
    return True
