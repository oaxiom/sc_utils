
'''

The C4 uses a double barcode system, I guess as it may help in barcode corection.

However, in practive barcode correction is around ~5% of total reads, and at best any
advantage from a Hamming distance on two barcodes is not much more than the Hamming distance from
a single barcode. So merging the two barcodes to make a whitelist seems reasonable
as it can now go through typical 10x-style analysis pipelines with no further modifications.

'''

oh = open('DNBelabC4barcodes_singles.txt', 'rt')
barcodes = []
for lin in oh:
    barcodes.append(lin.strip())
oh.close()

done = 0
oh = open('DNBelabC4barcodes_doubles.txt', 'wt')
for b1 in barcodes:
    for b2 in barcodes:
        oh.write(f"{b1}{b2}\n")
        done += 1

oh.close()
print(f"Did {done:,} barcodes")

