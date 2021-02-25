from FastaReader import FastaReader
from BamReader import BamReader

import numpy as np

import sys
import bisect
import time

class SequenceAlignment:
    _name = None
    _reads = None
    _ins = None
    _del = None
    _sequence = None
    _nreads = None
    _quality = None
    _quality_mean = None

    def __init__(self, name, sequence):
        length = len(sequence)
        self._name = name
        self._reads = [0] * length
        self._mis = [0] * length
        self._ins = [0] * length
        self._del = [0] * length
        self._sequence = sequence
        self._nreads = 0
        self._quality = list()
        self._qualitySum = [0] * length
        for _ in range(length):
            self._quality.append([])

    def setAlignment(self, cigar, qual, pos, seq, flag):
        if (flag & 0x04 or flag & 0x10 or flag & 0x100 or flag & 0x200
            or flag & 0x400 or flag & 0x800):
            return False
        self._nreads += 1
        refI = pos # indice on the reference
        alnI = 0   # indice on seq
        aa = 0
        for e in cigar:
            n = e[0]
            op = e[1]
            if op == 0: # M
                for _ in range(n):
                    if refI >= 0:
                        self._reads[refI] += 1
                        if qual[alnI] != 255:
                            bisect.insort(self._quality[refI], qual[alnI])
                            self._qualitySum[refI] += qual[alnI]
                        if seq[alnI] != self._sequence[refI]:
                            self._mis[refI] += 1
                    refI += 1
                    alnI += 1
            if op == 1: # I
                self._ins[refI] += 1
                if refI > 0:
                    self._ins[refI - 1] += 1
                alnI += n
            if op == 2: # D
                for _ in range(n):
                    if refI >= 0:
                        self._reads[refI] += 1
                        self._del[refI] += 1
                    refI += 1
            if op == 3: # N
                refI += n
            if op == 4: # S
                alnI += n
            if op == 5: # H
                pass
            if op == 6: # P
                pass
            if op == 7: # =
                for _ in range(n):
                    if refI >= 0:
                        self._reads[refI] += 1
                        self._quality[refI].append(qual[alnI])
                    refI += 1
                    alnI += 1
            if op == 8: # X
                for _ in range(n):
                    if refI >= 0:
                        self._reads[refI]
                        self._quality[refI].append(qual[alnI])
                        self._mis[refI]
                    refI += 1
                    alnI += 1
            aa += 1
        return True

    def write(self):
        t1 = time.clock_gettime(time.CLOCK_MONOTONIC)
        f = open("{}.csv".format(self._name), "w")
        for n in range(len(self._sequence)):
            if self._reads[n] != 0:
                if self._quality[n]:
                    q_mean = self._qualitySum[n] / len(self._quality[n])
                    mid = len(self._quality[n]) // 2
                    q_median = (self._quality[n][mid] + self._quality[n][~mid]) / 2
                    q_std = np.sqrt(np.mean((np.array(self._quality[n]) - q_mean)**2))
                else:
                    q_mean = 0
                    q_median = 0
                    q_std = 0
                f.write("{},{},{},{},{},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f}\n"
                        .format(self._name,
                                n + 1,
                                self._sequence[n],
                                '+',
                                self._reads[n], # cov
                                q_mean, # q_mean
                                q_median, # q_median
                                q_std, # q_std
                                self._mis[n] / self._reads[n],
                                self._ins[n] / self._reads[n],
                                self._del[n] / self._reads[n]))
        f.close()
        t2 = time.clock_gettime(time.CLOCK_MONOTONIC)
        print(t2 - t1)



f = FastaReader(sys.argv[2])
b = BamReader(sys.argv[1])

sequenceIndex = 0
sequenceName = b.getSequenceName(sequenceIndex)
sequenceLength = b.getSequenceLength(sequenceIndex)
try:
    sequence = f.getSequence(sequenceName)
    if sequenceLength != len(sequence):
        raise Exception("[Error] sequence length does not correspond between "
                        "BAM file and fasta")
except Exception as ex:
    print(str(ex))
    sys.exit()

print("Open sequence {}".format(sequenceName))
seqAlign = SequenceAlignment(sequenceName, sequence)

nReads = 0
while True:
    align = b.getNextAlignment()
    if align is None:
        print("Done. {} sequences and {} reads".format(sequenceIndex, nReads))
        break

    if align.getRefId() != sequenceIndex:
        print("Export alignments of sequence {} with {} reads"
              .format(sequenceName, nReads))
        seqAlign.write()
        sequenceIndex += 1
        sequenceName = b.getSequenceName(sequenceIndex)
        sequenceLength = b.getSequenceLength(sequenceIndex)
        if sequenceName is None:
            print("Done. {} sequences and {} reads"
                  .format(sequenceIndex, nReads))
            break
        try:
            sequence = f.getSequence(sequenceName)
            if sequenceLength != len(sequence):
                raise Exception("[Error] sequence length does not correspond "
                                "between BAM file and fasta")
        except Exception as ex:
            print(str(ex))
            break
        print("Open sequence {}".format(sequenceName))
        seqAlign = SequenceAlignment(sequenceName, sequence)
        nReads = 0

    done = seqAlign.setAlignment(align.getCigar(),
                                 align.getQual(),
                                 align.getPos(),
                                 align.getSeq(),
                                 align.getFlag())
    if done:
        nReads += 1
    if nReads % 100 == 0:
        print("    Add read {}".format(nReads))

