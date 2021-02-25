from FastaReader import FastaReader
from BamReader import BamReader

import numpy as np

import sys
import bisect
import time

class SequenceAlignment:
    _name = None
    _sequence = None
    _readsPlus = None
    _readsMinus = None
    _insPlus = None
    _insMinus = None
    _delPlus = None
    _delMinus = None
    _qualityPlus = None
    _qualityMinus = None
    _qualitySumPlus = None
    _qualitySumMinus = None

    def __init__(self, name, sequence):
        length = len(sequence)
        self._name = name
        self._sequence = sequence
        self._readsPlus = [0] * length
        self._readsMinus = [0] * length
        self._misPlus = [0] * length
        self._misMinus = [0] * length
        self._insPlus = [0] * length
        self._insMinus = [0] * length
        self._delPlus = [0] * length
        self._delMinus = [0] * length
        self._qualityPlus = list()
        self._qualityMinus = list()
        self._qualitySumPlus = [0] * length
        self._qualitySumMinus = [0] * length
        for _ in range(length):
            self._qualityPlus.append([])
            self._qualityMinus.append([])

    def setAlignment(self, cigar, qual, pos, seq, flag):
        if (flag & 0x004 or flag & 0x200 or flag & 0x400 or flag & 0x800):
            return False
        if (flag & 0x010):
            _reads = self._readsMinus
            _mis = self._misMinus
            _del = self._delMinus
            _ins = self._insMinus
            _quality = self._qualityMinus
            _qualitySum = self._qualitySumMinus
        else:
            _reads = self._readsPlus
            _mis = self._misPlus
            _del = self._delPlus
            _ins = self._insPlus
            _quality = self._qualityPlus
            _qualitySum = self._qualitySumPlus

        refI = pos # indice on the reference
        alnI = 0   # indice on seq
        aa = 0
        for e in cigar:
            n = e[0]
            op = e[1]
            if op == 0: # M
                for _ in range(n):
                    if refI >= 0:
                        _reads[refI] += 1
                        if qual[alnI] != 255:
                            bisect.insort(_quality[refI], qual[alnI])
                            _qualitySum[refI] += qual[alnI]
                        if seq[alnI] != self._sequence[refI]:
                            _mis[refI] += 1
                    refI += 1
                    alnI += 1
            if op == 1: # I
                _ins[refI] += 1
                if refI > 0:
                    _ins[refI - 1] += 1
                alnI += n
            if op == 2: # D
                for _ in range(n):
                    if refI >= 0:
                        _reads[refI] += 1
                        _del[refI] += 1
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
                        _reads[refI] += 1
                        _quality[refI].append(qual[alnI])
                    refI += 1
                    alnI += 1
            if op == 8: # X
                for _ in range(n):
                    if refI >= 0:
                        _reads[refI]
                        _quality[refI].append(qual[alnI])
                        _mis[refI]
                    refI += 1
                    alnI += 1
            aa += 1
        return True

    def write(self):
        t1 = time.clock_gettime(time.CLOCK_MONOTONIC)
        fPlus = open("{}_plus.csv".format(self._name), "w")
        fMinus = open("{}_minus.csv".format(self._name), "w")
        for n in range(len(self._sequence)):
            if self._readsPlus[n] != 0:
                if self._qualityPlus[n]:
                    qMeanPlus = self._qualitySumPlus[n] / len(self._qualityPlus[n])
                    mid = len(self._qualityPlus[n]) // 2
                    qMedianPlus = (self._qualityPlus[n][mid] + self._qualityPlus[n][~mid]) / 2
                    qStdPlus = np.sqrt(np.mean((np.array(self._qualityPlus[n]) - qMeanPlus)**2))
                else:
                    qMeanPlus = 0
                    qMedianPlus = 0
                    qStdPlus = 0
                fPlus.write("{},{},{},{},{},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f}\n"
                            .format(self._name,
                                    n + 1,
                                    self._sequence[n],
                                    '+',
                                    self._readsPlus[n],
                                    qMeanPlus,
                                    qMedianPlus,
                                    qStdPlus,
                                    self._misPlus[n] / self._readsPlus[n],
                                    self._insPlus[n] / self._readsPlus[n],
                                    self._delPlus[n] / self._readsPlus[n]))
            if self._readsMinus[n] != 0:
                if self._qualityMinus[n]:
                    qMeanMinus = self._qualitySumMinus[n] / len(self._qualityMinus[n])
                    mid = len(self._qualityMinus[n]) // 2
                    qMedianMinus = (self._qualityMinus[n][mid] + self._qualityMinus[n][~mid]) / 2
                    qStdMinus = np.sqrt(np.mean((np.array(self._qualityMinus[n]) - qMeanMinus)**2))
                else:
                    qMeanMinus = 0
                    qMedianMinus = 0
                    qStdMinus = 0
                fMinus.write("{},{},{},{},{},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f}\n"
                            .format(self._name,
                                    n + 1,
                                    self._sequence[n],
                                    '-',
                                    self._readsMinus[n],
                                    qMeanMinus,
                                    qMedianMinus,
                                    qStdMinus,
                                    self._misMinus[n] / self._readsMinus[n],
                                    self._insMinus[n] / self._readsMinus[n],
                                    self._delMinus[n] / self._readsMinus[n]))
        fPlus.close()
        fMinus.close()
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

