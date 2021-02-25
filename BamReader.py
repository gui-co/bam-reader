from BgzfReader import BgzfReader

import struct
import sys

class BamAlignment:
    _name = None
    _refId = None
    _pos = None
    _seq = None
    _cigar = None
    _qual = None

    def setName(self, name):
        self._name = name

    def getName(self):
        return self._name

    def setRefId(self, refId):
        self._refId = refId

    def getRefId(self):
        return self._refId

    def setPos(self, pos):
        self._pos = pos

    def getPos(self):
        return self._pos

    def setSeq(self, seq):
        self._seq = seq

    def getSeq(self):
        return self._seq

    def setCigar(self, cigar):
        self._cigar = cigar

    def getCigar(self):
        return self._cigar

    def setQual(self, qual):
        self._qual = qual

    def getQual(self):
        return self._qual


class BamReader:
    _reader = None
    _sequenceNames = None
    _sequenceLengths = None

    def __init__(self, filename):
        self._reader = BgzfReader(filename)
        self._sequenceNames = list()
        self._sequenceLengths = list()
        magic = self._reader.read(4)
        if len(magic) != 4 or struct.unpack("<4s", magic)[0] != b"BAM\x01":
            raise Exception("[Error] file {} badly formatted. This is not a "
                            "BAM file".format(filename))
        try:
            l_text = struct.unpack("<I", self._reader.read(4))[0]
            header = self._reader.read(l_text)

            n_ref = struct.unpack("<I", self._reader.read(4))[0]
            for _ in range(n_ref):
                l_name = struct.unpack("<I", self._reader.read(4))[0]
                name = struct.unpack("<{}s".format(l_name),
                                     self._reader.read(l_name))[0]
                l_ref = struct.unpack("<I", self._reader.read(4))[0]
                self._sequenceNames.append(name[0:-1].decode())
                self._sequenceLengths.append(l_ref)
        except:
            raise Exception("[Error] file {} badly formatted. Unable to parse "
                            "the header".format(filename))

    def getSequenceName(self, i):
        if i < len(self._sequenceNames):
            return self._sequenceNames[i]
        return None

    def getSequenceLength(self, i):
        if i < len(self._sequenceLengths):
            return self._sequenceLengths[i]
        return 0

    def getNextAlignment(self):
        blockSize = self._reader.read(4)
        if blockSize == 0:
            return None
        if len(blockSize) != 4:
            raise Exception("[Error] file {} badly formatted. Unable to read "
                            "next alignment".format(filename))
        blockSize = struct.unpack("<I", blockSize)[0]
        block = self._reader.read(blockSize)
        if len(block) != blockSize:
            raise Exception("[Error] file {} badly formatted. Unable to read "
                            "next alignment".format(filename))
        p = 0
        refId = struct.unpack_from("<I", block, p)[0]
        p += 4
        pos = struct.unpack_from("<i", block, p)[0]
        p += 4
        lReadName = struct.unpack_from("<B", block, p)[0]
        p += 1
        mapq = struct.unpack_from("<B", block, p)[0]
        p += 1
        binindex = struct.unpack_from("<H", block, p)[0]
        p += 2
        lCigar = struct.unpack_from("<H", block, p)[0]
        p += 2
        flag = struct.unpack_from("<H", block, p)[0]
        p += 2
        lSeq = struct.unpack_from("<I", block, p)[0]
        p += 4
        next_refID = struct.unpack_from("<i", block, p)[0]
        p += 4
        next_pos = struct.unpack_from("<i", block, p)[0]
        p += 4
        tlen = struct.unpack_from("<i", block, p)[0]
        p += 4
        readName = struct.unpack_from("<{}s".format(lReadName),
                                 block, p)[0]
        p += lReadName
        i = 0
        cigar = list()
        while i < lCigar:
            c = struct.unpack_from("<I", block, p)[0]
            n = int(c >> 4)
            op = int(0x0F & c)
            cigar.append((n, op))
            i += 1
            p += 4
        i = 0
        seq = list()
        while i < int((lSeq + 1) / 2):
            nuc = struct.unpack_from("<B", block, p)[0]
            n1 = int(nuc >> 4)
            n2 = int(nuc & 0x0F)
            for n in [n1, n2]:
                if n == 1:
                    seq.append('A')
                elif n == 2:
                    seq.append('C')
                elif n == 4:
                    seq.append('G')
                elif n == 8:
                    seq.append('T')
            i += 1
            p += 1
        qual = list()
        i = 0
        while i < lSeq:
            qual.append(struct.unpack_from("<b", block, p)[0])
            i += 1
            p += 1

        align = BamAlignment()
        align.setName(readName)
        align.setRefId(refId)
        align.setPos(pos)
        align.setSeq(seq)
        align.setCigar(cigar)
        align.setQual(qual)
        return align

