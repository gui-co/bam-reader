import struct
import zlib

class BgzfReader:
    _f = None
    _currentBlock = None
    _currentIndex = None

    def __init__(self, filename):
        try:
            self._f = open(filename, mode="rb")
        except:
            raise Exception("[Error] unable to open {}".format(filename))
        self._decompressNextBlock()

    def __del__(self):
        self._f.close()

    def _decompressNextBlock(self):
        self._currentIndex = 0
        self._currentBlock = None
        header = self._f.read(12)
        if not header:
            self._f.close()
            return
        if len(header) != 12:
            self._f.close()
            raise Exception("[Error] file {} is badly formated. Unable to read "
                            "BGZF block header".format(filename))

        id1, id2, cm, flg, mtime, xfl, os, xlen = \
                struct.unpack("<BBBBIBBH", header)

        if id1 != 31 or id2 != 139 or cm != 8 or flg != 4:
            raise Exception("[Error] file is badly formatted. Header does not "
                            "correspond to a BGZF block")

        subfields = self._f.read(xlen)
        if len(subfields) != xlen:
            self._f.close()
            raise Exception("[Error] file {} is badly formated. Unable to read "
                            "BGZF block header".format(filename))

        si1, si2, slen, bsize = struct.unpack("<BBHH", subfields)
        cdata = self._f.read(bsize - xlen - 19)
        crc32 = self._f.read(4)
        isize = self._f.read(4)
        try:
            data = zlib.decompress(cdata, wbits=-zlib.MAX_WBITS)
        except Exception as ex:
            self._f.close()
            raise Exception("[Error] file {} is badly formatted. Unable to "
                            "deflate the BGZF data: "
                            .format(filename) + str(ex))
        self._currentBlock = bytearray(data)

    def read(self, size):
        if (self._currentIndex + size) < len(self._currentBlock):
            beg = self._currentIndex
            end = self._currentIndex + size
            self._currentIndex = end
            return self._currentBlock[beg:end]
        else:
            buff = self._currentBlock[self._currentIndex:]
            self._decompressNextBlock()
            beg = 0
            end = size - len(buff)
            self._currentIndex = end
            return buff + self._currentBlock[0:end]
