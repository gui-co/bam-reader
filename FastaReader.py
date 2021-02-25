class FastaReader:
    _f = None

    def __init__(self, filename):
        try:
            self._f = open(filename, "r")
        except Exception as ex:
            raise Exception("[Error] unable to open fasta file: " + str(ex))

    def __del__(self):
        self._f.close()

    def getSequence(self, sequenceName):
        self._f.seek(0)
        while True:
            line = self._f.readline()
            if not line:
                raise Exception("[Error] sequence {} not found in fasta file"
                                .format(sequenceName))
            if line[0] == '>':
                name = line.split(' ')[0]
                if name[1:] == sequenceName or name == sequenceName:
                    break

        sequence = list()
        while True:
            line = self._f.readline()
            for i in line:
                if i in "ACTGN":
                    sequence.append(i)
                elif i in "\r\n\t ":
                    continue
                elif i == '>':
                    return sequence
                else:
                    raise Exception("[Error] sequence {} is corrupted")

