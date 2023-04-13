from sequence import Sequence
class MiRNA(Sequence):
    def __init__(self, h, s):
        Sequence.__init__(self, h, s)


    def getSeedSequence(self):
        return self._sequence[1:8]